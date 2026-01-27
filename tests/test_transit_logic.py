"""
Tests for transit detection: _wait_deg, _lst_deg, crossing logic, fire-time sanity.

Run: python -m pytest tests/ -v
Or:  python -m unittest tests.test_transit_logic -v
Or:  python tests/test_transit_logic.py  (from project root)
"""

import sys
from pathlib import Path

# Allow running as python tests/test_transit_logic.py (project root not on path)
_root = Path(__file__).resolve().parent.parent
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))

import unittest

from astropy.time import Time

# Import private helpers for testing
from celestial_musicbox.transit_scheduler import (
    _norm_360,
    _lst_deg,
    _diff_deg,
    _wait_deg,
    _lst_rate_deg_per_sec,
    _next_transit_unix,
    _stars_just_crossed,
    _time,
    SIDEREAL_DAY_S,
    LST_RATE_DEG_S,
    JUST_CROSSED_DEG,
)


class TestNorm360(unittest.TestCase):
    def test_zero(self):
        self.assertEqual(_norm_360(0.0), 0.0)

    def test_mid(self):
        self.assertAlmostEqual(_norm_360(180.0), 180.0)

    def test_wraps(self):
        self.assertAlmostEqual(_norm_360(360.0), 0.0)
        self.assertAlmostEqual(_norm_360(720.0), 0.0)
        self.assertAlmostEqual(_norm_360(-90.0 + 360.0), 270.0)


class TestWaitDeg(unittest.TestCase):
    """wait_deg(lst, ra) = degrees LST must advance to reach RA. [0,360). >180 = past."""

    def test_at_transit(self):
        self.assertAlmostEqual(_wait_deg(100.0, 100.0), 0.0)

    def test_before_transit(self):
        # LST=10, RA=20 → need to advance 10°
        self.assertAlmostEqual(_wait_deg(10.0, 20.0), 10.0)
        self.assertAlmostEqual(_wait_deg(0.0, 1.0), 1.0)

    def test_past_transit(self):
        # LST=20, RA=10 → past. Advance 350° to wrap to RA.
        w = _wait_deg(20.0, 10.0)
        self.assertGreater(w, 180.0)
        self.assertLess(w, 360.0)
        self.assertAlmostEqual(w, 350.0)

    def test_anti_transit(self):
        # LST = RA + 180 → exactly 180°
        self.assertAlmostEqual(_wait_deg(10.0, 190.0), 180.0)
        self.assertAlmostEqual(_wait_deg(0.0, 180.0), 180.0)


class TestDiffDeg(unittest.TestCase):
    def test_same(self):
        self.assertAlmostEqual(_diff_deg(100.0, 100.0), 0.0)

    def test_shortest_arc(self):
        self.assertAlmostEqual(_diff_deg(10.0, 350.0), 20.0)
        self.assertAlmostEqual(_diff_deg(350.0, 10.0), 20.0)

    def test_180(self):
        self.assertAlmostEqual(_diff_deg(0.0, 180.0), 180.0)


class TestLstAndRate(unittest.TestCase):
    def test_lst_range(self):
        t = _time(1700000000.0)
        lst = _lst_deg(t, -87.0)
        self.assertGreaterEqual(lst, 0.0)
        self.assertLess(lst, 360.0)

    def test_rate_positive(self):
        r = _lst_rate_deg_per_sec(1700000000.0, -87.0)
        self.assertGreater(r, 1e-6)
        self.assertLess(r, 0.01)

    def test_next_transit_at_transit(self):
        t = 1700000000.0
        lst = _lst_deg(_time(t), -87.0)
        n = _next_transit_unix(lst, -87.0, t, skip_immediate=False)
        self.assertAlmostEqual(n, t, delta=0.01)


class TestCrossingSimulation(unittest.TestCase):
    """
    Simulate the wait loop: we must only 'fire' when we've crossed (wait_d > 180)
    and |LST - RA| is small. We must NEVER fire at anti-transit (wait_d == 180).
    """

    def test_fire_only_when_just_past_not_anti_transit(self):
        # At anti-transit: wait_d = 180. We must not fire.
        lst = 100.0
        ra_anti = 280.0  # LST + 180
        w = _wait_deg(lst, ra_anti)
        self.assertAlmostEqual(w, 180.0)
        # Logic: break when wait_d > 180 (fire) vs >= 180 (would fire at anti-transit).
        self.assertFalse(w > 180.0)
        self.assertTrue(w >= 180.0)

    def test_just_past_has_wait_d_gt_180(self):
        lst = 100.0
        ra_just_past = 100.0 - 0.01  # LST - 0.01, we're just past
        w = _wait_deg(lst, ra_just_past)
        self.assertGreater(w, 180.0)
        self.assertLess(w, 360.0)
        self.assertAlmostEqual(_diff_deg(lst, ra_just_past), 0.01, delta=0.001)

    def test_before_transit_has_wait_d_lt_180(self):
        lst = 100.0
        ra_before = 100.5
        w = _wait_deg(lst, ra_before)
        self.assertLess(w, 180.0)
        self.assertGreater(w, 0.0)


class TestFireTimeVerification(unittest.TestCase):
    """At fire time we expect |LST - RA| small. If diff > 1° we should not fire."""

    def test_diff_at_just_past(self):
        lst = 100.0
        for offset_deg in (0.0, 0.005, 0.5, 1.0):
            ra = (lst - offset_deg + 360.0) % 360.0
            d = _diff_deg(lst, ra)
            self.assertAlmostEqual(d, offset_deg, delta=0.001)

    def test_diff_at_anti_transit(self):
        lst = 100.0
        ra = 280.0
        self.assertAlmostEqual(_diff_deg(lst, ra), 180.0)


class TestStarsJustCrossed(unittest.TestCase):
    """Poll-based 'just crossed' search: only stars with wait_d > 180 and >= 360 - JUST_CROSSED_DEG."""

    def test_just_crossed_returns_star_within_window(self):
        stars = [
            {"name": "A", "ra_deg": 100.0, "dec_deg": 0.0},
            {"name": "B", "ra_deg": 200.0, "dec_deg": 0.0},
        ]
        lst = 100.01  # 0.01° past A's meridian
        cooldown: dict[str, float] = {}
        got = _stars_just_crossed(stars, 1.0, lst, cooldown, 0.0)
        self.assertEqual(len(got), 1)
        self.assertEqual(got[0][0]["name"], "A")
        self.assertLess(got[0][1], 0.1)

    def test_before_meridian_not_returned(self):
        stars = [{"name": "A", "ra_deg": 100.0, "dec_deg": 0.0}]
        lst = 99.9  # 0.1° before A
        got = _stars_just_crossed(stars, 1.0, lst, {}, 0.0)
        self.assertEqual(len(got), 0)

    def test_cooldown_excludes_star(self):
        stars = [{"name": "A", "ra_deg": 100.0, "dec_deg": 0.0}]
        lst = 100.01
        cooldown = {"A": 0.0}  # just fired
        got = _stars_just_crossed(stars, 1.0, lst, cooldown, 1.0)  # 1s later, still on cooldown
        self.assertEqual(len(got), 0)

if __name__ == "__main__":
    unittest.main()
