"""
Map star properties to MIDI: one chord per star.
- Pitch 1: B–V (blue→high, red→low).
- Pitch 2: distance (near→high, far→low).
- Velocity (both): magnitude (brighter→louder).
"""

import math

NOTE_LO, NOTE_HI = 0, 60


def _clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def _lerp(t: float, a: float, b: float) -> float:
    return a + t * (b - a)


def _to_note(t: float) -> int:
    return int(round(_lerp(_clamp(t, 0, 1), NOTE_LO, NOTE_HI)))


def bv_to_note(bv: float) -> int:
    """B–V → pitch. Blue (≈-0.4) → high, red (≈+2) → low."""
    t = (2.0 - _clamp(bv, -0.4, 2.0)) / 2.4
    return _to_note(t)


def magnitude_to_velocity(mag: float, vlo: int = 40, vhi: int = 115) -> int:
    """Brighter (lower mag) → higher velocity."""
    t = (_clamp(mag, -1.5, 8) + 1.5) / 9.5
    return int(round(_lerp(1.0 - t, vlo, vhi)))


def distance_to_note(distance_ly: float) -> int:
    """Distance (ly) → pitch. Near → high, far → low. Log scale."""
    d = distance_ly if distance_ly > 0 else 10.0
    x = math.log10(_clamp(d, 1, 15000))
    t = 1.0 - x / 4.2
    return _to_note(t)


def star_to_dyads(rec: dict) -> list[tuple[int, int]]:
    """
    One chord per star: two notes, same velocity.
    rec must have bv, vmag. distance_ly/distance optional (default 10 ly).
    """
    pitch1 = bv_to_note(float(rec["bv"]))
    vel = magnitude_to_velocity(float(rec["vmag"]))
    d = rec.get("distance_ly") or rec.get("distance") or 10.0
    pitch2 = distance_to_note(max(1.0, float(d)))
    return [(pitch1, vel), (pitch2, vel)]
