"""
Meridian transit scheduler: computes LST = RA and fires at transit times.

Uses a precomputed star catalog (ra_deg, dec_deg, ...) and the observer's lon/lat
to compute when each star crosses the meridian (LST = RA). Fires at those times.
"""

import heapq
import itertools
import json
import math
import time
from pathlib import Path
from typing import Optional

from astropy import units as u
from astropy.time import Time

from .midi_sender import open_mido_output, send_dyads

# Minimum seconds between stars when we're "catching up" (transit was in the past)
_CATCH_UP_SPACING = 2.0
from .stellarium_client import slew_to
from .star_to_midi import star_to_dyads
from . import visualizer as viz

# Fallback sidereal day in UT seconds (mean) if we can't derive rate from astropy
_SIDEREAL_DAY_UTC = 86164.0905


def _lst_deg(t: Time, lon_deg: float) -> float:
    """Local sidereal time in degrees [0, 360). Longitude: East positive (e.g. -120 for 120°W).

    Uses UT1 when IERS is available so LST follows Earth rotation; falls back to t as-is (UTC)
    if UT1 conversion fails (e.g. no IERS). Stellarium and system clock should both be in sync
    with real UT for best accuracy.
    """
    try:
        t_use = t.ut1
    except Exception:
        t_use = t
    h = t_use.sidereal_time("apparent", longitude=lon_deg * u.deg).hour
    return (h * 15.0) % 360.0


def _lst_rate_deg_per_sec(t: Time, lon_deg: float) -> float:
    """Rate of change of LST (deg/s) at time t, from astropy's apparent sidereal time."""
    lst0 = _lst_deg(t, lon_deg)
    lst1 = _lst_deg(t + 1.0 * u.s, lon_deg)
    d = lst1 - lst0
    if d < 0:
        d += 360.0
    return d


def next_transit_utc(ra_deg: float, lon_deg: float, after: Time) -> Time:
    """Next UTC time when LST = RA (mod 360) at or after `after`."""
    lst = _lst_deg(after, lon_deg)
    wait_deg = (ra_deg - lst + 360.0) % 360.0
    if wait_deg < 0.05:
        wait_deg = 0.0  # at transit now; fire immediately (was 360, which pushed it 24h)
    rate = _lst_rate_deg_per_sec(after, lon_deg)
    if rate < 1e-9:
        rate = 360.0 / _SIDEREAL_DAY_UTC
    wait_sec = wait_deg / rate
    return after + (wait_sec * u.s)


def altitude_at_transit_deg(dec_deg: float, lat_deg: float) -> float:
    """Altitude (degrees) when star is on the meridian: sin(alt)=cos(dec-lat)."""
    x = math.radians(dec_deg - lat_deg)
    cos_x = max(-1.0, min(1.0, math.cos(x)))
    return math.degrees(math.asin(cos_x))


def run_scheduler(
    catalog_path: Path,
    supplement_path: Path,
    lon_deg: float,
    lat_deg: float,
    midi_port_name: Optional[str] = None,
    quiet: bool = False,
    stellarium_url: Optional[str] = None,
    note_duration: Optional[float] = None,
) -> None:
    supplement: dict = {}
    if supplement_path.is_file():
        try:
            supplement = json.loads(supplement_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            pass

    raw = json.loads(catalog_path.read_text(encoding="utf-8"))
    if not isinstance(raw, list):
        raise SystemExit("star_catalog.json must be a JSON array of star objects.")

    # Filter: dec such that star can be above horizon at transit (alt > 0 => |dec-lat| < 90)
    lo, hi = lat_deg - 90.0, lat_deg + 90.0
    stars = [
        s
        for s in raw
        if isinstance(s, dict)
        and s.get("name")
        and s.get("ra_deg") is not None
        and s.get("dec_deg") is not None
        and lo <= float(s["dec_deg"]) <= hi
    ]

    # HYG and our build output ra_deg in degrees (0–360). Old catalogs may have ra in HOURS (0–24).
    # If max ra_deg <= 24, assume hours and scale so LST=RA math is correct.
    ras = [float(s["ra_deg"]) for s in stars]
    ra_scale = 15.0 if (ras and max(ras) <= 24.0) else 1.0
    if ra_scale != 1.0:
        import sys as _sys
        print("Note: catalog ra_deg looks like hours (0–24); converting to degrees. Rebuild with build_star_catalog.py for correct ra_deg.\n", file=_sys.stderr)

    # Build heap using now taken immediately before the loop (so LST matches).
    now = Time.now()
    lst_deg = _lst_deg(now, lon_deg)
    tie = itertools.count()
    heap: list[tuple[float, int, dict]] = []
    # Each star: next_transit_utc(ra, lon, now). Min‑heap by t_unix → soonest transit.
    for s in stars:
        ra = float(s["ra_deg"]) * ra_scale
        t = next_transit_utc(ra, lon_deg, now)
        heapq.heappush(heap, (t.to_value("unix"), next(tie), s))

    port = open_mido_output(midi_port_name)

    if not quiet and heap:
        t0, _, s0 = heap[0]
        d = max(0.0, t0 - time.time())
        h, r = int(d) // 3600, int(d) % 3600
        m, s = r // 60, r % 60
        ra0 = float(s0["ra_deg"]) * ra_scale
        msg = f"Scheduled {len(stars)} stars. Next transit: {s0['name']} in {h}:{m:02d}:{s:02d} (soonest of all).\n"
        # If >30 min, show LST and star RA so user can sanity‑check (e.g. at 22h LST, next should be RA≈22h, not 0h).
        if d > 1800:
            msg += f"  (LST≈{lst_deg/15:.1f}h, this star RA≈{ra0/15:.1f}h — if your meridian is at 22h, stars near 22h RA should come first; check --lon/--lat and that catalog has ra_deg in degrees)\n"
        print(msg, flush=True)

    while True:
        t_unix, _, star = heapq.heappop(heap)
        now_unix = time.time()
        wait = t_unix - now_unix

        # Build dict for star_to_dyads and visualizer: name, vmag, spectral?, mass?, altitude?, distance-ly?, ra_deg, dec_deg
        ra = float(star["ra_deg"]) * ra_scale
        dec = float(star["dec_deg"])
        rec = {
            "name": star["name"],
            "vmag": star.get("vmag", 5.0),
            "ra_deg": ra,
            "dec_deg": dec,
        }
        if star.get("spectral"):
            rec["spectral"] = star["spectral"]
        if star.get("mass") is not None:
            rec["mass"] = star["mass"]
        else:
            rec["altitude"] = altitude_at_transit_deg(dec, lat_deg)
        if star.get("distance_ly") is not None:
            rec["distance_ly"] = star["distance_ly"]
        elif star.get("distance") is not None:
            rec["distance_ly"] = star["distance"]

        dyads = star_to_dyads(rec, supplement)

        if wait < 0:
            time.sleep(_CATCH_UP_SPACING)  # min spacing when catching up (avoids rapid-fire)

        if not quiet:
            if wait > 0:
                viz.countdown(wait, star["name"])  # wait for this transit
            # Current LST at display time for sanity-check (RA and LST in h should match at transit).
            rec["lst_deg"] = _lst_deg(Time(time.time(), format="unix"), lon_deg)
            # Upcoming = next 2 soonest with future transit times; skip past/"now" so we show real countdowns.
            now_t = time.time()
            upcoming = [(s["name"], t - now_t) for (t, _, s) in heapq.nsmallest(50, heap) if t > now_t][:2]
            viz.print_transit(rec, dyads, upcoming=upcoming)
        elif wait > 0:
            time.sleep(wait)

        if stellarium_url:
            slew_to(rec["ra_deg"], rec["dec_deg"], base_url=stellarium_url, target=rec["name"])

        send_dyads(port, dyads, note_duration=note_duration)

        # Schedule next transit for this star (~24h later)
        ra = float(star["ra_deg"]) * ra_scale
        t_next = next_transit_utc(ra, lon_deg, Time(t_unix, format="unix") + (1.0 * u.s))
        heapq.heappush(heap, (t_next.to_value("unix"), next(tie), star))
