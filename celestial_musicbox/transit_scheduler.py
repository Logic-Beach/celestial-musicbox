"""
Meridian transit scheduler.

Stars transit when LST = RA. We poll until we observe that, then fire (MIDI, optionally Stellarium).
"""

from __future__ import annotations

import heapq
import itertools
import json
import math
import sys
import time
from pathlib import Path
from typing import Optional

from astropy import units as u
from astropy.time import Time

from . import visualizer as viz
from .midi_sender import open_mido_output, send_dyads
from .stellarium_client import slew_to
from .star_to_midi import star_to_dyads

SIDEREAL_DAY_S = 86164.0905
LST_RATE_DEG_S = 360.0 / SIDEREAL_DAY_S
AT_TRANSIT_DEG = 0.01  # fire when |LST - RA| <= this (36 arcsec)


def _unix() -> float:
    return time.time()


def _time(t: float) -> Time:
    return Time(t, format="unix")


def _norm_360(x: float) -> float:
    """Angle in [0, 360). Clamp 360.0 -> 0.0 if float rounding produces it."""
    v = float(x) % 360.0
    return 0.0 if v >= 360.0 else v


def _lst_deg(t: Time, lon_deg: float) -> float:
    """Local sidereal time in degrees [0, 360). East positive. Uses UT1 when IERS available."""
    try:
        t_use = t.ut1
    except Exception:
        t_use = t
    h = t_use.sidereal_time("apparent", longitude=lon_deg * u.deg).hour
    return _norm_360(h * 15.0)


def _diff_deg(lst: float, ra: float) -> float:
    """Shortest angular distance |LST - RA| in [0, 180] degrees."""
    d = (lst - ra + 180.0) % 360.0 - 180.0
    return abs(d)


def _wait_deg(lst: float, ra: float) -> float:
    """Degrees LST must advance to reach RA. In [0, 360). >= 180 means past meridian."""
    return _norm_360(ra - lst + 360.0)


def _lst_rate_deg_per_sec(unix: float, lon_deg: float) -> float:
    """LST rate (deg/s) at unix, from astropy. Fallback: fixed sidereal rate."""
    try:
        t0 = _time(unix)
        t1 = _time(unix + 1.0)
        lst0 = _lst_deg(t0, lon_deg)
        lst1 = _lst_deg(t1, lon_deg)
        d = _norm_360(lst1 - lst0 + 360.0)
        if 1e-9 < d < 360.0 - 1e-9:
            return d
    except Exception:
        pass
    return LST_RATE_DEG_S


def _next_transit_unix(ra_deg: float, lon_deg: float, after_unix: float, *, skip_immediate: bool = False) -> float:
    """Unix time of next transit (LST = RA) at or after after_unix."""
    t = _time(after_unix)
    lst = _lst_deg(t, lon_deg)
    w = _wait_deg(lst, ra_deg)
    if w < AT_TRANSIT_DEG:
        if skip_immediate:
            w = 360.0
        else:
            return after_unix
    rate = _lst_rate_deg_per_sec(after_unix, lon_deg)
    return after_unix + w / rate


def _altitude_at_transit(dec_deg: float, lat_deg: float) -> float:
    x = math.radians(dec_deg - lat_deg)
    return math.degrees(math.asin(max(-1.0, min(1.0, math.cos(x)))))


def _star_rec(star: dict, ra_scale: float, lat_deg: float) -> dict:
    ra = float(star["ra_deg"]) * ra_scale
    dec = float(star["dec_deg"])
    rec = {"name": star["name"], "vmag": star.get("vmag", 5.0), "ra_deg": ra, "dec_deg": dec}
    if star.get("spectral"):
        rec["spectral"] = star["spectral"]
    if star.get("mass") is not None:
        rec["mass"] = star["mass"]
    else:
        rec["altitude"] = _altitude_at_transit(dec, lat_deg)
    rec["distance_ly"] = star.get("distance_ly") or star.get("distance")
    if star.get("bv") is not None:
        rec["bv"] = star["bv"]
    return rec


def _stellarium_candidates(rec: dict, star: dict) -> list[str]:
    out = [rec["name"]]
    for k, fmt in [("hip", "HIP {}"), ("hd", "HD {}"), ("hr", "HR {}")]:
        v = star.get(k)
        if v is not None:
            s = fmt.format(v)
            if s not in out:
                out.append(s)
    return out


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
        raise SystemExit("star_catalog.json must be a JSON array.")

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

    ras = [float(s["ra_deg"]) for s in stars]
    ra_scale = 15.0 if (ras and max(ras) <= 24.0) else 1.0
    if ra_scale != 1.0:
        print("Note: catalog ra_deg in hours; converting to degrees.\n", file=sys.stderr)

    now = _unix()
    tie = itertools.count()
    heap: list[tuple[float, int, dict]] = []
    for s in stars:
        ra = float(s["ra_deg"]) * ra_scale
        t = _next_transit_unix(ra, lon_deg, now, skip_immediate=False)
        heapq.heappush(heap, (t, next(tie), s))

    port = open_mido_output(midi_port_name)

    if not quiet and heap:
        t0, _, s0 = heap[0]
        d = max(0.0, t0 - _unix())
        h, r = int(d) // 3600, int(d) % 3600
        m, sec = r // 60, r % 60
        ra0 = float(s0["ra_deg"]) * ra_scale
        lst0 = _lst_deg(_time(now), lon_deg)
        msg = f"Scheduled {len(stars)} stars. Next: {s0['name']} in {h}:{m:02d}:{sec:02d}\n"
        if d > 1800:
            msg += f"  LST≈{lst0/15:.1f}h  next star RA≈{ra0/15:.1f}h\n"
        print(msg, flush=True)

    while True:
        t_unix, _, star = heapq.heappop(heap)
        ra = float(star["ra_deg"]) * ra_scale
        rec = _star_rec(star, ra_scale, lat_deg)
        dyads = star_to_dyads(rec, supplement)

        # Wait until the star has *crossed* the meridian (LST > RA). We only fire/slew
        # when we observe that crossing, not when we're still approaching.
        skipped = False
        waiting = False
        while True:
            u = _unix()
            lst = _lst_deg(_time(u), lon_deg)
            wait_d = _wait_deg(lst, ra)

            if wait_d >= 180.0:
                if not waiting:
                    t_next = _next_transit_unix(ra, lon_deg, u, skip_immediate=True)
                    heapq.heappush(heap, (t_next, next(tie), star))
                    skipped = True
                    break
                if wait_d <= 180.0:
                    t_next = _next_transit_unix(ra, lon_deg, u, skip_immediate=True)
                    heapq.heappush(heap, (t_next, next(tie), star))
                    skipped = True
                    break
                break
            waiting = True

            rate = _lst_rate_deg_per_sec(u, lon_deg)
            wait_s = wait_d / rate
            if wait_s > 600.0:
                t_next = _next_transit_unix(ra, lon_deg, u, skip_immediate=True)
                heapq.heappush(heap, (t_next, next(tie), star))
                skipped = True
                break

            chunk = min(0.5, max(0.05, wait_s))
            if not quiet and sys.stdout.isatty():
                print(viz.format_next(star["name"], wait_s), end="\r", flush=True)
            time.sleep(chunk)

        if skipped:
            continue
        if not quiet and sys.stdout.isatty():
            print(" " * 60, end="\r", flush=True)

        u_now = _unix()
        lst_now = _lst_deg(_time(u_now), lon_deg)
        rec["lst_deg"] = lst_now
        diff_deg = _diff_deg(lst_now, ra)
        if diff_deg > 1.0:
            sys.stderr.write(
                f"Transit skip: |LST−RA| = {diff_deg:.3f}° for {star['name']} — too far off meridian, rescheduling\n"
            )
            t_next = _next_transit_unix(ra, lon_deg, u_now, skip_immediate=True)
            heapq.heappush(heap, (t_next, next(tie), star))
            continue
        if not quiet:
            upcoming = [(s["name"], t - u_now) for t, _, s in heapq.nsmallest(50, heap) if t > u_now][:2]
            viz.print_transit(rec, dyads, upcoming=upcoming)

        if stellarium_url:
            slew_to(
                rec["ra_deg"],
                rec["dec_deg"],
                base_url=stellarium_url,
                target_candidates=_stellarium_candidates(rec, star),
            )

        send_dyads(port, dyads, note_duration=note_duration)

        t_next = _next_transit_unix(ra, lon_deg, _unix() + SIDEREAL_DAY_S, skip_immediate=True)
        heapq.heappush(heap, (t_next, next(tie), star))
