"""
Meridian transit trigger: poll-based "who's at meridian now?".

No heap scheduling. Every POLL_INTERVAL we get LST (from Stellarium if enabled, else
--lon + system time), search the catalog for stars that just crossed (wait_d in (180°, 180°+ε)),
pick the closest, fire (MIDI / Stellarium), cooldown that star for ~1 sidereal day. Repeat.

Using Stellarium's time and location when enabled keeps us aligned with Stellarium's sky.
"""

from __future__ import annotations

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
from .stellarium_client import slew_to, _get_status as get_status
from .star_to_midi import star_to_dyads

SIDEREAL_DAY_S = 86164.0905
LST_RATE_DEG_S = 360.0 / SIDEREAL_DAY_S
POLL_INTERVAL_S = 0.5
# Only consider stars that crossed in the last JUST_CROSSED_DEG (e.g. ~2 min of arc).
JUST_CROSSED_DEG = 0.5
COOLDOWN_FRAC = 0.98  # don't re-fire same star until 98% of sidereal day has passed


def _unix() -> float:
    return time.time()


def _time(t: float) -> Time:
    return Time(t, format="unix")


def _norm_360(x: float) -> float:
    v = float(x) % 360.0
    return 0.0 if v >= 360.0 else v


def _lst_deg(t: Time, lon_deg: float) -> float:
    try:
        t_use = t.ut1
    except Exception:
        t_use = t
    h = t_use.sidereal_time("apparent", longitude=lon_deg * u.deg).hour
    return _norm_360(h * 15.0)


def _lst_from_jd(jd: float, lon_deg: float) -> float:
    """LST in degrees from Julian day and longitude."""
    t = Time(jd, format="jd")
    return _lst_deg(t, lon_deg)


def _diff_deg(lst: float, ra: float) -> float:
    d = (lst - ra + 180.0) % 360.0 - 180.0
    return abs(d)


def _wait_deg(lst: float, ra: float) -> float:
    return _norm_360(ra - lst + 360.0)


def _lst_rate_deg_per_sec(unix: float, lon_deg: float) -> float:
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
    """Unix time of next transit (LST = RA) at or after after_unix. Kept for tests / upcoming."""
    t = _time(after_unix)
    lst = _lst_deg(t, lon_deg)
    w = _wait_deg(lst, ra_deg)
    if w < 0.01:
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


def _normalize_name(s: str) -> str:
    if not s or not isinstance(s, str):
        return ""
    return " ".join(str(s).split())


def _stellarium_candidates(rec: dict, star: dict) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()

    def add(c: str) -> None:
        if not c or c in seen:
            return
        seen.add(c)
        out.append(c)

    for k, fmt in [("hip", "HIP {}"), ("hd", "HD {}"), ("hr", "HR {}")]:
        v = star.get(k)
        if v is not None:
            add(fmt.format(v))
    name = rec.get("name")
    if name:
        add(_normalize_name(name))
    return out


def _log(msg: str, verbose: bool) -> None:
    if verbose:
        print(f"[cmb] {msg}", file=sys.stderr, flush=True)


def _get_lst_and_lon(
    stellarium_url: Optional[str],
    lon_deg: float,
    verbose: bool,
) -> tuple[float, float]:
    """Return (lst_deg, lon_deg). Always use caller's lon_deg (full precision).
    When stellarium_url set, use Stellarium's time (jday) only; never overwrite lon.
    """
    if stellarium_url:
        st = get_status(stellarium_url)
        if st:
            tinfo = st.get("time")
            if isinstance(tinfo, dict):
                jd = tinfo.get("jday")
                if jd is not None:
                    try:
                        lst = _lst_from_jd(float(jd), lon_deg)
                        _log(f"LST from Stellarium jd={jd} lon={lon_deg:.8f}° → lst={lst:.8f}°", verbose)
                        return (lst, lon_deg)
                    except (TypeError, ValueError):
                        pass
            _log("Stellarium status missing jday; using --lon and system time", verbose)
    lst = _lst_deg(_time(_unix()), lon_deg)
    return (lst, lon_deg)


def _stars_just_crossed(
    stars: list[dict],
    ra_scale: float,
    lst: float,
    cooldown: dict[str, float],
    now: float,
) -> list[tuple[dict, float]]:
    """Stars that just crossed: wait_d > 180 (past) and wait_d >= 360 - JUST_CROSSED_DEG (within ε of meridian).
    Not on cooldown. Return (star, diff) sorted by diff asc (closest first).
    """
    cooldown_s = COOLDOWN_FRAC * SIDEREAL_DAY_S
    lo = 360.0 - JUST_CROSSED_DEG  # e.g. 359.5: only 0–0.5° past
    candidates: list[tuple[dict, float]] = []
    for s in stars:
        name = s.get("name") or ""
        if name and (now - cooldown.get(name, -1e9)) < cooldown_s:
            continue
        ra = float(s["ra_deg"]) * ra_scale
        wait_d = _wait_deg(lst, ra)
        if wait_d <= 180.0 or wait_d < lo:
            continue
        diff = _diff_deg(lst, ra)
        candidates.append((s, diff))
    candidates.sort(key=lambda x: x[1])
    return candidates


def _upcoming(stars: list[dict], ra_scale: float, lst: float, lon_deg: float, now: float, n: int = 2):
    """Next n stars to cross (smallest wait_d in 0..180). Return [(name, seconds)]."""
    rate = _lst_rate_deg_per_sec(now, lon_deg)
    cand: list[tuple[dict, float, float]] = []
    for s in stars:
        ra = float(s["ra_deg"]) * ra_scale
        wait_d = _wait_deg(lst, ra)
        if wait_d >= 180.0:
            continue
        sec = wait_d / rate
        cand.append((s, wait_d, sec))
    cand.sort(key=lambda x: x[1])
    out = []
    for s, _, sec in cand[:n]:
        out.append((s.get("name") or "?", sec))
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
    verbose: bool = False,
) -> None:
    _log(f"loading catalog {catalog_path}", verbose)
    _log(f"lon={lon_deg:.8f}° lat={lat_deg:.8f}° (full precision preserved)", verbose)
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
    _log(f"catalog: {len(raw)} total, {len(stars)} visible (dec {lo:.0f}°..{hi:.0f}°)", verbose)

    ras = [float(s["ra_deg"]) for s in stars]
    ra_scale = 15.0 if (ras and max(ras) <= 24.0) else 1.0
    if ra_scale != 1.0:
        print("Note: catalog ra_deg in hours; converting to degrees.\n", file=sys.stderr)

    port = open_mido_output(midi_port_name, verbose=verbose)
    cooldown: dict[str, float] = {}
    last_upcoming = 0.0

    if not quiet:
        lst, lon = _get_lst_and_lon(stellarium_url, lon_deg, verbose)
        up = _upcoming(stars, ra_scale, lst, lon, _unix(), 2)
        msg = f"Poll-based meridian trigger. {len(stars)} stars. LST≈{lst/15:.2f}h.\n"
        if up:
            msg += "  Next: " + ", ".join(f"{n} in {s:.0f}s" for n, s in up) + "\n"
        print(msg, flush=True)

    while True:
        now = _unix()
        lst, lon = _get_lst_and_lon(stellarium_url, lon_deg, verbose)
        just_crossed = _stars_just_crossed(stars, ra_scale, lst, cooldown, now)

        if just_crossed:
            star = just_crossed[0][0]
            ra = float(star["ra_deg"]) * ra_scale
            rec = _star_rec(star, ra_scale, lat_deg)
            rec["lst_deg"] = lst
            dyads = star_to_dyads(rec)
            diff = just_crossed[0][1]
            _log(f"fire {star['name']} just crossed lst={lst:.4f}° ra={ra:.4f}° |LST-RA|={diff:.4f}°", verbose)

            if not quiet and sys.stdout.isatty():
                print(" " * 60, end="\r", flush=True)
            if not quiet:
                up = _upcoming(stars, ra_scale, lst, lon, now, 2)
                viz.print_transit(rec, dyads, upcoming=up)

            if stellarium_url:
                cand = _stellarium_candidates(rec, star)
                slew_to(
                    rec["ra_deg"],
                    rec["dec_deg"],
                    base_url=stellarium_url,
                    target_candidates=cand,
                    verbose=verbose,
                )

            send_dyads(port, dyads, note_duration=note_duration, verbose=verbose)
            cooldown[star["name"]] = now

        elif not quiet and sys.stdout.isatty() and (now - last_upcoming >= 30.0 or last_upcoming == 0.0):
            up = _upcoming(stars, ra_scale, lst, lon, now, 2)
            if up:
                print(viz.format_next(up[0][0], up[0][1]), end="\r", flush=True)
            last_upcoming = now

        time.sleep(POLL_INTERVAL_S)
