"""
Meridian transit trigger: Fire when LST = RA (star crosses meridian).

A star crosses the meridian when Local Sidereal Time equals the star's Right Ascension.
Poll every POLL_INTERVAL and check which stars have LST ≈ RA (within tolerance).
Cooldown each star for ~1 sidereal day after firing.
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path
from typing import Optional

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from . import visualizer as viz
from .midi_sender import open_mido_output, send_dyads
from .stellarium_client import slew_to, _get_status as get_status, get_object_info, get_objects_batch_info
from .star_to_midi import star_to_dyads

SIDEREAL_DAY_S = 86164.0905
POLL_INTERVAL_S = 0.5
# Trigger tolerance: account for query delay to Stellarium
# We trigger slightly BEFORE meridian to compensate for ~1s processing/query delay
# Range: -0.005° (1.2s before) to +0.015° (3.6s after)
# This centers the actual azimuth query at ~180°
TRANSIT_TOLERANCE_MIN_DEG = -0.005  # Allow triggering slightly before meridian
TRANSIT_TOLERANCE_MAX_DEG = 0.015   # Up to this far past meridian
COOLDOWN_FRAC = 0.98  # don't re-fire same star until 98% of sidereal day has passed


def _unix() -> float:
    return time.time()


def _time(t: float) -> Time:
    return Time(t, format="unix")


def _norm_360(x: float) -> float:
    """Normalize angle to [0, 360)."""
    v = float(x) % 360.0
    return 0.0 if v >= 360.0 else v


def _lst_deg(t: Time, lon_deg: float) -> float:
    """Calculate Local Sidereal Time in degrees."""
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


def _j2000_to_apparent_ra_fast(ra_j2000_deg: float, dec_j2000_deg: float, time: Time) -> float:
    """Convert J2000 RA/Dec to apparent RA for given time (accounts for precession/nutation).
    Returns apparent RA in degrees.
    
    Uses a simple linear approximation for precession:
    - Average precession rate is ~50 arcsec/year = 0.0139°/year
    - For RA, the rate depends on declination
    """
    from astropy.coordinates import CIRS
    # Quick calculation using Astropy's built-in precession
    # This is faster than full transformation but still accurate
    coord_j2000 = SkyCoord(ra=ra_j2000_deg * u.deg, dec=dec_j2000_deg * u.deg, frame="icrs", obstime=Time("J2000"))
    coord_apparent = coord_j2000.transform_to(CIRS(obstime=time))
    return float(coord_apparent.ra.deg)


def _angular_separation(angle1: float, angle2: float) -> float:
    """Calculate the minimum angular separation between two angles (0-360°).
    Returns value in [0, 180]."""
    diff = abs(angle1 - angle2)
    if diff > 180.0:
        diff = 360.0 - diff
    return diff


def _get_location_from_stellarium(stellarium_url: str, verbose: bool) -> Optional[tuple[float, float]]:
    """Get (lon_deg, lat_deg) from Stellarium. Returns None if not available."""
    st = get_status(stellarium_url)
    if not st:
        return None
    loc = st.get("location")
    if not isinstance(loc, dict):
        return None
    lon = loc.get("longitude")
    lat = loc.get("latitude")
    if lon is None or lat is None:
        return None
    try:
        lon_deg = float(lon)
        lat_deg = float(lat)
        if verbose:
            print(f"[cmb] Using Stellarium location: lon={lon_deg:.8f}° lat={lat_deg:.8f}°", 
                  file=sys.stderr, flush=True)
        return (lon_deg, lat_deg)
    except (TypeError, ValueError):
        return None


def _get_lst_and_time(
    stellarium_url: Optional[str],
    lon_deg: float,
    verbose: bool,
) -> tuple[float, float, Time]:
    """Return (lst_deg, lon_deg, time). 
    When stellarium_url set, use Stellarium's time and location for perfect sync.
    Returns the Time object for coordinate conversions.
    """
    if stellarium_url:
        st = get_status(stellarium_url)
        if st:
            # Try to get location from Stellarium
            loc = st.get("location")
            if isinstance(loc, dict):
                st_lon = loc.get("longitude")
                if st_lon is not None:
                    try:
                        lon_deg = float(st_lon)
                    except (TypeError, ValueError):
                        pass
            
            # Get time from Stellarium
            tinfo = st.get("time")
            if isinstance(tinfo, dict):
                jd = tinfo.get("jday")
                if jd is not None:
                    try:
                        t = Time(float(jd), format="jd")
                        lst = _lst_deg(t, lon_deg)
                        return (lst, lon_deg, t)
                    except (TypeError, ValueError):
                        pass
            if verbose:
                print("[cmb] Stellarium status missing jday; using --lon and system time", 
                      file=sys.stderr, flush=True)
    t = _time(_unix())
    lst = _lst_deg(t, lon_deg)
    return (lst, lon_deg, t)


def _altitude_at_transit(dec_deg: float, lat_deg: float) -> float:
    """Calculate altitude when star transits the meridian."""
    x = math.radians(dec_deg - lat_deg)
    return math.degrees(math.asin(max(-1.0, min(1.0, math.cos(x)))))


def _star_rec(star: dict, ra_scale: float, lat_deg: float) -> dict:
    """Build star record with computed fields."""
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
    """Generate candidate names for Stellarium search."""
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


def _format_upcoming_stars(
    stars: list[dict],
    ra_scale: float,
    lst_deg: float,
    stellarium_url: Optional[str],
    n: int = 5
) -> str:
    """Format a list of upcoming stars with their distance to meridian from Stellarium."""
    # Find stars approaching meridian
    candidates: list[tuple[str, float, float]] = []
    
    for s in stars:
        name = s.get("name") or ""
        if not name:
            continue
        ra_j2000_deg = float(s["ra_deg"]) * ra_scale
        
        # Quick check: is star approaching meridian? (within 5° ahead)
        diff_j2000 = ra_j2000_deg - lst_deg
        if diff_j2000 < 0:
            diff_j2000 += 360.0
        if diff_j2000 > 180.0:
            diff_j2000 -= 360.0
        
        if 0 < diff_j2000 <= 5.0:
            candidates.append((name, ra_j2000_deg, diff_j2000))
    
    if not candidates:
        return ""
    
    # Sort by distance (closest first)
    candidates.sort(key=lambda x: x[2])
    
    # Query Stellarium for accurate apparent RA of top candidates
    lines = []
    colors = ["\033[92m", "\033[93m", "\033[94m", "\033[95m", "\033[96m"]  # Green, yellow, blue, magenta, cyan
    reset = "\033[0m"
    bold = "\033[1m"
    dim = "\033[2m"
    
    for i, (name, ra_j2000, _) in enumerate(candidates[:n]):
        if stellarium_url:
            # Find this star in the catalog to get identifiers
            star = next((s for s in stars if s.get("name") == name), None)
            if star:
                cand_names = _stellarium_candidates({"name": name}, star)
                if cand_names:
                    obj_info = get_object_info(stellarium_url, cand_names[0])
                    if obj_info:
                        ra_apparent = obj_info.get("ra")
                        if ra_apparent is not None:
                            dist_to_meridian = float(ra_apparent) - lst_deg
                            if dist_to_meridian < 0:
                                dist_to_meridian += 360.0
                            if dist_to_meridian > 180.0:
                                dist_to_meridian -= 360.0
                            
                            color = colors[i % len(colors)]
                            # Calculate time to meridian (degrees / rate)
                            # LST advances at ~0.004178° per second
                            lst_rate = 360.0 / 86164.0905  # degrees per second
                            seconds_to_meridian = dist_to_meridian / lst_rate
                            
                            # Format time as MM:SS
                            mins = int(seconds_to_meridian // 60)
                            secs = int(seconds_to_meridian % 60)
                            time_str = f"{mins:2d}m {secs:02d}s"
                            
                            lines.append(f"{color}{bold}{i+1}.{reset} {color}{name:20s}{reset} {dim}│{reset} {color}{dist_to_meridian:5.2f}°{reset} {dim}│{reset} {color}{time_str}{reset}")
    
    if lines:
        header = f"\n{bold}\033[96m╔═══════════════════════════════════════════════════════════════╗{reset}"
        title = f"{bold}\033[96m║{reset}  {bold}⏳ APPROACHING MERIDIAN{reset}                                  {bold}\033[96m║{reset}"
        divider = f"{bold}\033[96m╚═══════════════════════════════════════════════════════════════╝{reset}"
        return header + "\n" + title + "\n" + divider + "\n" + "\n".join(lines)
    return ""


def _stars_transiting_now(
    stars: list[dict],
    ra_scale: float,
    lst_deg: float,
    current_time: Time,
    cooldown: dict[str, float],
    now: float,
    stellarium_url: Optional[str],
    verbose: bool,
) -> list[tuple[dict, float, float, float]]:
    """Find stars transiting NOW (LST just passed apparent RA, within tolerance).
    A star transits when LST catches up to and equals its APPARENT RA (not J2000).
    
    If stellarium_url is provided, we query Stellarium for the apparent RA (slow but accurate).
    Otherwise we calculate it ourselves (fast but may have small discrepancies).
    
    Returns list of (star, ra_j2000_deg, ra_apparent_deg, how_far_past) sorted by how_far_past."""
    cooldown_s = COOLDOWN_FRAC * SIDEREAL_DAY_S
    
    # First pass: find candidates based on J2000 RA (fast pre-filter)
    # We use a wider tolerance for this pass
    pre_candidates: list[tuple[dict, float]] = []
    for s in stars:
        name = s.get("name") or ""
        if name and (now - cooldown.get(name, -1e9)) < cooldown_s:
            continue
        ra_j2000_deg = float(s["ra_deg"]) * ra_scale
        # Quick check: is star roughly near meridian? (±2° tolerance for pre-filter)
        diff_j2000 = lst_deg - ra_j2000_deg
        if diff_j2000 < -180.0:
            diff_j2000 += 360.0
        elif diff_j2000 > 180.0:
            diff_j2000 -= 360.0
        if -2.0 <= diff_j2000 <= 2.0:
            pre_candidates.append((s, ra_j2000_deg))
    
    if not pre_candidates:
        return []
    
    # Second pass: get accurate apparent RA from Stellarium
    candidates: list[tuple[dict, float, float, float]] = []
    
    if stellarium_url:
        # Query Stellarium for apparent RA of each candidate
        for s, ra_j2000_deg in pre_candidates:
            name = s.get("name") or ""
            # Build candidate names for Stellarium
            cand_names = _stellarium_candidates({"name": name}, s)
            if not cand_names:
                continue
            
            # Query Stellarium for this star
            obj_info = get_object_info(stellarium_url, cand_names[0])
            if not obj_info:
                continue
            
            # Get apparent RA from Stellarium (in degrees)
            ra_apparent = obj_info.get("ra")
            if ra_apparent is None:
                continue
            
            ra_apparent_deg = float(ra_apparent)
            
            # Calculate difference from LST
            diff = lst_deg - ra_apparent_deg
            if diff < -180.0:
                diff += 360.0
            elif diff > 180.0:
                diff -= 360.0
            
            # Check if within tolerance
            if TRANSIT_TOLERANCE_MIN_DEG <= diff <= TRANSIT_TOLERANCE_MAX_DEG:
                candidates.append((s, ra_j2000_deg, ra_apparent_deg, diff))
                if verbose:
                    status = "approaching" if diff < 0 else "past"
                    _log(f"  candidate: {name} J2000_RA={ra_j2000_deg:.3f}° Stellarium_apparent_RA={ra_apparent_deg:.3f}° LST={lst_deg:.3f}° {status} by {abs(diff):.4f}°", verbose)
    else:
        # Fallback: calculate apparent RA ourselves
        for s, ra_j2000_deg in pre_candidates:
            dec_j2000_deg = float(s["dec_deg"])
            ra_apparent_deg = _j2000_to_apparent_ra_fast(ra_j2000_deg, dec_j2000_deg, current_time)
            
            diff = lst_deg - ra_apparent_deg
            if diff < -180.0:
                diff += 360.0
            elif diff > 180.0:
                diff -= 360.0
            
            if TRANSIT_TOLERANCE_MIN_DEG <= diff <= TRANSIT_TOLERANCE_MAX_DEG:
                candidates.append((s, ra_j2000_deg, ra_apparent_deg, diff))
    
    # Sort by diff (smallest = most recent transit)
    candidates.sort(key=lambda x: x[3])
    return candidates


def _wait_deg(lst: float, ra: float) -> float:
    """Calculate degrees of LST rotation until star transits (RA - LST, normalized)."""
    return _norm_360(ra - lst)


def _upcoming(
    stars: list[dict], 
    ra_scale: float, 
    lst: float, 
    lon_deg: float, 
    now: float, 
    n: int = 2
) -> list[tuple[str, float]]:
    """Next n stars to transit. Return [(name, seconds)]."""
    rate = 360.0 / SIDEREAL_DAY_S  # degrees per second
    cand: list[tuple[dict, float, float]] = []
    
    for s in stars:
        ra = float(s["ra_deg"]) * ra_scale
        wait_d = _wait_deg(lst, ra)
        # Only consider stars approaching transit (not already past)
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
    """Run the meridian transit scheduler."""
    _log(f"Loading catalog {catalog_path}", verbose)
    
    # If Stellarium URL provided, get location from Stellarium for perfect sync
    if stellarium_url:
        st = get_status(stellarium_url)
        if st:
            # Display Stellarium configuration
            print("=" * 60, file=sys.stderr)
            print("STELLARIUM CONFIGURATION:", file=sys.stderr)
            print("=" * 60, file=sys.stderr)
            
            # Location
            loc = st.get("location", {})
            if isinstance(loc, dict):
                st_lon = loc.get("longitude")
                st_lat = loc.get("latitude")
                st_alt = loc.get("altitude", 0)
                st_name = loc.get("name", "Unknown")
                st_planet = loc.get("planet", "Earth")
                print(f"Location: {st_name} ({st_planet})", file=sys.stderr)
                print(f"  Lon: {st_lon:.8f}°  Lat: {st_lat:.8f}°  Alt: {st_alt}m", file=sys.stderr)
            
            # Time settings
            time_info = st.get("time", {})
            if isinstance(time_info, dict):
                jd = time_info.get("jday")
                time_rate = time_info.get("timerate", 1.0)
                is_paused = time_info.get("timepaused", False)
                print(f"Time: JD={jd}", file=sys.stderr)
                print(f"  Rate: {time_rate}x  Paused: {is_paused}", file=sys.stderr)
                if abs(time_rate - 1.0) > 0.01:
                    print(f"  WARNING: Time rate is not 1.0! This will cause issues.", file=sys.stderr)
            
            # View settings
            view = st.get("view", {})
            if isinstance(view, dict):
                fov = view.get("fov")
                print(f"View FOV: {fov}°", file=sys.stderr)
            
            # Check for settings that might affect coordinates
            # Note: Stellarium API doesn't expose all settings, but we can check what we can
            print("\nNOTE: Check these Stellarium settings manually:", file=sys.stderr)
            print("  - F4 > Sky tab > Atmosphere: Refraction (should match real observations)", file=sys.stderr)
            print("  - F2 > Tools tab > Coordinates: Should be using J2000 epoch", file=sys.stderr)
            print("=" * 60, file=sys.stderr)
            print(file=sys.stderr)
        
        stellarium_loc = _get_location_from_stellarium(stellarium_url, verbose)
        if stellarium_loc:
            lon_deg, lat_deg = stellarium_loc
            print(f"Using location from Stellarium: lon={lon_deg:.6f}° lat={lat_deg:.6f}°\n", 
                  file=sys.stderr, flush=True)
        else:
            print(f"Could not get location from Stellarium, using --lon/--lat: lon={lon_deg:.6f}° lat={lat_deg:.6f}°\n", 
                  file=sys.stderr, flush=True)
    else:
        _log(f"Location: lon={lon_deg:.8f}° lat={lat_deg:.8f}°", verbose)
    
    raw = json.loads(catalog_path.read_text(encoding="utf-8"))
    if not isinstance(raw, list):
        raise SystemExit("star_catalog.json must be a JSON array.")

    # Filter stars that cross at azimuth 180° (south meridian)
    # Stars cross at az=180° when: -(90-lat) < dec < lat
    # For lat=36°N: -54° < dec < 36°
    # This excludes:
    #   - Circumpolar stars (dec > 36°) that cross at az=0° (north)
    #   - Stars below horizon (dec < -54°)
    min_dec = -(90.0 - abs(lat_deg))
    max_dec = lat_deg
    stars = [
        s
        for s in raw
        if isinstance(s, dict)
        and s.get("name")
        and s.get("ra_deg") is not None
        and s.get("dec_deg") is not None
        and min_dec < float(s["dec_deg"]) < max_dec
    ]
    _log(f"Catalog: {len(raw)} total, {len(stars)} cross at Az=180° (dec >{min_dec:.1f}° and <{max_dec:.1f}°)", verbose)

    # Detect if RA is in hours (<=24) or degrees (<=360)
    ras = [float(s["ra_deg"]) for s in stars]
    ra_scale = 15.0 if (ras and max(ras) <= 24.0) else 1.0
    if ra_scale != 1.0:
        print("Note: catalog ra_deg in hours; converting to degrees.\n", file=sys.stderr)

    # Setup
    port = open_mido_output(midi_port_name, verbose=verbose)
    cooldown: dict[str, float] = {}

    # Initialize cooldown: mark stars currently transiting as already fired
    now = _unix()
    lst_init, _, time_init = _get_lst_and_time(stellarium_url, lon_deg, verbose)
    print("Checking for stars currently at meridian (this may take a moment)...", file=sys.stderr, flush=True)
    transiting = _stars_transiting_now(stars, ra_scale, lst_init, time_init, {}, now, stellarium_url, verbose=False)
    for star, _, _, _ in transiting:
        cooldown[star.get("name", "")] = now
    if transiting:
        _log(f"Startup: marked {len(transiting)} stars currently transiting as already fired", verbose)

    if not quiet:
        lst, lon, current_time = _get_lst_and_time(stellarium_url, lon_deg, verbose)
        up = _upcoming(stars, ra_scale, lst, lon, _unix(), 2)
        
        # L33t styled startup banner
        cyan = "\033[96m"
        green = "\033[92m"
        yellow = "\033[93m"
        bold = "\033[1m"
        reset = "\033[0m"
        
        print(f"\n{bold}{cyan}╔═══════════════════════════════════════════════════════╗{reset}", flush=True)
        print(f"{bold}{cyan}║{reset}  {bold}🎵 CELESTIAL MUSICBOX - MERIDIAN TRANSIT TRIGGER{reset}  {bold}{cyan}║{reset}", flush=True)
        print(f"{bold}{cyan}╚═══════════════════════════════════════════════════════╝{reset}", flush=True)
        print(f"{green}★{reset} Tracking {yellow}{len(stars)}{reset} stars crossing at Az=180°", flush=True)
        print(f"{green}★{reset} LST: {yellow}{lst/15:.2f}h{reset} ({yellow}{lst:.2f}°{reset})", flush=True)
        print(f"{green}★{reset} Location: {yellow}{lat_deg:.4f}°N{reset} {yellow}{lon_deg:.4f}°W{reset}\n", flush=True)

    # Main loop
    while True:
        now = _unix()
        lst, lon, current_time = _get_lst_and_time(stellarium_url, lon_deg, verbose)
        
        # Find stars transiting now (LST ≈ apparent RA from Stellarium)
        transiting = _stars_transiting_now(stars, ra_scale, lst, current_time, cooldown, now, stellarium_url, verbose)

        if transiting:
            # Trigger the star closest to exact transit (smallest |LST - apparent RA|)
            star, ra_j2000_deg, ra_apparent_deg, how_far_past = transiting[0]
            name = star.get("name", "")
            
            # L33t trigger notification with colors
            if not quiet:
                status_text = f"{'past' if how_far_past >= 0 else 'before'}"
                # Cyan for star name, magenta for numbers, bold for TRANSIT
                cyan = "\033[96m"
                magenta = "\033[95m"
                bold = "\033[1m"
                reset = "\033[0m"
                print(f"\n🎵 {bold}[TRANSIT]{reset} {cyan}{name}{reset} │ {magenta}{abs(how_far_past):.3f}°{reset} {status_text} meridian", 
                      file=sys.stderr, flush=True)
            
            # Build star record (use J2000 for catalog compatibility)
            rec = _star_rec(star, ra_scale, lat_deg)
            rec["lst_deg"] = lst
            dyads = star_to_dyads(rec)
            
            # Display transit info box
            if not quiet:
                if sys.stdout.isatty():
                    print(" " * 60, end="\r", flush=True)
                up = _upcoming(stars, ra_scale, lst, lon, now, 2)
                # Add apparent RA to the record for display
                rec["ra_apparent_deg"] = ra_apparent_deg
                viz.print_transit(rec, dyads, upcoming=up)
            
            # Slew Stellarium and verify azimuth
            if stellarium_url:
                cand = _stellarium_candidates(rec, star)
                
                # Check azimuth BEFORE slewing
                if cand and not quiet:
                    obj_info = get_object_info(stellarium_url, cand[0])
                    if obj_info:
                        azimuth = obj_info.get("azimuth")
                        altitude = obj_info.get("altitude")
                        if azimuth is not None:
                            az_error = azimuth - 180.0
                            if az_error > 180:
                                az_error -= 360
                            elif az_error < -180:
                                az_error += 360
                            # Color code based on accuracy
                            if abs(az_error) < 0.1:
                                color = "\033[92m"  # Green - excellent
                            elif abs(az_error) < 0.3:
                                color = "\033[93m"  # Yellow - good
                            else:
                                color = "\033[91m"  # Red - needs adjustment
                            reset = "\033[0m"
                            print(f"  {color}Stellarium: Az={azimuth:.2f}° (Δ={az_error:+.2f}°) Alt={altitude:.1f}°{reset}", 
                                  file=sys.stderr, flush=True)
                
                slew_to(
                    rec["ra_deg"],
                    rec["dec_deg"],
                    base_url=stellarium_url,
                    target_candidates=cand,
                    verbose=verbose,
                )
            
            # Send MIDI
            send_dyads(port, dyads, note_duration=note_duration, verbose=verbose)
            
            # Add to cooldown
            cooldown[name] = now

        elif not quiet and sys.stdout.isatty():
            # Show upcoming stars with their distance to meridian
            # Update every 10 seconds to avoid too many Stellarium queries
            import time
            current_time_s = time.time()
            if not hasattr(_format_upcoming_stars, '_last_update'):
                _format_upcoming_stars._last_update = 0
            
            if current_time_s - _format_upcoming_stars._last_update >= 10:
                upcoming_display = _format_upcoming_stars(stars, ra_scale, lst, stellarium_url, n=5)
                if upcoming_display:
                    # Clear previous output and show new list
                    print("\033[2J\033[H", end="")  # Clear screen and move to top
                    print(upcoming_display, flush=True)
                _format_upcoming_stars._last_update = current_time_s

        time.sleep(POLL_INTERVAL_S)
