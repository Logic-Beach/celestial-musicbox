"""
Stellarium Remote Control client: slew view to J2000 RA/Dec and optionally select by name.

Requires Stellarium running with Remote Control plugin on port 8090 (or --stellarium-url).
"""

import json
import math
import sys
from typing import Iterable, Optional
from urllib import request
from urllib.error import URLError
from urllib.parse import urlencode

_STELLARIUM_WARNED = False


def _get_status(base_url: str) -> dict | None:
    """GET /api/main/status; returns parsed JSON or None. Includes view.fov."""
    url = f"{base_url.rstrip('/')}/api/main/status"
    try:
        req = request.Request(url, method="GET")
        with request.urlopen(req, timeout=5) as r:
            return json.loads(r.read().decode("utf-8"))
    except (URLError, OSError, json.JSONDecodeError):
        return None


def get_object_info(base_url: str, object_name: str) -> dict | None:
    """GET /api/objects/info?name=... returns object info including azimuth/altitude."""
    url = f"{base_url.rstrip('/')}/api/objects/info"
    try:
        from urllib.parse import urlencode
        req = request.Request(f"{url}?{urlencode({'name': object_name, 'format': 'json'})}", method="GET")
        with request.urlopen(req, timeout=5) as r:
            return json.loads(r.read().decode("utf-8"))
    except (URLError, OSError, json.JSONDecodeError):
        return None


def get_objects_batch_info(base_url: str, object_names: list[str]) -> dict[str, dict]:
    """Query multiple objects and return their info. Returns dict of {name: info}."""
    results = {}
    for name in object_names:
        info = get_object_info(base_url, name)
        if info:
            results[name] = info
    return results


def _set_fov(base_url: str, fov_deg: float) -> bool:
    """POST /api/main/fov with fov (degrees). Returns True if no HTTP error."""
    url = f"{base_url.rstrip('/')}/api/main/fov"
    body = urlencode({"fov": fov_deg}).encode("utf-8")
    try:
        req = request.Request(
            url,
            data=body,
            method="POST",
            headers={"Content-Type": "application/x-www-form-urlencoded"},
        )
        request.urlopen(req, timeout=5)
        return True
    except (URLError, OSError):
        return False


def _find_objects(base_url: str, query: str) -> list[str]:
    """GET /api/objects/find?str=query; returns list of matching object names."""
    url = f"{base_url.rstrip('/')}/api/objects/find"
    try:
        req = request.Request(f"{url}?{urlencode({'str': query})}", method="GET")
        with request.urlopen(req, timeout=5) as r:
            raw = r.read().decode("utf-8")
            out = json.loads(raw)
            return out if isinstance(out, list) else []
    except (URLError, OSError, json.JSONDecodeError):
        return []


def _focus(base_url: str, target: str, mode: str = "mark") -> bool:
    """POST /api/main/focus with target=name, mode. Returns True if no HTTP error."""
    url = f"{base_url.rstrip('/')}/api/main/focus"
    body = urlencode({"target": target, "mode": mode}).encode("utf-8")
    try:
        req = request.Request(
            url,
            data=body,
            method="POST",
            headers={"Content-Type": "application/x-www-form-urlencoded"},
        )
        request.urlopen(req, timeout=5)
        return True
    except (URLError, OSError):
        return False


def _variants(s: str) -> list[str]:
    """Yield query variants: original, then e.g. 'HIP 677' -> 'HIP677', 'HD 358' -> 'HD358'."""
    s = str(s).strip()
    if not s:
        return []
    out = [s]
    if s.startswith("HIP ") and len(s) > 4:
        out.append("HIP" + s[4:])
    elif s.startswith("HD ") and len(s) > 3:
        out.append("HD" + s[3:])
    elif s.startswith("HR ") and len(s) > 3:
        out.append("HR" + s[3:])
    return out


def slew_to(
    ra_deg: float,
    dec_deg: float,
    base_url: str = "http://localhost:8090",
    mode: str = "mark",
    target: Optional[str] = None,
    target_candidates: Optional[Iterable[str]] = None,
    preserve_fov: bool = True,
    verbose: bool = False,
) -> None:
    """Set Stellarium view to J2000 (ra_deg, dec_deg) and optionally select object.

    1. If preserve_fov: GET /api/main/status, remember view.fov.
    2. POST /api/main/view with j2000 to set view (no animation).
    3. If preserve_fov and we have fov: POST /api/main/fov to restore it.
    4. If target or target_candidates: try /api/objects/find + /api/main/focus for each.
       Prefer HIP/HD/HR; use variants (HIP 677 / HIP677). Log to stderr when no selection.
    """
    global _STELLARIUM_WARNED

    def _log(msg: str) -> None:
        if verbose:
            print(f"[stellarium] {msg}", file=sys.stderr, flush=True)

    fov_deg: Optional[float] = None
    if preserve_fov:
        st = _get_status(base_url)
        if st:
            v = st.get("view") if isinstance(st.get("view"), dict) else {}
            f = v.get("fov")
            if f is not None:
                try:
                    fov_deg = float(f)
                except (TypeError, ValueError):
                    pass
        if fov_deg is not None:
            _log(f"current FOV {fov_deg:.4f}° (will restore after view)")
        else:
            _log("could not read FOV; skipping restore")

    ra_rad = math.radians(ra_deg)
    dec_rad = math.radians(dec_deg)
    x = math.cos(dec_rad) * math.cos(ra_rad)
    y = math.cos(dec_rad) * math.sin(ra_rad)
    z = math.sin(dec_rad)
    j2000 = [x, y, z]

    url_view = f"{base_url.rstrip('/')}/api/main/view"
    body_view = urlencode({"j2000": json.dumps(j2000)}).encode("utf-8")

    _log(f"setting view → RA={ra_deg:.4f}° Dec={dec_deg:.4f}° ({base_url})")
    try:
        req_view = request.Request(
            url_view,
            data=body_view,
            method="POST",
            headers={"Content-Type": "application/x-www-form-urlencoded"},
        )
        request.urlopen(req_view, timeout=5)
    except (URLError, OSError) as e:
        if not _STELLARIUM_WARNED:
            _STELLARIUM_WARNED = True
            print(
                f"Stellarium view: {e}. Is Stellarium running with Remote Control on {base_url}?",
                file=sys.stderr,
            )
        return
    _log("view set")

    if fov_deg is not None and _set_fov(base_url, fov_deg):
        _log(f"FOV restored to {fov_deg:.4f}°")

    candidates: list[str] = []
    if target_candidates is not None:
        candidates = list(target_candidates)
    elif target:
        candidates = [target]

    if not candidates:
        return

    _log(f"candidates: {candidates}")
    tried: list[str] = []

    for q in candidates:
        if not q or not str(q).strip():
            continue
        for v in _variants(q):
            if v in tried:
                continue
            tried.append(v)
            matches = _find_objects(base_url, v)
            _log(f"  find({v!r}) → {len(matches)} match(es)")
            for m in matches:
                if m and _focus(base_url, m, mode):
                    _log(f"  focus({m!r}) → ok (selected)")
                    return
            if not matches and _focus(base_url, v, mode):
                _log(f"  focus({v!r}) direct → ok")
                return
            _log(f"  focus({v!r}) → no")

    display = candidates[0] if candidates else "(none)"
    print(
        f"Stellarium: could not select {display!r} (tried {len(tried)} variants). "
        "View updated; star may be missing from Stellarium catalogs.",
        file=sys.stderr,
    )
