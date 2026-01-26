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


def slew_to(
    ra_deg: float,
    dec_deg: float,
    base_url: str = "http://localhost:8090",
    mode: str = "mark",
    target: Optional[str] = None,
    target_candidates: Optional[Iterable[str]] = None,
) -> None:
    """Set Stellarium view to J2000 (ra_deg, dec_deg) and optionally select object.

    1. POST /api/main/view with j2000 to set view (no animation).
    2. If target or target_candidates: use /api/objects/find then /api/main/focus to select.
       Tries each candidate: find(c), then focus(first match). Uses target if candidates omitted.
    """
    global _STELLARIUM_WARNED

    ra_rad = math.radians(ra_deg)
    dec_rad = math.radians(dec_deg)
    x = math.cos(dec_rad) * math.cos(ra_rad)
    y = math.cos(dec_rad) * math.sin(ra_rad)
    z = math.sin(dec_rad)
    j2000 = [x, y, z]

    url_view = f"{base_url.rstrip('/')}/api/main/view"
    body_view = urlencode({"j2000": json.dumps(j2000)}).encode("utf-8")

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

    candidates: list[str] = []
    if target_candidates is not None:
        candidates = list(target_candidates)
    elif target:
        candidates = [target]

    for q in candidates:
        if not q or not str(q).strip():
            continue
        matches = _find_objects(base_url, str(q).strip())
        for m in matches:
            if m and _focus(base_url, m, mode):
                return
        if not matches and _focus(base_url, str(q).strip(), mode):
            return
