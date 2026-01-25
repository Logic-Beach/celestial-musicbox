"""
Stellarium Remote Control client: slew view to J2000 RA/Dec and optionally select by name.

Requires Stellarium running with Remote Control plugin on port 8090 (or --stellarium-url).
"""

import json
import math
import sys
from typing import Optional
from urllib import request
from urllib.error import URLError
from urllib.parse import urlencode

_STELLARIUM_WARNED = False


def slew_to(
    ra_deg: float,
    dec_deg: float,
    base_url: str = "http://localhost:8090",
    mode: str = "center",
    target: Optional[str] = None,
) -> None:
    """Slew Stellarium to J2000 (ra_deg, dec_deg) and optionally select the star by name.

    1. POST /api/main/focus with position (J2000 unit vector) to slew the view.
    2. If target is given, POST /api/main/focus with target=name to select that object
       (so it is marked, info can be shown, etc.). Stellarium looks up by localized then
       english name; catalog names like HIP 12345, HR 123, or Vega often match.

    On connection/HTTP errors for the slew, prints a one-time warning and returns. Failure
    of the select-by-name call is ignored (view is already at the right position).
    """
    global _STELLARIUM_WARNED

    ra_rad = math.radians(ra_deg)
    dec_rad = math.radians(dec_deg)
    x = math.cos(dec_rad) * math.cos(ra_rad)
    y = math.cos(dec_rad) * math.sin(ra_rad)
    z = math.sin(dec_rad)
    position = [x, y, z]

    url = f"{base_url.rstrip('/')}/api/main/focus"
    body = urlencode({"position": json.dumps(position), "mode": mode}).encode("utf-8")

    try:
        req = request.Request(
            url,
            data=body,
            method="POST",
            headers={"Content-Type": "application/x-www-form-urlencoded"},
        )
        with request.urlopen(req, timeout=5) as resp:
            if resp.status >= 400 and not _STELLARIUM_WARNED:
                _STELLARIUM_WARNED = True
                print(
                    f"Stellarium focus: HTTP {resp.status} from {url}",
                    file=sys.stderr,
                )
    except (URLError, OSError) as e:
        if not _STELLARIUM_WARNED:
            _STELLARIUM_WARNED = True
            print(
                f"Stellarium focus: {e}. Is Stellarium running with Remote Control on {base_url}?",
                file=sys.stderr,
            )
        return

    # Select the star by name so it is marked and info can be shown
    if target:
        body2 = urlencode({"target": target, "mode": mode}).encode("utf-8")
        try:
            req2 = request.Request(
                url,
                data=body2,
                method="POST",
                headers={"Content-Type": "application/x-www-form-urlencoded"},
            )
            request.urlopen(req2, timeout=5)
        except (URLError, OSError):
            pass  # Slew already done; skip select if name not in Stellarium's catalog
