#!/usr/bin/env python3
"""
Build data/star_catalog.json from the HYG catalog.

  python scripts/build_star_catalog.py --lat 36
  python scripts/build_star_catalog.py --hyg path/to/hygdata_v42.csv --max-mag 7 --lat 36

Output: data/star_catalog.json
"""

from __future__ import annotations

import csv
import gzip
import json
import sys
from pathlib import Path
from urllib.request import urlretrieve

HYG_URL = "https://astronexus.com/downloads/catalogs/hygdata_v42.csv.gz"
# Parsecs to light-years
PC_TO_LY = 1.0 / 0.306601

# HYG columns we use (v3/v4): id,hip,hd,hr,gl,bf,proper,ra,dec,dist,pmra,pmrv,mag,absmag,spect,...
# proper=name; bf=Bayer/Flamsteed; hr=HR number. dist in parsecs, mag=Vmag, spect=spectral.
# HYG "ra" is in HOURS (0-24). We convert to degrees (ra*15) for ra_deg.


def _project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def _download_hyg(dest: Path) -> Path:
    print("Downloading HYG catalog (may take a moment)...")
    urlretrieve(HYG_URL, dest)
    return dest


def _parse_float(v: str | None, default: float | None = None) -> float | None:
    if v is None or v == "" or (isinstance(v, str) and v.strip() == ""):
        return default
    try:
        return float(v)
    except (ValueError, TypeError):
        return default


def _name_from_hyg(row: dict) -> str | None:
    # Prefer proper name, then Bayer/Flamsteed, then "HR N", else "HIP N" or "HD N".
    # Many HYG stars have only hip/hd; without these we'd only get ~9k (BSC + named).
    n = (row.get("proper") or "").strip()
    if n:
        return n
    n = (row.get("bf") or "").strip()
    if n:
        return n
    hr = row.get("hr")
    if hr is not None:
        try:
            h = int(float(hr))
            if h > 0:
                return f"HR {h}"
        except (ValueError, TypeError):
            pass
    hip = row.get("hip")
    if hip is not None:
        try:
            h = int(float(hip))
            if h > 0:
                return f"HIP {h}"
        except (ValueError, TypeError):
            pass
    hd = row.get("hd")
    if hd is not None:
        try:
            h = int(float(hd))
            if h > 0:
                return f"HD {h}"
        except (ValueError, TypeError):
            pass
    return None


def build_from_hyg(
    hyg_path: Path | None,
    out_path: Path,
    supplement: dict,
    max_mag: float = 9.0,
    lat_deg: float | None = None,
) -> int:
    if hyg_path is None or not hyg_path.is_file():
        # Download to data/hygdata_v42.csv.gz
        data_dir = out_path.parent
        data_dir.mkdir(parents=True, exist_ok=True)
        gz = data_dir / "hygdata_v42.csv.gz"
        if not gz.is_file():
            _download_hyg(gz)
        hyg_path = gz

    # Support .gz
    if str(hyg_path).endswith(".gz"):
        f = gzip.open(hyg_path, "rt", encoding="utf-8", errors="replace")
    else:
        f = hyg_path.open("r", encoding="utf-8", errors="replace")

    rows: list[dict] = []
    with f:
        reader = csv.DictReader(f)
        for row in reader:
            mag = _parse_float(row.get("mag"))
            if mag is not None and mag > max_mag:
                continue
            name = _name_from_hyg(row)
            if not name or name == "Sol":
                continue
            ra_h = _parse_float(row.get("ra"))
            dec = _parse_float(row.get("dec"))
            if ra_h is None or dec is None:
                continue
            # HYG ra is in hours (0-24). Convert to degrees for ra_deg.
            ra = ra_h * 15.0
            # Only stars that can ever be above the horizon at this latitude.
            # dec in [lat-90, lat+90]; at 36°N that's dec >= -54°.
            if lat_deg is not None:
                if dec < lat_deg - 90.0 or dec > lat_deg + 90.0:
                    continue
            rec: dict = {
                "name": name,
                "ra_deg": ra,
                "dec_deg": dec,
                "vmag": mag if mag is not None else 5.0,
            }
            dist_pc = _parse_float(row.get("dist"))
            if dist_pc is not None and dist_pc > 0:
                rec["distance_ly"] = round(dist_pc * PC_TO_LY, 2)
            sp = (row.get("spect") or "").strip()
            if sp:
                rec["spectral"] = sp
            # Merge supplement for mass (and overwrite spectral/distance if we want; supplement usually better for named)
            sup = supplement.get(name, {})
            if isinstance(sup, dict):
                if sup.get("mass") is not None:
                    rec["mass"] = sup["mass"]
                if sup.get("spectral") is not None and not rec.get("spectral"):
                    rec["spectral"] = sup["spectral"]
                if sup.get("distance_ly") is not None:
                    rec["distance_ly"] = sup["distance_ly"]
            rows.append(rec)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out:
        json.dump(rows, out, indent=2)
    msg = f"Wrote {len(rows)} stars to {out_path}"
    if lat_deg is not None:
        dec_lo = max(-90, lat_deg - 90)
        dec_hi = min(90, lat_deg + 90)
        msg += f" (dec [{dec_lo:.0f}°, {dec_hi:.0f}°] for lat {lat_deg}°)"
    print(msg)
    return len(rows)


def main() -> int:
    import argparse

    ap = argparse.ArgumentParser(description="Build star_catalog.json from HYG.")
    ap.add_argument("--hyg", type=Path, default=None, help="Path to HYG CSV/CSV.GZ. If not set, download.")
    ap.add_argument("--max-mag", type=float, default=8.0, help="Max magnitude for HYG (default 8).")
    ap.add_argument("--lat", type=float, default=None, help="Observer latitude (deg). Keep only stars that can rise at this latitude (dec in [lat-90, lat+90]).")
    ap.add_argument("--supplement", type=Path, default=None, help="Supplement JSON (default: data/star_supplement.json).")
    ap.add_argument("-o", "--output", type=Path, default=None, help="Output path (default: data/star_catalog.json).")
    args = ap.parse_args()

    root = _project_root()
    out = args.output or root / "data" / "star_catalog.json"
    sup_path = args.supplement or root / "data" / "star_supplement.json"
    supplement = {}
    if sup_path.is_file():
        try:
            supplement = json.loads(sup_path.read_text(encoding="utf-8"))
        except Exception:
            pass

    n = build_from_hyg(args.hyg, out, supplement, max_mag=args.max_mag, lat_deg=args.lat)
    return 0 if n else 1


if __name__ == "__main__":
    sys.exit(main())
