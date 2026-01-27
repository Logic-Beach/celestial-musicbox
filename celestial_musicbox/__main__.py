import argparse
import os
from pathlib import Path

from .midi_sender import get_mido_output_names
from .transit_scheduler import run_scheduler


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    ap = argparse.ArgumentParser(
        description="Celestial Musicbox: meridian transits from a star catalog → MIDI.",
    )
    ap.add_argument(
        "--catalog",
        type=Path,
        default=root / "data" / "star_catalog.json",
        help="Path to star_catalog.json (default: data/star_catalog.json).",
    )
    ap.add_argument("--lon", type=float, default=None, help="Observer longitude (degrees, East positive). Full precision kept (e.g. -86.808250).")
    ap.add_argument("--lat", type=float, default=None, help="Observer latitude (degrees). Full precision kept.")
    ap.add_argument(
        "--supplement",
        type=Path,
        default=root / "data" / "star_supplement.json",
        help="Path to star_supplement.json (mass, spectral, distance).",
    )
    ap.add_argument(
        "--midi-port",
        type=str,
        default=None,
        help="MIDI output port name (substring match). If not set, use MIDI_PORT env or first available.",
    )
    ap.add_argument(
        "--list-ports",
        action="store_true",
        help="List MIDI output ports and exit.",
    )
    ap.add_argument(
        "--quiet",
        action="store_true",
        help="No terminal visualizer (transit blocks, countdown, MIDI).",
    )
    ap.add_argument(
        "--stellarium",
        action="store_true",
        help="Slew Stellarium to each transiting star (Remote Control on localhost:8090).",
    )
    ap.add_argument(
        "--stellarium-url",
        type=str,
        default=None,
        help="Stellarium Remote Control base URL (e.g. http://localhost:8090). Implied by --stellarium.",
    )
    ap.add_argument(
        "--note-duration",
        type=float,
        default=None,
        help="Seconds the chord (both notes) rings per star (default 0.6).",
    )
    ap.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Debug logging to stderr (catalog, heap, wait, fire, Stellarium, MIDI).",
    )
    args = ap.parse_args()

    if args.list_ports:
        for n in get_mido_output_names():
            print(n)
        return

    if args.lon is None or args.lat is None:
        ap.error("--lon and --lat are required.")
    if not args.catalog.is_file():
        ap.error(f"Catalog not found: {args.catalog}. Run: python scripts/build_star_catalog.py --lat <latitude>")

    midi_port = args.midi_port or os.environ.get("MIDI_PORT")
    stellarium_url = args.stellarium_url or ("http://localhost:8090" if args.stellarium else None)
    run_scheduler(
        catalog_path=args.catalog,
        supplement_path=args.supplement,
        lon_deg=args.lon,
        lat_deg=args.lat,
        midi_port_name=midi_port,
        quiet=args.quiet,
        stellarium_url=stellarium_url,
        note_duration=args.note_duration,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
