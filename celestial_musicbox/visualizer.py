"""
Terminal visualizer: transit blocks, ASCII stars by spectral type, MIDI dyads, countdown.
"""

import sys
import time
from typing import Optional

# Spectral type → symbol (hot→cool). Unicode for prettiness; fallback for plain ASCII.
SPECTRAL_SYMBOLS = {
    "O": "\u2726",   # ✦ (hot, blue)
    "B": "\u2605",   # ★
    "A": "\u2729",   # ⁑
    "F": "\u25cf",   # ●
    "G": "\u2022",   # • (Sun-like)
    "K": "\u25e6",   # ◦
    "M": "\u25cb",   # ○ (cool, red)
    "L": "\u25d0",   # ◐
    "T": "\u25d1",   # ◑
}
DEFAULT_SYMBOL = "\u2022"  # •

NOTES = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
# One chord: [color (B-V), distance]; velocity = magnitude
CHORD_LABELS = ["color", "dist"]


def spectral_to_symbol(spectral: Optional[str]) -> str:
    if not spectral or not isinstance(spectral, str):
        return DEFAULT_SYMBOL
    letter = spectral.strip().upper()[:1] if spectral else "G"
    return SPECTRAL_SYMBOLS.get(letter, DEFAULT_SYMBOL)


def symbol_repeat(symbol: str, vmag: float) -> str:
    """More copies for brighter stars (lower vmag)."""
    if vmag < 1:
        return symbol * 3
    if vmag < 3.5:
        return symbol * 2
    return symbol


def midi_note_to_name(n: int) -> str:
    octave = n // 12 - 1
    return f"{NOTES[n % 12]}{octave}"


def _format_dyads(dyads: list[tuple[int, int]]) -> str:
    """Format chord: spec C4, mag G4; velocity (from distance) on both, e.g. 'spec C4  mag G4  @90 (dist)'."""
    if not dyads:
        return ""
    parts = []
    for i, (n, v) in enumerate(dyads):
        label = CHORD_LABELS[i] if i < len(CHORD_LABELS) else ""
        parts.append(f"{label}: {midi_note_to_name(n)}")
    vel = dyads[0][1] if dyads else 0
    return "  ".join(parts) + f"  @{vel} (vel=mag)"


def _prop_line(rec: dict) -> str:
    """Stellar data: vmag, mass/alt, spectral, distance."""
    bits = []
    bits.append(f"vmag {rec.get('vmag', '?')}")
    if rec.get("mass") is not None:
        bits.append(f"mass {rec['mass']} M\u2609")
    if rec.get("altitude") is not None:
        bits.append(f"alt {rec['altitude']:.0f}\u00b0")
    if rec.get("spectral"):
        bits.append(rec["spectral"])
    d = rec.get("distance_ly") or rec.get("distance")
    if d is not None:
        bits.append(f"dist {d} ly")
    return "  ".join(bits)


def _shortest_arc_deg(a_deg: float, b_deg: float) -> float:
    """Shortest angular distance |a - b| in [0, 180] degrees. Preserves precision."""
    a, b = float(a_deg), float(b_deg)
    d = (a - b + 180.0) % 360.0 - 180.0
    return abs(d)


def _coord_line(rec: dict) -> tuple[str, str]:
    """LST, RA, Dec for sanity checking. Returns (main_line, warning_or_empty).
    At transit, LST and apparent RA should match.
    """
    ra_j2000 = rec.get("ra_deg")
    ra_apparent = rec.get("ra_apparent_deg")
    dec = rec.get("dec_deg")
    lst = rec.get("lst_deg")
    parts = []
    if lst is not None:
        parts.append(f"LST {float(lst)/15:.2f}h")
    
    # Show apparent RA if available (this should match LST at transit)
    if ra_apparent is not None:
        parts.append(f"RA(now) {float(ra_apparent):.2f}\u00b0 ({float(ra_apparent)/15:.2f}h)")
    elif ra_j2000 is not None:
        parts.append(f"RA {float(ra_j2000):.2f}\u00b0 ({float(ra_j2000)/15:.2f}h)")
    
    if dec is not None:
        parts.append(f"Dec {float(dec):+.2f}\u00b0")
    main = "  ".join(parts) if parts else ""

    # Check if LST matches apparent RA
    warning = ""
    ra_to_check = ra_apparent if ra_apparent is not None else ra_j2000
    if ra_to_check is not None and lst is not None:
        short_deg = _shortest_arc_deg(ra_to_check, lst)
        if short_deg > 1.0:  # More than 1° off
            diff_deg = short_deg
            warning = f"  \u2502  \u26a0 LST and RA differ by {diff_deg:.2f}\u00b0 \u2014 star not at meridian!\n"
    return main, warning


def format_transit_block(
    rec: dict,
    dyads: list[tuple[int, int]],
    upcoming: list[tuple[str, float]] | None = None,
) -> str:
    sym = spectral_to_symbol(rec.get("spectral"))
    sym_str = symbol_repeat(sym, float(rec.get("vmag", 5)))
    name = rec.get("name", "?")
    prop = _prop_line(rec)
    coord_main, coord_warning = _coord_line(rec)
    dyad_str = _format_dyads(dyads)
    lines = [
        "\n  \u256d\u2500 TRANSIT \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u256e\n",
        f"  \u2502  {sym_str}  {name}\n",
        f"  \u2502  Stellar: {prop}\n",
    ]
    if coord_main:
        lines.append(f"  \u2502  {coord_main}\n")
    if coord_warning:
        lines.append(coord_warning)
    lines.append(f"  \u2502  Notes: {dyad_str}\n")
    if upcoming:
        parts = [
            f"{n} in {_fmt_delta(d)}" if d > 0 else f"{n} now"
            for n, d in upcoming[:3]
        ]
        lines.append(f"  \u2502  Up next: {', '.join(parts)}\n")
    lines.append("  \u2570\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u256f\n")
    return "".join(lines)


def print_transit(
    rec: dict,
    dyads: list[tuple[int, int]],
    upcoming: list[tuple[str, float]] | None = None,
) -> None:
    print(format_transit_block(rec, dyads, upcoming), flush=True)


def _fmt_delta(sec: float) -> str:
    s = int(max(0, sec))
    h, s = s // 3600, s % 3600
    m, s = s // 60, s % 60
    return f"{h}:{m:02d}:{s:02d}"


def format_next(name: str, seconds: float) -> str:
    s = int(seconds)
    h, s = s // 3600, s % 3600
    m, s = s // 60, s % 60
    return f"\u23f3 Next: {name}  in  {h:02d}:{m:02d}:{s:02d}  "


def countdown(seconds: float, next_name: str) -> None:
    """Sleep for `seconds`, updating a one-line countdown when stdout is a TTY."""
    if seconds <= 0:
        return
    if sys.stdout.isatty():
        s = seconds
        while s > 0:
            print(format_next(next_name, s), end="\r", flush=True)
            time.sleep(min(1.0, s))
            s -= 1
        print(" " * 60, end="\r", flush=True)  # clear line
    else:
        time.sleep(seconds)
