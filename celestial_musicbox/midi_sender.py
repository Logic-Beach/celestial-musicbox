"""
Send MIDI dyads (note_on/note_off pairs) to an external synth.
"""

import sys
import time
from typing import Optional

import mido
from mido import Message

# How long the chord (both notes) rings per star
NOTE_DURATION = 0.6

_MIDI_BACKEND_HINT = """
mido needs a MIDI backend. On Linux, use one of:

  1. ALSA (recommended):
       sudo apt install libasound2-dev
       pip install python-rtmidi

  2. PortMidi:
       sudo apt install libportmidi-dev
       MIDO_BACKEND=mido.backends.portmidi python -m celestial_musicbox --lon ... --lat ...
"""


def get_mido_output_names() -> list[str]:
    """List MIDI output port names. Tries PortMidi if rtmidi is missing."""
    try:
        return mido.get_output_names()
    except ModuleNotFoundError as e:
        if getattr(e, "name", None) == "rtmidi" or "rtmidi" in str(e):
            try:
                mido.set_backend("mido.backends.portmidi")
                return mido.get_output_names()
            except Exception:
                raise SystemExit(_MIDI_BACKEND_HINT.strip()) from e
        raise


def open_mido_output(port_name: Optional[str] = None, verbose: bool = False):
    """Open MIDI output. If port_name given, substring-match; else first port."""
    names = get_mido_output_names()
    if not names:
        raise SystemExit("No MIDI output ports found. Connect a synth or create a virtual port.")
    if port_name:
        cand = [n for n in names if port_name.lower() in n.lower()]
        name = cand[0] if cand else names[0]
        if verbose:
            print(f"[midi] matching {port_name!r} → {name}", file=sys.stderr, flush=True)
    else:
        name = names[0]
        if verbose:
            print(f"[midi] using first port: {name}", file=sys.stderr, flush=True)
    return mido.open_output(name)


def send_dyads(
    port: mido.ports.BaseOutput,
    dyads: list[tuple[int, int]],
    *,
    note_duration: Optional[float] = None,
    verbose: bool = False,
) -> None:
    """Play all notes as one simultaneous chord: note_on each, hold, note_off each."""
    nd = note_duration if note_duration is not None else NOTE_DURATION
    if verbose:
        parts = [f"n{n}@v{v}" for n, v in dyads]
        print(f"[midi] send_dyads {parts} dur={nd}s", file=sys.stderr, flush=True)
    for note, vel in dyads:
        port.send(Message("note_on", note=note, velocity=vel))
    time.sleep(nd)
    for note, vel in dyads:
        port.send(Message("note_off", note=note, velocity=0))
