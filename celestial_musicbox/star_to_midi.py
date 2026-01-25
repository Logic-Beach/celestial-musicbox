"""
Map star properties to MIDI: one chord per star.
Two notes simultaneously: pitch = spectral type, pitch = magnitude; velocity = distance.
"""

import math
import re
from typing import Optional

# Spectral class to numeric: O=0, B=1, A=2, F=3, G=4, K=5, M=6, L=7, T=8
# Subclass 0–9 gives fractional part.
SPECTRAL_ORDER = {"O": 0, "B": 1, "A": 2, "F": 3, "G": 4, "K": 5, "M": 6, "L": 7, "T": 8, "W": 9}
# Roman numerals and extra chars are stripped; we only use the leading letter + digit


def parse_spectral(s: Optional[str]) -> float:
    """Parse spectral type to a number in [0, 10). O0=0, M9≈6.9."""
    if not s or not isinstance(s, str):
        return 4.0  # G2-like default
    s = s.strip().upper()
    m = re.match(r"^([OBAGKMLTW])\s*(\d)?", s)
    if not m:
        return 4.0
    cls = SPECTRAL_ORDER.get(m.group(1), 4)
    sub = int(m.group(2)) if m.group(2) else 5
    return cls + sub / 10.0


def clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def lerp(t: float, a: float, b: float) -> float:
    return a + t * (b - a)


# MIDI note range we use (C2=36 to C6=84)
NOTE_LO, NOTE_HI = 36, 84


def prop_to_note(t: float) -> int:
    """Map t in [0,1] to a note in [NOTE_LO, NOTE_HI]."""
    return int(round(lerp(clamp(t, 0, 1), NOTE_LO, NOTE_HI)))


def prop_to_velocity(t: float, vlo: int = 40, vhi: int = 115) -> int:
    """Map t in [0,1] to velocity. t=1 → vhi."""
    return int(round(lerp(clamp(t, 0, 1), vlo, vhi)))




def magnitude_to_midi(mag: float) -> tuple[int, int]:
    """Brighter (lower mag) → higher note and higher velocity. mag typically -1.5..7."""
    # 0 = very bright (-1), 1 = dim (7)
    t = (clamp(mag, -1.5, 8) + 1.5) / 9.5
    t = 1 - t  # brighter = higher
    return (prop_to_note(t), prop_to_velocity(t))


def mass_to_midi(mass: float) -> tuple[int, int]:
    """Log scale: 0.1–50 solar masses → [0,1]."""
    if mass <= 0:
        mass = 0.1
    x = math.log10(clamp(mass, 0.08, 60))
    # log10(0.1)≈-1, log10(50)≈1.7 → map to [0,1]
    t = (x + 1) / 2.7
    return (prop_to_note(t), prop_to_velocity(t))


def altitude_to_midi(alt_deg: float) -> tuple[int, int]:
    """Altitude in degrees 0–90 → [0,1]. Higher in sky = higher note."""
    t = clamp(alt_deg, 0, 90) / 90.0
    return (prop_to_note(t), prop_to_velocity(t))


def spectral_to_midi(spectral: float) -> tuple[int, int]:
    """O0≈0 to M9≈6.9 → [0,1]. Hot=high note."""
    t = clamp(spectral, 0, 7) / 7.0
    return (prop_to_note(t), prop_to_velocity(t))


def distance_to_midi(distance_ly: float) -> tuple[int, int]:
    """Log scale: 1–10000 ly → [0,1]. Near=high note. Kept for possible reuse."""
    if distance_ly <= 0:
        distance_ly = 10
    x = math.log10(clamp(distance_ly, 1, 15000))
    t = 1 - x / 4.2  # near = 1
    return (prop_to_note(t), prop_to_velocity(t))


def distance_to_velocity(distance_ly: float, vlo: int = 40, vhi: int = 115) -> int:
    """Map distance (ly) to MIDI velocity. Near = high velocity."""
    if distance_ly <= 0:
        distance_ly = 10.0
    x = math.log10(clamp(distance_ly, 1, 15000))
    t = 1 - x / 4.2  # near = 1
    return prop_to_velocity(t, vlo, vhi)


def star_to_dyads(star: dict, supplement: dict) -> list[tuple[int, int]]:
    """
    One chord per star: two notes played simultaneously.
    - Note 1: pitch = spectral type (O→low, M→high in our range)
    - Note 2: pitch = magnitude (dim→low, bright→high)
    - Velocity (both): distance (near→loud, far→quiet)

    Returns 2 (note, velocity) pairs. supplement supplies spectral/distance_ly when star lacks them.
    """
    name = star.get("name", "?")
    sup = supplement.get(name, {}) if isinstance(supplement, dict) else {}

    mag = star.get("vmag") or star.get("vmage")
    if mag is None or not isinstance(mag, (int, float)):
        mag = 5.0

    spectral_str = star.get("spectral") or sup.get("spectral")
    spectral = parse_spectral(spectral_str)

    dist = star.get("distance_ly") or star.get("distance") or sup.get("distance_ly")
    if dist is None or not isinstance(dist, (int, float)) or dist <= 0:
        dist = 50.0

    pitch_spectral = spectral_to_midi(spectral)[0]
    pitch_magnitude = magnitude_to_midi(float(mag))[0]
    vel = distance_to_velocity(float(dist))

    return [(pitch_spectral, vel), (pitch_magnitude, vel)]
