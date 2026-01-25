# celestial-musicbox

Meridian transits → MIDI. More stars = more notes.

## 1. Build catalog

Download HYG and build `data/star_catalog.json`:

```bash
python scripts/build_star_catalog.py --lat 36
# optional: --max-mag 7 (fewer, brighter), --hyg path/to/hyg.csv
```

`--lat 36` keeps only stars that can rise at 36°N (dec ≥ −54°).

## 2. Run

```bash
pip install -r requirements.txt
python -m celestial_musicbox --lon -122.4 --lat 36 --midi-port "part of your port name"
```

`--catalog` defaults to `data/star_catalog.json`. **Don’t know your MIDI port name?** Plug in the device and run `python -m celestial_musicbox --list-ports`; use any part of the line that’s your box with `--midi-port` (substring match). Or omit `--midi-port` to use the first port. `MIDI_PORT=…` in the environment works too.

The terminal shows each transit: ASCII star by spectral type (O→✦ B→★ A→⁑ F→● G→• K→◦ M→○), properties (vmag, mass/alt, spectral, dist), and the MIDI dyads (note@velocity) for mag, mass, spec, dist. A countdown runs until the next transit. Use `--quiet` to disable.

### Time sync (NTP) — important for Stellarium

**NTP** (Network Time Protocol) keeps your system clock in sync with atomic clocks on the internet. If the PC clock is off by even 30–60 seconds, transits will fire early or late and Stellarium will disagree with the real sky.

**Check status:**
```bash
timedatectl
```
Look for `System clock synchronized: yes` and `NTP service: active`. If not:

**Enable NTP (systemd, most Linux):**
```bash
sudo timedatectl set-ntp true
```
If your distro uses **chrony** instead:
```bash
sudo systemctl enable --now chronyd
```

**If you can't use NTP** (air‑gapped or no net): set the clock manually and avoid long runs—drift will build up. In Stellarium: **Time rate = 1** and time set to “real time” so it matches the system.

### MIDI on Linux

**Finding your MIDI port name (plug in the device first):**
```bash
python -m celestial_musicbox --list-ports
```
You’ll see one line per output, e.g. `Midi Through Port-0` or `audioUSB MIDI 1`. Use **any part** of the correct line with `--midi-port` (substring match), e.g. `--midi-port audioUSB` or `--midi-port "MIDI 1"`.

**If that fails**, try `aconnect -l` (from `alsa-utils`). Look for a client name that matches your hardware; mido’s names are usually similar. Then run `sudo modprobe snd-seq` and `--list-ports` again.

**Using the name:** `--midi-port` does a substring match, so you don’t need the exact string. Or set `export MIDI_PORT=audioUSB` and omit `--midi-port`. If you omit both, the first available port is used.

**If `open /dev/snd/seq failed`:** load the ALSA sequencer and ensure you’re in the `audio` group:
```bash
sudo modprobe snd-seq
groups  # should include audio; if not: sudo usermod -aG audio $USER then log out/in
```

**Backend / build:** If you get `No module named 'rtmidi'` or an ALSA build error:

1. **ALSA (recommended):** `sudo apt install libasound2-dev` then `pip install python-rtmidi`
2. **PortMidi:** `sudo apt install libportmidi-dev` and run with `MIDO_BACKEND=mido.backends.portmidi`
