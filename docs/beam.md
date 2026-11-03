# The beam in run 3062 (and runs 3060-3064): species, energy, rate

Sources: the Hall D run log (`context/runlog_and_map.pdf`), the TIC deck of 2026-09-15
(`context/HallD-BIC-09-15-2026.pdf`, slides 2-7), the September analysis decks, and the run-3062 data
(`sync/beam_rate.py`, `decoder_tests/kcu_event_anatomy`, `docs/method.md`). Items marked "inferred"
are not in the log and should be confirmed with the Hall D beam people (Sasha Somov).

## Where and what

- **Location:** the positron arm of the Hall D Pair Spectrometer (PS), parasitic to GlueX running.
  The CEBAF electron beam (12 GeV era, continuous-wave, 2 ns bunch spacing, so no spill structure)
  hits the GlueX radiator and the bremsstrahlung photon beam converts in the TPOL converter; the PS
  dipole sweeps the e+ e- pairs and the positron side exits into the test area where the setup stood.
- **Species:** positrons. Contamination is small by construction (a charged beam of one sign out of a
  photon beam) but was not measured here; muons and pions are negligible at a pair-spectrometer arm,
  and the only other particles reaching the paddles are low-energy background and showers off the
  collimators and the hall floor. The "paddle-only" triggers (39-46% of all triggers, no SFIL or Baby
  BCAL pulse, `method.md` 4.7) are particles crossing a paddle outside the calorimeter acceptance, not
  another species.
- **Energy:** the run-log settings for run 3062 are PS field 1.662 T (the "300 A" setting), radiator
  JD80-213 (80 um diamond), TPOL converter B (Be, about 75 um) at 70.000 mm, PS radiator out, electron
  beam current 191 nA. At this field the positron arm accepts about 3 to 6 GeV (the collaboration's
  "3-6 GeV e+" for the August run; at the nominal 1.8 T the PS arms cover 3.0-6.2 GeV, scaling with
  the field: inferred). The energy of a given positron is fixed by where it crosses the PS focal
  plane, i.e. by its position along the module: in May the four paddles (then calibrated against the PS
  counters) tagged 0.48 / 0.55 / 0.57 / 0.68 GeV at 0.19 T for P4 / P3 / P6 / P5, rising from North
  to South. **In August the paddles were not calibrated, so run 3062 has no per-event energy;** the
  paddle that fired is only an ordered tag (P4 < P3 < P6 < P5 in energy). The spread inside one paddle
  is set by the paddle width and the PS dispersion and is also unknown here. A per-event energy needs
  the PS hodoscope information (the PS-A/B counters), which was not read into RCDAQ.
- **Profile at the setup:** the beam band on the AstroPix module is at y = 6-11 mm (rows 12-22) and
  spans all nine chips along x, densest on chips 3, 5, 6 and 8 (`channel_map.md` 3.5); on the SFILs
  the beam row is cell 3 of every layer. The paddle pattern of every trigger places the positron along
  the module: P3 27.8%, P4 25.4%, P5 35.1%, P6 27.3% of triggers, with two paddles in 9.5% and
  three or four in 2.6% (`method.md` 4.11).

## Rate

| quantity | value | source |
| --- | --- | --- |
| trigger = paddle OR (R580 PMTs, NIM) | issued at 4.10 kHz, constant over the run | `TRIGGER_OUT` counter of both KCUs |
| issued triggers over the 1486 s of KCU0 data | 6.06 M | 4.10 kHz x 1479-1486 s |
| KCU0 recorded triggers | 3.77 M (2.53 kHz mean, 2.1-3.1 kHz 10-s bins); 88% complete waveforms | `kcu0_3062.root` |
| KCU1 recorded | 1.53 M (1.04 kHz) | `kcu1_3062.root` |
| CAEN recorded | 1.51 M (1.03 kHz, 1.00 kHz median) | `caen_3062.root` |
| AstroPix | self-triggered, 2.86 M hit records, about 1.9 kHz of records | `astropix_3062.root` |
| minimum spacing of recorded triggers | 100 us (trigger busy) | `method.md` 4.11 |
| triggers with a shower in the beam row | 65% of KCU0 triggers; 81% have any SFIL pulse | `method.md` 4.11 |
| pile-up | 2-2.5% of beam-row waveforms start on the tail of an earlier particle | `method.md` 4.10 |

Each board records what its own dead time allows: the KCUs lose triggers to UDP loss and to the
100 us busy, the CAEN to its 1 kHz readout limit, so the three recorded streams are different
subsets of the same 4.1 kHz trigger (see `completeness_plan.md`). The recorded rates versus time are in
`docs/beam_rate.png`: KCU0 shows a dip at 540-560 s and the end of the run at 1440 s, the CAEN and
KCU1 ride at their readout limits throughout. The 191 nA electron current and the 80 um radiator
set the photon flux and hence this positron rate; the May runs at 0.19 T and 0.2-0.3 nA gave a few
Hz to a few hundred Hz.

## What this means for the analysis

- Energy resolution and linearity need the beam energy per event, which run 3062 does not have; the
  options are (a) the May-style paddle calibration against the PS counters for a future run, (b)
  reading the PS hodoscope into the DAQ, or (c) using the AstroPix x position as the energy tag
  (dispersion along the module), calibrated once against the PS. Until then run 3062 gives shape and
  relative-response measurements, and the simulation can be compared at 3-6 GeV as a band.
- The 4.1 kHz rate with 100 us busy is far above the CAEN capability and above what the KCU links
  deliver loss-free; the DAQ side of `completeness_plan.md` applies to the next run.
- Pile-up at 2% and the per-event baseline shift are rate effects that the per-event pedestal handles.
