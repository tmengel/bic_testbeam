# decoding/rawcheck: raw-byte format checks for the run 3062 streams

These tools read the raw packet words (`Packet::copyMe`) next to the stock decoder output to
document the on-the-wire format. **They are diagnostics, not pipeline code**: the analysis
pipeline in `sync/` uses only the stock `online_distribution` decoders, by the analyst's
requirement. Results are written up in `docs/data_format.md`.

    source ../setup.sh && make

| Tool | What it does |
| --- | --- |
| `rawfmt <kcu\|caen\|astro> file.evt nevt` | field-by-field histograms of the UDP/line headers (KCU), event/group headers (CAEN) or record lengths and hit fields (AstroPix) |
| `kcu2 file.evt nevt` | KCU line-header study (ASIC header word, CMN, CALIB, Tc/Tp flags) plus TOA/TOT occupancy inside stock-decoded waveforms |
| `kcu3 file.evt nevt tot.txt hdr.txt` | full 10-bit TOT from the raw word next to the stock 6-bit TOT, with the number of saturated ADC samples; and a dump of (ts, link, BX, ev6, evctr, trig_out) per line |
| `hexpkt file.evt pktid skip nshow w0 w1 [bytes]` | hex dump of one packet's payload words or bytes |
| `pktbin file.evt nevt out.bin pktid` | writes raw payloads to a binary file for python inspection |

Key finding (2026-10-21): the stock KCU decoder masks TOT with `0x3f`; the field is 10 bits
(19:10) in the HGCROC3 compressed form, and the full value tracks pulse size (Spearman 0.74/0.89
against the number of saturated samples on KCU0/KCU1, versus 0.50/0.20 for the 6-bit value).
