# Code and data: sidelink relaying to direct-to-cell satellites in street canyons (IEEE T-ITS, 2026)

Companion material for M. Uzair, "Which Neighbor Should Forward? Sidelink Relaying of Vehicle Status to Direct-to-Cell Satellites in Urban Street Canyons," manuscript submitted to IEEE Transactions on Intelligent Transportation Systems, 2026.

The street-canyon relay example (`examples/ntn-v2x-canyon-relay.cc`), the urban sidelink channel and the forwarding rules are in this module. This folder holds the paper pipeline.

| Path | Contents |
|---|---|
| `ns3_patch*/` | Source of record of each campaign build (`ns3_patch/` is the base, `_r2` to `_r4` the later builds with their regression checks) |
| `run_ns3_canyon5.py` | Main campaign: 15 forwarding rules, three canyons, twenty starts, plus every sensitivity study (`run_ns3_canyon*.py` are the earlier builds) |
| `make_traces.py`, `traces/` | SUMO 1.18 stop-and-go vehicle traces |
| `realcity/` | Manhattan facade heights from NYC Open Data building footprints and the street stretches used |
| `p8_model.py`, `p8_aoi.py`, `p8_predict*.py`, `p8_theory_check.py`, `mc_aoi.py` | Age-of-information model, its predictions and Monte Carlo checks |
| `out/` | Every result file the paper uses (`ns3_canyon5.jsonl`, `holm5.json`, `model*_rep*.json`, regression checks) |
| `make_macros5.py`, `make_figs5.py` | Numbers and figures of the paper, generated from `out/` |

## Reproduce
Build ns-3.43 with this module in `contrib/`, run `ns3_patch_r4/build.sh`, then `python3 run_ns3_canyon5.py <jobs> primary` (add study names for the sensitivity runs). `python3 make_macros5.py` and `python3 make_figs5.py` regenerate the numbers and figures from `out/`.

License: as the repository. NYC building footprints are NYC Open Data; element sets are public CelesTrak data.
Contact: Muhammad Uzair, muhammaduzairr69@gmail.com, ORCID 0009-0002-4104-2680.
