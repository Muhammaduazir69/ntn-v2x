#!/usr/bin/env python3
"""Relay terms of the model, exact in the facade heights (p8_model.relay_none_exact), per replication:
writes out/relay_rep<r>.json with, per input/canyon/density, P(source blocked), P(all blocked), the share of
blocked messages a relay can reach, and the delivered-share bound 1 - P(all blocked)."""
import json, os, sys
import p8_model as P
from p8_predict import realized, HZ, LAW, DENS

HERE = os.path.dirname(os.path.abspath(__file__))
rep = int(sys.argv[1]); out = {}
V = P.views(P.sats(), 1792800000.0 + 977.0 * rep, 300.0, 2.0)
for inp, F in (("law", LAW), ("realized", realized(HZ[str(rep)]["z"]))):
    P.K.Street.F = F
    for c in P.CANYONS:
        for n in DENS:
            none, src = P.relay_none_exact(c, V, n, n_pos=48, seed=rep)
            out[f"{inp}_{c}_{n}_{rep}"] = dict(p_src_blocked=src, p_all_blocked=none,
                                               reach=(src - none) / src if src > 0 else 1.0, relay_bound=1 - none)
            print(inp, c, n, rep, round(src, 3), round(none, 3), flush=True)
            json.dump(out, open(os.path.join(HERE, "out", f"relay_rep{rep}.json"), "w"), indent=1)
