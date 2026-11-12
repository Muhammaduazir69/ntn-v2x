#!/usr/bin/env python3
"""Model ages of the announced family and of thinning at the primary geometry: out/model7_rep<rep>.json."""
import json, os, sys
import numpy as np
import p8_model as PM, p8_aoi as A
HERE = os.path.dirname(os.path.abspath(__file__))
rep = int(sys.argv[1]); res = {}
V = PM.views(PM.sats(), 1792800000.0 + 4320.0 * rep, 300.0, 2.0)
for c, W, H in (("urban", 20, 18), ("dense", 15, 30)):
    st = A.Street(W, H, 0.5, gap_every=240.0, gap_len=18.0)
    res[c] = A.evaluate(st, V, 2.0, 1.0, 0.0, A.TraceTraj(st, 60, [os.path.join(HERE, "traces", f"n60_r{rep}.csv")]), np.random.default_rng(1000 * rep + 7))
json.dump(res, open(os.path.join(HERE, "out", f"model7_rep{rep}.json"), "w"), indent=1)
