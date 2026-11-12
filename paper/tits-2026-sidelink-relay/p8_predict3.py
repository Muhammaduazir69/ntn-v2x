#!/usr/bin/env python3
"""Model predictions for the v3 campaign (ten starts over a day), one replication per call:
    python3 p8_predict3.py <rep>   ->  out/model3_rep<rep>.json
For the law and for the facade heights ns-3 drew in that replication: delivered share (serving satellite and
any satellite), AoI over the pass (Proposition, eq. aoit), and the relay terms exact in the heights
(Proposition relay reach) for 30, 60 and 120 vehicles."""
import json, os, sys
import numpy as np
import p8_model as P

HERE = os.path.dirname(os.path.abspath(__file__))
EPOCH0, SPACING = 1792800000.0, 8640.0
HZ = json.load(open(os.path.join(HERE, "out", "heights_ns3_v3.json")))
LAW = P.K.Street.F


def realized(z):
    z = np.sort(np.asarray(z))
    def F(self, x):
        return np.searchsorted(self.H50 * np.exp(z), np.asarray(x, float), side="left") / len(z)
    return F


rep = int(sys.argv[1]); res = {}
V = P.views(P.sats(), EPOCH0 + SPACING * rep, 300.0, 2.0)
for inp, F in (("law", LAW), ("realized", realized(HZ[str(rep)]["z"]))):
    P.K.Street.F = F
    for c in P.CANYONS:
        d, _ = P.predict(c, V)
        m = P.predict2(c, V)
        A = P.aoi_timeline(c, V)
        for n in ((30, 60, 120) if c == "urban" else (60,)):
            none, src = P.relay_none_exact(c, V, n, n_pos=48, seed=rep)
            res[f"{inp}_{c}_{n}_{rep}"] = dict(direct=d, multisat=m, aoi_direct=A["direct"], aoi_multisat=A["multisat"],
                                              p_src_blocked=src, p_all_blocked=none,
                                              reach=(src - none) / src if src > 0 else 1.0, relay_bound=1 - none)
            print(inp, c, n, rep, {k: round(v, 3) for k, v in res[f"{inp}_{c}_{n}_{rep}"].items()}, flush=True)
        json.dump(res, open(os.path.join(HERE, "out", f"model3_rep{rep}.json"), "w"), indent=1)
P.K.Street.F = LAW
