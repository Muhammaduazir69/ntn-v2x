#!/usr/bin/env python3
"""Share of a blocked source's messages with no clear neighbor in range: facades shared (Proposition reach) against neighbors
drawn independently, the assumption of relay coverage models that draw blockage per relay. Same geometry and SUMO trajectories,
height law. Writes out/reach_compare.json; ns-3's counterpart is 1 - relay_asked / blocked_gen of the announced rule."""
import json, os
import numpy as np
import p8_model as PM, p8_aoi as A
HERE = os.path.dirname(os.path.abspath(__file__))
EPOCH0, SPACING = 1792800000.0, 4320.0
CY = {"urban": (20, 18), "dense": (15, 30)}
S = PM.sats(); out = {}
for rep in range(10, 20):
    V = PM.views(S, EPOCH0 + SPACING * rep, 300.0, 10.0)
    for c, (W, H) in CY.items():
        st = A.Street(W, H, 0.5, gap_every=240.0, gap_len=18.0)
        tr = A.TraceTraj(st, 60, [os.path.join(HERE, "traces", f"n60_r{rep}.csv")])
        rng = np.random.default_rng(rep)
        sh, ind, wsum = 0.0, 0.0, 0.0
        for E in V:
            if len(E) == 0:
                continue
            for _ in range(8):
                veh = [tr(rng) for _ in range(60)]
                t = float(rng.uniform(10, 300))
                ps, ys = veh[0]
                R0 = A.rays(st, ps(t), ys, E)
                if R0 is None:
                    continue
                top = {}; A.absorb(top, R0); Ps = A.product(st, top)
                allt = dict(top); q_ind = 1.0; dead = False
                for pf, y in veh[1:]:
                    if abs(pf(t) - ps(t)) > 300.0:
                        continue
                    Rn = A.rays(st, pf(t), y, E)
                    if Rn is None:
                        dead = True; q_ind = 0.0; continue
                    A.absorb(allt, Rn)
                    tn = {}; A.absorb(tn, Rn); q_ind *= A.product(st, tn)
                Pall = 0.0 if dead else A.product(st, allt)
                sh += Pall; ind += Ps * q_ind; wsum += Ps
        out[f"{c}_{rep}"] = dict(shared=sh / max(wsum, 1e-12), indep=ind / max(wsum, 1e-12))
    print(rep, {k: {a: round(b, 5) for a, b in v.items()} for k, v in out.items() if k.endswith(str(rep))}, flush=True)
    json.dump(out, open(os.path.join(HERE, "out", "reach_compare.json"), "w"), indent=1)
