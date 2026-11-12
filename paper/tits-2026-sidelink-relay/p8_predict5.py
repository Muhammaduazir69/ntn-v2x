#!/usr/bin/env python3
"""Model ages for the ns-3 runs of the revision, per start: python3 p8_predict5.py <rep> -> out/model5_rep<rep>.json
Each case (canyon, period, mobility) with the height law and with the empirical CDF of the run's own facades."""
import json, os, sys, time
import numpy as np
import p8_model as PM
import p8_aoi as A
HERE = os.path.dirname(os.path.abspath(__file__))
EPOCH0, SPACING = 1792800000.0, 4320.0
HZ = json.load(open(os.path.join(HERE, "out", "heights_ns3_v4.json")))
CY = {"suburban": (30, 10), "urban": (20, 18), "dense": (15, 30)}
rep = int(sys.argv[1])
out = os.path.join(HERE, "out", f"model6_rep{rep}.json")
res = json.load(open(out)) if os.path.exists(out) else {}
V = PM.views(PM.sats(), EPOCH0 + SPACING * rep, 300.0, 2.0)
tr = [os.path.join(HERE, "traces", f"n60_r{rep}.csv")]
CASES = [("urban", 1.0, "sumo"), ("dense", 1.0, "sumo"), ("suburban", 1.0, "sumo"), ("urban", 1.0, "const"), ("dense", 1.0, "const"), ("urban", 1.0, "constgap"), ("dense", 1.0, "constgap")]
CASES += [(c, T, "sumo") for T in (0.1, 0.5, 2.0, 5.0) for c in ("urban", "dense")] + [("urban", 0.1, "const0")]
# held-out law from the other half of each real stretch, and the geometry study: model at the same rep
for c, T, mob in CASES:
    for inp in ("law", "realized"):
        k = f"{c}_{T}_{mob}_{inp}"
        if k in res:
            continue
        W, H = CY[c]
        gaps = dict(gap_every=240.0, gap_len=18.0) if mob in ("sumo", "constgap") else {}
        emp = None if inp == "law" else H * np.exp(np.array(HZ[str(rep)]["z"]))
        st = A.Street(W, H, 0.5, emp=emp, **gaps)
        traj = A.TraceTraj(st, 60, tr) if mob == "sumo" else A.ConstTraj(st, 60)
        rng = np.random.default_rng(1000 * rep + 7)
        t0 = time.time()
        res[k] = A.evaluate(st, V, 2.0, T, 0.0, traj, rng, pi_sl=1.0)
        print(k, {a: round(b, 3) for a, b in res[k].items()}, round(time.time() - t0, 1), flush=True)
        json.dump(res, open(out + ".tmp", "w"), indent=1); os.replace(out + ".tmp", out)

# Manhattan stretches: the law fitted to the whole stretch (gaps as buildings of height zero)
ST = json.load(open(os.path.join(HERE, "realcity", "stretches_ns3.json")))
for name, I in ST.items():
    for inp in ("streetlaw",):
        k = f"real_{name}_{inp}"
        if k in res:
            continue
        hs = []
        for line in open(I["file"]):
            sd, a, b, h = line.split(",")
            hs += [float(h)] * (int(b) - int(a))
        hs = np.array(hs)
        g = float(np.mean(hs == 0)); pos = np.log(hs[hs > 0])
        PM.AZ = I["az"]
        V2 = PM.views(PM.sats(), EPOCH0 + SPACING * rep, 300.0, 2.0)
        PM.AZ = 29.0
        st = A.Street(I["W"], float(np.exp(np.median(pos))), float(pos.std()), Lb=I["Lb"], L=I["L"], gap_p=g)
        traj = A.TraceTraj(st, 60, tr)
        res[k] = A.evaluate(st, V2, 2.0, 1.0, 0.0, traj, np.random.default_rng(1000 * rep + 11), pi_sl=1.0)
        print(k, {a: round(b, 3) for a, b in res[k].items()}, flush=True)
        json.dump(res, open(out + ".tmp", "w"), indent=1); os.replace(out + ".tmp", out)

# Geometry study of the urban canyon
for tag, kw in (("az119", dict(AZ=119.0)), ("az74", dict(AZ=74.0)), ("lat25", dict(LAT=25.77, LON=-80.19)), ("lat52", dict(LAT=52.5, LON=13.4)), ("mask10", dict(MASK=10.0)), ("mask30", dict(MASK=30.0))):
    k = f"geo_{tag}"
    if k in res:
        continue
    old = {q: getattr(PM, q) for q in kw}
    for q, v in kw.items():
        setattr(PM, q, v)
    V2 = PM.views(PM.sats(), EPOCH0 + SPACING * rep, 300.0, 2.0)
    for q, v in old.items():
        setattr(PM, q, v)
    st = A.Street(20, 18, 0.5, gap_every=240.0, gap_len=18.0)
    res[k] = A.evaluate(st, V2, 2.0, 1.0, 0.0, A.TraceTraj(st, 60, tr), np.random.default_rng(1000 * rep + 13), pi_sl=1.0)
    print(k, {a: round(b, 3) for a, b in res[k].items()}, flush=True)
    json.dump(res, open(out + ".tmp", "w"), indent=1); os.replace(out + ".tmp", out)
