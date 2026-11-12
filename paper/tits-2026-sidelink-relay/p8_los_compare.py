#!/usr/bin/env python3
"""LOS probability against elevation on real Manhattan streets (M3, M6). 'Measured' is the fraction of positions along the
stretch and lanes at which the ray to a satellite at (e, a) clears the real facade, over the satellite directions of twenty
starts. Compared with: TR 38.811 Table 6.6.1-1 (urban, dense urban), the ITU-R P.1410 building-statistics form used by
Al-Hourani et al. (urban, dense urban), our Lemma 1 with the study's urban law, and Lemma 1 with the law fitted to the other
half of the same stretch (held out). Writes out/los_compare.json."""
import json, math, os, sys
import numpy as np
import p8_model as P
HERE = os.path.dirname(os.path.abspath(__file__))
RC = os.path.join(HERE, "realcity")
TR38811 = {"urban": {10:.246,20:.386,30:.493,40:.613,50:.726,60:.805,70:.919,80:.968,90:.992},
           "dense": {10:.282,20:.331,30:.398,40:.468,50:.537,60:.612,70:.738,80:.820,90:.981}}
# Assumed values: in the literature these are attributed to P.1410, but editions -2, -4, -5 and -6 give only the ranges
# (alpha 0.1 to 0.8, beta 750 to 100) and one example (0.11, 750, 7.63). So they are varied below.
P1410 = {"urban": (0.3, 500, 15), "dense": (0.5, 300, 20)}
GRID = [(0.11, 750, 7.63)] + [(a, b, g) for a in (0.1, 0.3, 0.5, 0.8) for b in (100, 300, 500, 750) for g in (8, 15, 20, 30, 50)]


def p1410(e, al, be, ga, h=1.5):
    """Probability that no building crosses a ray from height h at elevation e (ITU-R P.1410, Holis and Pechac form)."""
    q = 1.0
    step = 1000.0 / math.sqrt(al * be)                # mean spacing of buildings along the path
    for n in range(0, 4000):
        z = h + (n + 0.5) * step * math.tan(e)
        q *= 1 - math.exp(-z * z / (2 * ga * ga))
        if z > 6 * ga:
            break
    return q


def lemma(e, a, W, lanes, F, h=1.5):
    v = []
    for y in lanes:
        sa = abs(math.sin(a))
        if sa < 1e-6:
            v.append(1.0); continue
        d = W - y if math.sin(a) > 0 else y
        v.append(F(h + d * math.tan(e) / sa))
    return float(np.mean(v))


def law_F(h50, s, g):
    def F(z):
        return g + (1 - g) * 0.5 * math.erfc(-math.log(max(z, 1e-9) / h50) / (s * math.sqrt(2))) if z > 0 else g
    return F


S = P.sats()
ST = json.load(open(os.path.join(RC, "stretches_ns3.json")))
EPOCH0, SPACING = 1792800000.0, 4320.0
dirs = []                      # (e, a from north) over the starts, every 20 s
for rep in range(20):
    P.AZ = 0.0
    V = P.views(S, EPOCH0 + SPACING * rep, 300.0, 20.0)
    for E in V:
        for e, a in E:
            dirs.append((e, a))
dirs = np.array(dirs)
print("directions", len(dirs), flush=True)
bins = np.arange(20, 91, 10)
res = {}
for name, I in ST.items():
    W, az, L = I["W"], math.radians(I["az"]), int(I["L"])
    hl, hr = np.zeros(L), np.zeros(L)
    for line in open(I["file"]):
        sd, a, b, h = line.split(",")
        (hr if int(sd) else hl)[int(a):int(b)] = float(h)
    lanes = np.array([4.0, 7.5, 12.5, 16.0]) / 20.0 * W
    half = L // 2
    def fit(lo, hi):
        z = np.r_[hl[lo:hi], hr[lo:hi]]
        # per-building law from the 1 m raster: gaps share, then the log law of the positive heights
        g = float(np.mean(z == 0)); pos = z[z > 0]
        return np.exp(np.median(np.log(pos))), float(np.log(pos).std()), g
    h50A, sA, gA = fit(0, half); h50B, sB, gB = fit(half, L)
    F_held = {0: law_F(h50B, sB, gB), 1: law_F(h50A, sA, gA)}   # half 0 is predicted from half 1
    F_std = law_F(18.0, 0.5, 0.0)
    rows = []
    for lo_b, hi_b in zip(bins[:-1], bins[1:]):
        m = (np.degrees(dirs[:, 0]) >= lo_b) & (np.degrees(dirs[:, 0]) < hi_b)
        D = dirs[m]
        if len(D) > 400:
            D = D[np.random.default_rng(1).choice(len(D), 400, replace=False)]
        meas, pred_held, pred_std = [], [], []
        for e, a0 in D:
            a = a0 - az
            sa = math.sin(a)
            for half_id, (x_lo, x_hi) in enumerate(((0, half), (half, L))):
                c = 0; t = 0
                for y in lanes:
                    if abs(sa) < 1e-6:
                        c += x_hi - x_lo; t += x_hi - x_lo; continue
                    d = W - y if sa > 0 else y
                    s = d / abs(sa)
                    xc = (np.arange(x_lo, x_hi) + s * math.cos(a)).astype(int) % L
                    H = hr[xc] if sa > 0 else hl[xc]
                    c += int(np.sum(H < 1.5 + s * math.tan(e))); t += x_hi - x_lo
                meas.append(c / t)
                pred_held.append(lemma(e, a, W, lanes, F_held[half_id]))
                pred_std.append(lemma(e, a, W, lanes, F_std))
        cen = 0.5 * (lo_b + hi_b)
        # TR 38.811 clause 6.6.1: the table value at the reference angle nearest to each direction's elevation
        near = lambda d: int(min(10 * round(d / 10.0), 90)) if d >= 5 else 10
        tr = lambda env: float(np.mean([TR38811[env][near(np.degrees(e))] for e, _ in D]))
        pgrid = lambda p: float(np.mean([p1410(e, *p) for e, _ in D]))
        rows.append(dict(e=float(cen), measured=float(np.mean(meas)), held=float(np.mean(pred_held)), std=float(np.mean(pred_std)),
                         tr_urban=tr("urban"), tr_dense=tr("dense"),
                         p1410_urban=pgrid(P1410["urban"]), p1410_dense=pgrid(P1410["dense"]), **{f"pg_{i}": pgrid(p) for i, p in enumerate(GRID)}))
    res[name] = dict(rows=rows, law=[h50A, sA, gA, h50B, sB, gB])
    err = {k: float(np.mean([abs(r[k] - r["measured"]) for r in rows])) for k in ("held", "std", "tr_urban", "tr_dense", "p1410_urban", "p1410_dense")}
    res[name]["mae"] = err
    gm = [float(np.mean([abs(r[f"pg_{i}"] - r["measured"]) for r in rows])) for i in range(len(GRID))]
    res[name]["grid_best"] = dict(mae=min(gm), params=GRID[int(np.argmin(gm))], median=float(np.median(gm)), worst=max(gm), example=gm[0])
    for r in rows:
        for i in range(len(GRID)):
            r.pop(f"pg_{i}")
    print(name, {k: round(v, 3) for k, v in err.items()}, flush=True)
json.dump(res, open(os.path.join(HERE, "out", "los_compare.json"), "w"), indent=1)
