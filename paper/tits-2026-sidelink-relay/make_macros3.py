#!/usr/bin/env python3
"""Write ../macros_p8.tex from the v3 campaign (out/ns3_canyon3.jsonl, ten starts over a day) and the model
(out/model3_rep*.json). Statistics: medians over the ten starts, min and max, paired Wilcoxon signed-rank
tests on log AoI between policies (paired by start, which fixes the facade draw and the satellite pass)."""
import glob, json, os
from collections import defaultdict
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
CAN = {"suburban": "Sub", "urban": "Urban", "dense": "Dense"}
POL = {"direct": "Direct", "multisat": "Multi", "relay": "Relay", "blind": "Blind", "blindfar": "BlindFar",
       "blindlb": "BlindLb", "random": "Random", "blindbest": "BlindBest"}
rows = [json.loads(l) for l in open(os.path.join(OUT, "ns3_canyon3.jsonl"))]
rows = [r for r in rows if not r.get("error")]
Mo = {}
for p in glob.glob(os.path.join(OUT, "model3_rep*.json")):
    Mo.update(json.load(open(p)))
M = {}


def f(x, d=1):
    return f"{x:.{d}f}"


def sel(c, n, p, **ex):
    out = {}
    for r in rows:
        if r["canyon"] == c and r["n"] == n and r["policy"] == p and all(r.get(k) == v for k, v in ex.items()) \
                and not any(k in r for k in ("hoMs", "periodMs") if k not in ex):
            out[r["rep"]] = r
    return out


M["NStarts"] = str(len({r["rep"] for r in rows}))
for c, tc in CAN.items():
    for p, tp in POL.items():
        R = sel(c, 60, p)
        if not R:
            continue
        d = np.array([r["delivered_share"] for r in R.values()]); a = np.array([r["aoi_mean_s"] for r in R.values()])
        M[f"Del{tp}{tc}"] = f(100 * np.median(d), 1); M[f"Aoi{tp}{tc}"] = f(np.median(a), 1)
        M[f"AoiLo{tp}{tc}"], M[f"AoiHi{tp}{tc}"] = f(a.min(), 1), f(a.max(), 1)
        M[f"DelLo{tp}{tc}"], M[f"DelHi{tp}{tc}"] = f(100 * d.min(), 0), f(100 * d.max(), 0)
# paired tests on log AoI
def paired(c, pa, pb, n=60):
    A, B = sel(c, n, pa), sel(c, n, pb); ks = sorted(set(A) & set(B))
    x = np.log(np.array([A[k]["aoi_mean_s"] for k in ks])); y = np.log(np.array([B[k]["aoi_mean_s"] for k in ks]))
    p = stats.wilcoxon(x, y).pvalue if len(ks) >= 5 and np.any(x != y) else float("nan")
    return float(np.exp(np.median(y - x))), float(p), len(ks)
for c, tc in CAN.items():
    for (pa, pb, nm) in (("multisat", "relay", "RelayVsMulti"), ("blind", "blindfar", "FarVsBlind"),
                         ("direct", "multisat", "MultiVsDirect"), ("blind", "relay", "RelayVsBlind"),
                         ("blindlb", "blindfar", "FarVsLb"), ("random", "blindfar", "FarVsRandom"),
                         ("blind", "blindbest", "BestVsBlind"), ("random", "blindbest", "BestVsRandom"),
                         ("relay", "blindbest", "BestVsRelay"), ("blind", "random", "RandomVsBlind"),
                         ("relay", "random", "RandomVsRelay")):
        ratio, p, k = paired(c, pa, pb)
        M[f"Ratio{nm}{tc}"] = f(ratio, 2); M[f"P{nm}{tc}"] = f(p, 3) if np.isfinite(p) else "--"
        A_, B_ = sel(c, 60, pa), sel(c, 60, pb); ks = sorted(set(A_) & set(B_))
        dlog = np.log([B_[k]["aoi_mean_s"] / A_[k]["aoi_mean_s"] for k in ks])
        rng = np.random.default_rng(3)
        bs = [np.exp(np.median(rng.choice(dlog, len(dlog)))) for _ in range(2000)]
        lo_, hi_ = np.percentile(bs, [2.5, 97.5])
        M[f"RatioCi{nm}{tc}"] = f"{lo_:.2f}--{hi_:.2f}"
# model against ns-3, per start, serving and any satellite at 60 vehicles
for inp, t in (("realized", "Real"), ("law", "Law")):
    for p, tp in (("direct", "Direct"), ("multisat", "Multi")):
        de, ae, ar = [], [], []
        for r in rows:
            if r["policy"] != p or r["n"] != 60 or "hoMs" in r or "periodMs" in r:
                continue
            k = f"{inp}_{r['canyon']}_60_{r['rep']}"
            if k not in Mo:
                continue
            de.append(abs(Mo[k][p] - r["delivered_share"])); a = Mo[k][f"aoi_{p}"]
            ae.append(abs(a - r["aoi_mean_s"])); ar.append(abs(a - r["aoi_mean_s"]) / r["aoi_mean_s"])
        if de:
            M[f"DelErrMax{t}{tp}"], M[f"DelErrMean{t}{tp}"] = f(100 * max(de), 1), f(100 * np.mean(de), 1)
            M[f"AoiAbsErrMed{t}{tp}"], M[f"AoiAbsErrMax{t}{tp}"] = f(np.median(ae), 1), f(max(ae), 1)
            M[f"AoiRelErrMed{t}{tp}"] = f(100 * np.median(ar), 0)
            M[f"AoiRelErrMax{t}{tp}"] = f(100 * max(ar), 0)
# relay reach: model (realized) against the share of blocked messages for which ns-3 named a relay, both from
# the relay run itself (blocked = generated - sent directly)
re_, re_inf, succ = [], [], []
for r in rows:
    if r["policy"] != "relay" or "hoMs" in r or "periodMs" in r:
        continue
    k = f"realized_{r['canyon']}_{r['n']}_{r['rep']}"
    if k not in Mo:
        continue
    direct = r["sent_direct"] if r.get("sent_direct", -1) >= 0 else r["delivered"] * (1 - r["relay_share"])
    blocked = r["generated"] - direct
    if blocked > 100:
        e = Mo[k]["reach"] - r["relay_asked"] / blocked
        re_.append(e)
        if Mo[k]["reach"] < 0.99:
            re_inf.append(e)
        succ.append(r["relayed"] / max(r["relay_asked"], 1))
if re_:
    M["ReachErrMed"] = f(100 * np.median(np.abs(re_)), 1); M["ReachErrMax"] = f(100 * np.max(np.abs(re_)), 0)
    M["ReachBias"] = f(100 * np.median(re_), 1)
    M["RelaySuccLo"], M["RelaySuccHi"] = f(100 * min(succ), 0), f(100 * max(succ), 0)
if re_inf:
    M["ReachErrInfMed"] = f(100 * np.median(np.abs(re_inf)), 1); M["ReachErrInfMax"] = f(100 * np.max(np.abs(re_inf)), 0)
    M["ReachInfN"] = str(len(re_inf))
# density, urban
for n, tn in ((30, "Thirty"), (60, "Sixty"), (120, "OneTwenty")):
    for p, tp in (("relay", "Relay"), ("blind", "Blind"), ("blindfar", "BlindFar"), ("blindbest", "BlindBest")):
        R = sel("urban", n, p)
        if R:
            M[f"Aoi{tp}Urban{tn}"] = f(np.median([r["aoi_mean_s"] for r in R.values()]), 2)
# sensitivity: handover interruption and one-second period, urban, 60 vehicles
for ex, te in (({"hoMs": 100}, "Ho"), ({"periodMs": 1000}, "Slow")):
    for p, tp in (("direct", "Direct"), ("multisat", "Multi"), ("relay", "Relay")):
        R = sel("urban", 60, p, **ex)
        if R:
            M[f"Aoi{te}{tp}"] = f(np.median([r["aoi_mean_s"] for r in R.values()]), 1)
            M[f"Del{te}{tp}"] = f(100 * np.median([r["delivered_share"] for r in R.values()]), 1)
prr = [r["prr_100"] for r in rows if r["policy"] in ("relay", "blind", "blindfar")]
M["PrrMin"] = f(100 * min(prr), 1)
M["AnyErrMax"] = M.get("DelErrMaxRealMulti", "--")
# AoI error above the floor delta + T/2, any satellite, per canyon (realized input)
FLOOR = 0.005 + 0.05
for c_, tc in CAN.items():
    ab = []
    for r in rows:
        if r["policy"] != "multisat" or r["n"] != 60 or r["canyon"] != c_ or "hoMs" in r or "periodMs" in r:
            continue
        k = f"realized_{c_}_60_{r['rep']}"
        if k in Mo and r["aoi_mean_s"] - FLOOR > 0.05:
            ab.append(abs(Mo[k]["aoi_multisat"] - r["aoi_mean_s"]) / (r["aoi_mean_s"] - FLOOR))
    if ab:
        M[f"AoiAboveFloorErr{tc}"] = f(100 * np.median(ab), 0)
# Holm over all paired tests
pv = sorted(float(v) for k, v in M.items() if k.startswith("P") and k[1:].startswith(("RelayVs", "FarVs", "MultiVs")) and v != "--")
holm_ok = all(p_ <= 0.05 / (len(pv) - i) for i, p_ in enumerate(pv))
M["HolmTests"] = str(len(pv)); M["HolmAll"] = "all" if holm_ok else "not all"
# shared facades against independent rays (any-satellite delivered share), law input, first start
import p8_model as PM2
_V0 = PM2.views(PM2.sats(), 1792800000.0, 300.0, 2.0)
for c_, tc in CAN.items():
    _, ind = PM2.predict(c_, _V0)
    M[f"IndepMinusShared{tc}"] = f(100 * (ind - Mo[f"law_{c_}_60_0"]["multisat"]), 1)
# Proposition (pair): conditional / unconditional LOS of a neighbor 5 m ahead or behind, any lane, urban canyon,
# serving satellite of the first start's pass
import p8_model as PM
_V = PM.views(PM.sats(), 1792800000.0, 300.0, 2.0)[::5]
_st = PM.K.Street(W=20, H50=18, sH=0.5, Lb=30.0, lanes=np.array([4.0, 7.5, 12.5, 16.0]))
_num = _unc = 0.0
for _E in _V:
    if len(_E) == 0:
        continue
    _e, _a = _E[np.argmax(_E[:, 0])]
    for _y in _st.lanes:
        _p1 = float(_st.p_los_y(_e, _a, _y))
        if _p1 >= 0.999:
            continue
        for _y2 in _st.lanes:
            _p2 = float(_st.p_los_y(_e, _a, _y2))
            for _dx in (5.0, -5.0):
                _num += _p2 - _st.p_both(_e, _a, _y, _y2, _dx); _unc += _p2 * (1 - _p1)
M["PairRatioFive"] = f(_num / _unc, 2)
with open(os.path.join(HERE, "..", "macros_p8.tex"), "w") as fh:
    fh.write("% generated by sim/make_macros3.py; do not edit\n")
    for k, v in sorted(M.items()):
        fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
print(len(M), "macros")
