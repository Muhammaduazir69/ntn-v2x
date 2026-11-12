#!/usr/bin/env python3
"""Write ../macros_p8.tex from out/ns3_canyon.jsonl (ns-3 campaign) and out/model.json (model)."""
import json, os
from collections import defaultdict
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
CAN = {"suburban": "Sub", "urban": "Urban", "dense": "Dense"}
POL = {"direct": "Direct", "multisat": "Multi", "relay": "Relay", "blind": "Blind", "blindfar": "BlindFar"}
DENS = {30: "Thirty", 60: "Sixty", 120: "OneTwenty"}
rows = [json.loads(l) for l in open(os.path.join(OUT, "ns3_canyon.jsonl"))]
rows = [r for r in rows if not r.get("error")]
Mo = json.load(open(os.path.join(OUT, "model.json")))
R = defaultdict(list)
for r in rows:
    R[(r["canyon"], r["n"], r["policy"])].append(r)
M = {}


def f(x, d=1):
    return f"{x:.{d}f}"


def ci(v):
    v = np.asarray(v, float)
    if len(v) < 2:
        return float("nan")
    return float(stats.t.ppf(0.975, len(v) - 1) * v.std(ddof=1) / np.sqrt(len(v)))


nrep = min(len(v) for v in R.values()) if R else 0
M["NRepsDone"] = str(nrep)
for c, tc in CAN.items():
    for n, tn in DENS.items():
        for p, tp in POL.items():
            v = R.get((c, n, p), [])
            if not v:
                continue
            M[f"Del{tp}{tc}{tn}"] = f(100 * np.mean([x["delivered_share"] for x in v]), 1)
            M[f"Aoi{tp}{tc}{tn}"] = f(np.mean([x["aoi_mean_s"] for x in v]), 1)
            M[f"AoiCi{tp}{tc}{tn}"] = f(ci([x["aoi_mean_s"] for x in v]), 1)
            M[f"Stale{tp}{tc}{tn}"] = f(100 * np.mean([x["stale5_share"] for x in v]), 0)
# headline aliases at 60 vehicles
for c, tc in CAN.items():
    for p, tp in POL.items():
        if f"Aoi{tp}{tc}Sixty" in M:
            M[f"Aoi{tp}{tc}"] = M[f"Aoi{tp}{tc}Sixty"]; M[f"Del{tp}{tc}"] = M[f"Del{tp}{tc}Sixty"]
# model against ns-3, same replication (same start time and, for "realized", the same facade draw)
err = defaultdict(list); aerr = defaultdict(list); aabs = defaultdict(list)
for r in rows:
    if r["policy"] not in ("direct", "multisat"):
        continue
    for inp in ("law", "realized"):
        k = f"{inp}_{r['canyon']}_{r['n']}_{r['rep']}"
        if k not in Mo:
            continue
        err[(inp, r["policy"])].append(abs(Mo[k][r["policy"]] - r["delivered_share"]))
        if r["policy"] in ("direct", "multisat"):
            m = Mo[k][f"aoi_{r['policy']}"]
            aerr[(inp, r["policy"])].append(abs(m - r["aoi_mean_s"]) / r["aoi_mean_s"]); aabs[(inp, r["policy"])].append(abs(m - r["aoi_mean_s"]))
for (inp, p), v in err.items():
    t = "Law" if inp == "law" else "Real"
    M[f"DelErrMax{t}{POL[p]}"] = f(100 * max(v), 1); M[f"DelErrMean{t}{POL[p]}"] = f(100 * np.mean(v), 1)
for (inp, p), v in aabs.items():
    t = "Law" if inp == "law" else "Real"
    M[f"AoiAbsErrMean{t}{POL[p]}"] = f(np.mean(v), 1); M[f"AoiAbsErrMax{t}{POL[p]}"] = f(max(v), 1)
for (inp, p), v in aerr.items():
    t = "Law" if inp == "law" else "Real"
    M[f"AoiRelErrMean{t}{POL[p]}"] = f(100 * np.mean(v), 0); M[f"AoiRelErrMax{t}{POL[p]}"] = f(100 * max(v), 0)
if ("realized", "multisat") in err:
    M["AnyErrMax"] = M["DelErrMaxRealMulti"]
# blindfar against blind: AoI cut, paired by canyon, density, replication
cut = []
for r in rows:
    if r["policy"] != "blindfar":
        continue
    b = [x for x in R[(r["canyon"], r["n"], "blind")] if x["rep"] == r["rep"]]
    if b:
        cut.append(1 - r["aoi_mean_s"] / b[0]["aoi_mean_s"])
if cut:
    M["BlindFarGain"] = f(100 * np.mean(cut), 0); M["BlindFarGainMin"] = f(100 * min(cut), 0); M["BlindFarGainMax"] = f(100 * max(cut), 0)
    M["BlindFarPairs"] = str(len(cut))
# relay: reach of Proposition (relay reach) against the share of blocked messages for which ns-3 named a relay,
# the bound 1 - P_none against the delivered share, and the measured sidelink-and-staleness success
RL = {}
for r_ in range(3):
    pth = os.path.join(OUT, f"relay_rep{r_}.json")
    if os.path.exists(pth):
        RL.update(json.load(open(pth)))
reach_err, bound_gap, succ, naive_err = [], [], [], []
for r in rows:
    if r["policy"] != "relay":
        continue
    k = f"realized_{r['canyon']}_{r['n']}_{r['rep']}"
    m = [x for x in R[(r["canyon"], r["n"], "multisat")] if x["rep"] == r["rep"]]
    if k not in RL or not m:
        continue
    blocked = r["generated"] * (1 - m[0]["delivered_share"])
    if blocked > 100:
        reach_err.append(abs(RL[k]["reach"] - r["relay_asked"] / blocked))
        if k in Mo:
            naive_err.append(Mo[k]["relay_avail"] - r["relay_asked"] / blocked)
        succ.append(r["relayed"] / max(r["relay_asked"], 1))
    bound_gap.append(RL[k]["relay_bound"] - r["delivered_share"])
if reach_err:
    M["ReachErrMax"] = f(100 * max(reach_err), 0); M["ReachErrMean"] = f(100 * np.mean(reach_err), 1)
    M["RelaySuccLo"], M["RelaySuccHi"] = f(100 * min(succ), 0), f(100 * max(succ), 0)
    M["BoundGapMin"], M["BoundGapMax"] = f(100 * min(bound_gap), 1), f(100 * max(bound_gap), 1)
    M["NaiveOverMax"] = f(100 * max(naive_err), 0); M["NaiveOverMean"] = f(100 * np.mean(naive_err), 0)
# Proposition (pair): conditional against unconditional LOS of a neighbor 5 m away, urban canyon, serving
# satellite of the first replication's pass
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
        for _y2 in _st.lanes:   # neighbor in any lane, 5 m ahead or behind
            _p2 = float(_st.p_los_y(_e, _a, _y2))
            for _dx in (5.0, -5.0):
                _num += _p2 - _st.p_both(_e, _a, _y, _y2, _dx); _unc += _p2 * (1 - _p1)
M["PairRatioFive"] = f(_num / _unc, 2)   # conditional / unconditional, averaged over lanes and both senses
with open(os.path.join(HERE, "..", "macros_p8.tex"), "w") as fh:
    fh.write("% generated by sim/make_macros.py; do not edit\n")
    for k, v in sorted(M.items()):
        fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
print(len(M), "macros;", "reps done", nrep)
