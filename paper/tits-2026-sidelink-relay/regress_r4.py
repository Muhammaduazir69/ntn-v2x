#!/usr/bin/env python3
"""Build 4 (names in the tag, grantRtt passed as 0) against the stored primary rows of builds 1 and 2: same AoI, same
delivered count, same forwards, run by run. Writes out/regress_r4.json."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
import importlib.util
spec = importlib.util.spec_from_file_location("m5", os.path.join(HERE, "make_figs5.py"))
BASE = dict(periodMs=1000, ulKbps=300, grantRtt=1, hoMs=100, slUrban=1, alignResel=1, nbAgeS=1.5, gapEveryM=240, gapLenM=18, mob="sumo")
RES = {"sha256", "wall_s", "error", "generated", "delivered", "rx_rows", "delivered_share", "relay_share", "delay_ms", "aoi_mean_s", "aoi_p90_s",
       "stale_t2", "relay_asked", "relayed", "sent_direct", "blocked_gen", "ul_sent", "ul_drop_q", "ul_drop_los", "ul_bytes", "ul_wait_ms",
       "far_fallback", "p_fail", "p_fail_given_fail", "study", "canyon", "n", "policy", "rep", "cold", "pos_err_m", "pos_p95_m", "pos_gt20",
       "cold_time", "stale_t5"}
conf = lambda r: {k: v for k, v in r.items() if k not in RES and not k.startswith(("prr_", "prrnohd_"))}
old = [json.loads(l) for l in open(os.path.join(HERE, "out", "ns3_canyon5.jsonl"))]
old = {(r["canyon"], r["n"], r["policy"], r["rep"], json.dumps(conf(r), sort_keys=True)): r for r in old if not r.get("error") and r["sha256"] != "b26933c9a1f99d22"}
new = [json.loads(l) for l in open(os.path.join(HERE, "out", "ns3_canyon6_pos.jsonl"))]
out = []
for r in new:
    if r.get("error"):
        continue
    c = conf(r)
    if c.get("kappa") == 1:
        c.pop("kappa")
    o = old.get((r["canyon"], r["n"], r["policy"], r["rep"], json.dumps(c, sort_keys=True)))
    if o is None:
        continue
    same = all(o[k] == r[k] for k in ("aoi_mean_s", "delivered", "relayed", "blocked_gen", "ul_sent"))
    out.append(dict(canyon=r["canyon"], policy=r["policy"], rep=r["rep"], relayK=r.get("relayK"), old_sha=o["sha256"], same=same))
json.dump(out, open(os.path.join(HERE, "out", "regress_r4.json"), "w"), indent=0)
print(len(out), "compared,", sum(x["same"] for x in out), "identical")
for x in out:
    if not x["same"]:
        print("DIFF", x)
