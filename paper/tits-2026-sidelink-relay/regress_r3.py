#!/usr/bin/env python3
"""Does the round 2 binary reproduce stored rows of the first binary? Same flags, compare the recorded results. Exact equality."""
import json, sys, random
sys.argv = ["x"]
import run_ns3_canyon5 as R
from concurrent.futures import ThreadPoolExecutor
rows = [json.loads(l) for l in open(R.OUT)]
rows = [r for r in rows if r["sha256"] == R.SHA[:16] and r["n"] == 60 and r["periodMs"] >= 1000]
random.seed(8)
pick = []
for pol in ("relay", "anyp", "anyfwd", "random", "rand3", "randpersist", "multisat", "blindbest"):
    pick.append(random.choice([r for r in rows if r["policy"] == pol]))
pick.append(random.choice([r for r in rows if r.get("mob") == "const"]))
pick.append(random.choice([r for r in rows if r.get("relaySetupMs")]))
pick.append(random.choice([r for r in rows if r["canyon"] == "real"]))
pick.append(random.choice([r for r in rows if r.get("ulKbps") == 0]))
def one(r):
    j = {k: r[k] for k in R.KEYS(r)}
    j["study"] = "r3regress"; j["grantRtt"] = 0
    n = R.run(j)
    same = all(n.get(k) == r.get(k) for k in ("generated", "delivered", "rx_rows", "relayed", "blocked_gen", "ul_sent", "ul_drop_q", "aoi_mean_s", "p_fail", "delay_ms"))
    print(r["study"], r["canyon"], r["policy"], r["rep"], "same" if same else "DIFF", n.get("sha256"), r["aoi_mean_s"], n.get("aoi_mean_s"), n.get("error", "")[:80], flush=True)
    return same
with ThreadPoolExecutor(3) as ex:
    res = list(ex.map(one, pick))
print("all same" if all(res) else "NOT ALL SAME", len(res))
json.dump([dict(study=p["study"], canyon=p["canyon"], policy=p["policy"], rep=p["rep"], same=x) for p, x in zip(pick, res)], open(R.os.path.join(R.HERE, "out", "regress_r3.json"), "w"), indent=1)
