#!/usr/bin/env python3
"""ns-3 campaign of the T-ITS study: ntn-v2x-canyon-relay over canyon types, densities and policies.
Resumable (one jsonl line per run). Per run: delivered share, share via relay, delay, and the age of
information per vehicle at the server (time average after a 10 s warm-up)."""
import csv, itertools, json, os, subprocess, sys, tempfile, threading
from concurrent.futures import ThreadPoolExecutor
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NS3 = "/home/uzair/6g_ntn_ns3/ns-3-dev"
EXE = os.path.join(NS3, "build/contrib/ntn-v2x/examples/ns3.43-ntn-v2x-canyon-relay-optimized")
OUT = os.path.join(HERE, "out", "ns3_canyon.jsonl")
CANYONS = {"suburban": dict(width=30, h50=10), "urban": dict(width=20, h50=18), "dense": dict(width=15, h50=30)}
DENS = (30, 60, 120)                       # vehicles on 1.2 km
POLICIES = ("direct", "multisat", "relay", "blind", "blindfar")
REPS, SIM, WARM = 3, 300.0, 10.0


def aoi(gen, rx, t0, t1, thresholds=(1.0, 5.0)):
    """Age of information of one vehicle over [t0, t1]: time average and share of time above each threshold.
    The age at t0 counts from the newest message received by then, or from the start of the run (generation
    time 0) if none has arrived, so a vehicle that is never served is counted too, the same way under every
    policy."""
    order = np.argsort(rx)
    gen, rx = gen[order], rx[order]
    before = rx <= t0
    cur = gen[before].max() if before.any() else 0.0
    segs, last = [], t0
    for g, r in zip(gen, rx):
        if r <= t0 or g <= cur:
            continue
        if r > t1:
            break
        segs.append((last - cur, r - cur))
        last, cur = r, g
    segs.append((last - cur, t1 - cur))
    span = t1 - t0
    mean = sum((b * b - a * a) / 2 for a, b in segs) / span
    above = [sum(max(0.0, b - max(a, x)) for a, b in segs) / span for x in thresholds]
    return mean, above


def run(job):
    canyon, n, pol, rep = job
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tf:
        out = tf.name
    a = dict(CANYONS[canyon], nVeh=n, policy=pol, simSeconds=SIM, out=out, RngRun=rep + 1,
             epochUnix=1792800000.0 + 977.0 * rep)
    cmd = ["nice", "-n", "5", EXE] + [f"--{k}={v}" for k, v in a.items()]
    r = subprocess.run(cmd, cwd=NS3, capture_output=True, text=True, timeout=7200)
    base = dict(canyon=canyon, n=n, policy=pol, rep=rep)
    if r.returncode != 0:
        return dict(base, error=f"rc={r.returncode} " + r.stderr[-300:])
    summ = dict(kv.split("=", 1) for kv in r.stdout.split() if "=" in kv)
    rows = list(csv.DictReader(open(out)))
    os.unlink(out)
    src = np.array([int(x["src"]) for x in rows], int)
    gen = np.array([float(x["gen_s"]) for x in rows])
    rx = np.array([float(x["rx_s"]) for x in rows])
    via = np.array([int(x["via"]) for x in rows], int)
    A = [aoi(gen[src == k], rx[src == k], WARM, SIM) for k in range(n)]
    ages = np.array([a for a, _ in A])
    above1 = np.array([s[0] for _, s in A]); above5 = np.array([s[1] for _, s in A])
    served = np.array([np.any(src == k) for k in range(n)])
    gen_n = int(summ["generated"])
    return dict(base, generated=gen_n, delivered=len(rows), delivered_share=len(rows) / gen_n,
                relay_share=float(np.mean(via == 1)) if len(via) else 0.0,
                delay_ms=float(1e3 * np.mean(rx - gen)) if len(rx) else float("nan"),
                aoi_mean_s=float(np.mean(ages)), aoi_p90_s=float(np.percentile(ages, 90)),
                stale1_share=float(np.mean(above1)), stale5_share=float(np.mean(above5)),
                unserved=int(np.sum(~served)), relay_asked=int(summ["relay_asked"]),
                relayed=int(summ["relayed"]),
                **{k: float(v) for k, v in summ.items() if k.startswith(("prr_", "prrnohd_"))}, error="")


def main(par):
    done = set()
    if os.path.exists(OUT):
        for l in open(OUT):
            r = json.loads(l)
            if not r.get("error"):
                done.add((r["canyon"], r["n"], r["policy"], r["rep"]))
    todo = [j for j in itertools.product(CANYONS, DENS, POLICIES, range(REPS)) if j not in done]
    todo.sort(key=lambda j: j[3])
    print(len(todo), "runs to do", flush=True)
    lock = threading.Lock()

    def work(j):
        r = run(j)
        with lock:
            open(OUT, "a").write(json.dumps(r) + "\n")
        print(*j, "ERR" if r.get("error") else "ok", flush=True)

    with ThreadPoolExecutor(par) as ex:
        list(ex.map(work, todo))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
