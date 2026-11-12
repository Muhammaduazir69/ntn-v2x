#!/usr/bin/env python3
"""Explicit-geometry Monte Carlo of the canyon, mirroring the ns-3 example (one height draw per run over a
1.2 km wrapped street, lanes, speeds, 0.1 s messages lost while blocked), to separate model error from
simulator behavior. Prints delivered share and mean AoI (from run start, averaged over [10, 300] s)."""
import sys
import numpy as np
import p8_model as P


def run(canyon, V, n=60, seed=0, dt=0.1, dtv=2.0, Lg=1200.0, Lb=30.0, sH=0.5, h=1.5, warm=10.0):
    rng = np.random.default_rng(seed)
    W, H50 = P.CANYONS[canyon]["W"], P.CANYONS[canyon]["H50"]
    nb = int(np.ceil(Lg / Lb))
    HL, HR = H50 * np.exp(sH * rng.standard_normal(nb)), H50 * np.exp(sH * rng.standard_normal(nb))
    lanes = np.array([4.0, 7.5, 12.5, 16.0]) * W / 20.0
    y = lanes[np.arange(n) % 4]; x0 = rng.uniform(0, Lg, n); v = rng.uniform(6, 14, n) * np.where(y < W / 2, 1, -1)
    tend = len(V) * dtv; ts = np.arange(0, tend, dt)
    res = {}
    for pol in ("direct", "multisat"):
        up = np.zeros((len(ts), n), bool)
        for i, t in enumerate(ts):
            E = V[min(int(t / dtv), len(V) - 1)]
            if len(E) == 0:
                continue
            x = np.mod(x0 + v * t, Lg)
            sel = [np.argmax(E[:, 0])] if pol == "direct" else range(len(E))
            for k in sel:
                e, a = E[k]; sa = np.sin(a)
                if abs(sa) < 1e-6:
                    up[i] = True; continue
                right = sa > 0
                s = (W - y if right else y) / abs(sa)
                xc = np.mod(x + s * np.cos(a), Lg)
                Hf = (HR if right else HL)[(xc // Lb).astype(int) % nb]
                up[i] |= Hf < h + s * np.tan(e)
        # AoI: age = t - generation time of newest delivered message, generation 0 if none yet
        last = np.zeros(n); ages = []
        for i, t in enumerate(ts):
            last = np.where(up[i], t, last)
            if t >= warm:
                ages.append(t - last + dt / 2)
        res[pol] = (float(up[ts >= warm].mean()), float(np.mean(ages)))
    return res


if __name__ == "__main__":
    S = P.sats(); V = P.views(S, 1792800000.0, 300.0, 2.0)
    for c in sys.argv[1:] or ("suburban", "urban", "dense"):
        r = [run(c, V, seed=s) for s in range(4)]
        print(c, {k: tuple(np.round(np.mean([x[k] for x in r], 0), 3)) for k in r[0]},
              "per seed AoI direct", [round(x["direct"][1], 1) for x in r])
