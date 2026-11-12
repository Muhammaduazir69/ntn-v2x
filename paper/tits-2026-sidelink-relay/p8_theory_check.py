#!/usr/bin/env python3
"""Exact-model Monte Carlo of the renewal facades against the closed forms of the analysis (Corollaries on rays and on the
fresh draw). One source, one ray per vehicle; each vehicle passes a building every ell seconds at its own phase, each clear with
probability p, independently. Messages every T seconds. Writes out/theory.json."""
import json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))


def G(tau, p, ell):
    """Blocked throughout a window of length tau."""
    m = np.floor(tau / ell); f = tau / ell - m
    return (1 - p) ** (m + 1) * (1 - f * p)


def Gn(n, T, p, ell):
    """Blocked at all of the n + 1 newest messages: messages at least a passage apart see independent buildings."""
    return (1 - p) ** (n + 1) if T >= ell else G(n * T, p, ell)


def excess_fixed(T, p, ell, K=1, N=200000):
    n = np.arange(N); g = Gn(n, T, p, ell)
    return T * np.sum(g * g ** K)


def excess_fresh(T, p, ell, f=None, pi=1.0, k=1, N=200000):
    f = p if f is None else f
    phi = (1 - pi * f) ** k
    n = np.arange(N); g = Gn(n, T, p, ell)
    return T * np.sum(g * phi ** (n + 1))


def closed_rays(p, ell, K):
    q = 1 - p
    return ell * q ** K / (1 - q ** K) * (1 - q ** (K + 1)) / ((K + 1) * p)


def mc(T, p, ell, rule, steps=400000, m=20, seed=1):
    """Time-average of the excess above T/2: number of newest failed messages times T, sampled at message times."""
    rng = np.random.default_rng(seed)
    ph = rng.random(m + 1) * ell                  # next renewal of source (0) and m neighbors
    st = rng.random(m + 1) < p
    fails = 0; tot = 0.0
    fixed = 1
    for _ in range(steps):
        ph -= T
        due = ph < 0
        while due.any():
            st[due] = rng.random(due.sum()) < p
            ph[due] += ell
            due = ph < 0
        if rule == "none":
            ok = st[0]
        elif rule == "fixed":
            ok = st[0] or st[fixed]
        else:
            ok = st[0] or st[rng.integers(1, m + 1)]
        fails = 0 if ok else fails + 1
        tot += fails
    return T * tot / steps


if __name__ == "__main__":
    out = {"rays": [], "rule": []}
    for K, p in ((1, .3), (1, .7), (2, .3), (2, .7), (3, .5), (3, .7)):
        q = 1 - p
        # K independent rays with the sum over a fine grid, the integral form, and the old single-grid form with p_any
        a = closed_rays(p, 3.0, K)
        n = np.arange(400000); gridval = 0.005 * np.sum(G(n * 0.005, p, 3.0) ** K)
        P = 1 - q ** K
        out["rays"].append(dict(K=K, p=p, closed=a, numeric=float(gridval), single_grid=(1 - P) * 3.0 * (2 - P) / (2 * P)))
    for p in (0.3, 0.5, 0.7):
        for ratio in (0.03, 0.1, 0.2, 1 / 3, 0.5, 0.75, 1.0, 1.5, 2.0):
            T = ratio * 3.0
            fx, fr = excess_fixed(T, p, 3.0, N=4000), excess_fresh(T, p, 3.0, N=4000)
            steps = int(min(400000, 60000 / T)) if T > 0.5 else 200000
            mf = mc(T, p, 3.0, "fixed", steps); mr = mc(T, p, 3.0, "fresh", steps)
            out["rule"].append(dict(p=p, T_over_ell=ratio, fixed=float(fx), fresh=float(fr), rho=float(fr / fx), mc_fixed=mf, mc_fresh=mr, mc_rho=mr / mf))
            print(p, round(ratio, 3), round(fr / fx, 3), round(mr / mf, 3), round(fx, 4), round(mf, 4), round(fr, 4), round(mr, 4), flush=True)
    for r in out["rays"]:
        print(r)
    json.dump(out, open(os.path.join(HERE, "out", "theory.json"), "w"), indent=1)
