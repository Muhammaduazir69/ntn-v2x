#!/usr/bin/env python3
"""Data figures of the T-ITS manuscript on satellite visibility and sidelink relaying in street canyons."""
import json, os, sys
from collections import defaultdict
import numpy as np
from figstyle import *  # noqa: F401,F403
import matplotlib.pyplot as plt

setup()
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out"); FIG = os.path.join(HERE, "..", "figures")
CAN = [("suburban", "Suburban"), ("urban", "Urban"), ("dense", "Dense urban")]
POL = [("direct", "Serving satellite", GREY), ("multisat", "Any satellite", BLUE), ("blind", "Blind, nearest", ORANGE),
       ("blindlb", "Blind, $\\geq L_b$", GOLD), ("blindfar", "Blind, apart", PURPLE), ("blindbest", "Blind, best lane", TEAL),
       ("randpersist", "Random, kept", "#a0522d"), ("random", "Random per message", RED), ("relay", "Announced relay", GREEN)]


def load():
    rows = [json.loads(l) for l in open(os.path.join(OUT, "ns3_canyon4.jsonl"))]
    R = defaultdict(list)
    for r in rows:
        if not r.get("error") and "hoMs" not in r and "periodMs" not in r:
            R[(r["canyon"], r["n"], r["policy"])].append(r)
    return R


def save(fig, name):
    fig.savefig(os.path.join(FIG, name)); plt.close(fig); print("wrote", name)


def fig_policies(n=60):
    """Delivered share and mean AoI per policy and canyon: one dot per start (ten over a day), bar at the median."""
    R = load()
    fig, axs = plt.subplots(1, 2, figsize=(PAGE, 2.3))
    w = 0.09
    for k, (p, lab, col) in enumerate(POL):
        for ax, key, sc in ((axs[0], "delivered_share", 100), (axs[1], "aoi_mean_s", 1)):
            for i, (c, _) in enumerate(CAN):
                v = np.array([sc * r[key] for r in R.get((c, n, p), [])])
                if len(v) == 0:
                    continue
                x = i + (k - 4) * w
                ax.bar([x], [np.median(v)], w, color=col, alpha=0.55, label=lab if i == 0 else None)
                ax.plot(np.full(len(v), x), v, ls="none", marker="o", ms=1.8, color=col)
    axs[0].set_ylabel("Messages delivered (%)"); axs[1].set_ylabel("Mean AoI at the server (s)"); axs[1].set_yscale("log")
    for ax in axs:
        ax.set_xticks(range(3)); ax.set_xticklabels([l for _, l in CAN])
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=9, bbox_to_anchor=(0.5, 1.08), fontsize=6)
    fig.tight_layout(); save(fig, "fig2_policies.pdf")


def fig_model(n=60):
    """Model (realized facade heights) against ns-3: delivered share and AoI, every replication."""
    R = load(); Mo = {}
    import glob
    for pth in glob.glob(os.path.join(OUT, "model4_rep*.json")):
        Mo.update(json.load(open(pth)))
    fig, axs = plt.subplots(1, 2, figsize=(PAGE, 2.3))
    mk = {"suburban": "o", "urban": "s", "dense": "^"}
    for p, lab, col in POL[:2]:
        for c, cl in CAN:
            for r in R.get((c, n, p), []):
                k = f"realized_{c}_{n}_{r['rep']}"; kl = f"law_{c}_{n}_{r['rep']}"
                if k in Mo:
                    axs[0].plot(100 * Mo[k][p], 100 * r["delivered_share"], mk[c], color=col, mfc="white", ms=4)
                    axs[1].plot(Mo[k][f"aoi_{p}"], r["aoi_mean_s"], mk[c], color=col, mfc="white", ms=4)
                if kl in Mo:
                    axs[1].plot(Mo[kl][f"aoi_{p}"], r["aoi_mean_s"], mk[c], color=col, alpha=0.35, ms=3)
    for ax, lim in ((axs[0], (15, 100)), (axs[1], (0.1, 60))):
        ax.plot(lim, lim, color=INK2, lw=0.8)
    axs[1].set_xscale("log"); axs[1].set_yscale("log")
    axs[0].set_xlabel("Model delivered share (%)"); axs[0].set_ylabel("ns-3 delivered share (%)")
    axs[1].set_xlabel("Model mean AoI (s)"); axs[1].set_ylabel("ns-3 mean AoI (s)")
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], ls="none", marker=mk[c], color=INK, mfc="white", label=l) for c, l in CAN]
    hs += [Line2D([], [], ls="none", marker="o", color=col, label=lab) for _, lab, col in POL[:2]]
    hs += [Line2D([], [], ls="none", marker="o", color=INK, alpha=0.35, label="lognormal law (AoI)")]
    fig.legend(handles=hs, loc="upper center", ncol=6, bbox_to_anchor=(0.5, 1.08), fontsize=6.5)
    fig.tight_layout(); save(fig, "fig3_model.pdf")


def fig_density():
    """Median AoI of the relaying policies against the number of vehicles, urban canyon, first ten starts."""
    R = load()
    fig, ax = plt.subplots(figsize=(COL, 2.6))
    for (c, cl), ls in zip(CAN[1:2], ("-",)):
        for p, lab, col in POL[2:]:
            if p == "randpersist":   # not run in the density sweep
                continue
            xs, ys = [], []
            for n in (30, 60, 120):
                v = [r["aoi_mean_s"] for r in R.get((c, n, p), []) if r["rep"] < 10]  # density ran on the first ten starts
                if v:
                    xs.append(n); ys.append(np.median(v))
            if xs:
                ax.plot(xs, ys, ls=ls, marker="o", color=col, ms=3.5, label=lab)
    ax.set_xscale("log"); ax.set_xticks([30, 60, 120]); ax.set_xticklabels(["30", "60", "120"])
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xlabel("Vehicles on 1.2 km"); ax.set_ylabel("Median of the mean AoI (s)")
    ax.legend(fontsize=6, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False)
    fig.tight_layout(); save(fig, "fig4_density.pdf")


if __name__ == "__main__":
    for w in (sys.argv[1:] or ["policies", "model", "density"]):
        globals()[f"fig_{w}"]()
