#!/usr/bin/env python3
"""Data figures of the revised T-ITS manuscript."""
import glob, json, os, sys
from collections import defaultdict
import numpy as np
from figstyle import *  # noqa: F401,F403
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

setup()
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out"); FIG = os.path.join(HERE, "..", "figures")
BASE = dict(periodMs=1000, ulKbps=300, grantRtt=1, hoMs=100, slUrban=1, alignResel=1, nbAgeS=1.5, gapEveryM=240, gapLenM=18, mob="sumo")
RES = {"sha256", "wall_s", "error", "generated", "delivered", "rx_rows", "delivered_share", "relay_share", "delay_ms", "aoi_mean_s", "aoi_p90_s",
       "stale_t2", "relay_asked", "relayed", "sent_direct", "blocked_gen", "ul_sent", "ul_drop_q", "ul_drop_los", "ul_bytes", "ul_wait_ms",
       "far_fallback", "p_fail", "p_fail_given_fail", "study", "canyon", "n", "policy", "rep"}
rows = [json.loads(l) for l in open(os.path.join(OUT, "ns3_canyon5.jsonl"))]
rows = [r for r in rows if not r.get("error")]
def conf(r): return {k: v for k, v in r.items() if k not in RES and not k.startswith(("prr_", "prrnohd_"))}
def sel(c, p, n=60, **ex):
    w = {**BASE, **ex}
    return {r["rep"]: r for r in rows if r["canyon"] == c and r["n"] == n and r["policy"] == p and conf(r) == w}
POLX = [("relayrand", dict(relayK=1), "Random announcer", "#2e8b57"), ("relayrand", dict(relayK=2), "2 random announcers", "#0b5d3b"),
        ("relayk", dict(relayK=2), "2 nearest announcers", "#4caf50"), ("anyp", dict(kappa=3), "Thinned, $\\kappa=3$", "#5a6bb8"),
        ("anyp", dict(kappa=5), "Thinned, $\\kappa=5$", "#26338a"), ("randk", dict(relayK=8), "8 fresh", "#5e1a10"),
        ("relayfar", {}, "Farthest announcer", "#6b8e23"), ("relaylb", {}, r"$\geq L_{\mathrm{b}}$ announcer", "#2f6f4f")]
POL = [("direct", "Serving satellite", GREY), ("multisat", "Any satellite", BLUE), ("blind", "Nearest", ORANGE), ("blindfar", "Apart", PURPLE),
       ("blindlb", r"$\geq L_{\mathrm{b}}$ away", GOLD), ("randpersist", "Random, kept", "#a0522d"), ("blindbest", "Best lane", TEAL), ("random", "Fresh draw", RED), ("roundrobin", "Round robin", "#b0628f"),
       ("hybrid", "Fresh, best lane", "#5f7f2f"), ("rand2", "2 fresh", "#c9503b"), ("rand3", "3 fresh", "#912f1f"), ("anyp", "Thinned, $\\kappa=1$", "#3b4a9b"), ("anyfwd", "Every decoder", INK), ("relay", "Nearest announcer", GREEN)]
MK = {"direct": "s", "multisat": "D", "blind": "v", "blindfar": "<", "blindlb": ">", "randpersist": "P", "blindbest": "X", "random": "o", "roundrobin": "h", "hybrid": "p", "rand2": "^", "rand3": "*", "anyp": "d", "anyfwd": "8", "relay": "o"}


def save(fig, name):
    fig.savefig(os.path.join(FIG, name)); plt.close(fig); print("wrote", name)


def fig_policies():
    fig, axs = plt.subplots(1, 2, figsize=(PAGE, 2.05), gridspec_kw=dict(width_ratios=[1.75, 1.05]))
    ax = axs[0]
    cs = [("urban", "Urban", -0.17, "o"), ("dense", "Dense", 0.17, "s")]
    for i, (p, lab, col) in enumerate(POL):
        for c, cl, dx, m in cs:
            S = sel(c, p)
            v = [r["aoi_mean_s"] for r in S.values()]
            if not v:
                continue
            ax.plot(np.full(len(v), i + dx), v, ls="none", marker=m, ms=1.6, color=col, alpha=0.55)
            ax.plot([i + dx - 0.15, i + dx + 0.15], [np.median(v)] * 2, color=col, lw=1.6)
    for i, (p, ex, lab, col) in enumerate(POLX):
        x = len(POL) + i
        for c, cl, dx, m in cs:
            v = [r["aoi_mean_s"] for r in sel(c, p, **ex).values()]
            if not v:
                continue
            ax.plot(np.full(len(v), x + dx), v, ls="none", marker=m, ms=1.6, color=col, alpha=0.55)
            ax.plot([x + dx - 0.15, x + dx + 0.15], [np.median(v)] * 2, color=col, lw=1.6)
    fl = 0.5 + np.median([r["delay_ms"] / 1e3 for c in ("urban", "dense") for r in sel(c, "relay").values()])
    ax.axhline(fl, color=INK2, lw=0.6, ls=":")
    ax.set_yscale("log"); ax.set_ylim(0.4, 14); ax.set_ylabel("Mean AoI at the server (s)")
    ax.set_xticks(range(len(POL) + len(POLX))); ax.set_xticklabels([l for _, l, _ in POL] + [l for _, _, l, _ in POLX], rotation=55, ha="right", fontsize=6)
    ax.grid(axis="x", visible=False)
    ax.legend(handles=[Line2D([], [], ls="none", marker="o", ms=3, color=INK2, label="Urban (left of each pair)"), Line2D([], [], ls="none", marker="s", ms=3, color=INK2, label="Dense urban (right)"), Line2D([], [], ls=":", lw=0.6, color=INK2, label=r"Floor $\delta+T/2$")], loc="upper right", fontsize=6)
    ax = axs[1]
    # paired ratio of each forwarding rule's age to the nearest announcer's, median over starts with a 95 % percentile
    # bootstrap interval, against its forwards per blocked message; linear age axis so that the 1 to 5 % differences show
    FR = [("relayrand", dict(relayK=1), "1R"), ("relayrand", dict(relayK=2), "2R"), ("relayk", dict(relayK=2), "2N"),
          ("relayk", dict(relayK=3), "3N"), ("relayfar", {}, "Far"), ("relaylb", {}, "Lb"), ("anyp", {}, r"$\kappa$1"),
          ("anyp", dict(kappa=3), r"$\kappa$3"), ("anyp", dict(kappa=5), r"$\kappa$5"), ("anyp", dict(kappa=8), r"$\kappa$8"),
          ("anyfwd", {}, "Dec"), ("rand3", {}, "3F")]
    rng = np.random.default_rng(13)
    OFF = {("urban", "1R"): (-12, 6), ("urban", "Far"): (4, 8), ("urban", "Lb"): (5, -1), ("dense", "1R"): (-12, -9), ("dense", "Far"): (5, -3),
           ("dense", "Lb"): (-4, -12), ("urban", "2R"): (-9, 4), ("urban", "2N"): (3, 3), ("dense", "2R"): (-10, -6), ("dense", "2N"): (3, -2)}
    for c, mfc, col, dy in (("urban", INK, INK, 4), ("dense", "white", BLUE, -8)):
        ref = sel(c, "relay")
        for p, ex, tag in FR:
            S = sel(c, p, **ex)
            k = sorted(set(S) & set(ref))
            if len(k) < 5:
                continue
            x = np.array([S[i]["aoi_mean_s"] / ref[i]["aoi_mean_s"] for i in k])
            bs = np.median(x[rng.integers(0, len(x), (4000, len(x)))], axis=1)
            fw = np.median([S[i]["relayed"] / max(1, S[i]["blocked_gen"]) for i in k])
            m = np.median(x)
            ax.errorbar(fw, m, yerr=[[m - np.percentile(bs, 2.5)], [np.percentile(bs, 97.5) - m]], fmt="o", ms=2.8, color=col,
                        mfc=mfc, mew=0.7, elinewidth=0.6, capsize=1.2)
            off = OFF.get((c, tag), (3, dy))
            ax.annotate(tag, (fw, m), xytext=off, textcoords="offset points", fontsize=5, color=col)
    ax.axhline(1.0, color=INK2, lw=0.6, ls=":")
    ax.set_xscale("log"); ax.set_ylim(0.86, 1.34); ax.set_xlim(0.6, 25)
    ax.set_xlabel("Forwards per blocked message"); ax.set_ylabel("Age / nearest announcer")
    keyh = [Line2D([], [], ls="none", marker="o", ms=3, color=INK, label="Urban"),
            Line2D([], [], ls="none", marker="o", ms=3, color=BLUE, mfc="white", label="Dense urban")]
    ax.legend(handles=keyh, loc="upper center", fontsize=5.5, ncol=2)
    fig.tight_layout(); save(fig, "fig2_policies.pdf")


RM = {"direct": "direct", "multisat": "multisat", "random": "fresh", "rand2": "fresh2", "rand3": "fresh3", "relay": "announced", "blindlb": "fixed_lb", "randpersist": "fixed_rand"}


def fig_model():
    Mo = {}
    for pth in glob.glob(os.path.join(OUT, "model6_rep*.json")):
        Mo[int(pth.split("rep")[-1].split(".")[0])] = json.load(open(pth))
    fig, axs = plt.subplots(1, 2, figsize=(PAGE, 1.95))
    pc = {p: c for p, _, c in POL}
    lab = {p: l for p, l, _ in POL}
    for ax, c, ttl in ((axs[0], "urban", "Urban canyon"), (axs[1], "dense", "Dense urban canyon")):
        for p, q in RM.items():
            for rep, r in sel(c, p).items():
                k = f"{c}_1.0_sumo_law"
                if rep in Mo and k in Mo[rep]:
                    ax.plot(Mo[rep][k][q] + r["delay_ms"] / 1e3, r["aoi_mean_s"], ls="none", marker=MK[p], ms=2.6, color=pc[p], mfc="none", mew=0.6)
        lim = (0.45, 25)
        ax.plot(lim, lim, color=INK2, lw=0.7)
        ax.set_xlim(lim); ax.set_ylim(lim); ax.set_xscale("log"); ax.set_yscale("log"); ax.set_title(ttl, fontsize=7.5)
        ax.set_xlabel("Model mean AoI, height law only (s)"); ax.set_ylabel("ns-3 mean AoI (s)")
    hs = [Line2D([], [], ls="none", marker=MK[p], color=pc[p], mfc="none", ms=3.2, label=lab[p]) for p in RM]
    fig.legend(handles=hs, loc="upper center", ncol=8, bbox_to_anchor=(0.5, 1.07), fontsize=6)
    fig.tight_layout(); save(fig, "fig3_model.pdf")


def per_ex(T):
    if T == 1000:
        return {}
    if T == 100:
        return dict(periodMs=100, nbAgeS=1.0, ulKbps=3000)
    return dict(periodMs=T, nbAgeS=max(1.0, 1.5 * T / 1000))


def fig_period():
    Th = json.load(open(os.path.join(OUT, "theory.json")))
    fig, axs = plt.subplots(1, 2, figsize=(PAGE, 1.95))
    ax = axs[0]
    for p, ls in ((0.3, "-"), (0.5, "--"), (0.7, ":")):
        Q = [r for r in Th["rule"] if r["p"] == p]
        ax.plot([r["T_over_ell"] for r in Q], [r["rho"] for r in Q], color=INK2, ls=ls, lw=0.9, label=f"Analysis, $p={p}$")
        ax.plot([r["T_over_ell"] for r in Q], [r["mc_rho"] for r in Q], ls="none", marker="x", ms=3, color=INK2)
    trace = np.array([[float(x) for x in l.split(",")[2:]] for f in glob.glob(os.path.join(HERE, "traces", "n60_r1?.csv")) for l in open(f).read().splitlines()[1:]])
    vbar = float(np.mean(trace[:, -1] / (trace.shape[1] - 1)))
    ell = 30.0 / vbar
    for c, col, m, lab in (("urban", RED, "o", "ns-3 urban"), ("dense", BLUE, "s", "ns-3 dense")):
        xs, ys = [], []
        for T in (100, 500, 1000, 2000, 5000):
            A, B = sel(c, "random", **per_ex(T)), sel(c, "randpersist", **per_ex(T))
            k = [i for i in set(A) & set(B) if i >= 10]
            if len(k) >= 5:
                ex = lambda r: r["aoi_mean_s"] - T / 2000 - r["delay_ms"] / 1e3  # age above the floor, as in Proposition 5
                xs.append(T / 1000 / ell); ys.append(np.median([ex(A[i]) / ex(B[i]) for i in k]))
        ax.plot(xs, ys, marker=m, color=col, mfc="white", ms=3.5, lw=0.8, label=lab)
        Mo = [json.load(open(p)) for p in glob.glob(os.path.join(OUT, "model6_rep1?.json"))]  # held-out starts 10 to 19
        xm, ym = [], []
        for T in (0.5, 1.0, 2.0, 5.0):
            k = f"{c}_{T}_sumo_law"
            xm.append(T / ell); ym.append(np.median([(q[k]["fresh"] - T / 2) / (q[k]["fixed_rand"] - T / 2) for q in Mo if k in q]))
        ax.plot(xm, ym, color=col, ls=":", lw=0.9, label=f"Model, {c}")
    ax.set_xscale("log"); ax.set_xlabel(r"Message period over passage time, $T/\ell$"); ax.set_ylabel(r"Excess-age ratio $\Gamma$, fresh/kept")
    ax.axhline(1, color=INK2, lw=0.5); ax.set_ylim(0, 1.75); ax.legend(fontsize=5.5, loc="upper left", ncol=3)
    ax = axs[1]
    # paired ratio to the nearest announcer at the same budget (median over the second-half starts), linear axis
    for p, ex0, lab, col, m in (("relayrand", dict(relayK=2), "2 random announcers", "#0b5d3b", "s"), ("anyp", dict(kappa=5), "Thinned, $\\kappa=5$", "#26338a", "d"),
                                ("anyfwd", {}, "Every decoder", INK, "8"), ("anyfwd", dict(ulPrio=1), "Every decoder, own first", GREY, "x")):
        for c, ls in (("urban", "-"), ("dense", "--")):
            xs, ys = [], []
            for C, x in ((100, 100), (150, 150), (200, 200), (250, 250), (300, 300), (350, 350), (400, 400), (450, 450), (500, 500), (600, 600), (1000, 1000), (0, 3000)):
                exC = {} if C == 300 else dict(ulKbps=C)
                S = {i: r for i, r in sel(c, p, **dict(ex0, **exC)).items() if i >= 10}
                R = {i: r for i, r in sel(c, "relay", **(dict(exC, ulPrio=1) if ex0.get("ulPrio") else exC)).items() if i >= 10}
                k = sorted(set(S) & set(R))
                if len(k) >= 5:
                    xs.append(x); ys.append(np.median([S[i]["aoi_mean_s"] / R[i]["aoi_mean_s"] for i in k]))
            ax.plot(xs, ys, marker=m, ms=3, color=col, ls=ls, lw=0.9, label=lab if c == "urban" else None)
    ax.axhline(1.0, color=INK2, lw=0.6, ls=":")
    ax.set_xscale("log"); ax.set_xticks([100, 200, 300, 600, 1000, 3000]); ax.set_xticklabels(["100", "200", "300", "600", "1000", "no limit"])
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_xlabel("Uplink capacity per satellite (kb/s)"); ax.set_ylabel("Age / nearest announcer")
    ax.set_ylim(0.85, 2.0)
    ax.legend(fontsize=5.5, loc="upper right", ncol=2, title="solid urban, dashed dense", title_fontsize=5.5)
    fig.tight_layout(); save(fig, "fig4_period.pdf")


def fig_los():
    L = json.load(open(os.path.join(OUT, "los_compare.json")))
    fig, axs = plt.subplots(1, 3, figsize=(PAGE, 1.8), sharey=True)
    for ax, s, t in zip(axs, ("first", "eighth", "w45"), ("First Avenue", "Eighth Avenue", "West 45th Street")):
        R = L[s]["rows"]; x = [r["e"] for r in R]
        ax.plot(x, [r["measured"] for r in R], color=INK, marker="o", ms=3, label="Facade data")
        ax.plot(x, [r["held"] for r in R], color=RED, marker="s", mfc="white", ms=3, label="Lemma 1, other half's law")
        ax.plot(x, [r["std"] for r in R], color=GOLD, ls=":", label="Lemma 1, study's urban law")
        ax.plot(x, [r["tr_urban"] for r in R], color=BLUE, ls="--", label="TR 38.811, urban")
        ax.plot(x, [r["tr_dense"] for r in R], color=PURPLE, ls="-.", label="TR 38.811, dense urban")
        ax.plot(x, [r["p1410_urban"] for r in R], color=GREY, ls="-", lw=0.8, label="P.1410 form, urban")
        ax.set_xlabel("Elevation (deg)"); ax.set_title(t, fontsize=7.5)
    axs[0].set_ylabel("LOS probability")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=6, bbox_to_anchor=(0.5, 1.09), fontsize=6)
    fig.tight_layout(); save(fig, "fig5_los.pdf")


if __name__ == "__main__":
    for w in (sys.argv[1:] or ["policies", "model", "period", "los"]):
        globals()[f"fig_{w}"]()
