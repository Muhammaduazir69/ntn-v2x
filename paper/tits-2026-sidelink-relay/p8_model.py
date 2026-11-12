"""Model predictions for the canyon-relay study, from SGP4 geometry and the canyon line-of-sight law.

For each second of a run, the visible satellites (elevation >= mask) of the direct-to-cell constellation
give (e_k, a_k) at the street. Per lane y:
  direct    P_LOS(e_best, a_best | y)                     (serving cell = highest elevation)
  multisat  1 - prod_k (1 - P_LOS(e_k, a_k | y))           (facades crossed by different rays taken
                                                           as independent)
Time and lane averages give the predicted delivered share of the direct and multisat policies."""
import json, os, sys
import numpy as np
from sgp4.api import Satrec, jday
import canyon as K

HERE = os.path.dirname(os.path.abspath(__file__))
TLE = "/home/uzair/6g_ntn_ns3/ns-3-dev/contrib/ntn-v2x/data/starlink_dtc_20261002.tle"
CANYONS = {"suburban": dict(W=30, H50=10), "urban": dict(W=20, H50=18), "dense": dict(W=15, H50=30)}
LAT, LON, AZ, MASK = 40.75, -73.99, 29.0, 20.0


def sats():
    L = [l.rstrip() for l in open(TLE) if l.strip()]
    return [Satrec.twoline2rv(L[k + 1], L[k + 2]) for k in range(0, len(L) - 2, 3)]


def gmst(jd, fr):
    d = (jd - 2451545.0) + fr
    return np.mod(np.radians(280.46061837 + 360.98564736629 * d), 2 * np.pi)


def views(S, epoch_unix, secs, step=1.0, with_ids=False):
    """List over time of arrays (e, a) of the satellites above the mask, a relative to the street axis."""
    la, lo = np.radians(LAT), np.radians(LON)
    a_e, e2 = 6378137.0, 6.69437999014e-3
    n = a_e / np.sqrt(1 - e2 * np.sin(la) ** 2)
    g = np.array([n * np.cos(la) * np.cos(lo), n * np.cos(la) * np.sin(lo), n * (1 - e2) * np.sin(la)])
    up = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    east = np.array([-np.sin(lo), np.cos(lo), 0])
    north = np.array([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)])
    jd0 = 2440587.5 + np.floor(epoch_unix / 86400.0)
    fr0 = (epoch_unix / 86400.0) % 1.0
    out = []
    for t in np.arange(0, secs, step):
        fr = fr0 + t / 86400.0
        th = gmst(jd0, fr)
        E = []
        for sid, s in enumerate(S):
            err, r, _ = s.sgp4(jd0, fr)
            if err:
                continue
            r = np.array(r) * 1e3
            x = np.array([np.cos(th) * r[0] + np.sin(th) * r[1], -np.sin(th) * r[0] + np.cos(th) * r[1], r[2]])
            d = x - g
            rr = np.linalg.norm(d)
            e = np.arcsin(d @ up / rr)
            if e >= np.radians(MASK):
                E.append((e, np.arctan2(d @ east, d @ north) - np.radians(AZ), sid))
        A = np.array(E).reshape(-1, 3)
        out.append(A if with_ids else A[:, :2])
    return out


def predict(canyon, V):
    st = K.Street(W=CANYONS[canyon]["W"], H50=CANYONS[canyon]["H50"], sH=0.5, Lb=30.0,
                  lanes=np.array([4.0, 7.5, 12.5, 16.0]) * CANYONS[canyon]["W"] / 20.0)
    pd, pm = [], []
    for E in V:
        if len(E) == 0:
            pd.append(0.0); pm.append(0.0); continue
        kb = np.argmax(E[:, 0])
        per_lane_d, per_lane_m = [], []
        for y in st.lanes:
            p = np.array([float(st.p_los_y(e, a, y)) for e, a in E])
            per_lane_d.append(p[kb])
            per_lane_m.append(1 - np.prod(1 - p))
        pd.append(np.mean(per_lane_d)); pm.append(np.mean(per_lane_m))
    return float(np.mean(pd)), float(np.mean(pm))


if __name__ == "__main__":
    S = sats()
    res = {}
    for rep in range(3):
        V = views(S, 1792800000.0 + 977.0 * rep, 300.0, 2.0)
        for c in CANYONS:
            d, m = predict(c, V)
            res[f"{c}_{rep}"] = dict(direct=d, multisat=m, nvis=float(np.mean([len(E) for E in V])))
            print(c, rep, "direct %.3f multisat %.3f visible %.1f" % (d, m, res[f"{c}_{rep}"]["nvis"]), flush=True)
    json.dump(res, open(os.path.join(HERE, "out", "model_direct.json"), "w"), indent=1)


def aoi_model(canyon, V, T=0.1, delay=0.005, Lb=30.0, vmin=6.0, vmax=14.0, cap=None, multisat_exact=True):
    """Mean AoI of the direct and multisat policies from the facade renewal law.

    Along the street the facade a ray crosses changes every building, Lb / v seconds; each clears it with
    probability p, independently, so a blocked run spans K buildings with K geometric. A message generated
    at a random time is in a blocked run with probability 1 - p and then waits the run's forward recurrence
    time (Lb/v) E[K^2] / (2 E[K]) = (Lb/v)(2 - p)/(2p). The satellite geometry changes on a slower scale, so
    the wait is capped at `cap` seconds: by default the forward recurrence time of the serving satellite's
    tenure (out/tenure.json), the mean time from a random instant until the serving geometry changes.
    Mean AoI = delay + T/2 + E[(1 - p) min((Lb/v)(2 - p)/(2p), cap)]."""
    if cap is None:
        cap = json.load(open(os.path.join(HERE, "out", "tenure.json")))["residual_s"]
    st = K.Street(W=CANYONS[canyon]["W"], H50=CANYONS[canyon]["H50"], sH=0.5, Lb=Lb,
                  lanes=np.array([4.0, 7.5, 12.5, 16.0]) * CANYONS[canyon]["W"] / 20.0)
    inv_v = np.log(vmax / vmin) / (vmax - vmin)          # E[1/v] for v uniform on [vmin, vmax]
    out = {"direct": [], "multisat": []}
    for E in V:
        for y in st.lanes:
            if len(E) == 0:
                for k in out:
                    out[k].append(cap)
                continue
            p = np.array([float(st.p_los_y(e, a, y)) for e, a in E])
            pd = p[np.argmax(E[:, 0])]
            pm = p_any(st, E, y) if multisat_exact else 1 - np.prod(1 - p)
            for k, q in (("direct", pd), ("multisat", pm)):
                w = cap if q <= 1e-6 else min(Lb * inv_v * (2 - q) / (2 * q), cap)
                out[k].append((1 - q) * w)
    return {k: delay + T / 2 + float(np.mean(v)) for k, v in out.items()}


def p_any(st, E, y, n_off=24):
    """P(at least one satellite in line of sight) for a vehicle at lane y, with the facades shared.

    Each ray k reaches its wall at along-street offset c_k = s_k cos a_k and needs the facade there to be
    below z_k. Rays that land on the same building of the same wall share its height H, and all of them are
    blocked iff H >= max z_k over those rays. With the building grid at a uniform offset u relative to the
    vehicle, P(all blocked | u) = prod over buildings b of (1 - F(max_{k in b} z_k)); P_any = 1 - E_u[.]."""
    if len(E) == 0:
        return 0.0
    e, a = E[:, 0], E[:, 1]
    sa = np.abs(np.sin(a))
    right = np.sin(a) >= 0
    d = np.where(right, st.W - y, y)
    s = d / np.maximum(sa, 1e-9)
    z = st.h + s * np.tan(e)
    c = s * np.cos(a)
    clear_always = sa < 1e-9
    if clear_always.any():
        return 1.0
    tot = 0.0
    for u in (np.arange(n_off) + 0.5) / n_off * st.Lb:
        blocked = 1.0
        idx = np.floor((c + u) / st.Lb).astype(int)
        for side in (True, False):
            for b in np.unique(idx[right == side]):
                m = (right == side) & (idx == b)
                blocked *= 1 - float(st.F(z[m].max()))
        tot += blocked
    return 1 - tot / n_off


def predict2(canyon, V):
    st = K.Street(W=CANYONS[canyon]["W"], H50=CANYONS[canyon]["H50"], sH=0.5, Lb=30.0,
                  lanes=np.array([4.0, 7.5, 12.5, 16.0]) * CANYONS[canyon]["W"] / 20.0)
    return float(np.mean([np.mean([p_any(st, E, y) for y in st.lanes]) for E in V]))


def relay_avail(canyon, V, n_veh, length=1200.0, R=300.0, n_dx=24):
    """Probability that a vehicle blocked toward every satellite finds a neighbor within R that sees one.

    For a neighbor at along-street offset dx in lane y2, Proposition 2 gives, per satellite k, the
    probability that its ray is blocked given that the source's ray is blocked:
        q_k = (1 - p_1k - p_2k + P(both clear)) / (1 - p_1k);
    satellites are taken as independent, so the neighbor has no satellite with probability prod_k q_k.
    Neighbors are uniform over the street, n_veh - 1 of them, a share 2R/length within range; the source has
    no relay with probability (E[prod_k q_k])^m averaged over the binomial count m of neighbors in range."""
    from scipy import stats as _st
    st = K.Street(W=CANYONS[canyon]["W"], H50=CANYONS[canyon]["H50"], sH=0.5, Lb=30.0,
                  lanes=np.array([4.0, 7.5, 12.5, 16.0]) * CANYONS[canyon]["W"] / 20.0)
    dxs = (np.arange(n_dx) + 0.5) / n_dx * R
    pin = min(1.0, 2 * R / length)
    m = np.arange(n_veh)
    pm = _st.binom.pmf(m, n_veh - 1, pin)
    out = []
    for E in V:
        if len(E) == 0:
            out.append(0.0); continue
        avail = []
        for y1 in st.lanes:
            p1 = np.array([float(st.p_los_y(e, a, y1)) for e, a in E])
            if np.all(p1 > 0.999):
                continue
            none = []
            for y2 in st.lanes:
                p2 = np.array([float(st.p_los_y(e, a, y2)) for e, a in E])
                for dx in dxs:
                    q = 1.0
                    for (e, a), a1, a2 in zip(E, p1, p2):
                        both = st.p_both(e, a, y1, y2, dx)
                        q *= max(0.0, (1 - a1 - a2 + both)) / max(1e-9, 1 - a1)
                    none.append(q)
            qbar = float(np.mean(none))
            avail.append(1 - float(np.sum(pm * qbar ** m)))
        out.append(float(np.mean(avail)) if avail else 1.0)
    return float(np.mean(out))


def aoi_timeline(canyon, V, dt_v=2.0, T=0.1, delay=0.005, Lb=30.0, vmin=6.0, vmax=14.0, warm=10.0, n_v=9, n_ph=4):
    """Mean AoI from the facade renewal process along the actual satellite timeline, without a tenure cap.

    A vehicle at speed v passes a new facade every Lb / v seconds; at each passage the link toward the
    satellite(s) is clear with probability p(t) of the geometry at that time (serving satellite: direct;
    any satellite with shared facades: multisat), independently across passages. The age at time t exceeds
    tau iff no passage in (t - tau, t] was clear, so
        E[age(t)] = int_0^t prod_{passages s in (t - tau, t]} (1 - p(s)) dtau,
    counted from the start of the run as the measurement does, averaged over [warm, end], lanes, speeds
    uniform on [vmin, vmax] and the phase of the facade grid. The link stays up for the whole building after
    a clear passage, so with the passages up to t numbered newest first (the building at t is number 1, the
    building at the start of the run is the oldest), age(t) > tau for tau between the (r-1)-th and the r-th
    newest passage iff the r newest buildings are all blocked."""
    st = K.Street(W=CANYONS[canyon]["W"], H50=CANYONS[canyon]["H50"], sH=0.5, Lb=Lb,
                  lanes=np.array([4.0, 7.5, 12.5, 16.0]) * CANYONS[canyon]["W"] / 20.0)
    nT = len(V); tend = nT * dt_v
    P = {"direct": np.zeros((len(st.lanes), nT)), "multisat": np.zeros((len(st.lanes), nT))}
    for i, E in enumerate(V):
        if len(E) == 0:
            continue
        kb = np.argmax(E[:, 0])
        for j, y in enumerate(st.lanes):
            P["direct"][j, i] = float(st.p_los_y(E[kb, 0], E[kb, 1], y))
            P["multisat"][j, i] = p_any(st, E, y)
    out = {}
    tg = np.arange(warm, tend, 0.5)
    for k, Pk in P.items():
        acc = []
        for v in np.linspace(vmin, vmax, n_v):
            step = Lb / v
            for ph in (np.arange(n_ph) + 0.5) / n_ph * step:
                ts = np.concatenate([[0.0], np.arange(ph, tend, step)])   # start of the run, then passages
                for j in range(len(st.lanes)):
                    ps = Pk[j, np.minimum((ts / dt_v).astype(int), nT - 1)]
                    cum = np.concatenate([[0.0], np.cumsum(np.log(np.maximum(1 - ps, 1e-300)))])
                    for t in tg:
                        m = np.searchsorted(ts, t, side="right")          # buildings 0..m-1 started by t
                        r = np.arange(1, m + 1)
                        surv = np.exp(cum[m] - cum[m - r])                  # r newest buildings all blocked
                        edges = np.concatenate([[0.0], t - ts[1:m][::-1], [t]])   # 0, t - newest passage, ..., t
                        acc.append(float(np.sum(surv * np.diff(edges))))
        out[k] = delay + T / 2 + float(np.mean(acc))
    return out


def relay_qbars(canyon, V, R=300.0, n_dx=24):
    """Vectorized core of relay_avail: for every time and source lane, the mean over neighbor lanes and offsets
    within R of the probability that the neighbor sees no satellite given that the source sees none (None when
    the source lane is always clear). Independent of the density, so computed once per canyon."""
    st = K.Street(W=CANYONS[canyon]["W"], H50=CANYONS[canyon]["H50"], sH=0.5, Lb=30.0,
                  lanes=np.array([4.0, 7.5, 12.5, 16.0]) * CANYONS[canyon]["W"] / 20.0)
    dxs = (np.arange(n_dx) + 0.5) / n_dx * R
    out = []
    for E in V:
        if len(E) == 0:
            out.append(None); continue
        e, a = E[:, 0], E[:, 1]
        sa = np.maximum(np.abs(np.sin(a)), 1e-9); ca = np.cos(a)
        right = np.sin(a) >= 0
        dp = {y: np.where(right, st.W - y, y) for y in st.lanes}
        z = {y: st.h + dp[y] * np.tan(e) / sa for y in st.lanes}
        xw = {y: dp[y] / sa * ca for y in st.lanes}
        F = {y: np.asarray(st.F(z[y]), float) for y in st.lanes}
        row = []
        for y1 in st.lanes:
            p1 = F[y1]
            if np.all(p1 > 0.999):
                row.append(None); continue
            qs = []
            for y2 in st.lanes:
                gap = np.abs(dxs[None, :] + xw[y2][:, None] - xw[y1][:, None])            # (k, n_dx)
                same = np.maximum(0.0, 1 - gap / st.Lb)
                both = same * np.asarray(st.F(np.minimum(z[y1], z[y2])), float)[:, None] + (1 - same) * (p1 * F[y2])[:, None]
                q = np.maximum(0.0, 1 - p1[:, None] - F[y2][:, None] + both) / np.maximum(1e-9, 1 - p1[:, None])
                qs.append(np.prod(q, axis=0))                                             # (n_dx,)
            row.append(float(np.mean(qs)))
        out.append(row)
    return out


def relay_avail_from(qbars, n_veh, length=1200.0, R=300.0, n_x=48):
    """Relay availability for n_veh vehicles on a street of the given length. The sidelink street is a line,
    not a ring, so a source at position x has neighbors within R with probability
    pin(x) = (min(x + R, length) - max(x - R, 0)) / length; the count in range is binomial and the result is
    averaged over x uniform on the street."""
    from scipy import stats as _st
    xs = (np.arange(n_x) + 0.5) / n_x * length
    pins = (np.minimum(xs + R, length) - np.maximum(xs - R, 0.0)) / length
    m = np.arange(n_veh)
    pms = [_st.binom.pmf(m, n_veh - 1, p_) for p_ in pins]
    vals = []
    for row in qbars:
        if row is None:
            vals.append(0.0); continue
        av = [np.mean([1 - float(np.sum(pm * q ** m)) for pm in pms]) for q in row if q is not None]
        vals.append(float(np.mean(av)) if av else 1.0)
    return float(np.mean(vals))


def relay_none_exact(canyon, V, n_veh, length=1200.0, R=300.0, Lb=30.0, n_pos=64, seed=5):
    """P(the source and every neighbor within R are blocked toward every satellite), exact in the facade heights.

    Facade heights are independent across buildings with law F. Given the positions and lanes of the source and
    its neighbors, every ray (vehicle, satellite) lands on one building of one wall at height z; all vehicles are
    blocked iff every building's height is at least the largest z landing on it, so
        P(all blocked | positions) = prod over buildings B of (1 - F(max z on B)).
    The facades wrap around the street (as in the simulation), the sidelink range does not. Positions are
    sampled: the source uniform on the street and lane, the other n_veh - 1 vehicles uniform, those within R
    kept. Returns the time average of P(all blocked) and of P(source blocked), whose difference is the share of
    messages a relay can carry."""
    rng = np.random.default_rng(seed)
    st = K.Street(W=CANYONS[canyon]["W"], H50=CANYONS[canyon]["H50"], sH=0.5, Lb=Lb,
                  lanes=np.array([4.0, 7.5, 12.5, 16.0]) * CANYONS[canyon]["W"] / 20.0)
    nb = int(np.ceil(length / Lb))
    none_t, src_t = [], []
    for E in V:
        if len(E) == 0:
            none_t.append(1.0); src_t.append(1.0); continue
        e, a = E[:, 0], E[:, 1]
        sa = np.maximum(np.abs(np.sin(a)), 1e-9); right = np.sin(a) >= 0
        acc_n, acc_s = 0.0, 0.0
        for _ in range(n_pos):
            x = rng.uniform(0, length, n_veh); y = st.lanes[rng.integers(0, 4, n_veh)]
            keep = np.abs(x - x[0]) <= R; keep[0] = True
            x, y = x[keep], y[keep]
            dp = np.where(right[None, :], st.W - y[:, None], y[:, None])          # (v, k)
            s = dp / sa[None, :]
            z = st.h + s * np.tan(e)[None, :]
            b = (np.floor(np.mod(x[:, None] + s * np.cos(a)[None, :], length) / Lb).astype(int) % nb) + np.where(right, nb, 0)[None, :]
            zmax = np.full(2 * nb, -np.inf)
            np.maximum.at(zmax, b.ravel(), z.ravel())
            used = np.isfinite(zmax)
            acc_n += float(np.prod(1 - np.asarray(st.F(zmax[used]), float)))
            zs = np.full(2 * nb, -np.inf); np.maximum.at(zs, b[0], z[0]); us = np.isfinite(zs)
            acc_s += float(np.prod(1 - np.asarray(st.F(zs[us]), float)))
        none_t.append(acc_n / n_pos); src_t.append(acc_s / n_pos)
    return float(np.mean(none_t)), float(np.mean(src_t))
