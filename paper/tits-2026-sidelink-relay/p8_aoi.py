"""Age of information of every relaying rule from the facades its rays have swept (Section IV of the paper).

Messages leave a vehicle every T seconds, so at a time t the excess of the age over delta + T/2 is T times the number of
newest messages that all failed, and in expectation E[Delta] = delta + T/2 + T sum_{n>=1} P_n, with P_n the probability that
the n newest messages all failed. A message of a vehicle at position x fails toward every satellite iff each building its
rays land on is at least as high as the ray. Facades are independent per building (A1), so given the positions,

    P_n^src = prod over buildings b swept by the source's rays of the n newest messages of [1 - F(max z on b)]        (*)

which holds for any number of rays, any phase and any trajectory (stops included). Rules differ only in a relay factor:
    fixed relay r (apart or not):  (*) over the union of the rays of source and relay, plus the sidelink term
    k fresh relays per message:    P_n^src [ a_n + (1 - a_n) (1 - pi fbar)^(k n) ]
    announced relay:               P_n^src [ a_n + (1 - a_n) eps^n ],   eps = 1 - pi (1 - sigma), sigma the staleness of the flag
    any decoder:                   the fresh form with k the number of neighbors in range
where a_n = P_n^all / P_n^src is the chance that no neighbor in range is clear at any of the n messages (the same union
product over all neighbors), and fbar is the mean over the neighbors of Pr(neighbor clear | source blocked), computed per
building with the conditional law of the facades the source shares with the neighbor."""
import json, math, os
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
LANE_FRAC = np.array([4.0, 7.5, 12.5, 16.0]) / 20.0


class Street:
    def __init__(self, W, H50, s, Lb=30.0, L=1200.0, h=1.5, gap_every=0.0, gap_len=0.0, emp=None, gap_p=0.0):
        self.W, self.H50, self.s, self.Lb, self.L, self.h = W, H50, s, Lb, L, h
        self.gap_every, self.gap_len = gap_every, gap_len
        self.emp = None if emp is None else np.sort(np.asarray(emp, float))
        self.gap_p = gap_p              # share of buildings of height zero (cross streets of a real street, law form)
        self.lanes = LANE_FRAC * W

    def F(self, z):
        """Probability that a facade is lower than z."""
        if self.emp is not None:
            return np.searchsorted(self.emp, z, side="left") / len(self.emp)
        if z <= 0:
            return self.gap_p
        return self.gap_p + (1 - self.gap_p) * 0.5 * math.erfc(-math.log(z / self.H50) / (self.s * math.sqrt(2.0)))


def rays(st, x, y, E):
    """(wall, building index, ray height) per satellite of E for a vehicle at (x, y); a hit on a cross street or a ray along
    the street yields None, meaning that the vehicle is clear whatever the facades."""
    out = []
    for e, a in E:
        sa = math.sin(a)
        if abs(sa) < 1e-6:
            return None
        right = sa > 0
        s = (st.W - y if right else y) / abs(sa)
        xc = (x + s * math.cos(a)) % st.L
        if st.gap_every > 0 and (xc % st.gap_every) < st.gap_len:
            return None
        out.append((right, int(xc // st.Lb), st.h + s * math.tan(e)))
    return out


def absorb(top, R):
    """Update the highest ray per building with the rays R. Returns the change of the log of the survival product."""
    for w, b, z in R:
        k = (w, b)
        if z > top.get(k, -1.0):
            top[k] = z


def product(st, top):
    p = 1.0
    for z in top.values():
        p *= 1.0 - st.F(z)
        if p <= 0.0:
            return 0.0
    return p


def walk(st, tgrid, E_at, pos_src, y_src, pos_nb, y_nb, relay=None, grid_dt=None, nmax=240, repick=False):
    """For one time t and one placement: lists P_src[n], P_all[n], P_fix[n] for n = 1..; the positions are functions of time.
    pos_nb: functions of the neighbors within range at t; relay: index of the fixed relay among them."""
    t = tgrid
    dt = grid_dt
    src, alln, fix = {}, {}, {}
    P_s, P_a, P_f, t_i = [], [], [], t
    dead_a = dead_f = False
    for n in range(1, nmax + 1):
        if t_i < 0:
            break
        E = E_at(t_i)
        if len(E) == 0:                # no satellite above the mask: the message fails whatever the facades
            P_s.append(product(st, src) if src or n > 1 else 1.0)
            P_a.append(0.0 if dead_a else product(st, alln) if alln or n > 1 else 1.0)
            P_f.append(0.0 if dead_f else product(st, fix) if fix or n > 1 else 1.0)
            t_i -= dt
            continue
        R = rays(st, pos_src(t_i), y_src, E)
        if R is None:
            break
        absorb(src, R); absorb(alln, R); absorb(fix, R)
        now_rel = relay
        if repick:     # the rule names, at each message, the nearest neighbor at least a building length away
            xs = pos_src(t_i)
            ds = [(abs(pf(t_i) - xs), j) for j, pf in enumerate(pos_nb)]
            far = [q for q in ds if q[0] >= st.Lb]
            now_rel = min(far)[1] if far else (min(ds)[1] if ds else None)
        for j, (pf, y) in enumerate(zip(pos_nb, y_nb)):
            Rn = rays(st, pf(t_i), y, E)
            if Rn is None:
                dead_a = True
                if j == now_rel:
                    dead_f = True
            else:
                if not dead_a:
                    absorb(alln, Rn)
                if j == now_rel and not dead_f:
                    absorb(fix, Rn)
        ps = product(st, src)
        P_s.append(ps); P_a.append(0.0 if dead_a else product(st, alln)); P_f.append(0.0 if dead_f else product(st, fix))
        if ps < 1e-4:
            break
        t_i -= dt
    return np.array(P_s), np.array(P_a), np.array(P_f)


def fbar_at(st, t, E_at, pos_src, y_src, pos_nb, y_nb):
    """Mean over the neighbors of Pr(neighbor sees a satellite | source sees none) at time t, with shared facades treated
    through the conditional law H >= z_src on the buildings the source's rays land on."""
    E = E_at(t)
    if len(E) == 0:
        return 0.0, 0
    Rs = rays(st, pos_src(t), y_src, E)
    if Rs is None:
        return 1.0, len(pos_nb)
    top = {}
    absorb(top, Rs)
    vals = []
    for pf, y in zip(pos_nb, y_nb):
        Rn = rays(st, pf(t), y, E)
        if Rn is None:
            vals.append(1.0)
            continue
        tn = {}
        absorb(tn, Rn)
        q = 1.0
        for k, z in tn.items():
            if k in top:
                den = max(1e-9, 1.0 - st.F(top[k]))
                q *= (1.0 - st.F(max(z, top[k]))) / den
            else:
                q *= 1.0 - st.F(z)
        vals.append(1.0 - q)
    return (float(np.mean(vals)) if vals else 0.0), len(vals)


def staleness(st, t, T, E_at, ps, ys, pos_nb, y_nb, d_nb, n_near=4, uniform=False):
    """Chance that the nearest neighbor that announced an uplink has lost it by the time the blocked source names it, exact in the
    facades given the positions. The flag is A old, A uniform on (0, T). For neighbor j, with S the event that the source is blocked
    now, C that j was clear at t - A, and N that j is blocked now:
        Pr(N | S, C) = [P(S, N) - P(S, N, not C)] / [P(S) - P(S, not C)],
    each joint blocked event being the product over buildings of the highest ray that must be blocked. The neighbor named is the
    nearest with a clear flag, so the neighbors are weighted by that event."""
    E = E_at(t)
    Rs = rays(st, ps(t), ys, E)
    if Rs is None or not pos_nb:
        return 0.0
    def prod(*Rl):
        top = {}
        for R in Rl:
            absorb(top, R)
        return product(st, top)
    near = np.argsort(d_nb)[:n_near] if not uniform else np.argsort(d_nb)[::3]
    P_s = prod(Rs)
    if P_s <= 0:
        return 0.0
    out, weight, left = 0.0, 0.0, 1.0
    for j in near:
        v, w = [], []
        for A in (0.25 * T, 0.75 * T):
            Rn, Ro = rays(st, pos_nb[j](t), y_nb[j], E), rays(st, pos_nb[j](t - A), y_nb[j], E_at(max(t - A, 0.0)))
            if Ro is None:                                     # clear for certain then
                den = P_s
                num = 0.0 if Rn is None else prod(Rs, Rn)
            else:
                den = P_s - prod(Rs, Ro)
                num = 0.0 if Rn is None else prod(Rs, Rn) - prod(Rs, Rn, Ro)
            v.append(num); w.append(den)
        c = float(np.mean(w)) / P_s                             # Pr(j clear at t - A | source blocked)
        lost = float(np.mean(v)) / max(float(np.mean(w)), 1e-12)
        if uniform:                                            # a uniform draw among those with a clear flag
            out += c * lost; weight += c
        else:
            out += left * c * lost; weight += left * c
            left *= 1.0 - c
    return out / weight if weight > 0 else 0.0


def evaluate(st, V, dtV, T, delta, traj, rng, pi_sl=1.0, n_t=32, n_place=4, warm=10.0, sim=300.0, R=300.0, direct=False, serve_idx=None):
    """Mean AoI per rule. V: list of arrays (e, a) every dtV seconds. traj: function(rng) returning (x(t) function, lane index);
    returns dict with the rule ages. The grid of messages in the sum is max(T, 0.5 s), and its Riemann weight goes with it."""
    def E_at(t):
        return V[min(int(t / dtV), len(V) - 1)]
    def E_serv(t):
        E = E_at(t)
        return E[np.argsort(E[:, 0])[-1:]] if len(E) else E
    dt = max(T, 0.5)
    ts = np.linspace(warm + 1, sim - 1, n_t)
    acc = {k: [] for k in ("direct", "multisat", "fixed_lb", "fixed_rand", "fresh", "fresh2", "fresh3", "announced", "announced2", "announced3", "any", "fresh5", "fresh8", "thin1", "thin2", "thin3", "thin5", "thin8", "fbar", "nnb", "reach", "blocked", "sigma", "sigmau", "announced_u", "announced_u2")}
    nvec = None
    for t in ts:
        for _ in range(n_place):
            veh = [traj(rng) for _ in range(traj.n)]
            ps, ys = veh[0]
            pos_nb, y_nb = [], []
            d_nb = []
            for pf, y in veh[1:]:
                d = abs(pf(t) - ps(t))
                if d <= R and not math.isclose(d, 0.0):
                    pos_nb.append(pf); y_nb.append(y); d_nb.append(d)
            # the direct rule: serving satellite only
            Pd, _, _ = walk(st, t, E_serv, ps, ys, [], [], None, dt)
            acc["direct"].append(dt * Pd.sum())
            # rays of all rules at the same geometry
            rd = rng.integers(0, len(pos_nb)) if pos_nb else None
            P_s, P_a, P_rand = walk(st, t, E_at, ps, ys, pos_nb, y_nb, rd, dt)
            acc["multisat"].append(dt * P_s.sum())
            acc["blocked"].append(float(P_s[0]) if len(P_s) else 0.0)
            if not pos_nb or len(P_s) == 0:
                for k in acc:
                    if k not in ("direct", "multisat", "fbar", "nnb", "reach", "blocked", "sigma", "sigmau"):
                        acc[k].append(dt * P_s.sum())
                continue
            # the nearest neighbor at least one building away
            _, _, P_lb = walk(st, t, E_at, ps, ys, pos_nb, y_nb, None, dt, repick=True)
            nv = np.arange(1, len(P_s) + 1) * (dt / T)
            side = lambda Pf: dt * (Pf + (P_s[:len(Pf)] - Pf) * (1 - pi_sl) ** nv[:len(Pf)]).sum()
            acc["fixed_lb"].append(side(P_lb)); acc["fixed_rand"].append(side(P_rand))
            f, m = fbar_at(st, t, E_at, ps, ys, pos_nb, y_nb)
            a = np.minimum(1.0, P_a / np.maximum(P_s, 1e-12))
            cond = f / np.maximum(1e-9, 1.0 - a[0]) if a[0] < 1 else 0.0
            cond = min(1.0, cond)
            w = float(P_s[0])
            acc["fbar"].append(f * w); acc["nnb"].append(m); acc["reach"].append((1.0 - a[0]) * w)
            for key, k in (("fresh", 1), ("fresh2", 2), ("fresh3", 3), ("fresh5", 5), ("fresh8", 8), ("any", m)):
                phi = (1.0 - pi_sl * cond) ** k
                acc[key].append(dt * (P_s * (a + (1.0 - a) * phi ** nv)).sum())
            # thinned: each of the u clear neighbors a forwarder knows forwards with probability min(1, kappa / (1 + u))
            un = max(1.0, m * cond)
            for kap in (1, 2, 3, 5, 8):
                phi = max(0.0, 1.0 - pi_sl * min(1.0, kap / (1.0 + un))) ** un
                acc[f"thin{kap}"].append(dt * (P_s * (a + (1.0 - a) * phi ** nv)).sum())
            # announced: fails while no neighbor is clear, or when the flag is stale or the sidelink loses the message; the stale-flag
            # chance sigma is eq. (sigma), evaluated exactly in the facades by staleness() below
            sig = staleness(st, t, T, E_at, ps, ys, pos_nb, y_nb, d_nb)
            acc["sigma"].append(sig * w)
            sigu = staleness(st, t, T, E_at, ps, ys, pos_nb, y_nb, d_nb, uniform=True)
            acc["sigmau"].append(sigu * w)
            epu = 1.0 - pi_sl * (1.0 - sigu)
            for key, k in (("announced_u", 1), ("announced_u2", 2)):
                acc[key].append(dt * (P_s * (a + (1.0 - a) * (epu ** k) ** nv)).sum())
            eps = 1.0 - pi_sl * (1.0 - sig)
            for key, k in (("announced", 1), ("announced2", 2), ("announced3", 3)):
                acc[key].append(dt * (P_s * (a + (1.0 - a) * (eps ** k) ** nv)).sum())
    b = float(np.mean(acc["blocked"])) if acc["blocked"] else 0.0
    res = {k: delta + T / 2 + float(np.mean(v)) if v else float("nan") for k, v in acc.items() if k not in ("fbar", "nnb", "reach", "blocked", "sigma", "sigmau")}
    # blocked-weighted mean clear share of the neighbors and reach; denominators over all samples
    res.update(sigmau=(float(np.sum(acc["sigmau"])) / (b * len(acc["blocked"]))) if b > 0 else 0.0, sigma=(float(np.sum(acc["sigma"])) / (b * len(acc["blocked"]))) if b > 0 else 0.0, blocked=b, fbar=(float(np.sum(acc["fbar"])) / (b * len(acc["blocked"]))) if b > 0 else 0.0,
               reach=(float(np.sum(acc["reach"])) / (b * len(acc["blocked"]))) if b > 0 else 1.0, nnb=float(np.mean(acc["nnb"])) if acc["nnb"] else 0.0)
    return res


class ConstTraj:
    """Vehicles at constant speed, as in (A3)."""
    def __init__(self, st, n, vmin=6.0, vmax=14.0):
        self.st, self.n, self.vmin, self.vmax = st, n, vmin, vmax
        self.vmean = (vmin + vmax) / 2
    def __call__(self, rng):
        x0, v = rng.uniform(0, self.st.L), rng.uniform(self.vmin, self.vmax)
        lane = int(rng.integers(0, 4))
        v = v if lane < 2 else -v
        return (lambda t, x0=x0, v=v: (x0 + v * t) % self.st.L), self.st.lanes[lane]


class TraceTraj:
    """Vehicles from the SUMO traces of the start: x0, dir, s(t) per second."""
    def __init__(self, st, n, files):
        self.st, self.n = st, n
        self.P = []
        for f in files:
            L = open(f).read().splitlines()[1:]
            for i, l in enumerate(L):
                q = l.split(",")
                self.P.append((float(q[0]), int(q[1]), np.array(q[2:], float), i % 4))
        self.vmean = float(np.mean([p[2][-1] / (len(p[2]) - 1) for p in self.P]))
    def __call__(self, rng):
        x0, d, s, lane = self.P[int(rng.integers(0, len(self.P)))]
        def pos(t, x0=x0, d=d, s=s):
            k = min(int(t), len(s) - 2)
            return (x0 + d * (s[k] + (t - k) * (s[k + 1] - s[k]))) % self.st.L
        return pos, self.st.lanes[lane]
