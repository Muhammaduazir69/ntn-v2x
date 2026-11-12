#!/usr/bin/env python3
"""Model predictions for every canyon, density and replication of the ns-3 campaign; writes out/model.json.
Satellite geometry: the same element sets and start times as the ns-3 runs, sampled every 2 s over 300 s.
Two facade-height inputs per replication: the lognormal law the street is drawn from ("law"), and the
empirical distribution of the 2 x 40 heights ns-3 actually drew for that replication ("realized",
out/heights_ns3.json), which separates the model's error from the sampling of one finite street."""
import json, os, sys
import numpy as np
import p8_model as P

HERE = os.path.dirname(os.path.abspath(__file__))
REPS, DENS = 3, (30, 60, 120)
S = P.sats()
HZ = json.load(open(os.path.join(HERE, "out", "heights_ns3.json")))
LAW = P.K.Street.F


def realized(z):
    z = np.sort(np.asarray(z))
    def F(self, x):
        return np.searchsorted(self.H50 * np.exp(z), np.asarray(x, float), side="left") / len(z)
    return F


if __name__ == "__main__":
    res = {}
    REP_LIST = [int(sys.argv[1])] if len(sys.argv) > 1 else list(range(REPS))
    OUTF = os.path.join(HERE, "out", f"model_rep{sys.argv[1]}.json" if len(sys.argv) > 1 else "model.json")
    for rep in REP_LIST:
        V = P.views(S, 1792800000.0 + 977.0 * rep, 300.0, 2.0)
        for inp, F in (("law", LAW), ("realized", realized(HZ[str(rep)]["z"]))):
            P.K.Street.F = F
            for c in P.CANYONS:
                d, _ = P.predict(c, V)
                m = P.predict2(c, V)
                A = P.aoi_timeline(c, V)
                QB = P.relay_qbars(c, V)
                for n in DENS:
                    av = P.relay_avail_from(QB, n)
                    res[f"{inp}_{c}_{n}_{rep}"] = dict(direct=d, multisat=m, relay=m + (1 - m) * av, relay_avail=av,
                                                      aoi_direct=A["direct"], aoi_multisat=A["multisat"])
                    print(inp, c, n, rep, {k: round(v, 3) for k, v in res[f"{inp}_{c}_{n}_{rep}"].items()}, flush=True)
                json.dump(res, open(OUTF, "w"), indent=1)
    P.K.Street.F = LAW
