"""How long the serving (highest-elevation) direct-to-cell satellite stays the serving one at the street;
writes out/tenure.json. This sets the time scale on which a blocked run can end through satellite motion."""
import json, numpy as np, p8_model as P
S = P.sats(); runs = []
for rep in range(3):
    V = P.views(S, 1792800000.0 + 977.0 * rep, 3600.0, 2.0, with_ids=True)
    ids = [int(E[np.argmax(E[:, 0]), 2]) if len(E) else -1 for E in V]
    cur, n = ids[0], 0
    for x in ids[1:]:
        n += 1
        if x != cur:
            runs.append(2.0 * n); cur, n = x, 0
r = np.array(runs[1:])
res = dict(mean_s=float(r.mean()), median_s=float(np.median(r)), residual_s=float((r ** 2).mean() / (2 * r.mean())))
json.dump(res, open("out/tenure.json", "w")); print(res)
