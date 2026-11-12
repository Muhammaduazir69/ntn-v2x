#!/usr/bin/env python3
"""Merge out/model_rep{0,1,2}.json into out/model.json."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__)); M = {}
for r in range(3):
    p = os.path.join(HERE, "out", f"model_rep{r}.json")
    if os.path.exists(p):
        M.update(json.load(open(p)))
json.dump(M, open(os.path.join(HERE, "out", "model.json"), "w"), indent=1); print(len(M), "entries")
