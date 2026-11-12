#!/usr/bin/env python3
"""Plain-text abstract for the submission portal, from sections/abstract.tex and macros_p8.tex."""
import os, re
H = os.path.dirname(os.path.abspath(__file__)); P = os.path.join(H, "..")
M = dict(re.findall(r"\\newcommand\{\\([A-Za-z]+)\}\{(.*)\}", open(os.path.join(P, "macros_p8.tex")).read()))
g = lambda s: re.sub(r"<<(\w+)>>", lambda m: M[m.group(1)], s)
abs_tex = open(os.path.join(P, "sections", "abstract.tex")).read().split(r"\begin{abstract}")[1].split(r"\end{abstract}")[0]
t = re.sub(r"\\([A-Za-z]+)\\? ?", lambda m: M.get(m.group(1), m.group(0)) + " ", abs_tex)
t = t.replace("~s", " s").replace("~", " ").replace("\\,\\%", " %").replace("\n", " ")
t = re.sub(r"\s+", " ", t).strip()
t = t.replace(" ,", ",").replace(" .", ".").replace("\\,\\%", " percent").replace(",\\%", " percent").replace(" %", " percent")
# The cover letter is kept by hand in ../cover_letter.md (its numbers checked against macros_p8.tex); this script
# writes only the plain-text abstract for the submission portal.
open(os.path.join(H, "out", "abstract_plain.txt"), "w").write(t + "\n")
print(len(t.split()), "words"); print(t)
