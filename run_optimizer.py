"""Optimize three clinical sites and build the phase map.

Usage: python run_optimizer.py results/radiotherapy.json
"""
import json, sys
import numpy as np
from lq import REF_BED_N, AB_N, optimize

SITES = [  # tumor a/b estimates from the clinical literature
    {"site": "Prostate", "abT": 1.5, "prolif": 0.0, "Tk": 0, "s": 0.8, "note": "Low a/b (Brenner & Hall 1999; CHHiP; PACE-B)"},
    {"site": "Breast", "abT": 3.5, "prolif": 0.0, "Tk": 0, "s": 0.9, "note": "a/b ~3.5 (START trials; FAST-Forward)"},
    {"site": "Head & neck", "abT": 10.0, "prolif": 0.9, "Tk": 21, "s": 0.9, "note": "a/b ~10, accelerated repopulation (Withers 1988)"},
]

def main():
    sites = []
    for st in SITES:
        best, curve = optimize(st["abT"], st["s"], st["prolif"], st["Tk"])
        sites.append({**st, "best": best, "curve": curve})
        print(st["site"], best)
    # phase map: optimal n over (a/b tumor, sparing factor), no repopulation
    abs_ = [round(x, 2) for x in np.linspace(1, 15, 57)]
    ss = [round(x, 3) for x in np.linspace(0.3, 1.0, 36)]
    phase = []
    for s in ss:
        row = []
        for abT in abs_:
            best, _ = optimize(abT, s, 0.0, 0)
            row.append(best["n"])
        phase.append(row)
    json.dump({"ref_bed_n": REF_BED_N, "abN": AB_N, "sites": sites, "phase": {"ab": abs_, "s": ss, "n": phase}},
              open(sys.argv[1], "w"), separators=(",", ":"))
    # sanity: theory says hypofractionate iff abT/abN < s (Mizuta et al. 2012)
    print("ref BED_N", REF_BED_N)
    for s in [0.5, 0.8, 1.0]:
        print(s, [optimize(ab, s, 0, 0)[0]["n"] for ab in [1, 1.5, 2.4, 3, 5, 10]])

if __name__ == "__main__":
    main()
