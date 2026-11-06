"""
How many fractions? Optimizing radiotherapy schedules with the linear-quadratic model.

For a schedule of n fractions of d Gy:
  tumor BED   = n d (1 + d/(a/b)_T) - D_prolif * max(0, T - T_k)
  normal BED  = n (s d) (1 + s d/(a/b)_N)    s = sparing factor (normal-tissue dose / tumor dose)
Constraint: normal BED <= BED of the reference schedule 35 x 2 Gy at s = 1 (116.7 Gy_3).
For each n, d is the largest dose that meets the constraint (closed form), and tumor control
probability uses a Poisson model averaged over a normal spread of radiosensitivity.

Standard parameter values: late-responding tissue a/b = 3 Gy; repopulation kickoff 21 days
and D_prolif 0.9 Gy/day for fast-growing (a/b ~10) tumors (Fowler 1989; Withers 1988).
"""
import math
import numpy as np

AB_N = 3.0
REF_BED_N = 35*2*(1+2/AB_N)

def overall_days(n):           # 5 fractions per week, Mon-Fri, start Monday
    if n <= 1: return 0
    last = n-1
    return (last//5)*7 + last % 5

def dose_per_fraction(n, s, abN=AB_N, limit=REF_BED_N):
    # n * s d (1 + s d/abN) = limit  ->  (n/abN) x^2 + n x - limit = 0, x = s d
    a = n/abN; b = n; c = -limit
    x = (-b + math.sqrt(b*b - 4*a*c))/(2*a)
    return x/s

def tumor_bed(n, d, abT, prolif, Tk):
    T = overall_days(n)
    return n*d*(1+d/abT) - prolif*max(0, T-Tk)

def tcp(bed, alpha=0.3, sd=0.06, N0=1e9):
    # BED_T is in units of Gy; survival = exp(-alpha * BED); average over alpha ~ N(mu, sd)
    al = np.linspace(max(0.01, alpha-4*sd), alpha+4*sd, 161)
    w = np.exp(-0.5*((al-alpha)/sd)**2); w /= w.sum()
    return float(np.sum(w*np.exp(-N0*np.exp(-al*bed))))

DMAX = 8.0   # LQ model is not trusted far above ~8 Gy per fraction; also the clinical SBRT range

def optimize(abT, s, prolif, Tk, nmax=40, dmax=DMAX):
    best = None; curve = []
    for n in range(1, nmax+1):
        d = dose_per_fraction(n, s)
        ok = d <= dmax
        b = tumor_bed(n, d, abT, prolif, Tk)
        curve.append({"n": n, "d": round(d, 3), "bed": round(b, 2), "ok": ok})
        if ok and (best is None or b > best["bed"] + 1e-9): best = curve[-1]
    ref = curve[34]
    best = dict(best, gain=round(best["bed"] - ref["bed"], 2), ref_bed=ref["bed"])
    return best, curve

