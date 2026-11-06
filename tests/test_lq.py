import math
from lq import dose_per_fraction, optimize, overall_days, REF_BED_N, AB_N

def test_reference_schedule_meets_limit_exactly():
    assert math.isclose(dose_per_fraction(35, 1.0), 2.0, rel_tol=1e-9)

def test_overall_time_skips_weekends():
    assert overall_days(5) == 4 and overall_days(6) == 7 and overall_days(16) == 21

def test_grid_matches_analytical_boundary():
    # without repopulation, fewer fractions win iff s * (a/b)_T < (a/b)_N
    for s in (0.5, 0.7, 0.9, 1.0):
        for abT in (1.5, 2.5, 4.0, 6.0, 8.0, 12.0):
            if abs(s * abT - AB_N) < 0.2: continue
            best, _ = optimize(abT, s, 0.0, 0)
            assert (best["n"] < 40) == (s * abT < AB_N)
