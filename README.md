# Radiotherapy fractionation optimizer

Given a fixed limit on late damage to healthy tissue (the BED of 35 x 2 Gy, alpha/beta = 3 Gy), find the number of fractions that maximizes tumor BED under the linear-quadratic model, with an optional tumor repopulation penalty and an 8 Gy cap on dose per fraction.

## Results

| Site | Tumor alpha/beta | Sparing | Optimum | Gain over 35 fractions |
|---|---:|---:|---|---:|
| Prostate | 1.5 Gy | 0.8 | 6 x 7.85 Gy | +60.6 Gy |
| Breast | 3.5 Gy | 0.9 | flat (5 to 40 fractions within 1.6%) | +0.1 Gy |
| Head and neck | 10 Gy | 0.9 | 16 x 3.79 Gy (ends day 21) | +11.1 Gy |

The grid search reproduces the analytical result (Mizuta et al. 2012) that fewer fractions win exactly when s x (alpha/beta)_T < (alpha/beta)_N. The test suite checks this.

## Run
```
pip install -r requirements.txt
pytest
python run_optimizer.py results/radiotherapy.json
```

## References
Fowler, *Br J Radiol* 1989; Withers et al., *Acta Oncol* 1988; Mizuta et al., *IJROBP* 2012; FAST-Forward, *Lancet* 2020; PACE-B, *Lancet Oncol* 2019.
