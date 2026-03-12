# Code Review Report: MainScript.R

**Date:** 2026-02-27
**File:** `2_Scripts/MainScript.R` (1,777 lines)
**Project:** Behavioral Input Price Risk Management (AgEnRes Deliverable 5.4)

---

## 1. Script Overview

The script implements a six-section analytical pipeline:

| Section | Lines | Description |
|---------|-------|-------------|
| 1. Nitrogen-Yield | 18-156 | Just & Pope production function via robust MM regression |
| 2. Price Dynamics | 158-761 | Stationarity tests, QVAR estimation, copula, stability analysis |
| 3. Price Simulation | 762-931 | Monte Carlo simulation of 1,000 price paths (24 periods) |
| 4. Contract Pricing | 933-1042 | Option premium calculation at 5 strike levels |
| 5. Broad Framing | 1044-1503 | EUT and CPT evaluation of full profit (yield + price risk) |
| 6. Narrow Framing | 1504-1777 | EUT and CPT evaluation of contract-only payoffs |

---

## 2. Strengths

- **Clear section organization** with consistent headers and numbered subsections.
- **Reproducibility**: `set.seed(1234)` is set at the start and before the simulation.
- **Well-chosen methods**: Robust MM regression (`lmrob`) for the production function handles outliers in simulation data well. QVAR is appropriate for capturing distributional dynamics.
- **Multiple stationarity tests** (ADF, KPSS, PP) provide a thorough check with cross-validation of results.
- **Dynamic stability analysis** (Section 2.9) with bootstrap eigenvalue checks ensures the QVAR is stable before simulation.
- **Systematic sensitivity analysis** in CPT sections varying alpha, lambda, and gamma individually around standard values.

---

## 3. Issues and Recommendations

### 3.1 Data Handling

| # | Issue | Location | Severity | Recommendation |
|---|-------|----------|----------|----------------|
| 1 | **Hardcoded grid location** `SimUID == 52338` with no justification | Line 31 | Medium | Document why this location was chosen or parameterize it |
| 2 | **Irregular time intervals** (14-day and 21-day gaps) are acknowledged but not handled | Lines 193-199 | Medium | Either interpolate to regular intervals or explicitly account for irregular spacing in the QVAR estimation |
| 3 | **PPI for 2025 is duplicated** from the last available year | Lines 414, 428 | Low | Document this assumption; consider nowcasting or flagging as provisional |
| 4 | **Date filter `>= 2009-01-01`** is hardcoded | Line 187 | Low | Add a comment explaining why pre-2009 data is excluded (monthly frequency) |

### 3.2 Econometric Modeling

| # | Issue | Location | Severity | Recommendation |
|---|-------|----------|----------|----------------|
| 5 | **Yield function omits several control variables** (weather scenarios, rotation) despite creating dummies for them | Lines 44-66 vs 110-121 | High | The yield function only includes `irr` and tillage dummies but not weather or rotation. Either include them or document why they are excluded |
| 6 | **Variation function uses only `sqrt(FTN)`** and ignores `FTN` and other controls used in the yield function | Lines 130-141 | Medium | Consider a consistent specification or document the rationale |
| 7 | **QVAR specification mismatch**: the formula `mf` uses `f_p ~ f_p1 + w_p1` but the `X` matrix for simulation uses `cbind(1, w_p1, f_p1)` — the column order is swapped | Lines 566 vs 654 | High | Verify the column order in `X` matches the coefficient order from `rq()`. Currently `coeff_f` row order is (intercept, f_p1, w_p1), but `X` has columns (1, w_p1, f_p1). **This may produce incorrect inverse distributions.** |
| 8 | **Simulation X matrix column order** also swapped for `Xfi` and `Xwi` | Lines 797-798 | High | `Xfi <- cbind(1, f_change1, w_change1)` but `mf = f_p ~ f_p1 + w_p1`, so coefficients are ordered (intercept, f_p1, w_p1). The matrix should be `cbind(1, f_change1, w_change1)` which is correct for `mf` but **not for `mw`** where `mw = w_p ~ f_p1 + w_p1`. Both equations use the same `Xfi`/`Xwi` construction which happens to be fine since they share the same RHS variables |
| 9 | **VAR serial correlation** is flagged but not resolved | Lines 488-504 | Medium | ARCH effects are detected. Consider a GARCH extension or document why QVAR sufficiently addresses this |
| 10 | **Copula model is linear** (`F_w ~ F_f`) | Lines 683-691 | Low | Consider whether a parametric copula (Gaussian, Clayton, Frank) might better capture tail dependence |
| 11 | **`cummax` rearrangement** to enforce monotonicity of quantile functions | Lines 658, 661 | Low | This is necessary but may mask crossing quantiles, which could indicate model misspecification. Consider reporting how many crossings occur |

### 3.3 Simulation Design

| # | Issue | Location | Severity | Recommendation |
|---|-------|----------|----------|----------------|
| 12 | **Single quantile draw per time step** (`qi <- sample(1:99, 1)`) applies the same quantile to both equations | Line 800 | High | This imposes perfect rank correlation between fertilizer and wheat price shocks. The copula estimated in Section 2.8 is never used in the simulation. Consider drawing correlated quantiles via the copula |
| 13 | **Starting price is the last observed value** | Lines 784-789 | Medium | This is noted in the TODO; consider using a forward-looking reference price |
| 14 | **Simulation horizon is 24 periods (1 year)** but code comments say "3 years" | CLAUDE.md vs Line 771 | Low | The CLAUDE.md says 3 years but `Nt=24` with `obs_per_year=24` is 1 year. Clarify |
| 15 | **`simulated_years = 1`** means all contract/profit analysis is based on a single exercise date | Line 938-940 | Medium | Consider extending to multiple years or documenting this as a single-period analysis |

### 3.4 EUT Implementation

| # | Issue | Location | Severity | Recommendation |
|---|-------|----------|----------|----------------|
| 16 | **Mean-Variance approximation** used for EUT instead of direct expected utility calculation | Lines 1269-1285 | Medium | The M-V approximation (`CE = E[pi] - (delta/2)*Var[pi]`) is only exact for quadratic utility or normal distributions. Since profits are likely skewed, consider computing `E[u(pi)]` directly and inverting |
| 17 | **`delta` parameter interpretation**: code uses absolute risk aversion (CARA) but calls it "power utility" in CLAUDE.md | CLAUDE.md vs Line 1257 | Low | Clarify whether this is CARA (`u = -exp(-delta*x)`) or CRRA (`u = x^(1-gamma)`) |
| 18 | **Different delta values** used for broad vs narrow framing with no explanation | Lines 1257 vs 1554 | Low | Document why the parameter ranges differ |

### 3.5 CPT Implementation

| # | Issue | Location | Severity | Recommendation |
|---|-------|----------|----------|----------------|
| 19 | **Reference point in broad framing** is `mean(profit_fc)` (mean forward contract profit) | Line 1333 | Medium | Document the theoretical justification. Common alternatives: status quo wealth, expected profit of no contract |
| 20 | **Reference point in narrow framing** is `0` | Line 1619 | Low | This is standard for narrow framing (pure contract evaluation) — consistent |
| 21 | **Weighting function `w(p) = exp(-(-log(p))^g)`** (Prelec, 1998) has `w(0) = 0` but `log(0)` is undefined | Line 1351 | Medium | The code avoids `p=0` because cumulative probabilities start at `1/n`, but add a guard clause or document this |
| 22 | **Decision weights computed from cumulative probabilities rather than decumulative** for gains | Lines 1386-1393 | High | In standard CPT, gains are weighted using decumulative probabilities (from best to worst), and losses using cumulative (from worst to best). The current code uses `cumsum` for both. Verify against Tversky & Kahneman (1992) formulation |
| 23 | **Growing data frames with `rbind` in loops** | Lines 1278, 1420, 1579, 1695 | Low | Slow for large iterations. Pre-allocate or use `lapply` + `bind_rows` |

### 3.6 Code Quality

| # | Issue | Location | Severity | Recommendation |
|---|-------|----------|----------|----------------|
| 24 | **Commented-out `write_csv` calls throughout** | Multiple | Low | Uncomment when generating final results, or use a flag variable to toggle output writing |
| 25 | **Variable shadowing**: `c` used as loop variable name | Lines 1270, 1358 | Low | Shadows base R function `c()`. Use `col` or `col_name` instead |
| 26 | **Duplicate CPT code** between broad framing (lines 1355-1424) and narrow framing (lines 1630-1699) | Sections 5.3, 6.2 | Medium | Extract into a reusable function `compute_cpt(profits, ref_point, a, l, g)` |
| 27 | **Duplicate EUT code** between broad framing (lines 1269-1285) and narrow framing (lines 1565-1581) | Sections 5.2, 6.1 | Medium | Extract into a reusable function |
| 28 | **No error handling** for edge cases (e.g., all gains, all losses, negative prices) | Throughout | Low | Add assertions or checks for simulation outputs |
| 29 | **Hardcoded label in plot**: `"Strike = 0.375"` | Line 1028 | Low | Should reference `strike_vec[4]` dynamically |
| 30 | **`v()` function depends on global `ref_point`** | Line 1346 | Low | Pass `ref_point` as a parameter for clarity and reusability |

---

## 4. Priority Summary

### Critical (verify correctness)
- ~~**#7/#8**: QVAR coefficient-matrix column order — may be producing incorrect price simulations~~ **RESOLVED 2026-03-03**
- ~~**#12**: Copula is estimated but never used in simulation; single quantile draw imposes perfect correlation~~ **RESOLVED 2026-03-03**
- ~~**#22**: CPT decision weight computation may not follow standard formulation for gains~~ **RESOLVED 2026-03-03**

### Important (affects results interpretation)
- ~~**#5**: Yield function omits weather/rotation controls despite creating dummies~~ **DOCUMENTED 2026-03-03** (excluded due to collinearity)
- ~~**#16**: M-V approximation for EUT may not hold for skewed profit distributions~~ **DOCUMENTED 2026-03-03**
- ~~**#19**: Reference point choice in broad framing needs theoretical justification~~ **DOCUMENTED 2026-03-03**

### Improvement (code quality / maintainability)
- ~~**#26/#27**: Extract duplicate CPT/EUT code into functions~~ **RESOLVED 2026-03-03**
- ~~**#24**: Toggle output writing with a flag~~ **RESOLVED 2026-03-03**
- ~~**#25**: Avoid shadowing `c()`~~ **RESOLVED 2026-03-03**

### Still open
- ~~**#1**: Hardcoded grid location~~ **DOCUMENTED 2026-03-03** (SimUID 52338 is the location of interest)
- **#2**: Irregular time intervals not handled (medium)
- **#3**: PPI 2025 duplicated from last year (low)
- **#9**: VAR serial correlation / ARCH effects (medium)
- **#10**: Linear copula model (low)
- **#13**: Starting price is last observed value (medium, in TODO)
- ~~**#14**: Simulation horizon mismatch in docs~~ **RESOLVED 2026-03-03** (CLAUDE.md already says "1 year")
- **#15**: Single exercise date (medium)
- **#17**: Delta parameter interpretation (low)
- **#18**: Different delta ranges undocumented (low)
- **#21**: Weighting function guard for p=0 (medium)
- **#28**: No edge-case error handling (low)
- **S1-F**: Yield residuals are skewed but simulated with `rnorm()` — understates downside yield risk. Consider a skewed distribution, bootstrapping residuals, or adding a third-moment function (medium)

---

## 4b. Section 1.0 Deep Review (2026-03-03)

Detailed analysis of the nitrogen-yield relationship (lines 20-170).

### Documented (comments added)

| Issue | Description | Location |
|-------|-------------|----------|
| **A** | Variation function uses absolute residuals (not squared). E[&#124;e&#124;] ≈ 0.798*sigma underestimates true SD by ~20%. Internally consistent with downstream `rnorm(sd=...)` usage. | Line ~158 |
| **B** | Only intercept + nitrogen coefficients used downstream, implicitly assuming rainfed + conventional BAU tillage (reference categories). Paper paragraph drafted. | Lines ~1349-1362 |
| **C** | N = 200 interpolates between FTN = 180 (BAU) and 250 in the EPIC data. | Line ~1342 |
| **D** | FTN distribution is heavily unbalanced — 67% at 180 kg/ha. | Line ~131 |
| **E** | Weather/rotation exclusion: effects captured in residual variance since downstream simulation uses fixed management scenario. | Line ~136 |
| **F** | Normal distribution assumption for stochastic yield errors. Skewness in residuals would warrant a different distribution. | Line ~1370 |

### No action needed

| Issue | Description | Reason |
|-------|-------------|--------|
| **G** | Yield can go negative in simulation. | P(yield < 0) ≈ 0.0001% at N = 200; negligible. |

---

## 5. Relationship to TODO Items

| TODO Item | Status | Related Issues |
|-----------|--------|---------------|
| Fix yield estimation under different scenarios | Open | #5, #6 |
| Deal with outliers in price time series | Open | Partially addressed by robust regression; price series needs explicit outlier treatment |
| Change starting price for the simulation | Open | #13 |
| Update graph with vertical bars and colors | Open | #29 (hardcoded label); strike bars already exist in `plot_original` |
| Integrate ambiguity theory | Open | New section needed; could extend CPT framework |
| ~~Create minimum price contract for output prices~~ | **DONE 2026-03-11** | Sections 4.2, 4.2.1 added; wc_* and both_* columns in Sections 5 and 6 |
