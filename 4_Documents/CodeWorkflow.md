# MainScript.R — Code Workflow

## Overview

The script evaluates **price risk management contracts** (fertilizer call option + wheat put option) for a German wheat farmer, comparing decisions under **Expected Utility Theory (EUT)** and **Cumulative Prospect Theory (CPT)** with broad vs. narrow framing.

---

## Pipeline Diagram

```
 INPUT DATA
 ==========

 EPIC.csv                FertilizerPrices.csv     WheatPrices.csv
 (crop simulation)       (CAN biweekly)           (bread wheat biweekly)
       |                        |                        |
       v                        +----------+-------------+
                                           |
 +-----------------+              +--------v---------+
 | 1. YIELD MODEL  |              | 2. PRICE MODEL   |
 | Just & Pope     |              | Stationarity,    |
 | production fn   |              | cointegration,   |
 | (lmrob)         |              | QVAR estimation, |
 |                 |              | copula fitting   |
 | OUT:            |              |                  |
 |  predicted_yield|              | OUT:             |
 |  predicted_sigma|              |  QVAR coeffs    |
 +---------+-------+              |  copula params   |
           |                      +--------+---------+
           |                               |
           |                      +--------v---------+
           |                      | 3. SIMULATION    |
           |                      | 1000 Monte Carlo |
           |                      | price paths      |
           |                      | (copula-linked)  |
           |                      |                  |
           |                      | OUT:             |
           |                      |  sim_f_p (fert)  |
           |                      |  sim_w_p (wheat) |
           |                      +--------+---------+
           |                               |
           |                      +--------v---------+
           |                      | 4. CONTRACT      |
           |                      |    PRICING       |
           |                      |                  |
           |                      | Fert call option |
           |                      |  5 strikes       |
           |                      |  (50th-90th pct) |
           |                      |                  |
           |                      | Wheat put option |
           |                      |  5 strikes       |
           |                      |  (10th-50th pct) |
           |                      |                  |
           |                      | OUT:             |
           |                      |  premium_results |
           |                      |  premium_results_w|
           |                      +--------+---------+
           |                               |
           +---------------+---------------+
                           |
                  +--------v---------+
                  | 5. HELPER FNS    |
                  |  v()             |  CPT value function
                  |  weighting_fn()  |  probability weighting
                  |  compute_eut()   |  mean-variance approx.
                  |  compute_cpt()   |  full CPT evaluation
                  +--------+---------+
                           |
              +------------+------------+
              |                         |
     +--------v---------+     +--------v---------+
     | 6. BROAD FRAMING |     | 7. NARROW FRAMING|
     |                  |     |                  |
     | Joint yield +    |     | Price risk only  |
     | price risk       |     | (fixed yield)    |
     |                  |     |                  |
     | 17 scenarios:    |     | 16 scenarios:    |
     |  nc  (no cont.)  |     |  fc (forward)    |
     |  fc  (forward)   |     |  ac_1..5 (fert)  |
     |  ac_1..5 (fert)  |     |  wc_1..5 (wheat) |
     |  wc_1..5 (wheat) |     |  both_1..5 (comb)|
     |  both_1..5 (comb)|     |                  |
     |                  |     | Ref point: 0     |
     | Ref point:       |     |                  |
     |  E[profit|fwd]   |     | EUT + CPT        |
     |                  |     | evaluation       |
     | EUT + CPT        |     +------------------+
     | evaluation       |
     +------------------+


 OUTPUTS (3_Outputs/)
 ====================

 Figures (PDF):                          Tables (CSV):
  Fig_Prices                              Table_PriceSummary
  Fig_CombinedPrices                      Table_ContractPayoff
  Fig_DifferencedPrices                   Table_ContractPayoff_Wheat
  Fig_PriceDensity_CAN                    Table_ResultsCostsStats
  Fig_PriceDensity_Wheat                  Table_ResultsProfitStats
  Fig_PriceDensity_deflated               Table_ResultsWheatPriceStats
  Fig_PriceDensities_SimulatedStrikes     Table_ResultsWheatRevenueStats
  Fig_PriceDensities_SimulatedStrikes_W   Table_Results_EUT_Broadframing
                                          Table_Results_EUT_NarrowFraming
                                          Table_Results_CPT_BroadFraming
                                          Table_Results_CPT_NarrowFraming
                                          Table_Results_Profits_NarrowFraming
```

---

## Section Details

### Section 1 — Nitrogen-Yield Relationship (lines 20–200)

- **Input:** `EPIC.csv` (EPIC crop simulation for winter wheat at grid cell 52338)
- **Method:** Just & Pope production function with robust MM regression (`lmrob`)
  - Mean function: `YLD ~ sqrt(FTN) + FTN` (nitrogen response)
  - Variance function: `|residuals| ~ sqrt(FTN) + FTN`
- **Output:** `predicted_yield` and `predicted_sigma_yield` at N=200 kg/ha
- **Used by:** Sections 6 (stochastic yield draws) and 7 (fixed yield for narrow framing)

### Section 2 — Price Dynamics (lines 205–850)

- **Input:** `FertilizerPrices.csv` (CAN), `WheatPrices.csv` (bread wheat), PPI indices
- **Steps:**
  1. Stationarity tests (ADF, KPSS) — both series non-stationary in levels
  2. First-differencing of log prices
  3. Cointegration test (weak evidence at 10%)
  4. QVAR model selection (lag order, exogenous variables)
  5. QVAR estimation at quantiles 0.25, 0.50, 0.75 (jointly for fert & wheat)
  6. Copula estimation on QVAR residuals (captures cross-commodity dependence)
  7. Dynamic stability analysis
- **Output:** QVAR coefficients, inverse CDF parameters, copula parameters

### Section 3 — Price Path Simulation (lines 850–1080)

- **Method:** Monte Carlo simulation of 1,000 copula-linked biweekly price paths over 24 periods (1 year)
- **Starting values:** Last observed (log) prices
- **Output:**
  - `sim_f_p`: 1000 simulated fertilizer prices at purchase time (t=5, March)
  - `sim_w_p`: 1000 simulated wheat prices at harvest time (t=13, July)

### Section 4 — Insurance Contract Pricing (lines 1082–1350)

- **Fertilizer call option** (max price / ceiling):
  - 5 strike prices at 50th, 60th, 70th, 80th, 90th percentile of `sim_f_p`
  - Premium = discounted expected payoff `E[max(P_f - K, 0)]`
  - Discount: 6 months (T=0.5 years)
- **Wheat put option** (min price / floor):
  - 5 strike prices at 10th, 20th, 30th, 40th, 50th percentile of `sim_w_p`
  - Premium = discounted expected payoff `E[max(K - P_w, 0)]`
  - Discount: 12 months (T=1.0 years)
- **Output:** `premium_results` (5x2), `premium_results_w` (5x2), density plots

### Section 5 — Helper Functions (lines 1355–1470)

| Function | Purpose |
|---|---|
| `v(x, alpha, lambda)` | CPT value function (Tversky & Kahneman 1992) |
| `weighting_function(p, gamma)` | Probability weighting (Prelec 1998) |
| `compute_eut(df, delta)` | EUT certainty equivalents via mean-variance approximation |
| `compute_cpt(df, ref, ...)` | Full CPT value with rank-dependent weighting |

### Section 6 — Broad Framing (lines 1473–1995)

- **Profit formula:** `Profit = WheatRevenue - FertilizerCost - OptionPremiums`
  - Revenue uses `simulated_yield` (stochastic, drawn from Just & Pope model)
  - Both yield risk and price risk are present
- **17 scenarios:** nc, fc, ac_1–5, wc_1–5, both_1–5
- **Reference point (CPT):** `mean(profit under forward contract)`
- **Evaluation:** EUT across 6 risk aversion levels; CPT across varying alpha, lambda, gamma

### Section 7 — Narrow Framing (lines 1997–2223)

- **Payoff formula:** Contract payoff only (not full farm profit)
  - Fertilizer: `CAN * (max(P_f - K_f, 0) - premium_f)`
  - Wheat: `predicted_yield * (max(K_w - P_w, 0) - premium_w)`
  - Uses `predicted_yield` (fixed/deterministic) to isolate price risk
- **16 scenarios:** fc, ac_1–5, wc_1–5, both_1–5 (no "nc" — no contract has zero payoff)
- **Reference point (CPT):** 0 (break-even on the contract itself)
- **Evaluation:** Same EUT and CPT as broad framing

---

## Key Design Choices

| Choice | Rationale |
|---|---|
| QVAR (not VAR) | Captures heterogeneous price dynamics across quantiles; handles non-normal distributions |
| Copula simulation | Preserves marginal distributions while modelling cross-commodity dependence |
| Mean-variance EUT | Analytical tractability; valid approximation for near-normal profit distributions |
| Full CPT evaluation | Rank-dependent probability weighting captures behavioral biases (loss aversion, probability distortion) |
| Broad vs. narrow framing | Tests whether farmers evaluate contracts as part of total farm profit (broad) or as standalone investments (narrow) |
| `both_i` diagonal pairing | Current simplification — pairs fert strike i with wheat strike i; planned replacement with full 5x5 grid |
