# Plan: Minimum Price Contract for Wheat (Output)

## Context

The existing analysis (Sections 4–6 of `MainScript.R`) evaluates a **maximum price contract** (call option) on fertilizer — hedging the farmer's input cost risk. This plan adds a symmetric **minimum price contract** (put option) on wheat — hedging the farmer's revenue risk — and evaluates both contracts together (diagonal paired strikes).

---

## Critical Files

- `2_Scripts/MainScript.R` — all substantive changes (Steps 1–8)
- `2_Scripts/ContractPayoff.R` — add wheat payoff diagram (Step 9)
- `TODO.md` — mark task complete after implementation

---

## Implementation Steps

### Step 1 — Move `sim_w_p` extraction to Section 4 (prerequisite)

`sim_w_p` is currently extracted in Section 5 (~line 1391) but is needed in Section 4 for pricing the put option. Move the block to Section 4, immediately after `sim_f_p` is built (~line 1101):

```r
# Wheat prices at harvest (t=13, early July)
output_selling_time <- 13 + (0:(simulated_years - 1)) * obs_per_year
sim_w_p <- array(dim = c(Ns, simulated_years))
for (t in 1:simulated_years) {
  sim_w_p[, t] <- exp(as.vector(Ys_levels[2, output_selling_time[t], ]))
}
```

Delete the duplicate block from Section 5, replace with a comment pointing to Section 4.

---

### Step 2 — Add wheat put option pricing (Section 4, new subsection 4.2)

After `premium_results` (fertilizer, ~line 1125), add:

```r
## 4.2 Pricing based on simulated wheat prices (minimum price / put option) ----

T_years_w <- 1.0                              # wheat sold ~12 months from contract initiation
discount_factor_w <- 1 / (1 + int)^T_years_w  # reuses int already defined

# Lower quantiles = downside protection (floor on revenue)
min_prices <- quantile(sim_w_p, probs = c(0.1, 0.2, 0.3, 0.4, 0.5))
premium_array_w <- numeric(length(min_prices))

for (i in seq_along(min_prices)) {
  strike_w <- min_prices[i]
  payoff_w <- pmax(strike_w - sim_w_p, 0)       # put option payoff per tonne
  premium_array_w[i] <- mean(payoff_w * discount_factor_w)
}

premium_results_w <- data.frame(
  strike_price   = min_prices,
  option_premium = premium_array_w
)
```

**Rationale for quantiles 0.1–0.5:** Symmetric to the fertilizer call's 0.5–0.9 (same 40-percentile span from the median, but lower tail — a put protects against low wheat prices).

**Rationale for T_years_w = 1.0:** Fertilizer purchased at t=5 (March ≈ 6 months from autumn), wheat sold at t=13 (July ≈ 12 months from autumn). The put matures at delivery.

---

### Step 3 — Add wheat density plot with strike levels (Section 4, subsection 4.2.1)

Mirror the fertilizer density plot (subsection 4.1.1) with:
- `sim_w_p` on x-axis
- `pmax(sim_w_p, k_w)` (effective floored price received) as the adjusted distribution
- Axis label: `"Simulated wheat prices and strike levels"`
- Output: `3_Outputs/Fig_PriceDensities_SimulatedStrikes_Wheat.pdf`

---

### Step 4 — Add wheat forward price benchmark (Section 5)

After `sim_f_p_forward` (~line 1406), add:

```r
# Last observed wheat price, discounted to harvest (eur/t, no /1000 unlike fertilizer)
sim_w_p_forward <- rep(tail(w_data$Avg_Price, 1) * discount_factor_w, Ns)
```

---

### Step 5 — Add wheat-only and combined profit columns to `df_simulated_profit` (Section 5)

After the existing fertilizer `ac_*` loop:

```r
# Wheat minimum-price contracts (wc_j): spot fertilizer + floored wheat revenue
for (j in 1:nrow(premium_results_w)) {
  k_w <- premium_results_w$strike_price[j]
  c_w <- premium_results_w$option_premium[j]
  profit_wc <- (pmax(sim_w_p, k_w) * simulated_yield) -
               (c_w * simulated_yield) -    # put premium per tonne × yield
               (sim_f_p * CAN)              # fertilizer at spot
  df_simulated_profit[[paste0("wc_", j)]] <- profit_wc
}

# Combined contracts (both_i): diagonal — fertilizer strike i paired with wheat strike i
for (i in 1:nrow(premium_results)) {
  k_f <- premium_results$strike_price[i];   c_f <- premium_results$option_premium[i]
  k_w <- premium_results_w$strike_price[i]; c_w <- premium_results_w$option_premium[i]
  profit_both <- (pmax(sim_w_p, k_w) * simulated_yield) - (c_w * simulated_yield) -
                 (pmin(sim_f_p, k_f) * CAN) - (c_f * CAN)
  df_simulated_profit[[paste0("both_", i)]] <- profit_both
}
```

**Premium units:** `c_f` (eur/kg) × `CAN` (kg/ha) = eur/ha. `c_w` (eur/t) × `simulated_yield` (t/ha) = eur/ha. Both consistent.

Final `df_simulated_profit` columns: `nc, fc, ac_1..5, wc_1..5, both_1..5` — **17 columns total**.

---

### Step 6 — EUT/CPT calls in Section 5: no changes needed

`compute_eut()` and `compute_cpt()` already iterate over all columns of `df_simulated_profit` automatically. The CPT reference point `mean(df_simulated_profit$fc)` remains correct (expected profit under fertilizer-forward status quo).

---

### Step 7 — Add narrow framing wheat payoffs (Section 6)

After the existing fertilizer `ac_*` loop in Section 6:

```r
# Wheat put option: narrow framing payoff per ha
# Uses predicted_yield (fixed) to isolate price risk — consistent with CAN being fixed
for (j in 1:nrow(premium_results_w)) {
  k_w <- premium_results_w$strike_price[j]
  c_w <- premium_results_w$option_premium[j]
  profit_wc_nf <- predicted_yield * (pmax(k_w - sim_w_p, 0) - c_w)
  df_profit_narrowframing[[paste0("wc_", j)]] <- profit_wc_nf
}

# Combined narrow framing: diagonal — fertilizer call i + wheat put i
for (i in 1:nrow(premium_results)) {
  k_f <- premium_results$strike_price[i];   c_f <- premium_results$option_premium[i]
  k_w <- premium_results_w$strike_price[i]; c_w <- premium_results_w$option_premium[i]
  profit_both_nf <- CAN * (pmax(sim_f_p - k_f, 0) - c_f) +
                    predicted_yield * (pmax(k_w - sim_w_p, 0) - c_w)
  df_profit_narrowframing[[paste0("both_", i)]] <- profit_both_nf
}
```

**Yield choice:** `predicted_yield` (deterministic) isolates pure price risk in the contract payoff, consistent with `CAN` being fixed on the fertilizer side. Reference point remains 0.

---

### Step 8 — EUT/CPT calls in Section 6: no changes needed

Same rationale as Step 6.

---

### Step 9 — Update `ContractPayoff.R`

Add a second plot block (after the existing fertilizer diagram) for the wheat minimum price contract:

- x-axis: market wheat price `s_w`
- y-axis: effective price **received** (not paid)
- Three lines: `Spot = s_w`, `Fixed = fixed_price_w`, `MinPrice = pmax(s_w, k_w)`
- Vertical dashed red line at `k_w` labeled "price floor"
- Legend positioned top-left (lines converge at low prices)

---

## Detailed TODO Checklist

### MainScript.R — Section 4

- [x] **S4.1** Locate the `sim_f_p` extraction block (line ~1095–1101) and insert the `sim_w_p` extraction block immediately after it (Step 1)
- [x] **S4.2** Locate the duplicate `sim_w_p` block in Section 5 (~line 1388–1397) and delete it; add comment `# sim_w_p extracted in Section 4`
- [x] **S4.3** After `premium_results` (~line 1125), add subsection header `## 4.2 Pricing based on simulated wheat prices (minimum price / put option) ----`
- [x] **S4.4** Define `T_years_w <- 1.0` and `discount_factor_w <- 1 / (1 + int)^T_years_w`
- [x] **S4.5** Define `min_prices <- quantile(sim_w_p, probs = c(0.1, 0.2, 0.3, 0.4, 0.5))`
- [x] **S4.6** Initialize `premium_array_w <- numeric(length(min_prices))`
- [x] **S4.7** Write the put premium loop: `payoff_w <- pmax(strike_w - sim_w_p, 0)` → `mean(payoff_w * discount_factor_w)`
- [x] **S4.8** Build `premium_results_w <- data.frame(strike_price = min_prices, option_premium = premium_array_w)`
- [x] **S4.9** Add `if (SAVE_OUTPUTS) write_csv(premium_results_w, "3_Outputs/Table_ContractPayoff_Wheat.csv")`
- [x] **S4.10** Add subsection header `### 4.2.1 Plot distribution density (wheat put) ----`
- [x] **S4.11** Define `strike_vec_w <- premium_results_w$strike_price`
- [x] **S4.12** Build `plot_original_w` density plot for `sim_w_p` with vertical lines at `strike_vec_w`
- [x] **S4.13** Build `plot_list_w` — list of 5 density plots comparing `sim_w_p` vs `pmax(sim_w_p, k_w)` for each strike
- [x] **S4.14** Combine into `final_plot_w` using `patchwork` (same 2×3 layout as fertilizer plot)
- [x] **S4.15** Add `if (SAVE_OUTPUTS) ggsave("3_Outputs/Fig_PriceDensities_SimulatedStrikes_Wheat.pdf", ...)`

### MainScript.R — Section 5

- [x] **S5.1** After `sim_f_p_forward` (~line 1406), add `sim_w_p_forward <- rep(tail(w_data$Avg_Price, 1) * discount_factor_w, Ns)`
- [x] **S5.2** After the existing fertilizer `ac_*` loop (after the last `df_simulated_profit[["ac_5"]]` assignment), add a comment `# Wheat minimum-price contracts (wc_j) ----`
- [x] **S5.3** Write the `wc_*` loop: `profit_wc <- (pmax(sim_w_p, k_w) * simulated_yield) - (c_w * simulated_yield) - (sim_f_p * CAN)`
- [x] **S5.4** Add a comment `# Combined fertilizer max-price + wheat min-price (both_i, diagonal) ----`
- [x] **S5.5** Write the `both_*` loop (diagonal only): pairs `premium_results[i, ]` with `premium_results_w[i, ]`
- [x] **S5.6** Verify `ncol(df_simulated_profit) == 17` (nc, fc, ac_1..5, wc_1..5, both_1..5)
- [x] **S5.7** Confirm `compute_eut(df_simulated_profit, delta)` call requires no changes
- [x] **S5.8** Confirm `compute_cpt(df_simulated_profit, ref_point_broad, ...)` call requires no changes
- [x] **S5.9** Confirm `ref_point_broad <- mean(df_simulated_profit$fc)` remains the correct reference point
- [x] **S5.10** Add subsections 5.3.3 (wheat price stats) and 5.3.4 (wheat revenue stats) summary tables, mirroring the fertilizer tables

### MainScript.R — Section 6

- [x] **S6.1** After the existing fertilizer `ac_*` loop in Section 6, add a comment `# Wheat put option: narrow framing (wc_j) ----`
- [x] **S6.2** Write the `wc_*` narrow framing loop: `profit_wc_nf <- predicted_yield * (pmax(k_w - sim_w_p, 0) - c_w)` — uses `predicted_yield` (fixed), not `simulated_yield`
- [x] **S6.3** Add a comment `# Combined narrow framing (both_i, diagonal) ----`
- [x] **S6.4** Write the `both_*` narrow framing loop: sum of fertilizer call payoff + wheat put payoff, both scaled by respective quantities (`CAN` and `predicted_yield`)
- [x] **S6.5** Verify `df_profit_narrowframing` column count (fc, ac_1..5, wc_1..5, both_1..5 = 16 columns)
- [x] **S6.6** Confirm `compute_eut(df_profit_narrowframing, delta_nf)` call requires no changes
- [x] **S6.7** Confirm `compute_cpt(df_profit_narrowframing, ref_point_narrow, ...)` call requires no changes
- [x] **S6.8** Confirm `ref_point_narrow <- 0` remains correct for wheat and combined payoffs

### ContractPayoff.R

- [x] **CP.1** After the existing fertilizer plot code, add section header comment `### Minimum price (put option) payoff graph for wheat`
- [x] **CP.2** Define parameters: `k_w <- 100`, `fixed_price_w <- 100`, `s_w <- seq(50, 150, length.out = 500)`
- [x] **CP.3** Build `df_w` with columns: `MarketPrice = s_w`, `Spot = s_w`, `Fixed = rep(fixed_price_w, ...)`, `MinPrice = pmax(s_w, k_w)`
- [x] **CP.4** Build ggplot with three lines (Spot, Fixed, MinPrice) using same color scheme as fertilizer plot
- [x] **CP.5** Add vertical dashed red line at `k_w` with annotation "price floor"
- [x] **CP.6** Position legend top-left (lines converge at low prices, unlike fertilizer plot where they converge at high prices)
- [x] **CP.7** Use same axis styling as fertilizer plot (no tick labels, arrow annotations for Low/High)

### Final

- [x] **F.1** Mark "Create a minimum price contract for output prices" as complete in `TODO.md`
- [ ] **F.2** Run Sections 4–6 end-to-end and verify all checks in the Verification section below pass

---

## Verification

1. `premium_results_w` is 5×2 with non-negative premiums; higher strikes → higher put premiums.
2. `df_simulated_profit` has 17 columns; mean of `wc_*` > mean of `nc` minus premium (floor adds value net of cost).
3. `df_profit_narrowframing` has correct column count; mean payoff of `wc_*` ≈ 0 (fair-valued option under risk-neutral pricing).
4. Spot-check `both_1`: verify `profit_both = max(P_w, k_w)*Y - c_w*Y - min(P_f, k_f)*CAN - c_f*CAN` and that it lies between `ac_1` and `wc_1` column means.
5. `ContractPayoff.R` wheat diagram: floor shape — flat at `k_w` for low prices, follows spot for high prices.
