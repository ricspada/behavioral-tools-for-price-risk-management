# TODO

## Open Items
- [ ] Integrate ambiguity theory
- [ ] Fix yield estimation under different scenarios
- [ ] Deal with outliers in price time series
- [ ] Change starting price for the simulation (Issue 13)
    - [ ] Use historical average or expected price
- [ ] Bootstrapping (Frederic)
- [ ] Strike prices based on historical or simulated prices?

## Completed
- [x] Fix X-matrix column order in QVAR (Issue #7)
- [x] Integrate copula into price simulation (Issue #12)
- [x] Fix CPT decision weights for gains (Issue #22)
- [x] Document yield function control variable exclusion (Issues #5, #6)
- [x] Document EUT Mean-Variance approximation (Issue #16)
- [x] Document reference point choice (Issue #19)
- [x] Add SAVE_OUTPUTS toggle flag (Issue #24)
- [x] Rename loop variable `c` to `col_name` (Issue #25)
- [x] Extract CPT/EUT into helper functions (Issues #26, #27, #30)
- [x] Fix hardcoded strike label (Issue #29)
- [x] Use list pre-allocation instead of rbind (Issue #23)
- [x] Create a minimum price contract for output prices
