# Methods note: fixed risk-aversion coefficient across framings

Draft text for the methods (and a companion sentence for results) explaining why
the CARA coefficient `delta` is held constant across broad and narrow framing.

---

## For the methods section

We model expected-utility preferences with a constant absolute risk aversion
(CARA) utility function, $u(\pi) = -e^{-\delta \pi}$, defined over per-hectare
profit $\pi$. CARA is chosen in preference to constant relative risk aversion
because the analysis evaluates profit changes and contract payoffs that may be
negative, for which a power utility is ill-defined. We calibrate the absolute
risk-aversion coefficient $\delta$ to span the relative-risk-aversion categories
of Hardaker et al. (2004), using the identity $R = \delta\,\bar{\pi}$, where
$\bar{\pi}$ is mean per-hectare profit; the resulting values
$\delta \in \{0.0005, 0.001, 0.002, 0.003, 0.004\}$ correspond to relative risk
aversion ranging from "hardly" (0.5) to "extremely" (4) risk-averse.

| Hardaker label     | RRA $R$ | $\delta = R/\bar{\pi}$ | $\delta$ (rounded) | RP $= (\delta/2)\sigma^2$ | RP as % of mean |
|--------------------|---------|------------------------|--------------------|---------------------------|-----------------|
| hardly risk-averse | 0.5     | 0.00045                | 0.0005             | EUR 48/ha                 | 4.4%            |
| somewhat           | 1       | 0.00091                | 0.001              | EUR 96/ha                 | 8.7%            |
| rather             | 2       | 0.00182                | 0.002              | EUR 192/ha                | 17.5%           |
| very               | 3       | 0.00272                | 0.003              | EUR 289/ha                | 26.2%           |
| extremely          | 4       | 0.00363                | 0.004              | EUR 385/ha                | 35.0%           |

*RRA-to-$\delta$ mapping uses mean broad-framing profit $\bar{\pi} \approx$ EUR
1101/ha and profit variance $\sigma^2 \approx 192{,}500$ (sd $\approx$ EUR
439/ha). RP is the mean-variance risk premium.*

The same coefficient $\delta$ is applied under both broad and narrow framing.
Because CARA exhibits no wealth effects — $e^{-\delta(W+X)} = e^{-\delta W}\,
e^{-\delta X}$ — any difference between the two framings reflects only the
statistical dependence between baseline farm profit $W$ and the contract payoff
$X$, and not a mechanical change in curvature. Holding $\delta$ fixed therefore
isolates the framing effect to the bracketing decision itself: broad framing
accounts for the covariance through which a contract reduces total-profit
variability, whereas narrow framing evaluates the contract payoff in isolation.

---

## Companion sentence for the results section

A direct consequence of local risk neutrality is that, under narrow framing, the
expected-utility risk premium for the fertilizer contracts is negligible: with a
single coefficient calibrated to whole-farm profit, the small isolated payoffs
(standard deviation of roughly EUR 6–12/ha) are valued in a near-risk-neutral
manner, as expected-utility theory requires for small stakes. This is the
intended rational benchmark against which the cumulative-prospect-theory results,
which display loss aversion over the same isolated payoffs, are contrasted.

---

## Reference to add

Hardaker, J.B., Huirne, R.B.M., Anderson, J.R., Lien, G. (2004). *Coping with
Risk in Agriculture*, 2nd ed. CABI Publishing.
