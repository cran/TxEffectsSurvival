---
output:
  html_document: default
  pdf_document: default
  word_document: default
---

# TxEffectsSurvival

The **TxEffectsSurvival** package implements inference procedures for treatment
effects on terminal and non-terminal events when both may be subject to competing
risks. An example is when the time to heart failure hospitalization (H) and the time
to cardiovascular death (CVD) are of interest, in which case H is subject to the
competing risk of all-cause mortality, and CVD is subject to the competing risk of
non-cardiovascular deaths. The setting is quite general and includes the often-studied
semi-competing risks model as a special case in which the only competing risk is that
of the terminal event to the non-terminal event. While the package is primarily designed
for joint analysis of terminal and non-terminal events, univariate survival analysis
with a single event is also available by assigning identical inputs to both event
types, as illustrated by the example below. To apply the package to a situation where more
than one non-terminal event is relevant, define the non-terminal event to be the
composite of the non-terminal events. 

The package provides inference procedures based on several different methods:

-   **Event-Specific Win Ratios (ESWR):** Event-specific win ratios are obtained separately for terminal and non-terminal events, and then used to form global testing procedures (linear combination, maximum, and $\chi^2$), as described in Yang et al. (2022), following Yang and Troendle (2021) which proposed the event-specific win ratios. Under the proportional hazards conditions (PH) as specified in Yang et al. (2022), the ESWRs consistently estimate relevant (cause-specific) hazard ratios.

-   **Log-Rank Generalized Rank Estimator (LRGRE):** This is a weighted ESWR and corresponds to the usual log-rank test. Under PH, the LRGREs consistently estimate relevant (cause-specific) hazard ratios. Global tests can also be obtained from the LRGREs for the terminal and non-terminal events.

-   **Ratio of Integrated Cumulative Hazard (RICH):** RICH can be interpreted as the ratio of the area under the Nelson–Aalen curve for a specific cause, treatment over control. It is a nonparametric measure, and reduces to the (cause-specific) hazard ratio under PH. Moreover, corresponding inference procedures remain valid and meaningful even if the proportional hazards conditions fail to hold. RICH related procedures are proposed by Yang (2025).

-   **Ratio of Integrated Transformed Cumulative Hazard (RITCH):** RITCH can be interpreted as the area under a transformed Nelson–Aalen estimator for a specific cause, treatment over control. It provides another nonparametric measure. Corresponding inference procedures remain valid and meaningful even if the proportional hazards conditions fail to hold. While other choices are possible, the package uses the transform $\Phi(t) = 1 - \exp(-t)$. When no competing risks are present, the corresponding RITCH is asymptotically equivalent to the ratio of restricted mean time lost. RITCH related procedures are proposed by Yang (2025).

-   **Linear Combination Tests:** Global tests are obtained by linear combination of event-specific estimates using pre-specified weights. The default weight is $l = (0.5, 0.5)$. They are obtained using ESWR, LRGRE, RICH, or RITCH. In addition, linear combinations with data-driven weights are also available using ESWR.

-   **Maximum Test:** A global test based on the maximum of event-specific test statistics. It ensures a confidence interval–testing coherency and tends to offer greater power when treatment effects are more differential across event types. They are also obtained using ESWR, LRGRE, RICH, or RITCH.

-   **Test of Proportional Hazards:** A test for assessing whether the proportional hazards assumption holds for the terminal event and the non-terminal event, as described in Yang et al. (2022).

-   **Test of Equal Hazard Ratios:** A test for evaluating whether the hazard ratios for the terminal and non-terminal events are equal, assuming proportional hazards hold for both. The test is also proposed by Yang et al. (2022).

------------------------------------------------------------------------

## Motivation

In many clinical trials:

-   Both a non-terminal event and a terminal event are clinically relevant.

-   The terminal event may censor the non-terminal event, but not vice versa.

-   Both the terminal and non-terminal events may be subject to competing risks.

-   Traditional time-to-first-event analyses (log-rank, Cox regression on composites) mix the terminal and non-terminal events in an opaque manner. More importantly, they do not utilize all data, resulting in reduced power.

The **TxEffectsSurvival** package explores the specific data structure to provide
valid and interpretable inference. For easy comparison, the package also provides
results for the win ratio (Finkelstein and Schoenfeld 1999, Pocock et al. 2012; Bebu
and Lachin 2016). 

## Details
The semiparametric results were developed in Yang et al. (2021), following the
introduction of event-specific win ratios by Yang and Troendle (2021), and the
nonparametric results were developed in Yang (2025). Compared with conventional
survival methods, these procedures utilize more information for improved
efficiency and accommodate competing risks within the semi-competing risks
framework. In this package, effect size estimates and confidence intervals are
obtained for each event type, and several testing procedures are implemented for
the global null hypothesis of no treatment effect on either terminal or
non-terminal events. Furthermore, tests of proportional hazards assumptions,
under which the event-specific win ratios converge to hazard ratios, and tests
of equal hazard ratios are provided. Asymptotic properties of these inference
procedures are discussed in Yang et al. (2022) and Yang (2025). Transformations
and finite-sample corrections described in the references are applied to achieve
better type I error control for moderately sized data sets.


------------------------------------------------------------------------

## Installation

``` r
install.packages("TxEffectsSurvival")
library(TxEffectsSurvival)
```
------------------------------------------------------------------------

## Example

The package includes the dataset `colon_wr`, which is derived from `survival::colon`. It is restricted to subjects in the `Obs` and `Lev+5FU` arms (excluding `Lev`), and organized for semi-competing risks analyses with a non-terminal event (recurrence) and a terminal event (death).\
Each row corresponds to a single subject.

``` r
## Example 1: joint analysis (non-terminal + terminal)
## Load example data
data("colon_wr")

## Map variables
yh   <- colon_wr$yh      # non-terminal event time
hcen <- colon_wr$hcen    # censoring indicator for non-terminal event
yd   <- colon_wr$yd      # terminal event time
dcen <- colon_wr$dcen    # censoring indicator for terminal event
z    <- colon_wr$z       # treatment (0 = control, 1 = treatment)

## Run the analysis
## RICH / RITCH (and other measures) with tests and CIs are included in the output
res <- tnt.analysis(yh, hcen, yd, dcen, z, 
  lin   = c(0.5, 0.5),   # optional linear combination weights
  alpha = 0.05
)

print(res)
# Example: pull selected results
# res$global.tests$lincomb
# res$ci.results$lincomb

## Example 2: univariate analysis (workaround)
## As shown in the example, you can handle a single-event setting by
## assigning identical inputs to both event types.
res_uni <- tnt.analysis(yh, hcen, yh, hcen, z)
res_uni$global.tests$lincomb
res_uni$ci.results$lincomb

```
------------------------------------------------------------------------
## References

- Yang, S. (2025). Testing and Estimation of Treatment Effects in Clinical Trials
  for Terminal and Non-Terminal Events Subject to Competing Risks.  
  *Statistics in Medicine.*

- Yang, S., Troendle, J., Pak, D., & Leifer, E. (2022). Event‐specific win ratios
  for inference with terminal and non‐terminal events.  
  *Statistics in Medicine*, 41(7), 1225–1241.

- Yang, S., & Troendle, J. (2021). Event-specific win ratios and testing with
  terminal and non-terminal events.  
  *Clinical Trials*, 18(2), 180–187.

- Bebu, I., & Lachin, J. M. (2016). Large sample inference for a win ratio
  analysis of a composite outcome based on prioritized components.  
  *Biostatistics*, 17(1), 178–187.

- Pocock, S. J., Ariti, C. A., Collier, T. J., & Wang, D. (2012). The win ratio:
  a new approach to the analysis of composite endpoints in clinical trials based on
  clinical priorities.  
  *European Heart Journal*, 33(2), 176–182.

- Finkelstein, D. M., & Schoenfeld, D. A. (1999). Combining mortality and
  longitudinal measures in clinical trials.  
  *Statistics in Medicine*, 18(11), 1341–1354.
