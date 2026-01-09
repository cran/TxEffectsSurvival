\name{tnt.analysis}
\alias{tnt.analysis}
\alias{tnt.analysis.default}
\alias{print.tnt.analysis}

\title{Inference for Treatment Effects on Terminal and Non-Terminal Events under Competing Risks}
\usage{
\method{tnt.analysis}{default}(nonterm_time, nonterm_event, term_time, term_event, z,
lin = c(0.5, 0.5), alpha = 0.05, repnum = 1e6, tau_r = 0.9, ...)
}

\description{
\code{tnt.analysis} provides several procedures for evaluating treatment effects on terminal and non-terminal events in the presence of competing risks. Based on the methodologies proposed in Yang et al. (2022) and Yang (2025), the function implements a variety of robust and flexible semiparametric and nonparametric estimation and hypothesis testing procedures that accommodate complex event-time structures with semi-competing risks.

Specifically, the function primarily provides the following inferential procedures for estimation and hypothesis testing:

\itemize{
  \item \strong{Event-Specific Win Ratios (ESWR)}: Event-specific win ratios obtained separately for terminal and non-terminal events,
  enabling powerful and interpretable global tests (linear combination, maximum,
  and \eqn{\chi^2}), as proposed by Yang and Troendle (2021).
  \item \strong{Log-Rank Generalized Rank Estimator (LRGRE)}: A generalized rank-based estimator extending the log-rank test to semi-competing risks settings, providing robust group comparisons for non-terminal and terminal events.
  \item \strong{Integrated Cumulative Hazard (ICH)}: A nonparametric estimator calculated as the ratio of integrated cause-specific hazard in the treatment group over that in the control group. Values less than 1 indicate a beneficial effect. Robust to violations of proportional hazards.
  \item \strong{Integrated Transformed Cumulative Hazard (ITCH)}: A stabilized version of ICH using a monotone transformation \eqn{\Phi(t)=1 - \exp(-t)} to improve performance under heavy censoring or non-PH conditions.
  \item \strong{Linear Combination Tests}: Global tests that combine event-specific estimates using either user-specified or data-adaptive weights. The default weight is \eqn{l = (0.5, 0.5)}. Applicable to ESWR, LRGRE, ICH, and ITCH.
  \item \strong{Maximum Test}: A global test based on the maximum of event-specific test statistics, which rejects the global null if and only if at least one event-specific estimate is significant. It ensures coherency and tends to offer greater power when treatment effects differ across event types.
  \item \strong{Test of Proportional Hazards}: A diagnostic test for assessing whether the proportional hazards assumption holds for the terminal event and the non-terminal event.
  \item \strong{Test of Equal Hazard Ratios}: A diagnostic test for evaluating whether the hazard ratios for the terminal and non-terminal events are equal, assuming proportional hazards hold for both.
}
}

\arguments{
  \item{nonterm_time}{A numeric vector indicating the time to the non-terminal event or censoring.}
  \item{nonterm_event}{A binary indicator for the non-terminal event (1 = event, 0 = censored).}
  \item{term_time}{A numeric vector indicating the time to the terminal event or censoring.}
  \item{term_event}{A binary indicator for the terminal event (1 = event, 0 = censored).}
  \item{z}{A numeric vector for the treatment group indicator (1 = treatment, 0 = control).}
  \item{lin}{A numeric vector of length 2 specifying the weights for the linear combination of event-specific estimates. The first element corresponds to the non-terminal event and the second to the terminal event. The weights must be non-negative and sum to one. Default is \code{lin=c(0.5, 0.5)}.}
  \item{alpha}{The significance level used for confidence intervals and hypothesis testing. Default is \code{alpha=0.05}.}
  \item{repnum}{The number of Monte Carlo replications used to simulate bivariate normal distributions for critical value and p-value estimation. Default is \code{repnum=1e6}.}
  \item{tau_r}{A numeric value between 0 and 1 specifying the truncation proportion for the integration range \eqn{[0, \tau]} used in computing ICH and ITCH, where \eqn{\tau = \texttt{tau\_r} \times \text{(total trial duration)}}. Default is \code{tau_r=0.9}. For trials with very long follow-up, a smaller value (e.g., 0.8) may improve numerical stability. See \code{details} for the definition of total trial duration.}
  \item{...}{Additional arguments (currently ignored). Included for future extensibility.}
}

\value{
An object of class \code{"tnt.analysis"}, which is a named list containing the following elements:
  \item{sum_table}{A contingency table summarizing the observed event patterns.}
  \item{wrlin}{Estimated win ratio from a linear combination of event-specific win ratios using user-specified weights \code{lin}.}
  \item{rho}{Estimated correlation between the event-specific test statistics (non-terminal and terminal).}
  \item{pvalinall}{Two-sided p-value for the linear combination test.}
  \item{lin_ar}{Data-driven linear combination weights.}
  \item{wrlin_ar}{Estimated win ratio from the data-adaptive linear combination \code{lin_ar}.}
  \item{lincileswr}{Confidence interval for the data-adaptive win ratio estimate, constructed on the log-transformed scale.}
  \item{peswrlog}{Two-sided p-value for the global test based on the log-transformed ESWR and data-adaptive weights.}
  \item{chi, pvachi}{Chi-square test statistic and its corresponding p-value.}
  \item{mx, pvalmx}{Maximum-type test statistic and corresponding Monte Carlo-based p-values.}
  \item{logrank}{Standard log-rank test result.}
  \item{wrtest0}{Unadjusted win ratio test result based on time to first event.}
  \item{zvaephl, pvaephl}{Z-statistic and p-value for testing the equality of hazard ratios between terminal and non-terminal events.}
  \item{pvalph}{P-value for testing the proportional hazards assumptions for both event types.}
  \item{ci1, ci2}{Pointwise confidence intervals for the event-specific win ratios for the non-terminal and terminal events, respectively.}
  \item{cil}{Confidence interval for the linear combination of win ratios using user-defined weights \code{lin}.}
}


\details{
These methods follow the statistical inference framework proposed by Yang (2025), without relying on parametric assumptions or restrictive model specifications. They are well-suited for semi-competing risks settings, where non-terminal events (e.g., hospitalizations) may be censored by terminal events (e.g., death), and both types are clinically important.

Global tests—including the linear combination, maximum, and \eqn{\chi^2} tests—offer flexible approaches to detecting treatment effects on either or both event types. Inference remains valid regardless of whether the proportional hazards assumption holds, with ICH and ITCH estimators providing additional robustness under non-proportional hazards.

In the implementation, the total trial duration is defined as \eqn{\tau = \max(\{yh, yd\}) + 0.001}. This upper bound is applied in the estimation of both ESWR and LRGRE. For ICH and ITCH, the analysis window is further truncated to \eqn{\tau_0 = \tau \times \texttt{tau\_r}}, where \code{tau_r} is a user-defined proportion (default = 0.9).
Log transformations are used for the data-driven test, while log-log transformations are used for confidence intervals for non-terminal and terminal events. To improve finite-sample performance, a small-sample correction factor (\eqn{cn}) is applied to standard errors and test statistics. Specifically, \eqn{cn = 1 + 15 / \max(n, 100)} is used for the chi-square test and the maximum test, while \eqn{cn = 1 + 10 / \max(n, 100)} is used for the data-driven test and the corresponding confidence intervals.

To apply the package for univariate survival data analysis, see Example 2 in the Examples section below.
}

\examples{
# Example 1: joint analysis of non-terminal and terminal events
data("colon_wr")

yh   <- colon_wr$yh     # non-terminal event time
hcen <- colon_wr$hcen   # censoring indicator for non-terminal event
yd   <- colon_wr$yd     # terminal event time
dcen <- colon_wr$dcen   # censoring indicator for terminal event
z    <- colon_wr$z      # treatment group (0 = control, 1 = treatment)

# Run the analysis
# The output includes test statistics and confidence intervals
# for RICH and RITCH measures (among others)
res <- tnt.analysis(
  yh, hcen, yd, dcen, z,
  lin   = c(0.5, 0.5),
  alpha = 0.05
)

# Display results
print(res)

# Example 2: univariate analysis (workaround)
# tnt.analysis is primarily designed for joint analysis of non-terminal and terminal events,
# but univariate analysis can be conducted by assigning identical inputs to both.
\donttest{
res <- tnt.analysis(yh, hcen, yh, hcen, z)
res$global.tests$lincomb
res$ci.results$lincomb
}
}

\references{
    Yang, S. (2025).
    Testing and Estimation of Treatment Effects in Clinical Trials for Terminal
    and Non-Terminal Events Subject to Competing Risks.
    \emph{Statistics in Medicine}.
    \doi{10.1002/sim.70205}

    Yang, S., Troendle, J., Pak, D., & Leifer, E. (2022).
    Event‐specific win ratios for inference with terminal and non‐terminal events.
    \emph{Statistics in Medicine}, 41(7), 1225–1241.
    \doi{10.1002/sim.9266}
}
