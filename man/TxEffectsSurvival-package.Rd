\name{TxEffectsSurvival-package}
\alias{TxEffectsSurvival-package}
\alias{TxEffectsSurvival}
\docType{package}

\title{
Inference for Treatment Effects on Terminal and Non-Terminal Events under Competing Risks
}

\description{
This package implements procedures for treatment effect inference on a
terminal event and a non-terminal event. It estimates the effect size using either
event-specific win ratios, RICH or RITCH, and the corresponding confidence intervals
are provided. For easy comparison, results for the win ratio (Finkelstein and Schoenfeld
1999, Pocock et al. 2012; Bebu and Lachin 2016) are also included. Moreover, testing
procedures are provided for the global null hypothesis of no treatment effect on either
event. Additionally, it offers confidence intervals for linear combinations of the
effect size estimates using fixed or data-driven weights. Finite-sample correction and
transformations are used to improve type I error control in moderate samples. For
further details, see Yang and Troendle (2021), Yang et al. (2022), and Yang (2025).
}

\details{
\tabular{ll}{
Package: \tab TxEffectsSurvival\cr
Type: \tab Package\cr
Version: \tab 0.0.1\cr
License: \tab GPL (>= 3)\cr
}
}

\author{
Daewoo Pak and Song Yang \\
Maintainer: Daewoo Pak \email{dpak@yonsei.ac.kr}
}

\value{
No return value, called for side effects.
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

    Yang, S., & Troendle, J. (2021).
    Event-specific win ratios and testing with terminal and non-terminal events.
    \emph{Clinical Trials}, 18(2), 180–187.
    \doi{10.1177/1740774520972408}

    Bebu, I., & Lachin, J. M. (2016).
    Large sample inference for a win ratio analysis of a composite outcome based
    on prioritized components.
    \emph{Biostatistics}, 17(1), 178–187.
    \doi{10.1093/biostatistics/kxv032}

    Pocock, S. J., Ariti, C. A., Collier, T. J., & Wang, D. (2012).
    The win ratio: a new approach to the analysis of composite endpoints in
    clinical trials based on clinical priorities.
    \emph{European Heart Journal}, 33(2), 176–182.
    \doi{10.1093/eurheartj/ehr352}

    Finkelstein, D. M., & Schoenfeld, D. A. (1999).
    Combining mortality and longitudinal measures in clinical trials.
    \emph{Statistics in Medicine}, 18(11), 1341–1354.
    \doi{10.1002/(SICI)1097-0258(19990615)18:11<1341::AID-SIM129>3.0.CO;2-7}
}


\keyword{package}

\seealso{
  \code{\link{tnt.analysis}}
}

\examples{
# Example: Treatment effect analysis using colon cancer trial data
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
}
