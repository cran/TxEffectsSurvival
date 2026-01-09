\name{colon_wr}
\alias{colon_wr}
\docType{data}

\title{Colon cancer trial example data (Obs vs Lev+5FU) for semi-competing risks analyses}

\description{
A derived dataset from \code{survival::colon}, restricted to subjects in the \code{Obs} and \code{Lev+5FU} arms (excluding \code{Lev}). It is organized for semi-competing risks analyses with a non-terminal event (recurrence) and a terminal event (death). Each row corresponds to one subject.
}

\usage{
data("colon_wr")
}

\format{
A data frame with one row per subject and the following variables:

\describe{
  \item{yh}{Time to \emph{non-terminal} event (recurrence).}
  \item{hcen}{Status for the non-terminal event: \code{1} = recurrence observed, \code{0} = censored.}
  \item{yd}{Time to \emph{terminal} event (death).}
  \item{dcen}{Status for the terminal event: \code{1} = death observed, \code{0} = censored.}
  \item{z}{Treatment indicator: \code{0} = \code{Obs}, \code{1} = \code{Lev+5FU}.}
}
}

\details{
The source data are \code{survival::colon}. We first restrict to \code{rx != "Lev"} and then, per subject, extract event time/status for recurrence (\code{etype = 1}) and death (\code{etype = 2}). The original treatment factor is recoded to a binary indicator \eqn{z \in \{0,1\}}.
}

\source{
\pkg{survival} package, dataset \code{colon}.
}

\examples{
# Load example dataset
data("colon_wr")

# Contents of colon_wr:
#   yh    : Time to non-terminal event (recurrence)
#   hcen  : Event indicator for recurrence (1 = recurrence, 0 = censored)
#   yd    : Time to terminal event (death)
#   dcen  : Event indicator for death (1 = death, 0 = censored)
#   z     : Treatment group (0 = Obs, 1 = Lev+5FU)

}
\keyword{datasets}
