\name{rootogram_zidw}
\title{Rootogram for a ZIDW fit}
\alias{rootogram_zidw}
\usage{
rootogram_zidw(object, type = c("hanging", "standing", "suspended"),
               sqrt = TRUE, ref_line = TRUE, warn_limits = TRUE, 
               fitted_colour = "steelblue", bar_colour = NA, 
               bar_fill = "grey", ref_line_colour = "black",
               warn_line_colour = "black", ylab = NULL, xlab = NULL, ...)
}
\description{
  A rootogram is a model diagnostic tool that assesses the goodness-of-fit of a statistical model. The rootogram is drawn using \code{\link[ggplot2:ggplot]{ggplot2::ggplot()}} graphics.
}

\arguments{
  \item{object}{A ZIDW model object.}
  \item{type}{Type of rootogram; default is \code{type = "hanging"}.}
  \item{sqrt}{Logical; show the observed and fitted frequencies (default is \code{sqrt = TRUE}).}
  \item{ref_line}{Logical; draw a reference line at zero (default is \code{ref_line = TRUE}).}
  \item{warn_limits}{Logical; draw Tukey's warning limit lines at +/- 1 (default is \code{warn_limits = TRUE}).}
  \item{fitted_colour, bar_colour, bar_fill, ref_line_colour, warn_line_colour}{Colors used to draw the respective element of the rootogram.}
  \item{xlab, ylab}{Character; labels for the x and y axis of the rootogram. May be missing (\code{NULL}), in which case suitable labels will be used.}
  \item{...}{Arguments passed to other methods.}
} 
  \value{
  A \code{ggplot} object.
  
} \seealso{
  \code{\link[gratia]{rootogram}}
} \references{
  Kleiber, C. and Zeileis, A. (2016). Visualizing Count Data
Regressions Using Rootograms. \emph{The American Statistician}, \bold{70}:3, 296--303.
}



\examples{ 
## data
data("rwm1984", package = 'COUNT')

zidw_out <- zidw_reg(docvis ~ outwork + female + age, 
                     ziformula = ~ age + female,
                     betaformula = ~ 1, data = rwm1984[1:100, ], beta = .1)
rootogram_zidw(zidw_out)
}

\keyword{file}


