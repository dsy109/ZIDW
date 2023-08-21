\name{rootogram_zidw}
\title{Rootogram for ZIDW fit}
\alias{rootogram_zidw}
\usage{
rootogram_zidw(object, type = c("hanging", "standing", "suspended"),
               sqrt = TRUE, ref_line = TRUE, warn_limits = TRUE, 
               fitted_colour = "steelblue", bar_colour = NA, 
               bar_fill = "grey", ref_line_colour = "black",
               warn_line_colour = "black", ylab = NULL, xlab = NULL, ...)
}
\description{
  A rootogram is a model diagnostic tool that assesses the goodness of fit of a statistical model. The rootogram is drawn using \code{\link[ggplot2:ggplot]{ggplot2::ggplot()}} graphics.
}

\arguments{
  \item{object}{zidw object}
  \item{type}{type of rootogram; default is hanging}
  \item{sqrt}{logical; show the observed and fitted frequencies}
  \item{ref_line}{logical; draw a reference line at zero?}
  \item{warn_limits}{logical; draw Tukey's warning limit lines at +/- 1?}
  \item{fitted_colour, bar_colour, bar_fill, ref_line_colour, warn_line_colour}{colours used to draw the respective element of the rootogram.}
  \item{xlab, ylab}{character; labels for the x and y axis of the rootogram. May be missing (NULL), in which case suitable labels will be used.}
  \item{...}{arguments passed to other methods.}
} 
  \value{
  A "ggplot" object.
  
} \seealso{
  \code{\link[gratia]{rootogram}}
} \references{
  Kleiber, C., Zeileis, A., (2016) Visualizing Count Data
Regressions Using Rootograms. \emph{Am. Stat.} \bold{70}, 296--303.
\doi{10.1080/00031305.2016.1173590}
}



\examples{ 
## data
data("rwm1984", package = 'COUNT')

zidw_out <- zidw_reg(docvis ~ outwork + female + age, ziformula = ~ age + female,
                     betaformula = ~ 1, data = rwm1984[1:100, ], beta = .1)
rootogram_zidw(zidw_out)
}

\keyword{file}


