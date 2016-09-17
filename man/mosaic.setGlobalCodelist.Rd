\name{mosaic.setGlobalCodelist}
\alias{mosaic.setGlobalCodelist}
\title{
setGlobalCodelist
}
\description{
set and parse a global code list for categorical data to be used in categorical plot descriptions
}
\usage{
mosaic.setGlobalCodelist(coding)
}

\arguments{
  \item{coding}{
list of code and value pairs, see example for details
}
}



\author{
The MOSAIC Project, Martin Bialke
}

\note{
Function call type: user
}
\examples{
mosaic.setGlobalCodelist(c('1=yes','2=no', '99996=no information'))
}
