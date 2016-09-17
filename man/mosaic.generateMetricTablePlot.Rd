\name{mosaic.generateMetricTablePlot}
\alias{mosaic.generateMetricTablePlot}

\title{
generateMetricTablePlot
}
\description{
Generate missing-ratio table for metric data (data, num of columns, column index, varname) 
}
\usage{
mosaic.generateMetricTablePlot(data, num_of_columns, index, varname)
}

\arguments{
  \item{data}{
 preprocessed data frame including 'valid value markers'
}
  \item{num_of_columns}{
absolute number of to be processed data columns
}
  \item{index}{
current column to be processed
}
  \item{varname}{
current name of variable to be used in table heading
}
}

\author{
The MOSAIC Project, Martin Bialke
}

\note{
Function call type: internal
}
