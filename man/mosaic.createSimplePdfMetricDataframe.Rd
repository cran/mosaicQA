\name{mosaic.createSimplePdfMetricDataframe}
\alias{mosaic.createSimplePdfMetricDataframe}
\title{
createSimplePdfMetricDataframe
}
\description{
Create simple PDF-file for metric data
}
\usage{
mosaic.createSimplePdfMetricDataframe(df, outputfolder)
}

\arguments{
  \item{df}{
path to input csv file
}
  \item{outputfolder}{
path to output folder 
}
}
\author{
The MOSAIC Project, Martin Bialke
}

\note{
Function call type: user
}


\examples{
# load MOSAIC package
library('mosaicQA')

# specify the metric dataframe with 1-n columns, here sample data is generated
metric_data=data.frame(matrix(rnorm(20), nrow=10))

#specify output folder
outputFolder="c:/mosaic/outputs/"

# set treshold to detect missings, default is 99900 (adjust this line to change this global value
# but be careful)
mosaic.setGlobalMissingTreshold(99900)

# create PDF-Files for vars, 
# uncomment to start report-generation
#mosaic.createSimplePdfMetricDataframe(metric_data, outputFolder)
}
\keyword{ pdf }
\keyword{ generate }
\keyword{ metric }
