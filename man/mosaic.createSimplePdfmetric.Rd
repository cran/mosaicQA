\name{mosaic.createSimplePdfmetric}
\alias{mosaic.createSimplePdfmetric}
\title{
createSimplePdfmetric
}
\description{
Create simple PDF-file for metric data
}
\usage{
mosaic.createSimplePdfmetric(inputfile, outputfolder)
}

\arguments{
  \item{inputfile}{
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

# specify the csv import file with metric data, use one column per variable
metric_datafile='c:/mosaic/metric_single_var.csv'

#specify output folder
outputFolder='c:/mosaic/output/'

#set missing threshold, optional, default is 99900
mosaic.setGlobalMissingTreshold(99900)

#set variable unit, optional
mosaic.setGlobalUnit('(cm)')

#set variable description, optional
mosaic.setGlobalDescription('Height')

#create PDF-report, uncomment to start report-generation
#mosaic.createSimplePdfmetric(metric_datafile, outputFolder)
}
\keyword{ pdf }
\keyword{ generate }
\keyword{ metric }
