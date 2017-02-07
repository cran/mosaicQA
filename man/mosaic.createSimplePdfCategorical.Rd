\name{mosaic.createSimplePdfCategorical}
\alias{mosaic.createSimplePdfCategorical}

\title{
createSimplePdfCategorical
}
\description{
Create simple PDF-file for categorical data 
}
\usage{
mosaic.createSimplePdfCategorical(inputfile, outputfolder)
}

\arguments{
  \item{inputfile}{
path to input csv-file
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

# specify the import file with categorial data
# first row has to contain variable names without special characters
categorial_datafile='c:/mosaic/cat_single_var_en.csv'

# specify output folder
outputFolder='c:/mosaic/outputs/'

# set treshold to detect missings, default is 99900 (adjust this line to change this global value, 
# but be careful)
mosaic.setGlobalMissingTreshold(99900)

# set description of var
mosaic.setGlobalCodelist(c('1=yes','2=no','99996=not specified','99997=not acquired'))

# create simple pdf file foreach variable column in categorial data file, uncomment to start 
# report-generation
# mosaic.createSimplePdfCategorical(categorial_datafile,outputFolder)
}

\keyword{ categorical }
\keyword{ pdf }
\keyword{ generate }
