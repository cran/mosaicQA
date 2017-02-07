\name{mosaicQA}
\alias{MOSAIC Basic Quality Assurance}
\docType{package}
\title{
\packageTitle{mosaicQA}
}
\description{
\packageDescription{mosaicQA}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{mosaicQA}
\packageIndices{mosaicQA}

The aim of the mosaicQA R-Package is to provide a basic assessment
of data quality and to generate a set of informative graphs. Especially, there should be no demand
for the potential researcher to master R. This R-package enables researchers to generate reports 
for various kinds of metric and categorical data. Additionally, general reports for multivariate 
input data and, if needed, detailed results for single-variable data can be produced.

CSV-files as well as dataframes can be used as input format to create a report. The results are instantly saved in an 
automatically generated PDF-file. For each study variable within the data input file a separate
PDF-file with standard or, if applicable, customized plots and tables is produced. These standard
reports enable the user to monitor and report the data integrity and completeness. However, for 
more specific reports the knowledge of metadata is necessary, including definition of units, 
variables, descriptions, code lists and categories of qualified missings.


Version 1.2
-----------
ADDED 	Support for metric and categorical dataframes
BUGFIX	Aborted report generation in case of non-existent missings in datacolumn

}

\author{
\packageAuthor{mosaicQA}

Maintainer: \packageMaintainer{mosaicQA}
}

\keyword{ package }
\keyword{ mosaicQA }
\keyword{ data assurance }
\keyword{ data quality }

\seealso{
mosaic-greifswald.de
}
\examples{

## Example 1: Generate pdf with graphs for a single metric data column, e.g. data of body height

# load mosaicQA package
library('mosaicQA')

# specify the csv import file with metric data, use one column per variable
metric_datafile='c:/mosaic/metric_single_var.csv'

#specify output folder
outputFolder='c:/mosaic/outputs/'

#set missing threshold, optional, default is 99900
mosaic.setGlobalMissingTreshold(99900)

#set variable unit, optional
mosaic.setGlobalUnit('(cm)')

#set variable description, optional, if not uses the name of the variable is displayed in 
#table heading
mosaic.setGlobalDescription('Height')

#create PDF-report, 
#uncomment to start report-generation
#mosaic.createSimplePdfmetric(metric_datafile, outputFolder)



## Example 2: Generate pdf with graphs for a single categorical data column

# load mosaicQA package
library('mosaicQA')

# specify the import file with Categorical data
# first row has to contain variable names without special characters
Categorical_datafile='c:/mosaic/cat_single_var_en.csv'

#specify output folder
outputFolder='c:/mosaic/outputs/'

#set treshold to detect missings, default is 99900 (adjust this line to change this global value,
#but be careful)
mosaic.setGlobalMissingTreshold(99900)

#set description of var
mosaic.setGlobalCodelist(c('1=yes','2=no','99996=not specified','99997=not acquired'))

# create simple pdf file foreach variable column in Categorical data file, 
# uncomment to start report-generation
# mosaic.createSimplePdfCategorical(Categorical_datafile,outputFolder)




## Example 3: Generate pdf with graphs for a multiple metric data columns, generates one pdf for 
# each column using the variable name for table headings

# load mosaicQA package
library('mosaicQA')

# specify the import file with metric data
# use one column per variable, first row should contain variable name, following rows should 
# contain data, csv Files with multiple rows are supported, decimal values should be formated 
# for example : 25.4
metric_datafile='c:/mosaic/metric_multi_var.csv'

#specify output folder
outputFolder="c:/mosaic/outputs/"

# set treshold to detect missings, default is 99900 (adjust this line to change this global value
# but be careful)
mosaic.setGlobalMissingTreshold(99900)

# create PDF-Files for vars, 
# uncomment to start report-generation
#mosaic.createSimplePdfmetric(metric_datafile, outputFolder)




## Example 4: Generate pdf with graphs for a multiple metric dataframe, generates one pdf for 
# each column using the variable name for table headings

# load mosaicQA package
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
