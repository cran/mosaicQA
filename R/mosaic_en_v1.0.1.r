
#   MOSAIC R Functions for BASIC Data Quality Assurance
#   Copyright (C) 2015  The MOSAIC-Project
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Affero General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Affero General Public License for more details.

#   You should have received a copy of the GNU Affero General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

#   For more information visit mosaic-greifswald.de or mail to mosaic@uni-greifswald.de

#   Author(s):  Martin Bialke, Thea Schwaneberg, Rene Walk
#   Date:       09. Sept. 2016
#   Version:    1.0.1
#

#internal environment
mosaicQA.env <- new.env()

# global labels
assign('label_normalverteilung', 'Distribution', envir=mosaicQA.env)
assign('labelPercentage', 'Percentage', envir=mosaicQA.env)
assign('labelCounts', 'Frequency', envir=mosaicQA.env)
assign('label_unit', '', envir=mosaicQA.env)
assign('label_description', '', envir=mosaicQA.env)
assign('label_boxplot', 'Box-Whisker-Plot', envir=mosaicQA.env)
assign('label_qnormplot', 'Q-Q-Plot', envir=mosaicQA.env)
assign('footnoteString', 'Created with mosaicQA R-package', envir=mosaicQA.env)

#codelist
assign('codelist', NULL, envir=mosaicQA.env)


#define interval of missings, every values above will be replaced with NA
assign('qualifiedMissingsTreshold', 99900, envir=mosaicQA.env)

#specify output file prefix
assign('outputPrefix', 'report', envir=mosaicQA.env)

###### About this Function-Library

mosaic.info = function() {
  paste("mosaicQA R package for BASIC Data Quality Assurance, More information from mosaic-greifswald.de",sep="")
}

###### Set Global Treshold for Missings

mosaic.setGlobalMissingTreshold = function (value) {
  assign("qualifiedMissingsTreshold", value, envir=mosaicQA.env)
}



###### Set Global Unit Label to used in graphs, e.g. "(cm)"

mosaic.setGlobalUnit = function (value) {
  assign("label_unit", value, envir=mosaicQA.env)
}



###### Set Global Description for  variable data, e.g. waist circumference

mosaic.setGlobalDescription = function (value) {
  assign("label_description", value, envir=mosaicQA.env)
}



###### Load CSV Data from file

mosaic.loadCsvData = function(filename){
  #load csv with n columns to data-variable
  as.data.frame(read.csv(filename,header=T,na.strings=NA,blank.lines.skip=T))

}

###### count occurence of searchvalue in datacolumn
mosaic.countValue=function(searchvalue, data_column){
  #read data as vector
  column_vector = data_column[,1] #col1 as vector
  summary_table<-table(column_vector)
  return(summary_table[names(summary_table)==searchvalue])
}

###### Preprocess metric data to allow missing-ratio table

mosaic.preProcessMetricData = function(data){
  num_of_columns=dim(data)[2]
  num_of_rows=dim(data)[1]

  # add column for each existing column filled with 1 (assume every value is a valid value)
  for (column_index in 1:num_of_columns){
    data[,num_of_columns+column_index]=rep(1,num_of_rows)
  }

  # process value columns, if value is a missing set "valid value marker" in previously added column to 0
  for (column_index in 1:num_of_columns){
    data[,num_of_columns+column_index][data[,column_index]>mosaicQA.env$qualifiedMissingsTreshold]=0
  }

  return (data)
}



###### identify unique values in data column, get abs, perc and cum stats
mosaic.preProcessCategorialData=function(data){
  num_of_columns=dim(data)[2]
  num_of_rows=dim(data)[1]

  #find unique values
  unique_values = unique(data[1])
  count_unique_values = dim(unique_values)[1]

  #sort unique values
  unique_values[,1]=sort(unique_values[,1])



  #matrix zum speichern der haeufigkeiten
  results=matrix(NA,count_unique_values,4)#abs,perc,abs_kum,perz_kum

  percent_cum=0
  abs_cum=0

  #count  unique values in data
  for (index in 1:count_unique_values){
    searchvalue=unique_values[index,1]
    #absolut
    results[index,1]=mosaic.countValue(searchvalue,data)
    #percent
    results[index,2]=round(results[index,1]*100/num_of_rows,2)
    #abs_cum
    abs_cum=abs_cum+results[index,1]
    results[index,3]=abs_cum
    #per_cum
    percent_cum=round(percent_cum+results[index,2],2)
    results[index,4]=percent_cum
  }

  sum_data = as.data.frame(cbind(results))

  #spalte benennen
  colnames(sum_data)[1]="Absolut"
  colnames(sum_data)[2]="%"
  colnames(sum_data)[3]="Abs. cum."
  colnames(sum_data)[4]="% cum."

  #reihen benennen, alles ueber treshold als missing benennen
  for (index in 1:count_unique_values){
    searchvalue=unique_values[index,1]

    #check if codelist exists
    if(!is.null(mosaicQA.env$codelist)){

      row.names(sum_data)[index]<-paste(mosaicQA.env$codelist[row.names(mosaicQA.env$codelist)[which(row.names(mosaicQA.env$codelist) == searchvalue)],1])

    } else {
      # use generic values as names
      if(searchvalue>mosaicQA.env$qualifiedMissingsTreshold){
        row.names(sum_data)[index]<-paste("Missing",searchvalue,sep=" ")
      }else{
        row.names(sum_data)[index]<-paste("Value",searchvalue,sep=" ")
      }
    }
  }

  #return summary
  return(sum_data)
}


###### generate missing-ratio table for metric data

mosaic.generateMetricTablePlot = function(data, num_of_columns,index, varname){

  # process previuously added columns containing "valid value markers"
  # calculate ratio of vaild values and misssings



  t=table(data[,num_of_columns+index])
  p=round(prop.table(table(data[,num_of_columns+index]))*100,2)
  t_kum=cumsum(t)
  p_kum=round(cumsum(p),2)
  tabelle1=cbind(t,p,t_kum,p_kum)

  # set header for temporary result table

  colnames(tabelle1)[1]<-paste(mosaicQA.env$labelCounts)
  colnames(tabelle1)[2]<-paste("Percentage")
  colnames(tabelle1)[3]<-paste("Frequency (cum.)")
  colnames(tabelle1)[4]<-paste("Percentage (cum.)")
  rownames(tabelle1)[1]<-paste("Missings")
  rownames(tabelle1)[2]<-paste("Valid Values")

  # plot result table
  textplot( tabelle1, valign="top"  )

  # set table heading
  title(main=paste('Valid Values and Missings',varname,sep=' '))
}



###### generate graphs for metric data

mosaic.generateMetricPlots = function(data_snippet, var_name){

  anzahl=dim(data_snippet)[1]

  #clean all qualified missings from column dataset and replace with NA
  tmpDataset=data_snippet[[var_name]]
  tmpDataset[tmpDataset>mosaicQA.env$qualifiedMissingsTreshold]=NA
  tmp_table=table(tmpDataset,useNA='always')

  #describe(tmpDataset)

  # calculate stats
  mue=round(mean(tmpDataset,na.rm=T),2)
  sigma=round(sd(tmpDataset,na.rm=T),2)
  max=max(tmpDataset,na.rm=T)
  min=min(tmpDataset,na.rm=T)
  quantile=quantile(tmpDataset,na.rm=T)

  #add first graph
  hist(tmpDataset,breaks=70,freq=F,col='cornflowerblue',xlab=mosaicQA.env$label_unit,ylab=mosaicQA.env$labelPercentage,main=paste(mosaicQA.env$label_normalverteilung,var_name,sep=" "))
  lines(density(tmpDataset,na.rm=T))

  #add legend for first graph  for stats
  legend("topright",
         legend = c(
                    paste("N",anzahl,sep=": "),
                    paste("Sigma",sigma,sep=": "),
                    paste("Mean Value",mue,sep=": "),
                    paste("Minimum",min,sep=": "),
                    paste("Maximum",max,sep=": "),
                    paste("Lower Quartile",quantile[1],sep=": "),
                    paste("Upper Quartile",quantile[4],sep=": ")
                    ),
                    cex=0.8)

  #add second graph
  boxplot(tmpDataset,col='cornflowerblue',xlab='',ylab=mosaicQA.env$label_unit,main=paste(mosaicQA.env$label_boxplot,var_name,sep=" "))

  #add third graph
  qqnorm(tmpDataset,col='cornflowerblue',xlab='Normal Theoretical Quantiles',ylab=mosaicQA.env$label_unit,main=paste(mosaicQA.env$label_qnormplot,var_name,sep=" "))
  qqline(tmpDataset)
}



###### begin plotting, generate PDF-File

mosaic.beginPlot= function(varname, outputfolder){
  #create outputpdf file name
  filename<-paste(mosaicQA.env$outputPrefix,"_",varname,"_",mosaic.getTimestamp(),'.pdf',sep="")

  #PDF with DIN A4 format, horizontal

  pdf(paste(outputfolder,filename,sep=""),width = 11.7, height = 8.3)
  #1/2 inch outer margins
  par(omi = rep(.5, 4))

  #place two graphs in two  row  s
  par(mfrow=c(2,2))

  #use varname or unit as label
  if(grepl(mosaicQA.env$label_unit, "") | grepl(mosaicQA.env$label_description, "")){
    #no unit specified, use varname
    assign("label_unit", varname, envir=mosaicQA.env)
  } else {
    #unit given, combine description and unit as label
    assign("label_unit", paste(mosaicQA.env$label_description,mosaicQA.env$label_unit,sep=" "), envir=mosaicQA.env)
  }
}




#### Add a Footnote

mosaic.addFootnote = function()

{
  footnote=c(paste(mosaicQA.env$footnoteString,", Date: ",format(Sys.time(), "%Y.%m.%d %H:%M:%S"),sep=""))

  size=.7
  color=grey(.5)

  pushViewport(viewport())
  grid.text(label=footnote,
            x = unit(1,"npc") - unit(2, "mm"),
            y= unit(2, "mm"),
            just=c("right", "bottom"),
            gp=gpar(cex= size, col=color))
  popViewport()
}


###### finish plotting, close PDF File

mosaic.finishPlot= function(){

  mosaic.addFootnote()

  dev.off()
}




###### set and parse a global codelist for categorial data

mosaic.setGlobalCodelist = function(coding){
  mosaicQA.env$codelist=matrix(coding)
  len=dim(mosaicQA.env$codelist)[1]

  #parse codelist
  if(len>0){
    parsedcodelist = as.data.frame(cbind(matrix(NA,len,1)))

    for (i in 1:len) {
    #use = as separator
    index_of_char = grep("=", strsplit(mosaicQA.env$codelist[i], "")[[1]])
    meaning = substr(mosaicQA.env$codelist[i],index_of_char+1,nchar(mosaicQA.env$codelist[i]))
    code = substr(mosaicQA.env$codelist[i],1,index_of_char-1)
    #names rows = code
    row.names(parsedcodelist)[i]<-paste(code)
    #cells = meaning
    parsedcodelist[i,1]<-paste(meaning)
    }

    assign("codelist", parsedcodelist, envir=mosaicQA.env)

  }
}


###### create plots for categorial data

mosaic.generateCategorialPlot= function(dataframe, varname){

  barLabel = paste(mosaicQA.env$label_normalverteilung, mosaicQA.env$labelCounts, varname, sep=" ")

  # plot result table
  textplot(dataframe, valign="top"  )

  # set table heading
  title(main=barLabel)

  #names
  list_of_names=c(row.names(dataframe))

  #boxplot

  barplot(dataframe[,1],main=barLabel,horiz=F,  ylab=mosaicQA.env$labelCounts,xlab='',names.arg=list_of_names,col='cornflowerblue',las=2,cex.axis=1, cex.names=0.8)
}


###### create simple pdf file for categorial data
###### (combines the above functionality)
mosaic.createSimplePdfcategorial = function(inputfile, outputfolder){

  # load csv with n columns to data-variable
  data = mosaic.loadCsvData(inputfile)
  num_of_columns=dim(data)[2]


  # create PDF Files  with missings statistics, histogram, boxplot and qqnorm  for each variable /column
  for (column_index in 1:num_of_columns){
    #use selected  column header as varname
    varname=names(data)[column_index]

    #calc abs, % and cums for selected column
    data_preprocessed = mosaic.preProcessCategorialData(data[column_index])

    #create PDF file with varname suffix

    mosaic.beginPlot(varname,outputfolder)

    #divide pdf page, 2 rows, 1 column
    par(mfrow=c(1,2))

    #generate missings-table based on results of preprocessing
    mosaic.generateCategorialPlot(data_preprocessed,varname)

    mosaic.finishPlot()
  }
}


###### create simple pdf file for metric data
###### (combines the above functionality)
mosaic.createSimplePdfmetric = function(inputfile, outputfolder){

    # load csv with n columns to data-variable
    metric_data = mosaic.loadCsvData(inputfile)

    num_of_columns=dim(metric_data)[2]

    # preprocess data to calculate ratio of valid values and missings
    data_preprocessed=mosaic.preProcessMetricData(metric_data)

    # create PDF Files  with missings statistics, histogram, boxplot and qqnorm  for each variable /column
    for (column_index in 1:num_of_columns){
      #use column header as varname
      varname=colnames(data_preprocessed)[column_index]

      #create PDF file with varname suffix
      mosaic.beginPlot(varname,outputfolder)

      #generate missings-table based on results of preprocessing
      mosaic.generateMetricTablePlot(data_preprocessed,num_of_columns,column_index,varname)

      # generate graphics
      mosaic.generateMetricPlots(data_preprocessed, varname)
      mosaic.finishPlot()
    }
}

##### get formatted timestamp
mosaic.getTimestamp = function(){
  return (format(Sys.time(), "%Y_%m_%d_%H%M%S"))
}
