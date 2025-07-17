#multiXdateR - Beta 2 - **paper version 17/07/2025**
#this sets the working directory to the source file location - nice
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#=================================****** call the libraries needed - ensure they are installed
library(readxl); library(zoo); library(tidyverse); library(ggpubr); library(dplR)
library(qpcR); library(dplyr); library(ggpubr); library(scales)
#=================================******
#=================================****** User Settings - lines 8-37 ******
runtitle <- "Test" #if you want a header for the graphics output
#=================================****** Choice of Graphical Output ******
CorrelationPlots <- "yes"   # "no" "yes" = sliding correlation plots
TvaluePlots <- "yes"        # "no" "yes" = Sliding T-value plots, T-value density histograms and P-value plots
#=================================
reftruncate <-  "no"   # "no" "yes". Truncation option for reference chronology. DANGER THERE BE DRAGONS!!!
refstart <- "xxx" # start year of reference chronologies if truncated
refend <- "xxx"   # end year of reference chronologies if truncated
#=================================
#read in Tucson file of reference chronologies - they must cover the same period
referenceRW <- read.rwl ("ExampleData/LGLrwTRUNC", header = FALSE) # RW reference
referenceEWBI <- read.rwl ("ExampleData/LGLewbiTRUNC", header = FALSE) # EWBI reference
referenceLWBI <- read.rwl ("ExampleData/LGLlwbiTRUNC", header = FALSE) # LWBI reference

#read in Tucson file of undated chronologies/series
undatedRW <- read.rwl ("ExampleData/LCLrwTRUNC", header = FALSE) # RW undated data
undatedEWBI <- read.rwl ("ExampleData/LCLewbiTRUNC", header = FALSE) # EWBI undated data
undatedLWBI <- read.rwl ("ExampleData/LCLlwbiTRUNC", header = FALSE) # LWBI undated data
#===============================
#minimum overlap between undated and reference chrons
minoverlap <- 30 #cannot be longer than the length of your undated series
#==================================
#detrending method is a fitted fixed spline removed via division using a biweight robust mean with prewhitening
#user can change spline window size, post detrend unfiltered or 1st difference transform and correlation type
#degrees of freedom are later adjusted related to the AC1 of the resultant chronologies
splinewindow <- 31
datatransform <- "pw" #residual chronology - options of "pw" = pre-whitened or "fd" = 1st differenced
correlationtype <- "spearman" # use pearson or spearman
#=================================******
#=================================******
#detrending the reference data
#RW
RWref.rwi <- detrend(rwl = referenceRW, method = c("Spline"), nyrs = splinewindow, f = 0.5, 
                    pos.slope = TRUE, difference = FALSE) 

#chronology build without variance stabilisation
RWrefgrow.crn <- chron(x = RWref.rwi, prefix = "ref", biweight = TRUE, prewhiten = TRUE)
RWrefyears <- row.names(RWrefgrow.crn)

#EWBI
EWBIref.EWBIi <- detrend(rwl = referenceEWBI, method = c("Spline"), nyrs = splinewindow, f = 0.5, 
                     pos.slope = TRUE, difference = FALSE) 

#chronology build without variance stabilisation
EWBIrefgrow.crn <- chron(x = EWBIref.EWBIi, prefix = "ref", biweight = TRUE, prewhiten = TRUE)
EWBIrefyears <- row.names(EWBIrefgrow.crn)

#LWBI
LWBIref.LWBIi <- detrend(rwl = referenceLWBI, method = c("Spline"), nyrs = splinewindow, f = 0.5, 
                         pos.slope = TRUE, difference = FALSE) 

#chronology build without variance stabilisation
LWBIrefgrow.crn <- chron(x = LWBIref.LWBIi, prefix = "ref", biweight = TRUE, prewhiten = TRUE)
LWBIrefyears <- row.names(LWBIrefgrow.crn)

#==================================
#detrending the undated data
#RW
RWund.rwi <- detrend(rwl = undatedRW, method = c("Spline"), nyrs = splinewindow, f = 0.5, 
                   pos.slope = TRUE, difference = FALSE)

#to derive length of the full undated chronology
ORIGundgrow.crn <- chron(x = RWund.rwi, prefix = "ref", biweight = TRUE, prewhiten = FALSE)
ORIGundyears <- row.names(ORIGundgrow.crn)

#chronology build without variance stabilisation
RWundgrow.crn <- chron(x = RWund.rwi, prefix = "ref", biweight = TRUE, prewhiten = TRUE)
RWundyears <- row.names(RWundgrow.crn)

#EWBI
EWBIund.EWBIi <- detrend(rwl = undatedEWBI, method = c("Spline"), nyrs = splinewindow, f = 0.5, 
                     pos.slope = TRUE, difference = FALSE) 

#chronology build without variance stabilisation
EWBIundgrow.crn <- chron(x = EWBIund.EWBIi, prefix = "ref", biweight = TRUE, prewhiten = TRUE)
EWBIundyears <- row.names(EWBIundgrow.crn)

#LWBI
LWBIund.LWBIi <- detrend(rwl = undatedLWBI, method = c("Spline"), nyrs = splinewindow, f = 0.5, 
                         pos.slope = TRUE, difference = FALSE) 

#chronology build without variance stabilisation
LWBIundgrow.crn <- chron(x = LWBIund.LWBIi, prefix = "ref", biweight = TRUE, prewhiten = TRUE)
LWBIundyears <- row.names(LWBIundgrow.crn)

#==================================
#then create dataframes for analysis based on data transform - pw or df
#RW
if (datatransform == "pw") {RWreference <- as.data.frame(cbind(RWrefyears, RWrefgrow.crn$res))} #unfiltered - use DF as is
if (datatransform == "pw") {RWundated <- as.data.frame(cbind(RWundyears, RWundgrow.crn$res))}
if (datatransform == "fd") {RWref1st <- diff(RWrefgrow.crn$res)
RWreference <- qpcR:::cbind.na(RWrefyears, RWref1st)
RWreference <- as.data.frame(RWreference)
}
if (datatransform == "fd") {RWund1st <- diff(RWundgrow.crn$res)
RWundated <- qpcR:::cbind.na(RWundyears, RWund1st)
RWundated <- as.data.frame(RWundated)
} 
colnames(RWreference) [1] <- "Year"
colnames(RWreference) [2] <- "RWReferenceChron"
colnames(RWundated) [1] <- "Year"
colnames(RWundated) [2] <- "RWUndatedChron"
RWreference <- na.omit(RWreference) #removes the NA values from the Residual chronology
if (reftruncate == "yes") {RWreference <- subset(RWreference, Year >= refstart & Year <= refend)} #truncate reference chron
RWundated <- na.omit(RWundated) #removes the NA values from the Residual chronology

#EWBI
if (datatransform == "pw") {EWBIreference <- as.data.frame(cbind(EWBIrefyears, EWBIrefgrow.crn$res))} #unfiltered - use DF as is
if (datatransform == "pw") {EWBIundated <- as.data.frame(cbind(EWBIundyears, EWBIundgrow.crn$res))}
if (datatransform == "fd") {EWBIref1st <- diff(EWBIrefgrow.crn$res)
EWBIreference <- qpcR:::cbind.na(EWBIrefyears, EWBIref1st)
EWBIreference <- as.data.frame(EWBIreference)
}
if (datatransform == "fd") {EWBIund1st <- diff(EWBIundgrow.crn$res)
EWBIundated <- qpcR:::cbind.na(EWBIundyears, EWBIund1st)
EWBIundated <- as.data.frame(EWBIundated)
} 
colnames(EWBIreference) [1] <- "Year"
colnames(EWBIreference) [2] <- "EWBIReferenceChron"
colnames(EWBIundated) [1] <- "Year"
colnames(EWBIundated) [2] <- "EWBIUndatedChron"
EWBIreference <- na.omit(EWBIreference) #removes the NA values from the Residual chronology
if (reftruncate == "yes") {EWBIreference <- subset(EWBIreference, Year >= refstart & Year <= refend)} #truncate reference chron
EWBIundated <- na.omit(EWBIundated) #removes the NA values from the Residual chronology

#LWBI
if (datatransform == "pw") {LWBIreference <- as.data.frame(cbind(LWBIrefyears, LWBIrefgrow.crn$res))} #unfiltered - use DF as is
if (datatransform == "pw") {LWBIundated <- as.data.frame(cbind(LWBIundyears, LWBIundgrow.crn$res))}
if (datatransform == "fd") {LWBIref1st <- diff(LWBIrefgrow.crn$res)
LWBIreference <- qpcR:::cbind.na(LWBIrefyears, LWBIref1st)
LWBIreference <- as.data.frame(LWBIreference)
}
if (datatransform == "fd") {LWBIund1st <- diff(LWBIundgrow.crn$res)
LWBIundated <- qpcR:::cbind.na(LWBIundyears, LWBIund1st)
LWBIundated <- as.data.frame(LWBIundated)
} 
colnames(LWBIreference) [1] <- "Year"
colnames(LWBIreference) [2] <- "LWBIReferenceChron"
colnames(LWBIundated) [1] <- "Year"
colnames(LWBIundated) [2] <- "LWBIUndatedChron"
LWBIreference <- na.omit(LWBIreference) #removes the NA values from the Residual chronology
if (reftruncate == "yes") {LWBIreference <- subset(LWBIreference, Year >= refstart & Year <= refend)} #truncate reference chron
LWBIundated <- na.omit(LWBIundated) #removes the NA values from the Residual chronology

#==================================
#series may not have a common first year after PW - need to truncate the series
#reference chrons

#original first year
RWfirstyearREF <- first(as.numeric(RWreference$Year))
EWBIfirstyearREF <- first(as.numeric(EWBIreference$Year))
LWBIfirstyearREF <- first(as.numeric(LWBIreference$Year))
combofirstyearREF <- as.numeric(max(RWfirstyearREF, EWBIfirstyearREF, LWBIfirstyearREF))
lastyearREF <- last(as.numeric(RWreference$Year))

#ensure 1st common year
RWreference$Year <- as.numeric(RWreference$Year)
RWreference <- RWreference[which(RWreference$Year >= combofirstyearREF),]
EWBIreference$Year <- as.numeric(EWBIreference$Year)
EWBIreference <- EWBIreference[which(EWBIreference$Year >= combofirstyearREF),]
LWBIreference$Year <- as.numeric(LWBIreference$Year)
LWBIreference <- LWBIreference[which(LWBIreference$Year >= combofirstyearREF),]

#common first year after correction
#original first year
RWfirstyearREF <- first(as.numeric(RWreference$Year))
EWBIfirstyearREF <- first(as.numeric(EWBIreference$Year))
LWBIfirstyearREF <- first(as.numeric(LWBIreference$Year))
combofirstyearREF <- as.numeric(max(RWfirstyearREF, EWBIfirstyearREF, LWBIfirstyearREF))

#undated chrons
RWfirstyearUND <- first(as.numeric(RWundated$Year))
EWBIfirstyearUND <- first(as.numeric(EWBIundated$Year))
LWBIfirstyearUND <- first(as.numeric(LWBIundated$Year))
combofirstyearUND <- as.numeric(max(RWfirstyearUND, EWBIfirstyearUND, LWBIfirstyearUND))
lastyearUND <- last(as.numeric(RWundated$Year))

#ensure 1st common year
RWundated$Year <- as.numeric(RWundated$Year)
RWundated <- RWundated[which(RWundated$Year >= combofirstyearUND),]
EWBIundated$Year <- as.numeric(EWBIundated$Year)
EWBIundated <- EWBIundated[which(EWBIundated$Year >= combofirstyearUND),]
LWBIundated$Year <- as.numeric(LWBIundated$Year)
LWBIundated <- LWBIundated[which(LWBIundated$Year >= combofirstyearUND),]

#=================================
#metrics needed for calculations
#RW
RWFULLlength <- length(as.numeric(RWundated$RWUndatedChron))
RWwindow_size <- RWFULLlength
RWreferencefirstyear <- min(as.numeric(RWreference$Year))

if (minoverlap > RWFULLlength) {
  print("DUDE! - your minimum overlap is longer than your undated chronology")
  print("Code has been stopped - Please re-run with lower minimum overlap value!!!")
  stop()
}

#EWBI
EWBIFULLlength <- length(as.numeric(EWBIundated$EWBIUndatedChron))
EWBIwindow_size <- EWBIFULLlength
EWBIreferencefirstyear <- min(as.numeric(EWBIreference$Year))
if (datatransform == "fd") {EWBIreferencefirstyear <- EWBIreferencefirstyear + 1} #currently in as a fidge fix

if (minoverlap > EWBIFULLlength) {
  print("DUDE! - your minimum overlap is longer than your undated chronology")
  print("Code has been stopped - Please re-run with lower minimum overlap value!!!")
  stop()
}

#LWBI
LWBIFULLlength <- length(as.numeric(LWBIundated$LWBIUndatedChron))
LWBIwindow_size <- LWBIFULLlength
LWBIreferencefirstyear <- min(as.numeric(LWBIreference$Year))
if (datatransform == "fd") {LWBIreferencefirstyear <- LWBIreferencefirstyear + 1} #currently in as a fidge fix

if (minoverlap > LWBIFULLlength) {
  print("DUDE! - your minimum overlap is longer than your undated chronology")
  print("Code has been stopped - Please re-run with lower minimum overlap value!!!")
  stop()
}

#=================================
#add padding to reference dataframe to allow extension beyond the reference chron range using overlap
#RW
RWpadding <- RWFULLlength - minoverlap #no. of NA values to be added

RWreference <- data.frame(Value = c(rep(NA, RWpadding), RWreference$RWReferenceChron))
#this adds NA values to beginning of DF, but deletes the year column
colnames(RWreference) [1] <- "RWReferenceChron"

#this codes adds NA to end - but need years to be added back
RWreference <- RWreference %>%
  slice(rep(2, RWpadding)) %>%
  mutate(across(everything(), ~ NA)) %>%
  bind_rows(RWreference, .)

RWreflength <- length(RWreference$RWReferenceChron)
RWnewstart <- RWreferencefirstyear - RWpadding

#create a new column of specific annual values
RWnewdates <- seq(from = RWnewstart, to = RWnewstart + RWreflength, by = 1)
RWreference <- qpcR:::cbind.na(RWnewdates, RWreference)
colnames(RWreference) [1] <- "Year"
colnames(RWreference) [2] <- "RWReferenceChron"

#EWBI
EWBIpadding <- EWBIFULLlength - minoverlap #no. of NA values to be added

EWBIreference <- data.frame(Value = c(rep(NA, EWBIpadding), EWBIreference$EWBIReferenceChron))
#this adds NA values to beginning of DF, but deletes the year column
colnames(EWBIreference) [1] <- "EWBIReferenceChron"

#this codes adds NA to end - but need years to be added back
EWBIreference <- EWBIreference %>%
  slice(rep(2, EWBIpadding)) %>%
  mutate(across(everything(), ~ NA)) %>%
  bind_rows(EWBIreference, .)

EWBIreflength <- length(EWBIreference$EWBIReferenceChron)
EWBInewstart <- EWBIreferencefirstyear - EWBIpadding

#create a new column of specific annual values
EWBInewdates <- seq(from = EWBInewstart, to = EWBInewstart + EWBIreflength, by = 1)
EWBIreference <- qpcR:::cbind.na(EWBInewdates, EWBIreference)
colnames(EWBIreference) [1] <- "Year"
colnames(EWBIreference) [2] <- "EWBIReferenceChron"

#LWBI
LWBIpadding <- LWBIFULLlength - minoverlap #no. of NA values to be added

LWBIreference <- data.frame(Value = c(rep(NA, LWBIpadding), LWBIreference$LWBIReferenceChron))
#this adds NA values to beginning of DF, but deletes the year column
colnames(LWBIreference) [1] <- "LWBIReferenceChron"

#this codes adds NA to end - but need years to be added back
LWBIreference <- LWBIreference %>%
  slice(rep(2, LWBIpadding)) %>%
  mutate(across(everything(), ~ NA)) %>%
  bind_rows(LWBIreference, .)

LWBIreflength <- length(LWBIreference$LWBIReferenceChron)
LWBInewstart <- LWBIreferencefirstyear - LWBIpadding

#create a new column of specific annual values
LWBInewdates <- seq(from = LWBInewstart, to = LWBInewstart + LWBIreflength, by = 1)
LWBIreference <- qpcR:::cbind.na(LWBInewdates, LWBIreference)
colnames(LWBIreference) [1] <- "Year"
colnames(LWBIreference) [2] <- "LWBIReferenceChron"

#=================================
#sliding correlation
#Function to calculate correlation coefficient
calculate_correlation <- function(x, y) {
  cor(x, y, method = correlationtype, use = "pairwise.complete.obs")
  #pairwise.complete.obs ignores NA values but still calculates r value.
}

#RW
# Calculate correlation coefficients and N using rollapply
RWcorrelation_result <- rollapply(as.numeric(RWreference$RWReferenceChron), width = RWwindow_size,  
                                FUN = calculate_correlation, align = "right", by = 1, y = as.numeric(RWundated$RWUndatedChron))

#count number of correlations per window
RWcount_non_na <- function(x) sum(!is.na(x))
RWnon_na_counts <- rollapply(RWreference$RWReferenceChron, width = RWwindow_size, FUN = RWcount_non_na, align = "right")

#EWBI
# Calculate correlation coefficients and N using rollapply
EWBIcorrelation_result <- rollapply(as.numeric(EWBIreference$EWBIReferenceChron), width = EWBIwindow_size,  
                                  FUN = calculate_correlation, align = "right", by = 1, y = as.numeric(EWBIundated$EWBIUndatedChron))

#count number of correlations per window
EWBIcount_non_na <- function(x) sum(!is.na(x))
EWBInon_na_counts <- rollapply(EWBIreference$EWBIReferenceChron, width = EWBIwindow_size, FUN = EWBIcount_non_na, align = "right")

#LWBI
# Calculate correlation coefficients and N using rollapply
LWBIcorrelation_result <- rollapply(as.numeric(LWBIreference$LWBIReferenceChron), width = LWBIwindow_size,  
                                  FUN = calculate_correlation, align = "right", by = 1, y = as.numeric(LWBIundated$LWBIUndatedChron))

#count number of correlations per window
LWBIcount_non_na <- function(x) sum(!is.na(x))
LWBInon_na_counts <- rollapply(LWBIreference$LWBIReferenceChron, width = LWBIwindow_size, FUN = LWBIcount_non_na, align = "right")

#=================================
#create mean correlation of 3 parameter correlation results
correlation_results <- cbind(RWcorrelation_result, EWBIcorrelation_result, LWBIcorrelation_result)
MEANcorrelation_result <- rowMeans(correlation_results)
ALLcorrelation_results <- as.data.frame(cbind(RWcorrelation_result, EWBIcorrelation_result, LWBIcorrelation_result, MEANcorrelation_result))

RWcorrSTDEV <- sd(RWcorrelation_result) #STDEV for RW correlations
EWBIcorrSTDEV <- sd(EWBIcorrelation_result) #STDEV for EWBI correlations
LWBIcorrSTDEV <- sd(LWBIcorrelation_result) #STDEV for LWBI correlations
MEANcorrSTDEV <- sd(MEANcorrelation_result) #STDEV for ALL correlations

#highest correlation values
RWmaxcorr <- max(RWcorrelation_result)
EWBImaxcorr <- max(EWBIcorrelation_result)
LWBImaxcorr <- max(LWBIcorrelation_result)
COMBOmaxcorr <- max(MEANcorrelation_result)


#2nd highest correlation values
NEXTRWmaxcorr <- max(RWcorrelation_result[RWcorrelation_result != RWmaxcorr])
NEXTEWBImaxcorr <- max(EWBIcorrelation_result[EWBIcorrelation_result != EWBImaxcorr])
NEXTLWBImaxcorr <- max(LWBIcorrelation_result[LWBIcorrelation_result != LWBImaxcorr])
NEXTCOMBOmaxcorr <- max(MEANcorrelation_result[MEANcorrelation_result != COMBOmaxcorr])

#=========================
# Find the position of the maximum correlation coefficient in COMBO MEANcorrelation timeseries
# if significant, this should be the "correct" date of the final year if the undated series

#RW
RWmax_corr_position <- which (RWcorrelation_result == RWmaxcorr, arr.ind=TRUE)
RWsec_corr_position <- which (RWcorrelation_result == NEXTRWmaxcorr, arr.ind=TRUE)

#RW top correlation position
RWouteryear <- RWfirstyearREF + RWmax_corr_position + RWwindow_size - RWpadding - 2 #should depend on 1st year of reference chron
if (datatransform == "pw") {
  RWouteryear <- RWouteryear
} else {
  RWouteryear <- RWouteryear +1
}

if (datatransform == "pw") {
  RWreferencefirstyear <- RWreferencefirstyear
} else {
  RWreferencefirstyear <- RWreferencefirstyear +1
}

#RW 2nd correlation position
RWnextyear <- RWfirstyearREF + RWsec_corr_position + RWwindow_size - RWpadding - 2 #should depend on 1st year of reference chron
if (datatransform == "pw") {
  RWnextyear <- RWnextyear
} else {
  RWnextyear <- RWnextyear +1
}

if (datatransform == "pw") {
  RWreferencefirstyear <- RWreferencefirstyear
} else {
  RWreferencefirstyear <- RWreferencefirstyear +1
}

#EWBI
EWBImax_corr_position <- which (EWBIcorrelation_result == EWBImaxcorr, arr.ind=TRUE)
EWBIsec_corr_position <- which (EWBIcorrelation_result == NEXTEWBImaxcorr, arr.ind=TRUE)

#EWBI top correlation position
EWBIouteryear <- EWBIfirstyearREF + EWBImax_corr_position + EWBIwindow_size - EWBIpadding - 2 #should depend on 1st year of reference chron
if (datatransform == "pw") {
  EWBIouteryear <- EWBIouteryear
} else {
  EWBIouteryear <- EWBIouteryear +1
}

if (datatransform == "pw") {
  EWBIreferencefirstyear <- EWBIreferencefirstyear
} else {
  EWBIreferencefirstyear <- EWBIreferencefirstyear +1
}

#EWBI 2nd correlation position
EWBInextyear <- EWBIfirstyearREF + EWBIsec_corr_position + EWBIwindow_size - EWBIpadding - 2 #should depend on 1st year of reference chron
if (datatransform == "pw") {
  EWBInextyear <- EWBInextyear
} else {
  EWBInextyear <- EWBInextyear +1
}

if (datatransform == "pw") {
  EWBIreferencefirstyear <- EWBIreferencefirstyear
} else {
  EWBIreferencefirstyear <- EWBIreferencefirstyear +1
}

#LWBI
LWBImax_corr_position <- which (LWBIcorrelation_result == LWBImaxcorr, arr.ind=TRUE)
LWBIsec_corr_position <- which (LWBIcorrelation_result == NEXTLWBImaxcorr, arr.ind=TRUE)

#LWBI top correlation position
LWBIouteryear <- LWBIfirstyearREF + LWBImax_corr_position + LWBIwindow_size - LWBIpadding - 2 #should depend on 1st year of reference chron
if (datatransform == "pw") {
  LWBIouteryear <- LWBIouteryear
} else {
  LWBIouteryear <- LWBIouteryear +1
}

if (datatransform == "pw") {
  LWBIreferencefirstyear <- LWBIreferencefirstyear
} else {
  LWBIreferencefirstyear <- LWBIreferencefirstyear +1
}

#LWBI 2nd correlation position
LWBInextyear <- LWBIfirstyearREF + LWBIsec_corr_position + LWBIwindow_size - LWBIpadding - 2 #should depend on 1st year of reference chron
if (datatransform == "pw") {
  LWBInextyear <- LWBInextyear
} else {
  LWBInextyear <- LWBInextyear +1
}

if (datatransform == "pw") {
  LWBIreferencefirstyear <- LWBIreferencefirstyear
} else {
  LWBIreferencefirstyear <- LWBIreferencefirstyear +1
}

#COMBO
COMBOmax_corr_position <- which (MEANcorrelation_result == COMBOmaxcorr, arr.ind=TRUE)
COMBOsec_corr_position <- which (MEANcorrelation_result == NEXTCOMBOmaxcorr, arr.ind=TRUE)

#COMBO top correlation position
COMBOouteryear <- combofirstyearREF + COMBOmax_corr_position + RWwindow_size - RWpadding - 2 #should depend on 1st year of reference chron
if (datatransform == "pw") {
  COMBOouteryear <- COMBOouteryear
} else {
  COMBOouteryear <- COMBOouteryear +1
}

if (datatransform == "pw") {
  RWreferencefirstyear <- RWreferencefirstyear
} else {
  RWreferencefirstyear <- RWreferencefirstyear +1
}

#COMBO 2nd correlation position
COMBOnextyear <- combofirstyearREF + COMBOsec_corr_position + RWwindow_size - RWpadding - 2 #should depend on 1st year of reference chron
if (datatransform == "pw") {
  COMBOnextyear <- COMBOnextyear
} else {
  COMBOnextyear <- COMBOnextyear +1
}

if (datatransform == "pw") {
  RWreferencefirstyear <- RWreferencefirstyear
} else {
  RWreferencefirstyear <- RWreferencefirstyear +1
}

#if statement for correlation graphics output
if (CorrelationPlots == "yes") {

#========================the plots code
# Plot the time series of the correlation coefficients
#RW
RWcorrelationplot <- ggplot() +
  annotate(geom = "point", x = RWouteryear, y = RWmaxcorr, shape = 25, colour = "black", fill = "green", size = 4) +
  annotate(geom = "point", x = RWnextyear, y = NEXTRWmaxcorr, shape = 25, colour = "black", fill = "purple", size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0+4*RWcorrSTDEV, linetype = "dashed", color = "grey", linewidth = 0.3) +
  {if(datatransform == "pw") geom_line(aes(x = (1:length(RWcorrelation_result)) + (RWreferencefirstyear + (RWwindow_size - RWpadding - 2)), y = RWcorrelation_result), 
                                       color = "red", linewidth = 0.3) } +
  {if(datatransform == "fd") geom_line(aes(x = (1:length(RWcorrelation_result)) + (RWreferencefirstyear + (RWwindow_size - RWpadding - 5)), y = RWcorrelation_result), 
                                       color = "red", linewidth = 0.3) } +
  geom_vline(xintercept = COMBOouteryear, linetype = "dashed", color = "blue", linewidth = 0.5) +
  annotate("text", x=(min(as.numeric(RWreference$Year)) + RWFULLlength), y=4*RWcorrSTDEV, 
           label= "4 STDEVs", color="grey", vjust = 1.5, fontface = "bold", size = 2.5) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(title = "RW Sliding correlation values", x = "Calendar Years CE", y = "Correlation") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        plot.title=element_text(face = "bold"),
        plot.margin = margin(0,0.5,0,0, "cm"))

#EWBI
EWBIcorrelationplot <- ggplot() +
  annotate(geom = "point", x = EWBIouteryear, y = EWBImaxcorr, shape = 25, colour = "black", fill = "green", size = 4) +
  annotate(geom = "point", x = EWBInextyear, y = NEXTEWBImaxcorr, shape = 25, colour = "black", fill = "purple", size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0+4*EWBIcorrSTDEV, linetype = "dashed", color = "grey", linewidth = 0.3) +
  {if(datatransform == "pw") geom_line(aes(x = (1:length(EWBIcorrelation_result)) + (RWreferencefirstyear + (EWBIwindow_size - EWBIpadding - 2)), y = EWBIcorrelation_result), 
                                       color = "red", linewidth = 0.3) } +
  {if(datatransform == "fd") geom_line(aes(x = (1:length(EWBIcorrelation_result)) + (RWreferencefirstyear + (EWBIwindow_size - EWBIpadding - 5)), y = EWBIcorrelation_result), 
                                       color = "red", linewidth = 0.3) } +
  geom_vline(xintercept = COMBOouteryear, linetype = "dashed", color = "blue", linewidth = 0.5) +
  annotate("text", x=(min(as.numeric(EWBIreference$Year)) + EWBIFULLlength), y=4*EWBIcorrSTDEV, 
           label= "4 STDEVs", color="grey", vjust = 1.5, fontface = "bold", size = 2.5) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(title = "EWBI Sliding correlation values", x = "Calendar Years CE", y = "Correlation") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        plot.title=element_text(face = "bold"),
        plot.margin = margin(0,0.5,0,0, "cm"))

#LWBI
LWBIcorrelationplot <- ggplot() +
  annotate(geom = "point", x = LWBIouteryear, y = LWBImaxcorr, shape = 25, colour = "black", fill = "green", size = 4) +
  annotate(geom = "point", x = LWBInextyear, y = NEXTLWBImaxcorr, shape = 25, colour = "black", fill = "purple", size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0+4*LWBIcorrSTDEV, linetype = "dashed", color = "grey", linewidth = 0.3) +
  {if(datatransform == "pw") geom_line(aes(x = (1:length(LWBIcorrelation_result)) + (RWreferencefirstyear + (LWBIwindow_size - LWBIpadding - 2)), y = LWBIcorrelation_result), 
                                       color = "red", linewidth = 0.3) } +
  {if(datatransform == "fd") geom_line(aes(x = (1:length(LWBIcorrelation_result)) + (RWreferencefirstyear + (LWBIwindow_size - LWBIpadding - 5)), y = LWBIcorrelation_result), 
                                       color = "red", linewidth = 0.3) } +
  geom_vline(xintercept = COMBOouteryear, linetype = "dashed", color = "blue", linewidth = 0.5) +
  annotate("text", x=(min(as.numeric(LWBIreference$Year)) + LWBIFULLlength), y=4*LWBIcorrSTDEV, 
           label= "4 STDEVs", color="grey", vjust = 1.5, fontface = "bold", size = 2.5) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(title = "LWBI Sliding correlation values", x = "Calendar Years CE", y = "Correlation") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        plot.title=element_text(face = "bold"),
        plot.margin = margin(0,0.5,0,0, "cm"))

#MEAN combo correlations
MEANcorrelationplot <- ggplot() +
  annotate(geom = "point", x = COMBOouteryear, y = COMBOmaxcorr, shape = 25, colour = "black", fill = "green", size = 4) +
  annotate(geom = "point", x = COMBOnextyear, y = NEXTCOMBOmaxcorr, shape = 25, colour = "black", fill = "purple", size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0+4*MEANcorrSTDEV, linetype = "dashed", color = "grey", linewidth = 0.3) +
  {if(datatransform == "pw") geom_line(aes(x = (1:length(MEANcorrelation_result)) + (RWreferencefirstyear + (LWBIwindow_size - LWBIpadding - 2)), y = MEANcorrelation_result), 
                                       color = "red", linewidth = 0.3) } +
  {if(datatransform == "fd") geom_line(aes(x = (1:length(MEANcorrelation_result)) + (RWreferencefirstyear + (LWBIwindow_size - LWBIpadding - 5)), y = MEANcorrelation_result), 
                                       color = "red", linewidth = 0.3) } +
  geom_vline(xintercept = COMBOouteryear, linetype = "dashed", color = "blue", linewidth = 0.5) +
  annotate("text", x=(min(as.numeric(LWBIreference$Year)) + LWBIFULLlength), y=4*MEANcorrSTDEV, 
           label= "4 STDEVs", color="grey", vjust = 1.5, fontface = "bold", size = 2.5) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(title = "All Parameter COMBO Sliding correlation values", x = "Calendar Years CE", y = "Correlation") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        plot.title=element_text(face = "bold"),
        plot.margin = margin(0,0.5,0,0, "cm"))

Correlationplots <- ggarrange(RWcorrelationplot, EWBIcorrelationplot, LWBIcorrelationplot, MEANcorrelationplot, 
                              ncol = 1, nrow = 4, align = "h")

Correlationplots<-annotate_figure(Correlationplots, top = text_grob("Vertical dashed lines represent year of strongest mean COMBO correlation\nSymbols denote 1st and 2nd highest correlations for each time-series\n ", 
                                                          color = "blue", face = "bold", size = 14))
ggsave("correlation_figure.tiff", width = 30, height = 20, units = "cm", dpi = 600, compression = "lzw", bg="white")
}

#=================================
#calculate 1st order AC to adjust N (FULLlength of reference)
#RW
RWrefAC1 <- acf(as.numeric(RWreference$RWReferenceChron), lag = 1 , pl=FALSE, na.action = na.pass)
RWrefAC1 <- RWrefAC1$acf[2]
RWundatedAC1 <- acf(as.numeric(RWundated$RWUndatedChron), lag = 1 , pl=FALSE)
RWundatedAC1 <- RWundatedAC1$acf[2]
RWmeanAC1 <- round(mean((RWrefAC1 + RWundatedAC1)/2), 2)

#if AC1 is negative, the following lines inverts the values for DF adjustment later
if (RWrefAC1 <0) {
  RWrefAC1 <- RWrefAC1*-1
} else {
  RWrefAC1 <- RWrefAC1
}

if (RWundatedAC1 <0) {
  RWundatedAC1 <- RWundatedAC1*-1
} else {
  RWundatedAC1 <- RWundatedAC1
}

#EWBI
EWBIrefAC1 <- acf(as.numeric(EWBIreference$EWBIReferenceChron), lag = 1 , pl=FALSE, na.action = na.pass)
EWBIrefAC1 <- EWBIrefAC1$acf[2]
EWBIundatedAC1 <- acf(as.numeric(EWBIundated$EWBIUndatedChron), lag = 1 , pl=FALSE)
EWBIundatedAC1 <- EWBIundatedAC1$acf[2]
EWBImeanAC1 <- round(mean((EWBIrefAC1 + EWBIundatedAC1)/2), 2)

#if AC1 is negative, the following lines inverts the values for DF adjustment later
if (EWBIrefAC1 <0) {
  EWBIrefAC1 <- EWBIrefAC1*-1
} else {
  EWBIrefAC1 <- EWBIrefAC1
}

if (EWBIundatedAC1 <0) {
  EWBIundatedAC1 <- EWBIundatedAC1*-1
} else {
  EWBIundatedAC1 <- EWBIundatedAC1
}

#LWBI
LWBIrefAC1 <- acf(as.numeric(LWBIreference$LWBIReferenceChron), lag = 1 , pl=FALSE, na.action = na.pass)
LWBIrefAC1 <- LWBIrefAC1$acf[2]
LWBIundatedAC1 <- acf(as.numeric(LWBIundated$LWBIUndatedChron), lag = 1 , pl=FALSE)
LWBIundatedAC1 <- LWBIundatedAC1$acf[2]
LWBImeanAC1 <- round(mean((LWBIrefAC1 + LWBIundatedAC1)/2), 2)

#if AC1 is negative, the following lines inverts the values for DF adjustment later
if (LWBIrefAC1 <0) {
  LWBIrefAC1 <- LWBIrefAC1*-1
} else {
  LWBIrefAC1 <- LWBIrefAC1
}

if (LWBIundatedAC1 <0) {
  LWBIundatedAC1 <- LWBIundatedAC1*-1
} else {
  LWBIundatedAC1 <- LWBIundatedAC1
}

#=================================
#adjust N(for each moving window) for autocorrelation - needed for individual Tvalue plots later
RWNadj <- RWnon_na_counts *((1-(RWrefAC1*RWundatedAC1)/(1+RWrefAC1*RWundatedAC1)))
EWBINadj <- EWBInon_na_counts *((1-(EWBIrefAC1*EWBIundatedAC1)/(1+EWBIrefAC1*EWBIundatedAC1)))
LWBINadj <- LWBInon_na_counts *((1-(LWBIrefAC1*LWBIundatedAC1)/(1+LWBIrefAC1*LWBIundatedAC1)))

#sum the adjusted N series for each parameter
TotalNadj <- cbind(RWNadj, EWBINadj, LWBINadj)
TotalNadj <- rowSums(TotalNadj)

#this Nadj needs to be further adjusted for the average intercorrelation between the reference and undated series
refcorr1 <- cor(as.numeric(RWreference$RWReferenceChron), as.numeric(EWBIreference$EWBIReferenceChron), use="complete.obs")
refcorr2 <- cor(as.numeric(RWreference$RWReferenceChron), as.numeric(LWBIreference$LWBIReferenceChron), use="complete.obs")
refcorr3 <- cor(as.numeric(EWBIreference$EWBIReferenceChron), as.numeric(LWBIreference$LWBIReferenceChron), use="complete.obs")

undcorr1 <- cor(as.numeric(RWundated$RWUndatedChron), as.numeric(EWBIundated$EWBIUndatedChron), use="complete.obs")
undcorr2 <- cor(as.numeric(RWundated$RWUndatedChron), as.numeric(LWBIundated$LWBIUndatedChron), use="complete.obs")
undcorr3 <- cor(as.numeric(EWBIundated$EWBIUndatedChron), as.numeric(LWBIundated$LWBIUndatedChron), use="complete.obs")

#if any correlations are inverse, then the signal is changed - a conservative tweak to the Nadj calculation
if (refcorr1 < 0) {refcorr1 <- refcorr1 * -1}
if (refcorr2 < 0) {refcorr2 <- refcorr2 * -1}
if (refcorr3 < 0) {refcorr3 <- refcorr3 * -1}
if (undcorr1 < 0) {undcorr1 <- undcorr1 * -1}
if (undcorr2 < 0) {undcorr2 <- undcorr2 * -1}
if (undcorr3 < 0) {undcorr3 <- undcorr3 * -1}

#now to calculate the final adjusted degrees of freedom for the combo T calculation
betweenParameterR <- mean((refcorr1 + refcorr2 + refcorr3 + undcorr1 + undcorr2 + undcorr3)/6)
betweenParameterR <- round((betweenParameterR), 2)
adjTotalNadj <- (((TotalNadj)*((1-betweenParameterR)/(1+betweenParameterR))) + ((TotalNadj)/3 * betweenParameterR)) 
#this is an extension of the Dawdy and Matalas 1964 correction of autocorrelation
#additional constant is the scaled (x RBAR) of the min N value for each year

#=================================
#calculate running T values
#=================================

#RW
RWcorrelation_n <- as.data.frame(cbind(RWcorrelation_result, RWNadj))
colnames(RWcorrelation_n) [1] <- "RWcorrelation"
colnames(RWcorrelation_n) [2] <- "RWNadj"

#RW Tvalues using Nadj - not a dataframe yet
RW_Tvalue <- (RWcorrelation_n$RWcorrelation*sqrt(RWcorrelation_n$RWNadj-2))/sqrt(1-RWcorrelation_n$RWcorrelation^2)

#=================================
#create new numerical vector of RW T values of >= 0. This is for the P value calculation
posRW_Tvalue <- ifelse(RW_Tvalue < 0, 0, RW_Tvalue)  
#=================================

#P-values, Bonferroni adjustment, Isolation Factor
RWNoCorrelations <- as.numeric(count(RWcorrelation_n))
log_RWpval <- pt(posRW_Tvalue, df = RWcorrelation_n$RWNadj, lower.tail = FALSE, log.p = TRUE)
RWrawPvalue <- exp(log_RWpval) #1 tailed p-value
RWBonADJPvalue <- RWrawPvalue*RWNoCorrelations # this is the Bonferroni adjustment of the raw P value: P x number of correlations
minRWadjP <- min(RWBonADJPvalue) 
nextminRWadjP <- min(RWBonADJPvalue[RWBonADJPvalue != minRWadjP]) #2nd highest 1/p
RW_IF <- round((nextminRWadjP/minRWadjP), 0) #isolation factor
#however, many adjP values might be > 1 - these need to be adjusted to they are never > 1
RWminadjP <- min(RWBonADJPvalue)
RWBonADJPvalue <- ifelse(RWBonADJPvalue >1 ,1, RWBonADJPvalue) #anything large than 1 gets reset to 1
RWreciprocalP <- as.data.frame(1/RWBonADJPvalue) # this present the adjust P value as 1/p
RWmaxRP <- max(RWreciprocalP) #max 1/p
RWmaxRP <- round(RWmaxRP, 0)

#adjustedP value for the RW data
PRW <- round(min(RWBonADJPvalue), 3)
if (PRW < 0.0001) {PRW = "<0.0001"}

if (RWmaxRP > 1000000) {
  RWmaxRP <- ">1 million"
}

if (RW_IF > 1000) {
  RW_IF <- ">1000"
}

#STDEV from all Tvalues generated - to highlight extreme 4x STDEV threshold
RWTvalSTDEV <- sd(RW_Tvalue, na.rm=TRUE) 

# Find the position of the maximum T-value
RWmax_t_position <- which.max(RW_Tvalue)
if (datatransform == "pw") {RW_Tmaxdate <- (combofirstyearREF + RWmax_t_position + (minoverlap-2))} 
if (datatransform == "fd") {RW_Tmaxdate <- (combofirstyearREF + RWmax_t_position + (minoverlap)-1)} 

#calculate maximum and Minimum T value
RWmaxT <- round(max(RW_Tvalue, na.rm=TRUE), 2)
RWminT <- min(RW_Tvalue, na.rm=TRUE)

#EWBI
EWBIcorrelation_n <- as.data.frame(cbind(EWBIcorrelation_result, EWBINadj))
colnames(EWBIcorrelation_n) [1] <- "EWBIcorrelation"
colnames(EWBIcorrelation_n) [2] <- "EWBINadj"

#EWBI Tvalues using Nadj - not a dataframe yet
EWBI_Tvalue <- (EWBIcorrelation_n$EWBIcorrelation*sqrt(EWBIcorrelation_n$EWBINadj-2))/sqrt(1-EWBIcorrelation_n$EWBIcorrelation^2)

#=================================
#create new numerical vector of EWBI T values of >= 0. This is for the P value calculation
posEWBI_Tvalue <- ifelse(EWBI_Tvalue < 0, 0, EWBI_Tvalue)  
#=================================

#P-values, Bonferroni adjustment, Isolation Factor
EWBINoCorrelations <- as.numeric(count(EWBIcorrelation_n))
log_EWBIpval <- pt(posEWBI_Tvalue, df = EWBIcorrelation_n$EWBINadj, lower.tail = FALSE, log.p = TRUE)
EWBIrawPvalue <- exp(log_EWBIpval) #1 tailed p-value
EWBIBonADJPvalue <- EWBIrawPvalue*EWBINoCorrelations # this is the Bonferroni adjustment of the raw P value: P x number of correlations
minEWBIadjP <- min(EWBIBonADJPvalue) 
nextminEWBIadjP <- min(EWBIBonADJPvalue[EWBIBonADJPvalue != minEWBIadjP]) #2nd highest 1/p
EWBI_IF <- round((nextminEWBIadjP/minEWBIadjP), 0) #isolation factor
#however, many adjP values might be > 1 - these need to be adjusted to they are never > 1
EWBIminadjP <- min(EWBIBonADJPvalue)
EWBIBonADJPvalue <- ifelse(EWBIBonADJPvalue >1 ,1, EWBIBonADJPvalue) #anything large than 1 gets reset to 1
EWBIreciprocalP <- as.data.frame(1/EWBIBonADJPvalue) # this present the adjust P value as 1/p
EWBImaxRP <- max(EWBIreciprocalP) #max 1/p
EWBImaxRP <- round(EWBImaxRP, 0)

#adjsutedP value for the EWBI data
PEWBI <- round(min(EWBIBonADJPvalue), 3)
if (PEWBI < 0.0001) {PEWBI = "<0.0001"}

if (EWBImaxRP > 1000000) {
  EWBImaxRP <- ">1 million"
}

if (EWBI_IF > 1000) {
  EWBI_IF <- ">1000"
}

#STDEV from all Tvalues generated - to highlight extreme 4x STDEV threshold
EWBITvalSTDEV <- sd(EWBI_Tvalue, na.rm=TRUE) 

# Find the position of the maximum T-value
EWBImax_t_position <- which.max(EWBI_Tvalue)
if (datatransform == "pw") {EWBI_Tmaxdate <- (combofirstyearREF + EWBImax_t_position + (minoverlap-2))} 
if (datatransform == "fd") {EWBI_Tmaxdate <- (combofirstyearREF + EWBImax_t_position + (minoverlap)-1)} 

#calculate maximum and Mininum T value
EWBImaxT <- round(max(EWBI_Tvalue, na.rm=TRUE), 2)
EWBIminT <- min(EWBI_Tvalue, na.rm=TRUE)

#LWBI
LWBIcorrelation_n <- as.data.frame(cbind(LWBIcorrelation_result, LWBINadj))
colnames(LWBIcorrelation_n) [1] <- "LWBIcorrelation"
colnames(LWBIcorrelation_n) [2] <- "LWBINadj"

#LWBI Tvalues using Nadj - not a dataframe yet
LWBI_Tvalue <- (LWBIcorrelation_n$LWBIcorrelation*sqrt(LWBIcorrelation_n$LWBINadj-2))/sqrt(1-LWBIcorrelation_n$LWBIcorrelation^2)

#=================================
#create new numerical vector of LWBI T values of >= 0. This is for the P value calculation
posLWBI_Tvalue <- ifelse(LWBI_Tvalue < 0, 0, LWBI_Tvalue)  
#=================================

#P-values, Bonferroni adjustment, Isolation Factor
LWBINoCorrelations <- as.numeric(count(LWBIcorrelation_n))
log_LWBIpval <- pt(posLWBI_Tvalue, df = LWBIcorrelation_n$LWBINadj, lower.tail = FALSE, log.p = TRUE)
LWBIrawPvalue <- exp(log_LWBIpval) #1 tailed p-value
LWBIBonADJPvalue <- LWBIrawPvalue*LWBINoCorrelations # this is the Bonferroni adjustment of the raw P value: P x number of correlations
minLWBIadjP <- min(LWBIBonADJPvalue) 
nextminLWBIadjP <- min(LWBIBonADJPvalue[LWBIBonADJPvalue != minLWBIadjP]) #2nd highest 1/p
LWBI_IF <- round((nextminLWBIadjP/minLWBIadjP), 0) #isolation factor
#however, many adjP values might be > 1 - these need to be adjusted to they are never > 1
LWBIminadjP <- min(LWBIBonADJPvalue)
LWBIBonADJPvalue <- ifelse(LWBIBonADJPvalue >1 ,1, LWBIBonADJPvalue) #anything large than 1 gets reset to 1
LWBIreciprocalP <- as.data.frame(1/LWBIBonADJPvalue) # this present the adjust P value as 1/p
LWBImaxRP <- max(LWBIreciprocalP) #max 1/p
LWBImaxRP <- round(LWBImaxRP, 0)

#adjsutedP value for the LWBI data
PLWBI <- round(min(LWBIBonADJPvalue), 3)
if (PLWBI < 0.0001) {PLWBI = "<0.0001"}

if (LWBImaxRP > 1000000) {
  LWBImaxRP <- ">1 million"
}

if (LWBI_IF > 1000) {
  LWBI_IF <- ">1000"
}

#STDEV from all Tvalues generated - to highlight extreme 4x STDEV threshold
LWBITvalSTDEV <- sd(LWBI_Tvalue, na.rm=TRUE) 

# Find the position of the maximum T-value
LWBImax_t_position <- which.max(LWBI_Tvalue)
if (datatransform == "pw") {LWBI_Tmaxdate <- (combofirstyearREF + LWBImax_t_position + (minoverlap-2))} 
if (datatransform == "fd") {LWBI_Tmaxdate <- (combofirstyearREF + LWBImax_t_position + (minoverlap)-1)}

#calculate maximum and Minmum T value
LWBImaxT <- round(max(LWBI_Tvalue, na.rm=TRUE), 2)
LWBIminT <- min(LWBI_Tvalue, na.rm=TRUE)

#MEAN
MEANcorrelation_n <- as.data.frame(cbind(MEANcorrelation_result, adjTotalNadj))
colnames(MEANcorrelation_n) [1] <- "MEANcorrelation"
colnames(MEANcorrelation_n) [2] <- "adjTotalNadj"

#MEAN Tvalues using Nadj - not a dataframe yet
MEAN_Tvalue <- (MEANcorrelation_n$MEANcorrelation*sqrt(MEANcorrelation_n$adjTotalNadj-2))/sqrt(1-MEANcorrelation_n$MEANcorrelation^2)

#=================================
#create new numerical vector of mean COMBO T values of >= 0. This is for the P value calculation
posMEAN_Tvalue <- ifelse(MEAN_Tvalue < 0, 0, MEAN_Tvalue)  
#=================================

#P-values, Bonferroni adjustment, Isolation Factor
MEANNoCorrelations <- as.numeric(count(MEANcorrelation_n))
log_MEANpval <- pt(posMEAN_Tvalue, df = MEANcorrelation_n$adjTotalNadj, lower.tail = FALSE, log.p = TRUE)
COMBOrawPvalue <- exp(log_MEANpval) #1 tailed p-value
COMBOBonADJPvalue <- COMBOrawPvalue*MEANNoCorrelations # this is the Bonferroni adjustment of the raw P value: P x number of correlations
minCOMBOadjP <- min(COMBOBonADJPvalue) 
nextminCOMBOadjP <- min(COMBOBonADJPvalue[COMBOBonADJPvalue != minCOMBOadjP]) #2nd highest 1/p
COMBO_IF <- round((nextminCOMBOadjP/minCOMBOadjP), 0) #isolation factor
#however, many adjP values might be > 1 - these need to be adjusted to they are never > 1
COMBOminadjP <- min(COMBOBonADJPvalue)
COMBOBonADJPvalue <- ifelse(COMBOBonADJPvalue >1 ,1, COMBOBonADJPvalue) #anything large than 1 gets reset to 1
COMBOreciprocalP <- as.data.frame(1/COMBOBonADJPvalue) # this present the adjust P value as 1/p
COMBOmaxRP <- max(COMBOreciprocalP) #max 1/p
COMBOmaxRP <- round(COMBOmaxRP, 0)

#adjsutedP value for the COMBO data
PCOMBO <- round(min(COMBOBonADJPvalue), 3)
if (PCOMBO < 0.0001) {PCOMBO = "<0.0001"}

if (COMBOmaxRP > 1000000) {
  COMBOmaxRP <- ">1 million"
}

if (COMBO_IF > 1000) {
  COMBO_IF <- ">1000"
}

#STDEV from all Tvalues generated - to highlight extreme 4x STDEV threshold
MEANTvalSTDEV <- sd(MEAN_Tvalue, na.rm=TRUE) 

# Find the position of the maximum T-value
MEANmax_t_position <- which.max(MEAN_Tvalue)
if (datatransform == "pw") {COMBO_Tmaxdate <- (combofirstyearREF + MEANmax_t_position + (minoverlap-2))} 
if (datatransform == "fd") {COMBO_Tmaxdate <- (combofirstyearREF + MEANmax_t_position + (minoverlap)-1)}

#calculate maximum and minimum T value
MEANmaxT <- round(max(MEAN_Tvalue, na.rm=TRUE), 2)
MEANminT <- min(MEAN_Tvalue, na.rm=TRUE)

#create value for degrees of freedom information
ORIGundatedN <- length(ORIGundyears) #length of original undated series in years
summarystats <- summary(undatedRW)
MSL <- round(mean(summarystats$year), 1)
PW_RWundated_df <- length(na.omit(RWundgrow.crn$res)) -2 #length of pre-whitened RW undated series in years
PW_EWBIundate_df <- length(na.omit(EWBIundgrow.crn$res)) -2 #length of pre-whitened EWBI undated series in years
PW_LWBIundated_df <- length(na.omit(LWBIundgrow.crn$res)) -2 #length of pre-whitened LWBI undated series in years

#AC1 adjust max N - for table at end
adjRWmax_df <- round (max(RWNadj), digits = 0) -2
adjEWBImax_df <- round (max(EWBINadj), digits = 0) -2
adjLWBIMax_df <- round (max(LWBINadj), digits = 0) -2

#final total adjusted N for combo series
adjTotalNadj_df <- round (max(adjTotalNadj), digits = 0) -2

#for x-axis range in histograms
histmin <- min(RWminT, EWBIminT, LWBIminT, MEANminT)
histmax <- max(RWmaxT, EWBImaxT, LWBImaxT, MEANmaxT)

#transform for -log10 # needed for the log10 scale for P values
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

#how many to shift out year of undated series to ensure correct date - based on COMBO
shiftyear <- as.numeric(COMBO_Tmaxdate) - as.numeric(lastyearUND)

if (datatransform == "fd") {shiftyear <- shiftyear - 1} #this ensures the correct date when using the FD transform

#if statement for graphical output
if (TvaluePlots == "yes") {
# Plot the Tvalue time series and T Distributions
#p-value minimum for Y axis value
GrandminP <- min(RWBonADJPvalue, EWBIBonADJPvalue, LWBIBonADJPvalue, COMBOBonADJPvalue)

#RW
RWtimeseries <- ggplot() +
  {if(RWmaxT >= 4) annotate(geom = "point", x = RW_Tmaxdate, y = RWmaxT, colour = "blue", size = 3)} + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0+4*RWTvalSTDEV, linetype = "dashed", color = "grey", linewidth = 0.3) +
  {if(datatransform == "pw") geom_line(aes(x = (1:length(RW_Tvalue) - RWpadding) + (RWreferencefirstyear + (RWwindow_size - 2)), y = RW_Tvalue), 
                                       color = "aquamarine4", linewidth = 0.3) } +
  {if(datatransform == "fd") geom_line(aes(x = (1:length(RW_Tvalue) - RWpadding) + (RWreferencefirstyear + (RWwindow_size - 4)), y = RW_Tvalue), 
                                       color = "aquamarine4", linewidth = 0.3) } +
  annotate("text", x=(min(as.numeric(RWreference$Year)) + RWFULLlength), y=4*RWTvalSTDEV, 
           label= "4 STDEVs", color="grey", vjust = 1.5, hjust = 0, fontface = "bold", size = 2.5) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), expand = c(0.01, 0)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  {if(RWmaxT >= 4) labs(title = "RW Sliding T values", x = "Calendar Years CE", y = "T-Value",
       subtitle = paste("Strongest Candidate Outer Year Date = ", RW_Tmaxdate, " CE", "   T = ", RWmaxT))} +
  {if(RWmaxT < 4) labs(title = "RW Sliding T values", x = "Calendar Years CE", y = "T-Value",
                       subtitle = paste("No T Value > 4.0"))} +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        plot.title=element_text(face = "bold"),
        plot.subtitle=element_text(size=12, face = "bold", color="blue"),
        plot.margin = margin(0,0.5,0,0, "cm"))

# p Value plot
#add in na values for non sig 1/p and IF
if (PRW > 0.05) {
  RWmaxRP <- "na"
  RW_IF <- "na"
}

RWPtimeseries <- ggplot() +
  {if(datatransform == "pw") geom_line(aes(x = (1:length(RWBonADJPvalue) - RWpadding) + (RWreferencefirstyear + (RWwindow_size - 2)), y = RWBonADJPvalue), 
                                       color = "black", linewidth = 0.3) } +
  {if(datatransform == "fd") geom_line(aes(x = (1:length(RWBonADJPvalue) - RWpadding) + (RWreferencefirstyear + (RWwindow_size - 4)), y = RWBonADJPvalue), 
                                       color = "black", linewidth = 0.3) } +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), expand = c(0.01, 0)) +
  scale_y_continuous(trans=reverselog_trans(10),
                     breaks = scales::trans_breaks("log10", function(x) 10^x), 
                     labels = scales::trans_format("log10", math_format(10^.x)),
                     limits = c(1, GrandminP)) +
  geom_hline(yintercept=0.01, linetype="dashed", color = "red", linewidth = 0.2) +
  annotate("text", x=(min(as.numeric(RWreference$Year)) + RWFULLlength), y=0.01, 
           label= "p = 0.01", color="red", vjust = 1.5, hjust = 0, fontface = "bold", size = 2.5) +
  geom_hline(yintercept=0.0001, linetype="dashed", color = "red", linewidth = 0.2) +
  annotate("text", x=(min(as.numeric(RWreference$Year)) + RWFULLlength), y=0.0001, 
           label= "p = 0.0001", color="red", vjust = 1.5, hjust = 0, fontface = "bold", size = 2.5) +
  labs(title = "RW Sliding Bonferroni adjusted p values", x = "Calendar Years CE", y = "adjusted P values",
       subtitle = paste("p = ", PRW, "   1/p = ", RWmaxRP, "   IF = ", RW_IF)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        plot.title=element_text(face = "bold"),
        plot.subtitle=element_text(size=12, face = "bold", color="blue"))

#plot the histogram
RW_TvalueDF <- as.data.frame(RW_Tvalue)
RWhistogram <- ggplot(RW_TvalueDF, aes(RW_Tvalue)) +
  geom_vline(xintercept = RWmaxT, linetype = "dashed", color = "blue", linewidth = 0.75) +
  geom_histogram(aes(y = after_stat(density)), colour = "black", fill = "aquamarine4") +
  xlim(histmin, histmax*1.1) +
  {if(RWmaxT >= 4) labs(title = "RW T-value density distribution", x = "T-value", y = "Density",
       subtitle = paste(RW_Tmaxdate, " CE", "   T = ", RWmaxT, "   p = ", PRW))} +
  {if(RWmaxT < 4) labs(title = "RW T-value density distribution", x = "T-value", y = "Density",
                        subtitle = paste("No T Value > 4.0"))} +
  ylim(0, 0.7) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        plot.title=element_text(face = "bold"),
        plot.subtitle=element_text(size=12, face = "bold", color="blue", hjust = 0.8))

RW <- ggarrange(RWtimeseries, RWhistogram, RWPtimeseries, ncol = 3, nrow = 1, align = "h",  widths = c(1, 0.5, 1))

#EWBI
EWBItimeseries <- ggplot() +
  {if(EWBImaxT >= 4) annotate(geom = "point", x = EWBI_Tmaxdate, y = EWBImaxT, colour = "blue", size = 3)} + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0+4*EWBITvalSTDEV, linetype = "dashed", color = "grey", linewidth = 0.3) +
  {if(datatransform == "pw") geom_line(aes(x = (1:length(EWBI_Tvalue) - EWBIpadding) + (EWBIreferencefirstyear + (EWBIwindow_size - 2)), y = EWBI_Tvalue), 
                                       color = "aquamarine4", linewidth = 0.3) } +
  {if(datatransform == "fd") geom_line(aes(x = (1:length(EWBI_Tvalue) - EWBIpadding) + (EWBIreferencefirstyear + (EWBIwindow_size - 4)), y = EWBI_Tvalue), 
                                       color = "aquamarine4", linewidth = 0.3) } +
  annotate("text", x=(min(as.numeric(EWBIreference$Year)) + EWBIFULLlength), y=4*EWBITvalSTDEV, 
           label= "4 STDEVs", color="grey", vjust = 1.5, hjust = 0, fontface = "bold", size = 2.5) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), expand = c(0.01, 0)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  {if(EWBImaxT >= 4) labs(title = "EWBI Sliding T values", x = "Calendar Years CE", y = "T-Value",
                        subtitle = paste("Strongest Candidate Outer Year Date = ", EWBI_Tmaxdate, " CE", "   T = ", EWBImaxT))} +
  {if(EWBImaxT < 4) labs(title = "EWBI Sliding T values", x = "Calendar Years CE", y = "T-Value",
                       subtitle = paste("No T Value > 4.0"))} +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        plot.title=element_text(face = "bold"),
        plot.subtitle=element_text(size=12, face = "bold", color="blue"),
        plot.margin = margin(0,0.5,0,0, "cm"))

# p Value plot
#add in na values for non sig 1/p and IF
if (PEWBI > 0.05) {
  EWBImaxRP <- "na"
  EWBI_IF <- "na"
}

EWBIPtimeseries <- ggplot() +
  {if(datatransform == "pw") geom_line(aes(x = (1:length(EWBIBonADJPvalue) - EWBIpadding) + (EWBIreferencefirstyear + (EWBIwindow_size - 2)), y = EWBIBonADJPvalue), 
                                       color = "black", linewidth = 0.3) } +
  {if(datatransform == "fd") geom_line(aes(x = (1:length(EWBIBonADJPvalue) - EWBIpadding) + (EWBIreferencefirstyear + (EWBIwindow_size - 4)), y = EWBIBonADJPvalue), 
                                       color = "black", linewidth = 0.3) } +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), expand = c(0.01, 0)) +
  scale_y_continuous(trans=reverselog_trans(10),
                     breaks = scales::trans_breaks("log10", function(x) 10^x), 
                     labels = scales::trans_format("log10", math_format(10^.x)),
                     limits = c(1, GrandminP)) +
  geom_hline(yintercept=0.01, linetype="dashed", color = "red", linewidth = 0.2) +
  annotate("text", x=(min(as.numeric(RWreference$Year)) + RWFULLlength), y=0.01, 
           label= "p = 0.01", color="red", vjust = 1.5, hjust = 0, fontface = "bold", size = 2.5) +
  geom_hline(yintercept=0.0001, linetype="dashed", color = "red", linewidth = 0.2) +
  annotate("text", x=(min(as.numeric(RWreference$Year)) + RWFULLlength), y=0.0001, 
           label= "p = 0.0001", color="red", vjust = 1.5, hjust = 0, fontface = "bold", size = 2.5) +
  labs(title = "EWBI Sliding Bonferroni adjusted P values", x = "Calendar Years CE", y = "adjusted P values",
       subtitle = paste("p = ", PEWBI, "   1/p = ", EWBImaxRP, "   IF = ", EWBI_IF)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        plot.title=element_text(face = "bold"),
        plot.subtitle=element_text(size=12, face = "bold", color="blue"))

#plot the histogram
EWBI_TvalueDF <- as.data.frame(EWBI_Tvalue)
EWBIhistogram <- ggplot(EWBI_TvalueDF, aes(EWBI_Tvalue)) +
  geom_vline(xintercept = EWBImaxT, linetype = "dashed", color = "blue", linewidth = 0.75) +
  geom_histogram(aes(y = after_stat(density)), colour = "black", fill = "aquamarine4") +
  xlim(histmin, histmax*1.1) +
  {if(EWBImaxT >= 4) labs(title = "EWBI T-value density distribution", x = "T-value", y = "Density",
                        subtitle = paste(EWBI_Tmaxdate, " CE", "   T = ", EWBImaxT, "   p = ",PEWBI))} +
  {if(EWBImaxT < 4) labs(title = "EWBI T-value density distribution", x = "T-value", y = "Density",
                       subtitle = paste("No T Value > 4.0"))} +
  ylim(0, 0.7) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        plot.title=element_text(face = "bold"),
        plot.subtitle=element_text(size=12, face = "bold", color="blue", hjust = 0.8))

EWBI <- ggarrange(EWBItimeseries, EWBIhistogram, EWBIPtimeseries, ncol = 3, nrow = 1, align = "h",  widths = c(1, 0.5, 1))

#LWBI
LWBItimeseries <- ggplot() +
  {if(LWBImaxT >= 4) annotate(geom = "point", x = LWBI_Tmaxdate, y = LWBImaxT, colour = "blue", size = 3)} + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0+4*LWBITvalSTDEV, linetype = "dashed", color = "grey", linewidth = 0.3) +
  {if(datatransform == "pw") geom_line(aes(x = (1:length(LWBI_Tvalue) - LWBIpadding) + (LWBIreferencefirstyear + (LWBIwindow_size - 2)), y = LWBI_Tvalue), 
                                       color = "aquamarine4", linewidth = 0.3) } +
  {if(datatransform == "fd") geom_line(aes(x = (1:length(LWBI_Tvalue) - LWBIpadding) + (LWBIreferencefirstyear + (LWBIwindow_size - 4)), y = LWBI_Tvalue), 
                                       color = "aquamarine4", linewidth = 0.3) } +
  annotate("text", x=(min(as.numeric(LWBIreference$Year)) + RWFULLlength), y=4*LWBITvalSTDEV, 
           label= "4 STDEVs", color="grey", vjust = 1.5, hjust = 0, fontface = "bold", size = 2.5) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), expand = c(0.01, 0)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  {if(LWBImaxT >= 4) labs(title = "LWBI Sliding T values", x = "Calendar Years CE", y = "T-Value",
                        subtitle = paste("Strongest Candidate Outer Year Date = ", LWBI_Tmaxdate, " CE", "   T = ", LWBImaxT))} +
  {if(LWBImaxT < 4) labs(title = "LWBI Sliding T values", x = "Calendar Years CE", y = "T-Value",
                       subtitle = paste("No T Value > 4.0"))} +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        plot.title=element_text(face = "bold"),
        plot.subtitle=element_text(size=12, face = "bold", color="blue"),
        plot.margin = margin(0,0.5,0,0, "cm"))

# p Value plot
#add in na values for non sig 1/p and IF
if (PLWBI > 0.05) {
  LWBImaxRP <- "na"
  LWBI_IF <- "na"
}

LWBIPtimeseries <- ggplot() +
  {if(datatransform == "pw") geom_line(aes(x = (1:length(LWBIBonADJPvalue) - LWBIpadding) + (LWBIreferencefirstyear + (LWBIwindow_size - 2)), y = LWBIBonADJPvalue), 
                                       color = "black", linewidth = 0.3) } +
  {if(datatransform == "fd") geom_line(aes(x = (1:length(LWBIBonADJPvalue) - LWBIpadding) + (LWBIreferencefirstyear + (LWBIwindow_size - 4)), y = LWBIBonADJPvalue), 
                                       color = "black", linewidth = 0.3) } +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), expand = c(0.01, 0)) +
  scale_y_continuous(trans=reverselog_trans(10),
                     breaks = scales::trans_breaks("log10", function(x) 10^x), 
                     labels = scales::trans_format("log10", math_format(10^.x)),
                     limits = c(1, GrandminP)) +
  geom_hline(yintercept=0.01, linetype="dashed", color = "red", linewidth = 0.2) +
  annotate("text", x=(min(as.numeric(RWreference$Year)) + RWFULLlength), y=0.01, 
           label= "p = 0.01", color="red", vjust = 1.5, hjust = 0, fontface = "bold", size = 2.5) +
  geom_hline(yintercept=0.0001, linetype="dashed", color = "red", linewidth = 0.2) +
  annotate("text", x=(min(as.numeric(RWreference$Year)) + RWFULLlength), y=0.0001, 
           label= "p = 0.0001", color="red", vjust = 1.5, hjust = 0, fontface = "bold", size = 2.5) +
  labs(title = "LWBI Sliding Bonferroni adjusted P values", x = "Calendar Years CE", y = "adjusted P values",
       subtitle = paste("p = ", PLWBI, "   1/p = ", LWBImaxRP, "   IF = ", LWBI_IF)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        plot.title=element_text(face = "bold"),
        plot.subtitle=element_text(size=12, face = "bold", color="blue"))

#plot the histogram
LWBI_TvalueDF <- as.data.frame(LWBI_Tvalue)
LWBIhistogram <- ggplot(LWBI_TvalueDF, aes(LWBI_Tvalue)) +
  geom_vline(xintercept = LWBImaxT, linetype = "dashed", color = "blue", linewidth = 0.75) +
  geom_histogram(aes(y = after_stat(density)), colour = "black", fill = "aquamarine4") +
  xlim(histmin, histmax*1.1) +
  {if(LWBImaxT >= 4) labs(title = "LWBI T-value density distribution", x = "T-value", y = "Density",
                        subtitle = paste(LWBI_Tmaxdate, " CE", "   T = ", LWBImaxT, "   p = ",PLWBI))} +
  {if(LWBImaxT < 4) labs(title = "LWBI T-value density distribution", x = "T-value", y = "Density",
                       subtitle = paste("No T Value > 4.0"))} +
  ylim(0, 0.7) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        plot.title=element_text(face = "bold"),
        plot.subtitle=element_text(size=12, face = "bold", color="blue", hjust = 0.8))

LWBI <- ggarrange(LWBItimeseries, LWBIhistogram, LWBIPtimeseries, ncol = 3, nrow = 1, align = "h",  widths = c(1, 0.5, 1))


#Combo
COMBOtimeseries <- ggplot() +
  {if(MEANmaxT >= 4) annotate(geom = "point", x = COMBO_Tmaxdate, y = MEANmaxT, colour = "blue", size = 3)} + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0+4*MEANTvalSTDEV, linetype = "dashed", color = "grey", linewidth = 0.3) +
  {if(datatransform == "pw") geom_line(aes(x = (1:length(MEAN_Tvalue) - RWpadding) + (RWreferencefirstyear + (RWwindow_size - 2)), y = MEAN_Tvalue), 
                                       color = "aquamarine4", linewidth = 0.3) } +
  {if(datatransform == "fd") geom_line(aes(x = (1:length(MEAN_Tvalue) - RWpadding) + (RWreferencefirstyear + (RWwindow_size - 4)), y = MEAN_Tvalue), 
                                       color = "aquamarine4", linewidth = 0.3) } +
  annotate("text", x=(min(as.numeric(RWreference$Year)) + RWFULLlength), y=4*MEANTvalSTDEV, 
           label= "4 STDEVs", color="grey", vjust = 1.5, hjust = 0, fontface = "bold", size = 2.5) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), expand = c(0.01, 0)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  {if(MEANmaxT >= 4) labs(title = "COMBO Sliding T values", x = "Calendar Years CE", y = "T-Value",
                        subtitle = paste("Strongest Candidate Outer Year Date = ", COMBO_Tmaxdate, " CE", "   T = ", MEANmaxT))} +
  {if(MEANmaxT < 4) labs(title = "COMBINED Sliding T values", x = "Calendar Years CE", y = "T-Value",
                       subtitle = paste("No T Value > 4.0"))} +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        plot.title=element_text(face = "bold"),
        plot.subtitle=element_text(size=12, face = "bold", color="blue"),
        plot.margin = margin(0,0.5,0,0, "cm"))

# p Value plot
#add in na values for non sig 1/p and IF
if (PCOMBO > 0.05) {
  COMBOmaxRP <- "na"
  COMBO_IF <- "na"
}

COMBOPtimeseries <- ggplot() +
  {if(datatransform == "pw") geom_line(aes(x = (1:length(COMBOBonADJPvalue) - RWpadding) + (RWreferencefirstyear + (RWwindow_size - 2)), y = COMBOBonADJPvalue), 
                                       color = "black", linewidth = 0.3) } +
  {if(datatransform == "fd") geom_line(aes(x = (1:length(COMBOBonADJPvalue) - RWpadding) + (RWreferencefirstyear + (RWwindow_size - 4)), y = COMBOBonADJPvalue), 
                                       color = "black", linewidth = 0.3) } +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), expand = c(0.01, 0)) +
  scale_y_continuous(trans=reverselog_trans(10),
                     breaks = scales::trans_breaks("log10", function(x) 10^x), 
                     labels = scales::trans_format("log10", math_format(10^.x)),
                     limits = c(1, GrandminP)) +
  geom_hline(yintercept=0.01, linetype="dashed", color = "red", linewidth = 0.2) +
  annotate("text", x=(min(as.numeric(RWreference$Year)) + RWFULLlength), y=0.01, 
           label= "p = 0.01", color="red", vjust = 1.5, hjust = 0, fontface = "bold", size = 2.5) +
  geom_hline(yintercept=0.0001, linetype="dashed", color = "red", linewidth = 0.2) +
  annotate("text", x=(min(as.numeric(RWreference$Year)) + RWFULLlength), y=0.0001, 
           label= "p = 0.0001", color="red", vjust = 1.5, hjust = 0, fontface = "bold", size = 2.5) +
  labs(title = "COMBO Sliding Bonferroni adjusted P values", x = "Calendar Years CE", y = "adjusted P values",
       subtitle = paste("p = ", PCOMBO, "   1/p = ", COMBOmaxRP, "   IF = ", COMBO_IF)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        plot.title=element_text(face = "bold"),
        plot.subtitle=element_text(size=12, face = "bold", color="blue"))

#=================================

#plot the histogram
MEAN_TvalueDF <- as.data.frame(MEAN_Tvalue)
COMBOhistogram <- ggplot(MEAN_TvalueDF, aes(MEAN_Tvalue)) +
  geom_vline(xintercept = MEANmaxT, linetype = "dashed", color = "blue", linewidth = 0.75) +
  geom_histogram(aes(y = after_stat(density)), colour = "black", fill = "aquamarine4") +
  xlim(histmin, histmax*1.1) +
  {if(MEANmaxT >= 4) labs(title = "COMBO T-value density distribution", x = "T-value", y = "Density",
                        subtitle = paste(COMBO_Tmaxdate, " CE", "   T = ", MEANmaxT, "   p = ",PCOMBO))} +
  {if(MEANmaxT < 4) labs(title = "COMBO T-value density distribution", x = "T-value", y = "Density",
                       subtitle = paste("No T Value > 4.0"))} +
  ylim(0, 0.7) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        plot.title=element_text(face = "bold"),
        plot.subtitle=element_text(size=12, face = "bold", color="blue", hjust = 0.8))

COMBO <- ggarrange(COMBOtimeseries, COMBOhistogram, COMBOPtimeseries, ncol = 3, nrow = 1, align = "h",  widths = c(1, 0.5, 1))

#dataframes for graphics figure table
TRvartiables<- as.data.frame(c("RW", "EWBI", "LWBI"))
AC1 <- as.data.frame(c(RWmeanAC1, EWBImeanAC1, LWBImeanAC1)) #see earlier for values
AC1adjDF <-as.data.frame(c(adjRWmax_df, adjEWBImax_df, adjLWBIMax_df))

sumtable<-cbind(TRvartiables, AC1, AC1adjDF)
names(sumtable) = c("TR Parameter", "mean AC1", "AC1adjDF") #need to rename the columns

#create formatted table
tab <- ggtexttable(sumtable, theme = ttheme("blank"), rows = NULL)

tab <- table_cell_font(tab, row = 2:4, column = 1, color = "black", face = "bold")

summarytable <- tab %>%
  tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(4), row.side = "bottom", linewidth = 3, linetype = 1) %>%
  tab_add_vline(at.column = 2:tab_ncol(tab), column.side = "left", from.row = 1, linetype = 2)

text1 <- text_grob(label = paste("Undated series length = ", ORIGundatedN, "years", "\nMSL = ", MSL), 
                   color = "red", size = 16, face = "bold")

text2 <- text_grob(label = paste("Between TR parameter \nRBAR = ", betweenParameterR), 
                   color = "red", size = 16, face = "bold")

if(MEANmaxT >= 4) {text3 <- text_grob(label = paste("Combo Multi-Parameter adjDF = ", adjTotalNadj_df,
                                                    "\n\nShift series by", shiftyear, "years to attain date"),
                         hjust = 0, color = "red", size = 16, face = "bold",
                         x = unit(30, "pt"))

} else {text3 <- text_grob(label = paste("NO SIGNIFICANT CROSSDATE"),
                         color = "red", size = 16, face = "bold")
}


lower <- ggarrange(text1, text2, summarytable, text3, ncol = 4, nrow = 1,  
                   widths = c(0.7, 0.55, 1, 1.2))

finalfigure <- ggarrange(RW, EWBI, LWBI, COMBO, lower, ncol = 1, nrow = 5, align = "v", 
                         heights = c(1, 1, 1, 1, 0.6))
finalfigure<-annotate_figure(finalfigure, top = text_grob(runtitle, 
                                                                color = "red", face = "bold", size = 28))

ggsave("T_value_figure.tiff", width = 50, height = 25, units = "cm", dpi = 600, compression = "lzw", bg="white")
}

#Final summary table - generated regardless of graphical output options
#correlation values at final COMBO max T value position
RW_Rcorr <- round(RWcorrelation_result[MEANmax_t_position], 2)
EWBI_Rcorr <- round(EWBIcorrelation_result[MEANmax_t_position], 2)
LWBI_Rcorr <- round(LWBIcorrelation_result[MEANmax_t_position], 2)
COMBO_Rcorr <- round(MEANcorrelation_result[MEANmax_t_position], 2)

if (PRW > 0.05) {RWmaxRP = "na"}
if (PEWBI > 0.05) {EWBImaxRP = "na"}
if (PLWBI > 0.05) {LWBImaxRP = "na"}
if (PCOMBO > 0.05) {COMBOmaxRP = "na"}

if (PRW > 0.05) {RW_IF = "na"}
if (PEWBI > 0.05) {EWBI_IF = "na"}
if (PLWBI > 0.05) {LWBI_IF = "na"}
if (PCOMBO > 0.05) {COMBO_IF = "na"}

if (PRW > 0.05) {RW_Tmaxdate = "no sig date"}
if (PEWBI > 0.05) {EWBI_Tmaxdate = "no sig date"}
if (PLWBI > 0.05) {LWBI_Tmaxdate = "no sig date"}
if (PCOMBO > 0.05) {COMBO_Tmaxdate = "no sig date"}

if (PRW < 0.05) {RWshiftyear <- shiftyear} else {RWshiftyear <- "no sig date"}
if (PEWBI < 0.05) {EWBIshiftyear <- shiftyear} else {EWBIshiftyear <- "no sig date"}
if (PLWBI < 0.05) {LWBIshiftyear <- shiftyear} else {LWBIshiftyear <- "no sig date"}
if (PCOMBO < 0.05) {COMBOshiftyear <- shiftyear} else {COMBOshiftyear <- "no sig date"}
  
#dataframes for graphics figure table
TRvartiables<- as.data.frame(c("RW", "EWBI", "LWBI", "COMBO"))
Rvalues <- as.data.frame(c(RW_Rcorr, EWBI_Rcorr, LWBI_Rcorr, COMBO_Rcorr)) #T values for each parameter
Tvalues <- as.data.frame(c(RWmaxT, EWBImaxT, LWBImaxT, MEANmaxT)) #T values for each parameter
Pvalues <- as.data.frame(c(PRW, PEWBI, PLWBI, PCOMBO)) #p values for each parameter
ReciprocalP <- as.data.frame(c(RWmaxRP, EWBImaxRP, LWBImaxRP, COMBOmaxRP)) #1/p values for each parameter
IFvalues <- as.data.frame(c(RW_IF, EWBI_IF, LWBI_IF, COMBO_IF)) #IF values for each parameter
proposeddate <- as.data.frame(c(RW_Tmaxdate, EWBI_Tmaxdate, LWBI_Tmaxdate, COMBO_Tmaxdate)) #proposed outer date
shiftyears <- as.data.frame(c(RWshiftyear, EWBIshiftyear, LWBIshiftyear, COMBOshiftyear))

sumtable<-cbind(TRvartiables, Rvalues, Tvalues, Pvalues, ReciprocalP, IFvalues, proposeddate, shiftyears)
names(sumtable) = c("TR Parameter", "r value", "T value", "p value", "1/p", "IF", "Proposed Outer Date", "No. of Years to shift") #table headers

#create formatted table
tab <- ggtexttable(sumtable, theme = ttheme("blank"), rows = NULL)

tab <- table_cell_font(tab, row = 2:5, column = 1, color = "black", face = "bold")

UBERsummarytable <- tab %>%
  tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(5), row.side = "bottom", linewidth = 3, linetype = 1) %>%
  tab_add_vline(at.column = 2:tab_ncol(tab), column.side = "left", from.row = 1, linetype = 2)

UBERsummarytable

#if you want to export the final table for Word or Excel - use this line
#write.csv (sumtable, "UBERsummarytable.csv")
