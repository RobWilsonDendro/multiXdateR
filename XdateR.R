#this sets the working directory to the source file location - nice
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#=================================
library(readxl); library(zoo); library(tidyverse); library(ggpubr)
library(dplR); library(qpcR); library(dplyr); library(scales)
#=================================
runtitle <- "test run"
#=================================
#Tucson file input and detrending
referenceTR <- read.rwl ("NCairnEWBI", header = FALSE) #read in reference file data
undatedTR <- read.rwl ("GOSewbi.rwl", header = FALSE) #read in reference file data

#=================================
#minimum overlap between undated and reference chrons
minoverlap <- 50 #ensure this is not longer than the length of your undated series
#I would also not make this value too small. 50 is already too small to be honest.
#However N is factored into the Tvalue calculation so no real harm - just be cautious.

#==================================
#detrending method is a fixed spline using a biweight robust mean with prewhitening
#power transform of the raw data is now an option.
#user can change spline window size, unfiltered or 1st difference transform and correlation type
#degrees of freedom are adjusted related to the AC1 of the resultant chronologies
powertransform <- "yes" # "yes" or "no" if you want to implement the powertransform (Cook & Peters 1997)
splinewindow <- 31
datatransform <- "pw" #residual chronology - options of "pw" = unfiltered or "fd" = 1st differenced
correlationtype <- "spearman" # use pearson or spearman
Tthreshold <- 4 #threshold of T-value acceptance

#=================================******
#=================================******
#detrends the reference data
if (powertransform == "yes") {referenceTR <- powt(referenceTR, rescale = FALSE) #this is the power transform of the data
}

ref.rwi <- detrend(rwl = referenceTR, method = c("Spline"), nyrs = splinewindow, f = 0.5, 
                    pos.slope = TRUE, difference = FALSE) 

#chronology build without variance stabilisation
refgrow.crn <- chron(x = ref.rwi, prefix = "ref", biweight = TRUE, prewhiten = TRUE)
refyears <- row.names(refgrow.crn)

#==================================
#detrends the undated data
if (powertransform == "yes") {undatedTR <- powt(undatedTR, rescale = FALSE) #this is the power transform of the data
}

und.rwi <- detrend(rwl = undatedTR, method = c("Spline"), nyrs = splinewindow, f = 0.5, 
                   pos.slope = TRUE, difference = FALSE) 

#chronology build without variance stabilisation
undgrow.crn <- chron(x = und.rwi, prefix = "ref", biweight = TRUE, prewhiten = TRUE)
undyears <- row.names(undgrow.crn)

#==================================
#then create dataframes for analysis based on data transform - pw or df
if (datatransform == "pw") {reference <- as.data.frame(cbind(refyears, refgrow.crn$res))} #unfiltered - use DF as is
if (datatransform == "pw") {undated <- as.data.frame(cbind(undyears, undgrow.crn$res))}
if (datatransform == "fd") {ref1st <- diff(refgrow.crn$res)
reference <- qpcR:::cbind.na(refyears, ref1st)
reference <- as.data.frame(reference)
}
if (datatransform == "fd") {und1st <- diff(undgrow.crn$res)
undated <- qpcR:::cbind.na(undyears, und1st)
undated <- as.data.frame(undated)
} 
colnames(reference) [1] <- "Year"
colnames(reference) [2] <- "ReferenceChron"
colnames(undated) [1] <- "Year"
colnames(undated) [2] <- "UndatedChron"
reference <- na.omit(reference) #removes the NA values from the Residual chronology
undated <- na.omit(undated) #removes the NA values from the Residual chronology

#=================================
#metrics needed for calculations
FULLlength <- length(as.numeric(undated$UndatedChron))
window_size <- FULLlength
referencefirstyear <- min(as.numeric(reference$Year))

if (minoverlap > FULLlength) {
  print("DUDE! - your minimum overlap is longer than your undated chronology")
  print("Code has been stopped!!!")
  stop()
}

#=================================
#add padding to reference dataframe to allow extension beyond the reference chron range using overlap
padding <- FULLlength - minoverlap #no. of NA values to be added

reference <- data.frame(Value = c(rep(NA, padding), reference$ReferenceChron))
#this adds NA values to beginning of DF, but deletes the year column
colnames(reference) [1] <- "ReferenceChron"

#this codes adds NA to end - but need years to be added back
reference <- reference %>%
  slice(rep(2, padding)) %>%
  mutate(across(everything(), ~ NA)) %>%
  bind_rows(reference, .)

reflength <- length(reference$ReferenceChron)
newstart <- referencefirstyear - padding

#create a new column of specific annual values
newdates <- seq(from = newstart, to = newstart + reflength, by = 1)
reference <- qpcR:::cbind.na(newdates, reference)
colnames(reference) [1] <- "Year"
colnames(reference) [2] <- "ReferenceChron"

#=================================
#sliding correlation
# Function to calculate correlation coefficient
calculate_correlation <- function(x, y) {
  cor(x, y, method = correlationtype, use = "pairwise.complete.obs")
  #pairwise.complete.obs ignores NA values but still calculates r value.
}

# Calculate correlation coefficients and N using rollapply
correlation_result <- rollapply(as.numeric(reference$ReferenceChron), width = window_size,  
                                FUN = calculate_correlation, align = "right", by = 1, y = as.numeric(undated$UndatedChron))

#count number of correlations per window
count_non_na <- function(x) sum(!is.na(x))
non_na_counts <- rollapply(reference$ReferenceChron, width = window_size, FUN = count_non_na, align = "right")

#=================================
#calculate 1st order AC to adjust N (FULLlength of reference)
refAC1 <- acf(as.numeric(reference$ReferenceChron), lag = 1 , pl=FALSE, na.action = na.pass)
refAC1 <- refAC1$acf[2]
undatedAC1 <- acf(as.numeric(undated$UndatedChron), lag = 1 , pl=FALSE)
undatedAC1 <- undatedAC1$acf[2]
if (refAC1 <0) {
  refAC1 <- refAC1*-1
} else {
  refAC1 <- refAC1
}

if (undatedAC1 <0) {
  undatedAC1 <- undatedAC1*-1
} else {
  undatedAC1 <- undatedAC1
}

#=================================
#adjust N(for each moving window) for autocorrelation
Nadj <- non_na_counts *((1-(refAC1*undatedAC1)/(1+refAC1*undatedAC1)))

#calculate running T values
correlation_n <- as.data.frame(cbind(correlation_result, Nadj))
colnames(correlation_n) [1] <- "correlation"
colnames(correlation_n) [2] <- "Nadj"

#Tvalues using Nadj - not a dataframe yet
Tvalue <- (correlation_n$correlation*sqrt(correlation_n$Nadj-2))/sqrt(1-correlation_n$correlation^2)

#calculate the Probability of the T-values
NadjfromT <- as.data.frame(cbind(Tvalue, correlation_n$Nadj))
colnames(NadjfromT) [1] <- "Tvalue"
colnames(NadjfromT) [2] <- "Nadj"

#P-values, Bonferroni adjustment, Isolation Factor
NoCorrelations <- as.numeric(count(correlation_n))
rawPvalue <- (2 * pt(-abs(Tvalue), NadjfromT$Nadj))
BonADJPvalue <- rawPvalue*NoCorrelations # this is the Bonferroni adjustment of the raw P value: P x number of correlations
#however, many adjP values might be > 1 - these need to be adjusted to they are never > 1
BonADJPvalue <- ifelse(BonADJPvalue >1 ,1, BonADJPvalue)
reciprocalP <- as.data.frame(1/BonADJPvalue) # this present the adjust P value as 1/p
maxRP <- max(reciprocalP) #max 1/p
NEXTmaxRP <- max(reciprocalP[reciprocalP != maxRP]) #2nd highest 1/p
IF <- round((maxRP - NEXTmaxRP), 0) #isolation factor
maxRP <- round(maxRP, 0)

if (maxRP > 1000000) {
  maxRP <- ">1 million"
}

if (IF > 1000) {
  IF <- ">1000"
}

#use this code to extract Nadj for max T value
NadjforTmax <- NadjfromT$Nadj[which.max(NadjfromT$Tvalue)]

#Pvalue for Max T with Bonferroni adjustment
PofT <- round(2 * (pt(-abs(max(Tvalue)), NadjforTmax)*NoCorrelations), 3)
if (PofT < 0.0001) {PofT = "<0.0001"}

#STDEV from all Tvalues generated - to highlight extreme 4x STDEV threshold
TvalSTDEV <- sd(Tvalue, na.rm=TRUE)

# Find the position of the maximum T-value
max_t_position <- which.max(Tvalue)
outeryear <- (referencefirstyear + max_t_position + window_size - 2) - padding #should depend on 1st year of reference chron
if (datatransform == "pw") {
  outeryear <- outeryear
} else {
  outeryear <- outeryear +1
}

if (datatransform == "pw") {
  referencefirstyear <- referencefirstyear
} else {
  referencefirstyear <- referencefirstyear +1
}

#calculate maximum T value
maxT <- round(max(Tvalue, na.rm=TRUE), 2)

# Plot the Tvalue time series
Ttimeseries <- ggplot() +
  annotate(geom = "point", x = outeryear, y = maxT, colour = "blue", size = 3) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0+4*TvalSTDEV, linetype = "dashed", color = "grey", linewidth = 0.3) +
  geom_line(aes(x = (1:length(Tvalue) - padding) + (referencefirstyear + (window_size - 2)), y = Tvalue), 
            color = "aquamarine4", linewidth = 0.3) +
  annotate("text", x=(min(as.numeric(reference$Year)) + FULLlength), y=4*TvalSTDEV, 
           label= "4 STDEVs", color="grey", vjust = 1.5, fontface = "bold", size = 3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), expand = c(0.01, 0)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  labs(title = "Sliding T values", x = "Calendar Years CE", y = "T-Values",
       subtitle = paste("Strongest Candidate Outer Year Date = ", outeryear, " CE", "   T = ", maxT)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        plot.title=element_text(face = "bold"),
        plot.subtitle=element_text(size=12, face = "bold", color="blue"))

#transform for -log10
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

# p Value plot
Ptimeseries <- ggplot() +
  geom_line(aes(x = (1:length(BonADJPvalue) - padding) + (referencefirstyear + (window_size - 2)), y = BonADJPvalue), 
            color = "black", linewidth = 0.3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), expand = c(0.01, 0)) +
  scale_y_continuous(trans=reverselog_trans(10),
                     breaks = scales::trans_breaks("log10", function(x) 10^x), 
                     labels = scales::trans_format("log10", math_format(10^.x))) +
  labs(title = "Sliding Bonferroni adjusted P values", x = "Calendar Years CE", y = "adjusted P values",
       subtitle = paste("P = ", PofT, "   1/p = ", maxRP, "   IF = ", IF)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        plot.title=element_text(face = "bold"),
        plot.subtitle=element_text(size=12, face = "bold", color="blue"))

#plot the histogram
TvalueDF <- as.data.frame(Tvalue)
histogram <- ggplot(TvalueDF, aes(Tvalue)) +
  geom_vline(xintercept = max(TvalueDF$Tvalue), linetype = "dashed", color = "blue", linewidth = 0.75) +
  geom_segment(aes(x = maxT-2, y = 0.4, xend = maxT-0.1, yend = 0.4),
               arrow = arrow(length = unit(0.33, "cm")), colour = "blue", linewidth = 0.75) +
  geom_histogram(aes(y = after_stat(density)), colour = "black", fill = "aquamarine4") +
  labs(title = "T-value density distribution", x = "T-value", y = "Density",
       subtitle = paste(outeryear, " CE", "   T = ", maxT)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        plot.title=element_text(face = "bold"),
        plot.subtitle=element_text(size=12, face = "bold", color="blue", hjust = 1))

#time-series comparison plot
undated <- as.data.frame(undated)
undatedlastyear <- max(undated$Year)
originalfinalyear <- tail(undated$Year, n=1)

if (datatransform == "pw") {
  referencefirstyear <- referencefirstyear
}
difference <- as.numeric(outeryear) - as.numeric(originalfinalyear)

newyears <- as.numeric(undated$Year) + difference
dated <- as.data.frame(cbind(newyears, undated$UndatedChron))
colnames(dated) <- c("Year", "Dated")

if (datatransform == "fd") {
  reference <- reference %>%
    mutate(ReferenceChron=lag(ReferenceChron)) %>%
    na.omit()
}

comparisonDF <- (merge(reference, dated, by = "Year"))
temporaryDF <- subset(comparisonDF, select = -1)
temporaryDF <- sapply(temporaryDF, as.numeric)
temporaryDF <- scale(temporaryDF)
comparisonDF <- cbind(comparisonDF$Year, temporaryDF)
colnames(comparisonDF)[1] <- "Year"
comparisonDF <- as.data.frame(comparisonDF)

library(reshape2)
comparisonDFlong <- melt(comparisonDF, id.vars = "Year")
comparisonDFlong$Year <- as.numeric(as.character(comparisonDFlong$Year)) 
comparisonDFlong$value <- as.numeric(as.character(comparisonDFlong$value)) 
comparisonDF$Dated <- as.numeric(as.character(comparisonDF$Dated)) 

#correlation value for comaprison plot
correlation_result <- as.data.frame(correlation_result)
Rcorr <- round(correlation_result$correlation_result[max_t_position], 2)

#sneak this in here so the shifting year is correct
if (datatransform == "fd") {difference <- difference - 1}

if (maxT >= Tthreshold) {
comparisonplot <- ggplot(comparisonDFlong, aes(x = Year, y = value, colour = variable)) +
  geom_line(linewidth = 0.4) +
  scale_color_manual(name = "", values = c("ReferenceChron" = "black", "Dated" = "red")) +
  geom_vline(xintercept = outeryear, linetype = "dashed", color = "blue", linewidth = 0.75) +
  annotate("text", x=outeryear, y=mean(comparisonDF$Dated) + max(comparisonDFlong$value*1.35), 
           label= paste(outeryear, " CE", "\nr = ", Rcorr, "   P = ", PofT), color="blue", hjust = 1.1, fontface = "bold") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(limits = c(min(comparisonDFlong$value*1.5), max(comparisonDFlong$value*1.5))) +
  labs(title = "Dated Chronology Comparison", x = "Calendar Years", y = "TR Index",
       subtitle = paste("shift undated data by", difference, "years to attain correct date")) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        legend.direction="horizontal",
        legend.position = c(0.5, 0.1),
        plot.title=element_text(face = "bold"),
        plot.subtitle=element_text(size=12, face = "bold", color="blue", hjust = 1),
        legend.background=element_rect(fill = alpha("white", 0.5)))
} else {

  comparisonplot <- ggplot(comparisonDFlong, aes(x = Year, y = value, colour = variable)) +
    geom_line(linewidth = 0.4) +
    scale_color_manual(name = "", values = c("ReferenceChron" = "white", "Dated" = "white")) +
    annotate("label", x = mean(comparisonDFlong$Year), y = 0, 
             label = paste("WARNING: Tvalue < ", Tthreshold), 
             size = 10, colour = "red", fontface =2, fill = "yellow", alpha = 0.5) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_y_continuous(limits = c(min(comparisonDFlong$value*1.5), max(comparisonDFlong$value*1.5))) +
    labs(title = "No Chronology Comparison", x = "Calendar Years", y = "TR Index",
         subtitle = paste("No signficant date identified")) +
    theme_classic() +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.title.x=element_blank(),
          legend.direction="horizontal",
          legend.position="none",
          plot.title=element_text(face = "bold"),
          plot.subtitle=element_text(size=12, face = "bold", color="blue", hjust = 1))
}

lower <- ggarrange(histogram, comparisonplot, ncol = 2, nrow = 1, align = "h",  widths = c(1.5,2))
finalfigure <- ggarrange(Ttimeseries, Ptimeseries, lower, ncol = 1, nrow = 3, heights = c(1.5, 1.5 ,2))
finalfigure<-annotate_figure(finalfigure, top = text_grob(runtitle, 
                                                          color = "red", face = "bold", size = 20))
finalfigure
ggsave("dating_result.tiff", width = 28, height = 20, units = "cm", dpi = 600, compression = "lzw", bg="white")
