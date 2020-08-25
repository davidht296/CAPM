# File:    CAPM Estimation - ASX LICs (5yr, Monthly)
# Author:  David Harris 
# Date:    15-04-2020

# INSTALL AND LOAD PACKAGES ################################

# Installs pacman ("package manager") if needed
if (!require("pacman")) install.packages("pacman")

# Use pacman to load add-on packages as desired
pacman::p_load(pacman, aplot, boot, car, caret, expss, GGally,
               ggthemes, ggvis, httr, huxtable, jtools, lars, 
               lmtest, lubridate, 
               MASS, NCmisc, olsrr, PerformanceAnalytics, plotly, plyr, 
               psych, quantmod, remotes, rio, rmarkdown, sandwich, shiny, sur, 
               tidyquant, tidyverse, xts)  

# LOAD t-TEST FUNCTION #####################################
ttest <- function(reg, coefnum, val){
  co <- coef(summary(reg))
  tstat <- (co[coefnum,1]-val)/co[coefnum,2]
  2 * pt(abs(tstat), reg$df.residual, lower.tail = FALSE)
}

# SET WORKING DIRECTORY ####################################
getwd()
setwd("C:/Users/dtunks/OneDrive - KPMG/Desktop/PD & Personal/R data analysis/CAPM") #*remember to change \ to /

# IMPORT AND TIDY LIC LIST #################################
LICs <- read.csv("https://www.asxlics.com/uploads/csv/20200401-lics.csv", header = TRUE)
n <- nrow(LICs)
LICs <- LICs[c(2:n), 1:3]
lic.colnames <- c("Code", "Company", "Market Cap")
names(LICs) <- lic.colnames
ticker <- as.character(LICs[,1])
row.names(LICs) <- ticker

# IMPORT AND TIDY RISK-FREE RATE DATA ######################
Rf <- import("f2.1-data.csv") ## need to have manually formatted dates to YYYY-MM-DD in Excel
n <- nrow(Rf)
Rf <- Rf[c(12:n), c(1, 4)]
Rf <- Rf[!apply(Rf == "", 1, all),]
Rf$V1 <- as.Date(Rf$V1)
Rf$V4 <- as.numeric(Rf$V4)
Rf$V4 <- ((1+(Rf$V4/100))^(1/12)-1)
Rf <- xts(Rf$V4, order.by = Rf$V1)
names(Rf) <- c("Rf")

# IMPORT AND TIDY BENCHMARK DATA ###########################
Rb <- read.csv("https://query1.finance.yahoo.com/v7/finance/download/VAS.AX?period1=1430265600&period2=1588118400&interval=1mo&events=history")
n <- nrow(Rb)
Rb <- Rb[c(1:n-1), c(1,6)]
Rb$Date <- as.Date(Rb[, 1])
Rb <- xts(Rb$`Adj.Close`, order.by = Rb$Date)
names(Rb) <- c("Rb")
Rb$Rb <- Return.calculate(Rb$Rb, method = "log")

# IMPORT AND TIDY LIC DATA #################################
url_f <- "https://query1.finance.yahoo.com/v7/finance/download/"
url_e <- ".AX?period1=1430265600&period2=1588118400&interval=1mo&events=history"
n <- nrow(LICs)
data <- merge(Rf, Rb)
k <- nrow(data)
for(i in 1:n){
  url_temp_ch <- as.character(LICs[i,1])
  url_temp <- paste(url_f, url_temp_ch, url_e, sep = "")
  Ra_temp <- data.frame(rep(NA, k))
  try(Ra_temp <- read.csv(url_temp, na.strings = c("null")), silent = T)
  n_temp <- nrow(Ra_temp)
  try(Ra_temp <- Ra_temp[c(1:n_temp-1), c(1,6)], silent = T)
  
  if(is.na(Ra_temp[1, 1]) != TRUE){
    Ra_temp$Date <- as.Date(Ra_temp[, 1])
    Ra_temp <- xts(Ra_temp$`Adj.Close`, order.by = Ra_temp$Date)
    header <- as.character(LICs[i,1])
    names(Ra_temp) <- header
    Ra_temp[, 1] <- Return.calculate(Ra_temp[, 1], method = "log")
    data <- merge(data, Ra_temp)
    rm(Ra_temp)
  }
  else if(is.na(Ra_temp[1, 1]) == TRUE){
    data_temp <- data
    data_temp$Rf <- rep(data_temp[1, 2], k)
    data_temp <- data_temp$Rf
    header <- as.character(LICs[i,1])
    names(data_temp) <- header
    data <- merge(data, data_temp)
    rm(data_temp)
  }
}
n <- nrow(data)
data <- data[complete.cases(data[1:n, c(1, 2)]),]
LIC.list <- names(data)
names(data) <- LIC.list
n <- ncol(data)
LIC.list <- LIC.list[3:n]

# GENERATE CAPM VARIABLES ##################################
n <- ncol(data)
capm.data <- as.xts(merge(data$Rf, data$Rb-data$Rf))
names(capm.data) <- c("Rf", "mrp")
for(i in 3:n){
  Ra.Er_temp <- as.xts(data[, i]-data$Rf)
  header <- as.character(names(data))
  header <- paste(header, ".Er", sep = "")
  names(Ra.Er_temp) <- header[i]
  capm.data <- merge(capm.data, Ra.Er_temp)
}
n <- ncol(capm.data)
LICs$Code <- LIC.list

# CALCULATE CAPM PARAMETERS ################################
n <- ncol(capm.data)
capm.para <- data.frame()
for(i in 3:n){
  try(
    capm <- lm(capm.data[, i] ~ capm.data$mrp)
    , silent = T)
  para.temp <- data.frame(rep(0, 4))
  try(para.temp <- capm$coefficients, silent = T)
  para.temp <- as.data.frame(para.temp)
  para.temp <- as.data.frame(transpose(para.temp))
  try(para.temp[1, 3] <- ttest(capm, 1, 0), silent = T)
  try(para.temp[1, 4] <- ttest(capm, 2, 1), silent = T)
  names(para.temp) <- c("alpha", "beta", "alpha(0) ~ Pr(>|t|)", "beta(1) ~ Pr(>|t|)")
  row.names(para.temp) <- as.character(LICs[i-2,1])
  capm.para[i-2, 1:4] <- para.temp[1, 1:4]
  try(rm(capm), silent = T)
  rm(para.temp)
}
row.names(LICs) <- LICs$Code
LICs <- merge(LICs, capm.para, by.x = 1, by.y = 0, all.x = TRUE, all.y = TRUE)

# CHECK FOR POSITIVE AND SIGNIFICANT ALPHA RETURNS #########
LICs$alpha.rtns <- ifelse(LICs$`alpha(0) ~ Pr(>|t|)`<= 0.05 & LICs$alpha > 0.0, "TRUE", "FALSE")

# CALCULATE SD(x) FOR LICs #################################
k <- ncol(data)-2
sd.temp <- as.numeric(vector())
er.list <- names(capm.data)
n <- nrow(er.list)
er.list <- er.list[3:(k+2)]
for(i in 1:k){
  sd.temp[i] <- STDEV(capm.data[, er.list[i]])
}
sd.temp <- as.data.frame(as.numeric(sd.temp))
row.names(sd.temp) <- LIC.list
names(sd.temp) <- c("SD(ER_at)")
LICs <- merge(LICs, sd.temp, by.x = 1, by.y = 0, all.x = TRUE, all.y = TRUE)

# CALCULATE SD(x) FOR REP. PORT. ###########################
k <- nrow(data)
j <- nrow(LICs)
sd.temp <- as.numeric(vector())
for(i in 1:j){
  beta.temp <- as.data.frame(rep(LICs[i, 5], k))
  rep.port.temp <- beta.temp
  Rf.temp <- as.numeric(data$Rf)
  rep.port.temp <- add_column(rep.port.temp, Rf.temp, .after = 100)
  rtn.temp <- as.data.frame(data[, 2])
  rep.port.temp <- add_column(rep.port.temp, rtn.temp, .after = 100)
  names(rep.port.temp) <- c("Beta", "Rf", "Rtn")
  port.temp <- (1-rep.port.temp$Beta)*rep.port.temp$Rf+rep.port.temp$Beta*rep.port.temp$Rtn
  rep.port.temp <- add_column(rep.port.temp, port.temp, .after = 100)
  names(rep.port.temp) <- c("Beta", "Rf", "Rtn", "Port. Rtn")
  rep.port.temp$`Port. Rtn` <- as.numeric(unlist(rep.port.temp$`Port. Rtn`))
  rep.port.temp$Rtn <- as.numeric(unlist(rep.port.temp$Rtn))
  rep.port.temp$Exc.Port.Rtn <- as.numeric(unlist(rep.port.temp$`Port. Rtn`-rep.port.temp$Rf))
  sd.temp[i] <- STDEV(rep.port.temp[, 5])
}
LICs$"SD(ER_pt)" <- sd.temp

# COMPARE SD(x) PERFORMANCE ################################
LICs$'Lower Rep. Port. Risk?' <- ifelse(LICs$`SD(ER_pt)` <= LICs$`SD(ER_at)`, "TRUE", "FALSE")

# SAVE DATA SETS ###########################################
write.table(LICs, file = "LICs_CAPM.csv", row.names = F, sep = ",")
write.table(capm.data, file = "LICs_CAPMdata.csv", row.names = F, sep = ",")

# LIST FUNCTIONS IN FILE ###################################
list.functions.in.file("C:/Users/dtunks/OneDrive - KPMG/Desktop/PD & Personal/R data analysis/CAPM/CAPM_LICs.R", alphabetic = T) 

# FURTHER IDEAS ############################################
#CUSUM Test on alpha and beta to see structural breaks?
#Run on All Ords companies 300-500??

# CLEAN UP #################################################

# Clear environment
rm(list = ls()) 

# Clear packages
p_unload(all)  # Remove all add-ons
detach("package:datasets", unload = TRUE)  # For base

# Clear plots
dev.off()  # But only if there IS a plot

# Clear console
cat("\014")  # ctrl+L

# Clear mind :)