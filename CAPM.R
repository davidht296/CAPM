# File:    CAPM Estimation - ASX (5yr, Monthly)
# Author:  David Harris 
# Date:    15-04-2020

# INSTALL AND LOAD PACKAGES ################################

# Installs pacman ("package manager") if needed
if (!require("pacman")) install.packages("pacman")

# Use pacman to load add-on packages as desired
pacman::p_load(pacman, aplot, boot, car, caret, expss, GGally,
               ggthemes, ggvis, httr, huxtable, jtools, lars, 
               lmtest, lubridate, 
               MASS, olsrr, PerformanceAnalytics, plotly, plyr, 
               psych, quantmod, remotes, rio, rmarkdown, sandwich, shiny, sur, 
               tidyquant, tidyverse, xts)  

# LOAD t-TEST FUNCTION #####################################
ttest <- function(reg, coefnum, val){
  co <- coef(summary(reg))
  tstat <- (co[coefnum,1]-val)/co[coefnum,2]
  2 * pt(abs(tstat), reg$df.residual, lower.tail = FALSE)
}

# Set Working Directory ####################################
getwd()
setwd("C:/Users/dtunks/OneDrive - KPMG/Desktop/PD & Personal/R data analysis/CAPM") #*remember to change \ to /

# ASSEST LIST ##############################################
asx200 <- read.csv("https://www.asx200list.com/uploads/csv/20200401-asx200.csv", header = TRUE) ## from ASX200 list - https://www.asx200list.com/
n <- nrow(asx200)
asx200 <- asx200[c(2:n), 1:4]
asx.colnames <- c("Code", "Company", "Sector", "Market Cap")
names(asx200) <- asx.colnames
ticker <- as.character(asx200[,1])
row.names(asx200) <- ticker

# RISK-FREE RATE DATA ######################################
Rf <- import("f2.1-data.csv")            ## from RBA - https://www.rba.gov.au/statistics/tables/#interest-rates
## need to manually format dates to YYYY-MM-DD in Excel
n <- nrow(Rf)
Rf <- Rf[c(12:n), c(1, 4)]
Rf <- Rf[!apply(Rf == "", 1, all),]
Rf$V1 <- as.Date(Rf$V1)
Rf$V4 <- as.numeric(Rf$V4)
Rf$V4 <- 100*((1+(Rf$V4/100))^(1/12)-1)
Rf <- xts(Rf$V4, order.by = Rf$V1)
names(Rf) <- c("Rf")

# BENCHMARK DATA ###########################################
Rb <- read.csv("https://query1.finance.yahoo.com/v7/finance/download/VAS.AX?period1=1430092800&period2=1587945600&interval=1mo&events=history") ## from Yahoo
n <- nrow(Rb)
Rb <- Rb[c(1:n-1), c(1,6)]
Rb$Date <- as.Date(Rb[, 1])
Rb <- xts(Rb$`Adj.Close`, order.by = Rb$Date)
names(Rb) <- c("Rb")
Rb$Rb <- Return.calculate(Rb$Rb, method = "log")

# GET ASSET DATA ###########################################
url_f <- "https://query1.finance.yahoo.com/v7/finance/download/"
url_e <- ".AX?period1=1430092800&period2=1587945600&interval=1mo&events=history"
n <- nrow(asx200)
data <- merge(Rf, Rb)
for(i in 1:n){
  url_temp_ch <- as.character(asx200[i,1])
  url_temp <- paste(url_f, url_temp_ch, url_e, sep = "")
  Ra_temp <- read.csv(url_temp)
  n_temp <- nrow(Ra_temp)
  Ra_temp <- Ra_temp[c(1:n_temp-1), c(1,6)]
  Ra_temp$Date <- as.Date(Ra_temp[, 1])
  Ra_temp <- xts(Ra_temp$`Adj.Close`, order.by = Ra_temp$Date)
  header <- as.character(asx200[i,1])
  header <- paste(header, ".Ra", sep = "")
  names(Ra_temp) <- header
  Ra_temp[, 1] <- Return.calculate(Ra_temp[, 1], method = "log")
  data <- merge(data, Ra_temp)
}

# TRIM DATA ################################################
data <- data[complete.cases(data[, c(1, 2)]),]

# GENERATE CAPM VARIABLES ##################################
n <- ncol(data)
capm.data <- as.xts(data$Rb-data$Rf)
names(capm.data) <- "mrp"
for(i in 3:n){
  Ra.Er_temp <- as.xts(data[, i]-data$Rf)
  header <- as.character(asx200[i-2,1])
  header <- paste(header, ".Er", sep = "")
  names(Ra.Er_temp) <- header
  capm.data <- merge(capm.data, Ra.Er_temp)
}

# CALCULATE PARAMETERS #####################################
n <- ncol(capm.data)
capm.para <- data.frame()
for(i in 2:n){
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
  row.names(para.temp) <- as.character(asx200[i-1,1])
  capm.para[i-1, 1:4] <- para.temp[1, 1:4]
  try(rm(capm), silent = T)
  rm(para.temp)
}

asx200 <- merge(asx200, capm.para, by = 0, all.x = TRUE, all.y = TRUE)
n <- ncol(asx200)
asx200 <- asx200[, 2:n]

write.table(asx200, file = "ASX200_CAPM.csv", row.names = F, sep = ",")

# FURTHER IDEAS ############################################
CUSUM Test on alpha and beta to see structural breaks?
Run on All Ords companies 300-500??


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
