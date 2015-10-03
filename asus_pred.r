#clear all
rm(list=ls())


library(plyr)
library(forecast)
library(reshape)
options(stringsAsFactors = FALSE)
#setwd("F:\\Dropbox\\AcademicsD\\ADM\\project2\\mycode")
setwd("C:\\Users\\Thusitha\\Dropbox\\AcademicsD\\ADM\\project2\\mycode")

DEFAULT_DECAY <- -.05
MONTH_CUTOFF <- 07
YEAR_CUTOFF <- 2009
TIME_CONSIDERATION <- (YEAR_CUTOFF*12+MONTH_CUTOFF)

###########################################################
# Functions

# Function for fitting exponential Moving average
emamod <- function(df){
  EMA(df[,"number_repair"], 3)
}

# data processing helper 
filtval <-function(df){
  t = 0:5
  block_id_val = df[1,"block_id"]
  block_id = rep(block_id_val,6)
  dft = data.frame(t,block_id)
  subset(merge(df, dft, all.y = TRUE ), select=c("t", "block_id", "number_repair")) 
}

#redefine log to be zero at x=0
log <- function(x) ifelse(x <= 0, 0, base::log(x))

#exponential decay function
expdecay <-function(df){
  mylm = lm(log(number_repair) ~ t, data = df)
  if (mylm$coefficients["t"]>0){
    mylm$coefficients["t"] = DEFAULT_DECAY
    mylm$coefficients["Intercept"] <- log(mean(df[,"number_repair"])) - 2.5*DEFAULT_DECAY
  }
  
  d <- round(exp(predict(mylm, data.frame(t=c(6:24)))),digits=0) 
  f <- ifelse(d < 0, 0, d)
  
}

#HoltWinters model
hwdecay <- function(df){
  myts <- ts(df[,"number_repair"])
  fit <- HoltWinters(myts, gamma=FALSE)
  
  if(fit$coefficients["a"]>0){
    fit$coefficients["a"] <= -0.05
    fit$coefficients["b"] <= -0.05
  }
  d = round(data.frame(forecast(fit, 19))[,"Point.Forecast"],digits=0)
  f = ifelse(d < 0, 0, d)
  
}


#non linear model: a/(b+t) 
nlsdecay <- function(df){
  nr <- df[,"number_repair"]
  df[,"number_repair"] <- jitter(nr, factor = 1, amount = NULL)
  nls1 <- nls(number_repair~a/(b+t),start=list(a=0.1,b=1),data=df, nls.control(warnOnly=TRUE))
  nd <- data.frame(t=c(6:24))
  d <- round(predict(nls1,newdata=nd),digits=0)
}

#non linear model: (a+ct/(b+t) 
nlsdecay2 <- function(df){
  nr <- df[,"number_repair"]
  df[,"number_repair"] <- jitter(nr, factor = 1, amount = NULL)
  nls1 <- nls(number_repair~(a+c*t)/(b+t),start=list(a=0.1,b=1, c=0.01),data=df, nls.control(warnOnly=TRUE))
  nd <- data.frame(t=c(6:24))
  d <- round(predict(nls1,newdata=nd),digits=0)
}

#########################################################
#Data processing 

SaleTrain <- read.csv("data\\SaleTrain.csv")
RepairTrain <- read.csv("data\\RepairTrain.csv")
mapping <- read.csv("data\\Output_TargetID_Mapping.csv")

mapping$id <- 1:nrow(mapping)

# Change the column names
names(SaleTrain)[3] <- "year_month_sale"
names(RepairTrain)[3] <- "year_month_sale"
names(RepairTrain)[4] <- "year_month_repair"

# process the strings to obtain year and month separately, and then 
RepairTrain <- transform(RepairTrain,
                         year_repair  = as.integer(substr(year_month_repair, 1, 4)),
                         month_repair = as.integer(substr(year_month_repair, 6, 7)))



RepairTrain <- transform(RepairTrain, 
                          year_month_repair = year_repair * 12 + month_repair,
                          number_repair = pmax(number_repair, 0))

# Right now just projecting off the last six months in the experience period
repair_train <- subset(RepairTrain, year_month_repair >= TIME_CONSIDERATION)

# repair_train is at the individual repair level, roll it up to make predictions
repair_summary <- aggregate(number_repair ~ module_category + component_category +
                                    year_month_repair, repair_train, sum)
repair_summary$t <- repair_summary$year_month_repair - TIME_CONSIDERATION

# Create a block_id for each module/component combination
data_mapping <- unique(mapping[ , c("module_category", "component_category")])
data_mapping$block_id <- 1:nrow(data_mapping)
repair_summary <- merge(repair_summary, data_mapping)

#End of initial data processing work
##############################################################



#prep data for fitting the models
repair_summary <- repair_summary[with(repair_summary, order(block_id, t)), ]
repair_sum_fill <- ddply(repair_summary, .(block_id), filtval)
repair_sum_fill[is.na(repair_sum_fill)] <- 0

#write.csv(repair_summary, "rep.csv", row.names=F)


#exp decay model
pred_vals <- ddply(repair_sum_fill, .(block_id), expdecay)
block_id <- 1:224
tmp_df =data.frame(block_id)
pred_vals <- merge(pred_vals, tmp_df, all.y = TRUE )
pred_vals[is.na(pred_vals)] <- 0
pred_vals <- melt.data.frame(pred_vals, id.vars="block_id")
pred_vals <- pred_vals[with(pred_vals, order(block_id)), ]
pred_vals$id <- 1:nrow(pred_vals)
pred_vals <- subset(pred_vals, select=c("id", "value")) 
sub <- pred_vals[, c("id", "value")]
colnames(sub) <- c("id", "target")
sub <- arrange(sub, id)
write.csv(sub, "submission_exp.csv", row.names=F)

#nls model
pred_vals <- ddply(repair_sum_fill, .(block_id), nlsdecay)
block_id <- 1:224
tmp_df =data.frame(block_id)
pred_vals <- merge(pred_vals, tmp_df, all.y = TRUE )
pred_vals[is.na(pred_vals)] <- 0
pred_vals <- melt.data.frame(pred_vals, id.vars="block_id")
pred_vals <- pred_vals[with(pred_vals, order(block_id)), ]
pred_vals$id <- 1:nrow(pred_vals)
pred_vals <- subset(pred_vals, select=c("id", "value")) 
sub <- pred_vals[, c("id", "value")]
colnames(sub) <- c("id", "target")
sub <- arrange(sub, id)
write.csv(sub, "submission_nls.csv", row.names=F)

#nls model2
pred_vals <- ddply(repair_sum_fill, .(block_id), nlsdecay2)
block_id <- 1:224
tmp_df =data.frame(block_id)
pred_vals <- merge(pred_vals, tmp_df, all.y = TRUE )
pred_vals[is.na(pred_vals)] <- 0
pred_vals <- melt.data.frame(pred_vals, id.vars="block_id")
pred_vals <- pred_vals[with(pred_vals, order(block_id)), ]
pred_vals$id <- 1:nrow(pred_vals)
pred_vals <- subset(pred_vals, select=c("id", "value")) 
sub <- pred_vals[, c("id", "value")]
colnames(sub) <- c("id", "target")
sub <- arrange(sub, id)
write.csv(sub, "submission_nls2.csv", row.names=F)


#HoltWinters model
pred_vals <- ddply(repair_sum_fill, .(block_id), nlsdecay)
block_id <- 1:224
tmp_df =data.frame(block_id)
pred_vals <- merge(pred_vals, tmp_df, all.y = TRUE )
pred_vals[is.na(pred_vals)] <- 0
pred_vals <- melt.data.frame(pred_vals, id.vars="block_id")
pred_vals <- pred_vals[with(pred_vals, order(block_id)), ]
pred_vals$id <- 1:nrow(pred_vals)
pred_vals <- subset(pred_vals, select=c("id", "value")) 
sub <- pred_vals[, c("id", "value")]
colnames(sub) <- c("id", "target")
sub <- arrange(sub, id)
write.csv(sub, "submission_hw.csv", row.names=F)
