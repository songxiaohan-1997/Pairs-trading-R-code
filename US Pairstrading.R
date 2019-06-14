############################################################################
install.packages("Hmisc")
install.packages("tseries")
install.packages("xts")
install.packages("zoo")
install.packages("readr")
install.packages("corrplot")
library(Hmisc)#to get rcorr
library(tseries)
library(xts)
library(zoo)
library(readr)
library(corrplot)
############################################################################
##read data
dir_dta<-"C:\\Users\\Lenovo X240\\Desktop\\Finance\\Data\\AMEX\\AMEX_2009"
file_list <- list.files(path=dir_dta,full.names=T)
#tradeday number in this year 
ytd_n2009<-length(file_list)
AMEX2009<- read.table("C:\\Users\\Lenovo X240\\Desktop\\Finance\\Data\\AMEX\\AMEX_2009\\AMEX_20090101.txt",sep=",",head=TRUE)
for(i in 2:ytd_n2009)
{
  a<-read.table(file_list[i],sep=",",head=TRUE)
  AMEX2009<-rbind(AMEX2009,a)
}
####change the date format
AMEX2009$X.date.<- as.Date(as.character(AMEX2009$X.date.), "%Y%m%d")
write.csv(AMEX2009,"AMEX2009.csv")
#########################################################
dir_dta<-"C:\\Users\\Lenovo X240\\Desktop\\Finance\\Data\\AMEX\\AMEX_2010"
file_list <- list.files(path=dir_dta,full.names=T)
#tradeday number in this year 
ytd_n2010<-length(file_list)
AMEX2010<- read.table("C:\\Users\\Lenovo X240\\Desktop\\Finance\\Data\\AMEX\\AMEX_2010\\AMEX_20100101.txt",sep=",",head=TRUE)
for(i in 2:ytd_n2010)
{
  a<-read.table(file_list[i],sep=",",head=TRUE)
  AMEX2010<-rbind(AMEX2010,a)
}
####change the date format
AMEX2010$X.date.<- as.Date(as.character(AMEX2010$X.date.), "%Y%m%d")
write.csv(AMEX2010,"AMEX2010.csv")
#########################################################
dir_dta<-"C:\\Users\\Lenovo X240\\Desktop\\Finance\\Data\\AMEX\\AMEX_2011"
file_list <- list.files(path=dir_dta,full.names=T)
#tradeday number in this year 
ytd_n2011<-length(file_list)
AMEX2011<- read.table("C:\\Users\\Lenovo X240\\Desktop\\Finance\\Data\\AMEX\\AMEX_2011\\AMEX_20110103.txt",sep=",",head=TRUE)
for(i in 2:ytd_n2011)
{
  a<-read.table(file_list[i],sep=",",head=TRUE)
  AMEX2011<-rbind(AMEX2011,a)
}
####change the date format
AMEX2011$X.date.<- as.Date(as.character(AMEX2011$X.date.), "%Y%m%d")
write.csv(AMEX2011,"AMEX2011.csv")
#########################################################
dir_dta<-"C:\\Users\\Lenovo X240\\Desktop\\Finance\\Data\\AMEX\\AMEX_2012"
file_list <- list.files(path=dir_dta,full.names=T)
#tradeday number in this year 
ytd_n2012<-length(file_list)
AMEX2012<- read.table("C:\\Users\\Lenovo X240\\Desktop\\Finance\\Data\\AMEX\\AMEX_2012\\AMEX_20120102.txt",sep=",",head=TRUE)
for(i in 2:ytd_n2012)
{
  a<-read.table(file_list[i],sep=",",head=TRUE)
  AMEX2012<-rbind(AMEX2012,a)
}
####change the date format
AMEX2012$X.date.<- as.Date(as.character(AMEX2012$X.date.), "%Y%m%d")
write.csv(AMEX2012,"AMEX2012.csv")
######################################################### 
dir_dta<-"C:\\Users\\Lenovo X240\\Desktop\\Finance\\Data\\AMEX\\AMEX_2013"
file_list <- list.files(path=dir_dta,full.names=T)
#tradeday number in this year 
ytd_n2013<-length(file_list)
AMEX2013<- read.table("C:\\Users\\Lenovo X240\\Desktop\\Finance\\Data\\AMEX\\AMEX_2013\\AMEX_20130101.txt",sep=",",head=TRUE)
for(i in 2:ytd_n2013)
{
  a<-read.table(file_list[i],sep=",",head=TRUE)
  AMEX2013<-rbind(AMEX2013,a)
}
####change the date format
AMEX2013$X.date.<- as.Date(as.character(AMEX2013$X.date.), "%Y%m%d")
write.csv(AMEX2013,"AMEX2013.csv")
######################################################### 
dir_dta<-"C:\\Users\\Lenovo X240\\Desktop\\Finance\\Data\\AMEX\\AMEX_2014"
file_list <- list.files(path=dir_dta,full.names=T)
#tradeday number in this year 
ytd_n2014<-length(file_list)
AMEX2014<- read.table("C:\\Users\\Lenovo X240\\Desktop\\Finance\\Data\\AMEX\\AMEX_2014\\AMEX_20140101.txt",sep=",",head=TRUE)
for(i in 2:ytd_n2014)
{
  a<-read.table(file_list[i],sep=",",head=TRUE)
  AMEX2014<-rbind(AMEX2014,a)
}
####change the date format
AMEX2014$X.date.<- as.Date(as.character(AMEX2014$X.date.), "%Y%m%d")
write.csv(AMEX2014,"AMEX2014.csv")
######################################################### 
dir_dta<-"C:\\Users\\Lenovo X240\\Desktop\\Finance\\Data\\AMEX\\AMEX_2015"
file_list <- list.files(path=dir_dta,full.names=T)
#tradeday number in this year 
ytd_n2015<-length(file_list)
AMEX2015<- read.table("C:\\Users\\Lenovo X240\\Desktop\\Finance\\Data\\AMEX\\AMEX_2015\\AMEX_20150101.txt",sep=",",head=TRUE)
for(i in 2:ytd_n2015)
{
  a<-read.table(file_list[i],sep=",",head=TRUE)
  AMEX2015<-rbind(AMEX2015,a)
}
####change the date format
AMEX2015$X.date.<- as.Date(as.character(AMEX2015$X.date.), "%Y%m%d")
write.csv(AMEX2015,"AMEX2015.csv")
######################################################### 
dir_dta<-"C:\\Users\\Lenovo X240\\Desktop\\Finance\\Data\\AMEX\\AMEX_2016"
file_list <- list.files(path=dir_dta,full.names=T)
#tradeday number in this year 
ytd_n2016<-length(file_list)
AMEX2016<- read.table("C:\\Users\\Lenovo X240\\Desktop\\Finance\\Data\\AMEX\\AMEX_2016\\AMEX_20160101.txt",sep=",",head=TRUE)
for(i in 2:ytd_n2016)
{
  a<-read.table(file_list[i],sep=",",head=TRUE)
  AMEX2016<-rbind(AMEX2016,a)
}
####change the date format
AMEX2016$X.date.<- as.Date(as.character(AMEX2016$X.date.), "%Y%m%d")
write.csv(AMEX2016,"AMEX2016.csv")
################################################################################
################################################################################
#Pairs trading
###############################################################################
#pair
#############################################################
#correlation mathod to find pairs
data_re[is.na(data_re)]=0#target dataset:2013+2014
datapair<-rbind(AMEX2013,AMEX2014)
#record the stock code
stockcode_c<-datapair$X.ticker.
stockcode_c_unique<- stockcode_c[!duplicated(stockcode_c)]
stock_no<-length(stockcode_c_unique)
#record the trading date
tradedate_c<-datapair$X.date.
tradedate_c_unique<- tradedate_c[!duplicated(tradedate_c)]
date_no<-length(tradedate_c_unique)
###modify the data structure
#calculate the start and end date of each stock 
S_E_time<-c(stockcode=NA,starttime=NA,endtime=NA)
for(i in 1:stock_no)
{
cname<-as.character(stockcode_c_unique[i])
pointdata<-datapair[which(datapair$X.ticker.==cname),]
start<-as.Date(pointdata$X.date.[1])
end<-as.Date(pointdata$X.date.[length(pointdata$X.date.)])
unit<-c(stockcode=cname,starttime=as.character(start),endtime=as.character(end))
S_E_time<-rbind(S_E_time,unit)
}
S_E_time<-as.data.frame(S_E_time)
S_E_time<-S_E_time[-c(1),]
write.csv(S_E_time,"start_end_time.csv")
###################################################################
#find out that not all the stocks start and end on the same day
#Since we donnot need to take those stocks who end before 2014/12/31 into consideration
#So we delete them
delete<-which(S_E_time$endtime!=S_E_time$endtime[1])
stockcode_c_unique_new<-stockcode_c_unique[-delete]
datanew<-datapair[which(datapair$X.ticker.%in%stockcode_c_unique_new),]
####################################################################
#now we have a new dataset "datanew"
#rearange the table 
stock_no_new<-length(stockcode_c_unique_new)
data_re<-rep(NA,date_no)
for (i in 1:stock_no_new)
{
  pointdata<-datanew[which(datanew$X.ticker.==stockcode_c_unique_new[i]),]
  pointclose<-pointdata$X.close.
closelength<-length(pointclose)
NAlength<-date_no-closelength
  a<-rep(NA,NAlength)
    a<-c(a,pointclose)
    data_re<-cbind(data_re,a)
}  
data_re<-data_re[,-c(1)]
names(data_re)<-as.character(stockcode_c_unique_new)
row.names(data_re)<-as.character(tradedate_c_unique)
data_re_f<-as.data.frame(data_re)
data_re_f_csv<-rbind(as.character(stockcode_c_unique_new),data_re_f)
write.csv(data_re_f_csv,"modify price.csv")
################################################################################
###################
cor<- cor(data_re)
cor_ex<-cor[1:50,1:50]
#corrplot(cor_ex, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
#rcorr <- rcorr(data_re)
cor[is.na(cor)]=0
#round(cor, 6)
pair_cor<-c(stock1_no=NA,stock2_no=NA,stock1_name=NA,stock2_name=NA,corr=NA,p=NA)
for(i in 1:stock_no_new)
{
  for(j in 1:stock_no_new)
  {
    if(cor[i,j]>0.98 & i<j)
      {
      unit<-c(stock1_no=i,stock2_no=j,stock1_name=as.character(stockcode_c_unique_new[i]),stock2_name=as.character(stockcode_c_unique_new[j]),corr=cor[i,j],p=sprintf("%.10f", rcorr$P[i,j]))
      pair_cor<-rbind(pair_cor,unit)
    }
  else
  {
    pair_cor<-pair_cor
  }
  }
}
pair_cor<-as.data.frame(pair_cor[-c(1),])
##################################
#test for same distribution
pair_cor_l=length(pair_cor$stock1_no)
ks_p<-rep(NA,pair_cor_l)
for(i in 1:pair_cor_l)
{
  ks_1<-data_re[, as.numeric(as.character(pair_cor$stock1_no[i]))]
  ks_2<-data_re[,as.numeric(as.character(pair_cor$stock2_no[i]))]
  #chi_new<-chisq.test(chi_1,chi_2)we donnot use chi-square test cause the distribution maynot be normal
  ks_new<-ks.test(ks_1,ks_2)
  #p[n]<-chi_new$p.value
  ks_p[i]<-sprintf("%.10f", ks_new$p.value) 
}
pair_cor<-cbind(pair_cor,ks_p)
#####################################
#take return into consideration
return<-rep(NA,pair_cor_l)
for (i in 1:pair_cor_l)
{
  stock1<-datanew[which(datanew$X.ticker.==as.character(pair_cor$stock1_name[i])),]
  stock1close<-stock1$X.close.
  stock2<-datanew[which(datanew$X.ticker.==as.character(pair_cor$stock2_name[i])),]
  stock2close<-stock2$X.close.
return1<-(stock1close[length(stock1close)]-stock1close[1])/stock1close[1]
return2<-(stock1close[length(stock1close)]-stock1close[2])/stock1close[2]
return[i]<-(return1+return2)/2
}  
pair_cor<-cbind(pair_cor,return)
write.csv(pair_cor,"pair_cor.csv")
#################################################################
#cointegration
############################3
######################################################################
# The code below is structured in the following manner:
# Part A -
# Helper Functions viz. Signal, ParameterEstimates,
# ParameterEstimatesHistorical, Stationarity, Weights, Returns, and
# .return
# Part B -
# Load the data and perform analysis based on the defined functions
######################################################################
# Part A
######################################################################

# Decision Rule enabling trade execution
Signal <- function(spread, threshold)
{
  signal <- ifelse(spread >=  threshold, -1, NA)
  signal <- ifelse(spread <= -threshold, 1, signal)
  return(na.locf(signal))
}

# Test that log(price) follows I(1) process and
# Calculate the spread between two log stock prices.
ParameterEstimates <- function(price.pair, method = lm)
{
  x <- log(price.pair)
  
  reg <- method(x[, 2] ~ x[, 1])
  hedgeRatio <- as.numeric(reg$coef[2])
  premium     <- as.numeric(reg$coef[1])
  spread      <- x[, 2] - (hedgeRatio * x[, 1] + premium)
  list(spread = spread, hedgeRatio = hedgeRatio, premium = premium)
}
ParameterEstimatesHistorical <- function(price.pair, period, method = lm)
{
  Apply <- function(price.pair){
    reg <- ParameterEstimates(price.pair, method)
    c(spread = as.numeric(last(reg$spread)), hedgeRatio = reg$hedgeRatio, premium = reg$premium)
  }
  as.xts(rollapplyr(price.pair, period, Apply, by.column = FALSE))
}

#Returns weather spread is stationary or not
Stationarity <- function(spread, threshold)
{
  Is.passed.PP.test  <- PP.test(as.numeric(spread))$p.value <= threshold
  Is.passed.adf.test <- adf.test(as.numeric(spread))$p.value <= threshold
  c(PP.test = Is.passed.PP.test, adf.test = Is.passed.adf.test)
}

Weights <- function(hedgeRatio)
{
  hedgeRatio <- abs(hedgeRatio) * (-1)
  normalization.factor <- 1 / (1 + abs(hedgeRatio))
  return(cbind(1 * normalization.factor, hedgeRatio * normalization.factor))
}

Returns <- function(price.pair, signal.lagged, hedgeRatio.lagged)
{
  signal      <- as.xts(na.omit(cbind(signal.lagged, -1*(signal.lagged))))
  return.pair <- as.xts(na.omit(.return(price.pair, type = "discrete")))
  weight.pair <- as.xts(na.omit(Weights(hedgeRatio.lagged)))
  x <-          as.xts(apply(merge(signal[, 1], weight.pair[, 1], return.pair[, 1], all = FALSE), 1, prod))
  x <- merge(x, as.xts(apply(merge(signal[, 2], weight.pair[, 2], return.pair[, 2], all = FALSE), 1, prod)))
  
  if(!length(dim(x))){
    xts(rep(NA, nrow(price.pair)), order.by = index(price.pair))
  }else{
    xts(rowSums(x), order.by = index(x))
  }
}
.return <- function(x, type = c("continuous", "discrete"), na.pad = TRUE) 
{
  type <- match.arg(type)
  if (type == "discrete") {
    result <- x/lag(x, na.pad = na.pad) - 1
  }else if (type == "continuous") {
    result <- diff(log(x), na.pad = na.pad)
  }
  return(result)
}
##############################################
###############################################################################
#test profit
IVV_no<-as.numeric(which(stockcode_c_unique_new=="IVV"))
MGC_no<-as.numeric(which(stockcode_c_unique_new=="MGC"))
SCHX_no<-as.numeric(which(stockcode_c_unique_new=="SCHX"))
GSPC<- read_csv("C:\\Users\\Lenovo X240\\Desktop\\Finance\\1st&2nd meeting pairs trading\\R code & raw forms &plots from R\\GSPC.csv")
GSPC<-as.data.frame(lapply(GSPC, as.numeric))
ivv_prices <- as.numeric(data_re[,IVV_no])
mgc_prices <- as.numeric(data_re[,MGC_no])
schx_prices <- as.numeric(data_re[,SCHX_no])
dates <- as.Date(tradedate_c_unique)
table<-cbind(mgc_prices,ivv_prices,schx_prices)
cor_3<-as.data.frame(cor(table))
write.csv(cor_3,"cor_3.csv")
rcorr(mgc_prices,ivv_prices)
rcorr(schx_prices ,ivv_prices)
rcorr(mgc_prices,schx_prices)
##########################################################3
# Store data in data frames for prices, log prices and log returns
#*****************************
prices <- data.frame(Dates = dates,MGC = mgc_prices,IVV = ivv_prices, SCHX = schx_prices,GSPC=GSPC$GSPC)
MGCPrices <- with(prices, zoo(MGC, order.by = Dates))
IVVPrices <- with(prices, zoo(IVV, order.by = Dates))
SCHXPrices <- with(prices, zoo(SCHX, order.by = Dates))
GSPCPrices <- with(prices, zoo(GSPC, order.by = Dates))
MGCscale<-scale(MGCPrices)
IVVscale<-scale(IVVPrices)
SCHXscale<-scale(SCHXPrices)
GSPCscale<-scale(GSPCPrices)
#######sign
MGC_return_sign<-sign(lag(MGCPrices)-MGCPrices)
IVV_return_sign<-sign(lag(IVVPrices)-IVVPrices)
SCHX_return_sign<-sign(lag(SCHXPrices)-SCHXPrices)
GSPC_return_sign<-sign(lag(GSPCPrices)-GSPCPrices)
sign<-cbind(MGC_return_sign,IVV_return_sign,SCHX_return_sign,GSPC_return_sign)
write.csv(sign,"sign.csv")
cor(sign)
tot<-(MGCscale+IVVscale+SCHXscale)/3
rcorr(GSPCscale,tot)
#plot(dates,GSPCscale,type="l",main="MGC(black) & SCHX(blue) & IVV(red)")
col <- rainbow(3)
plot(dates,MGCscale,type="l",lwd=1,col=col[1],main="MGC & SCHX & IVV Scaled Price",ylab = "Price")
  lines(dates,SCHXscale,lwd=1,col=col[2])
  lines(dates,IVVscale,lwd=1,col=col[3])
  legend("bottomright", c("scaled_MGC", "scaled_SCHX","scaled_IVV"), lty = 1, col = col,cex=0.5)
  
  plot(dates,MGCPrices,type="l",lwd=1,col=col[1],main="MGC& SCHX & IVV Raw Price",ylim=c(0,230),ylab = "Price")
  lines(dates,SCHXPrices,lwd=1,,col=col[2])
  lines(dates,IVVPrices,lwd=1,col=col[3])
  legend("bottomright", c("MGC", "SCHX","IVV"), lty = 1, col = col,cex=0.5)
  
#*****************************************