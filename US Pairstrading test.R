#Pairs trading
###############################################################################
#pair
#target dataset:2015+2016
datapair_test<-rbind(AMEX2015,AMEX2016)
#record the stock code
stockcode_c_test<-datapair_test$X.ticker.
stockcode_c_test_unique<- stockcode_c_test[!duplicated(stockcode_c_test)]
stock_no_test<-length(stockcode_c_test_unique)
#record the trading date
tradedate_c_test<-datapair_test$X.date.
tradedate_c_test_unique<- tradedate_c_test[!duplicated(tradedate_c_test)]
date_no_test<-length(tradedate_c_test_unique)
###modify the data structure
#calculate the start and end date of each stock 
S_E_time_test<-c(stockcode=NA,starttime=NA,endtime=NA)
for(i in 1:stock_no_test)
{
  cname<-as.character(stockcode_c_test_unique[i])
  pointdata<-datapair_test[which(datapair_test$X.ticker.==cname),]
  start<-as.Date(pointdata$X.date.[1])
  end<-as.Date(pointdata$X.date.[length(pointdata$X.date.)])
  unit<-c(stockcode=cname,starttime=as.character(start),endtime=as.character(end))
  S_E_time_test<-rbind(S_E_time_test,unit)
}
S_E_time_test<-as.data.frame(S_E_time_test)
S_E_time_test<-S_E_time_test[-c(1),]
write.csv(S_E_time_test,"start_end_test_time.csv") 
###################################################################
#find out that not all the stocks start and end on the same day
#Since we donnot need to take those stocks who end before 2014/12/31 into consideration
#So we delete them
delete<-which(S_E_time_test$endtime!=S_E_time_test$endtime[1])
stockcode_c_test_unique_new<-stockcode_c_test_unique[-delete]
datanew_test<-datapair_test[which(datapair_test$X.ticker.%in%stockcode_c_test_unique_new),]
####################################################################
#now we have a new dataset "datanew_test"
#rearange the table 
stock_no_test_new<-length(stockcode_c_test_unique_new)
data_re_test<-rep(NA,date_no_test)
for (i in 1:stock_no_test_new)
{
  pointdata<-datanew_test[which(datanew_test$X.ticker.==stockcode_c_test_unique_new[i]),]
  pointclose<-pointdata$X.close.
  closelength<-length(pointclose)
  NAlength<-date_no_test-closelength
  a<-rep(NA,NAlength)
  a<-c(a,pointclose)
  data_re_test<-cbind(data_re_test,a)
}  
data_re_test<-data_re_test[,-c(1)]
names(data_re_test)<-as.character(stockcode_c_test_unique_new)
row.names(data_re_test)<-as.character(tradedate_c_test_unique)
data_re_test_f<-as.data.frame(data_re_test)
data_re_test_f_csv<-rbind(as.character(stockcode_c_test_unique_new),data_re_test_f)
write.csv(data_re_test_f_csv,"modify test price.csv")
###############################################################################
#test profit
IVV_no<-as.numeric(which(stockcode_c_test_unique_new=="IVV"))
MGC_no<-as.numeric(which(stockcode_c_test_unique_new=="MGC"))
SCHX_no<-as.numeric(which(stockcode_c_test_unique_new=="SCHX"))
ivv_prices <- as.numeric(data_re_test[,IVV_no])
mgc_prices <- as.numeric(data_re_test[,MGC_no])
schx_prices <- as.numeric(data_re_test[,SCHX_no])
dates <- as.Date(tradedate_c_test_unique)
rcorr(mgc_prices,ivv_prices)
rcorr(schx_prices ,ivv_prices)
rcorr(mgc_prices,schx_prices)
cor(cbind(mgc_prices,schx_prices,ivv_prices))
##########################################################3
# Store data in data frames for prices, log prices and log returns

prices <- data.frame(Dates = dates,MGC = mgc_prices,IVV = ivv_prices, SCHX = schx_prices)


# Store the data in zoo objects
MGCPrices <- with(prices, zoo(MGC, order.by = Dates))
IVVPrices <- with(prices, zoo(IVV, order.by = Dates))
SCHXPrices <- with(prices, zoo(SCHX, order.by = Dates))
# Plot Prices
plot(MGCPrices, col = 'blue', main = 'MGC Prices ')
plot(IVVPrices, col = 'brown', main = 'IVV Prices ')
plot(SCHXPrices, col = 'purple', main = 'SCHX Prices ')

# MGC-IVV pair ######################################################

pair <- cbind(MGCPrices, IVVPrices)
regress <- ParameterEstimates(pair, method = lm)
str(regress)

plot(regress$spread, main = 'MGC-IVV Spread', col = 'brown')

#check stationarity
Stationarity(regress$spread, 0.05)

#estimate parameters for back test
params <- ParameterEstimatesHistorical(pair, period = 180)

#create & plot trading signals
signal <- Signal(params$spread, 0.004)
mycol <- rgb(20, 100, 255, max = 270, alpha = 125, names = "blue50")
plot(params$spread, main = 'MGC-IVV Spread and Signal',height=600,width=500)
par(new=TRUE)
barplot(signal,col= mycol,space = 0, border = mycol,xaxt="n",yaxt="n",xlab="",ylab="")

#Performance of pair trading
plot(params$hedgeRatio, main = 'MGC-IVV 180-day Moving-Window Hedge Ratio')
returns <- Returns(pair, lag(signal), lag(params$hedgeRatio))
test_start<-which(tradedate_c_test_unique=="2015-09-22")
MGC_return<-cumsum(lag(MGCPrices)-MGCPrices)/as.numeric(MGCPrices[test_start])
IVV_return<-cumsum(lag(IVVPrices)-IVVPrices)/as.numeric(IVVPrices[test_start])
SCHX_return<-cumsum(lag(SCHXPrices)-SCHXPrices)/as.numeric(SCHXPrices[test_start])
Date<-tradedate_c_test_unique[test_start:521]
MGC_return<-MGC_return[test_start:521]
IVV_return<-IVV_return[test_start:521]
SCHX_return<-SCHX_return[test_start:521]
returncompare<-as.data.frame(cbind(cumprod(1 + returns)-1,MGC_return,IVV_return),SCHX_return)
plot(100 * cumprod(1 + returns), main = 'MGC-IVV cumulative returns (%)')
plot(Date,100 * (1+MGC_return),type="l",main = 'cumulative returns of MGC(%)')
plot(Date,100 * (1+IVV_return),type="l",main = 'cumulative returns of IVV(%)')
plot(Date,100 * (1+SCHX_return),type="l",main = 'cumulative returns of SCHX(%)')
as.numeric(cumprod(1 + returns)[length(returns)])
############################comparison plot
gspc<-read.csv("C:\\Users\\Lenovo X240\\Desktop\\stock project\\1st&2nd meeting pairs trading\\R code & raw forms &plots from R\\GSPCZ_n.csv",header = F)
gspc<-as.data.frame(lapply(gspc, as.numeric))
gspc_cr<-as.numeric(100*gspc$V1)
pair_cr<-as.numeric(100 * cumprod(1 + returns))
MGC_cr<-as.numeric(100 * (1+MGC_return))
IVV_cr<-as.numeric(100 * (1+IVV_return))
SCHX_cr<-as.numeric(100 * (1+SCHX_return))
col <- rainbow(4)
plot(Date,gspc_cr,type="l",col=col[1],main = 'Comparison between the cumulative returns(%)',ylim = c(80,130))
lines(Date,pair_cr,col=col[2],type="l")
lines(Date,MGC_cr,col=col[3],type="l")
lines(Date,IVV_cr,col=col[4],,type="l")
legend("bottomright", c("S&P index", "MGC-IVV", "MGC","IVV"), lty = 1, col = col,cex=0.5)
# SCHX-IVV pair ######################################################

pair <- cbind(SCHXPrices, IVVPrices)
regress <- ParameterEstimates(pair, method = lm)
str(regress)

plot(regress$spread, main = 'SCHX-IVV Spread', col = 'brown')

#check stationarity
Stationarity(regress$spread, 0.10)

#estimate parameters for back test
params <- ParameterEstimatesHistorical(pair, period = 180)

#create & plot trading signals
signal <- Signal(params$spread, 0.004)
mycol <- rgb(0, 0, 255, max = 255, alpha = 125, names = "blue50")
barplot(signal,col= mycol,space = 0, border = mycol,xaxt="n",yaxt="n",xlab="",ylab="")
par(new=TRUE)
plot(params$spread, main = 'SCHX-IVV Spread and Signal')

#Performance of pair trading
plot(params$hedgeRatio, main = 'SCHX-IVV 180-day Moving-Window Hedge Ratio')
returns2 <- Returns(pair, lag(signal), lag(params$hedgeRatio))
plot(100 * cumprod(1 + returns2), main = 'SCHX-IVV cumulative returns (%)')
as.numeric(cumprod(1 + returns2)[length(returns)])

# SCHX-MGC pair #######################################################

pair <- cbind(SCHXPrices, MGCPrices)
regress <- ParameterEstimates(pair, method = lm)
str(regress)

plot(regress$spread, main = 'SCHX-MGC Spread', col = 'brown')

#check stationarity
Stationarity(regress$spread, 0.1)

#estimate parameters for back test
params <- ParameterEstimatesHistorical(pair, period = 180)

#create & plot trading signals
signal <- Signal(params$spread, 0.004)
mycol <- rgb(0, 0, 255, max = 255, alpha = 125, names = "blue50")
barplot(signal,col= mycol,space = 0, border = mycol,xaxt="n",yaxt="n",xlab="",ylab="")
par(new=TRUE)
plot(params$spread, main = 'Spread and Signal')

#Performance of pair trading
plot(params$hedgeRatio, main = 'SCHX-MGC 180-day Moving-Window Hedge Ratio')
returns3 <- Returns(pair, lag(signal), lag(params$hedgeRatio))
plot(100 * cumprod(1 + returns3), main = 'SCHX-MGC cumulative returns (%)')
as.numeric(cumprod(1 + returns3)[length(returns)])

cumulative_returns <- 100 * cumprod(1 + returns3)
############################################################
gspc<-read.csv("C:\\Users\\Lenovo X240\\Desktop\\Finance\\1st&2nd meeting pairs trading\\R code & raw forms &plots from R\\GSPCZ_n.csv",header = F)
gspc<-as.data.frame(lapply(gspc, as.numeric))
gspc_cr<-as.numeric(100*gspc$V1)
pair_cr<-as.numeric(100 * cumprod(1 + returns))
MGC_cr<-as.numeric(100 * (1+MGC_return))
IVV_cr<-as.numeric(100 * (1+IVV_return))
SCHX_cr<-as.numeric(100 * (1+SCHX_return))
Returns<-100 * cumprod(1 + returns)
Returns2<-100 * cumprod(1 + returns2)
Returns3<-100 * cumprod(1 + returns3)
Returns<-Returns[23:333]
Returns2<-Returns2[23:333]
date_new<-Date[23:333]
col <- rainbow(3)
plot(date_new,Returns,type="l",col=col[1],main = 'Comparison between the cumulative returns of different pairs(%)',ylim = c(90,130))
lines(date_new,Returns2,col=col[2],type="l")
lines(date_new,Returns3,col=col[3],type="l")
legend("bottomright", c("MGC-IVV", "IVV-SCHX", "SCHX-MGC"), lty = 1, col = col,cex=0.5)




col <- rainbow(4)
plot(Date,gspc_cr,type="l",col=col[1],main = 'Comparison between the cumulative returns(%)',ylim = c(80,130))
lines(Date,pair_cr,col=col[2],type="l")
lines(Date,MGC_cr,col=col[3],type="l")
lines(Date,IVV_cr,col=col[4],,type="l")
legend("bottomright", c("S&P index", "MGC-IVV", "MGC","IVV"), lty = 1, col = col,cex=0.5)
# SCHX-IVV pair ######################################################
