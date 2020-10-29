source('load_yasso15.r')
source('load_yasso15_month.r')
file.temp.2015.sinusoid<-"sinusoid.txt"
#----The sine (sinusoid) function to plot temperatures----  
Ta<-function(Tm.max,Tm.min){
    return((1.0/2.0)*(Tm.max-Tm.min))
}
temp.intra.annual<-function(t,Tm,Tm.max,Tm.min,tp){
    return (Tm + Ta(Tm.max,Tm.min)*sin(2.0*pi*t/tp))
}
#-------and the file for the sine function----------------------------
v<-c()
for (i in 1:12){
    v[i]<-temp.intra.annual(i,5.642288333,23.7,-18.5,12)
}
write(t(as.matrix(v)),file.temp.2015.sinusoid,ncolumns=12)
#-----------------------------------------------------------
#Yasso init.sf from GAF/KHK to run Yasso to steady state
init.sf <- c(8.2359634,0.8931869,0.9298725,9.1277685,44.4027225)
#Yasso07 parameters
#a.skan <- c( -0.5172509,-3.551512,-0.3458914,-0.2660175,0.044852223,0.0029265443,0.9779027,0.6373951,0.3124745,0.018712098,	0.022490378,0.011738963,0.00099046889,0.3361765,0.041966144,0.089885026,0.089501545,-0.0022709155,0.17,-0.0015,0.17,-0.0015,	0.17,-0.0015,0,-2.935411,0,101.8253, 260,-0.080983594,-0.315179,-0.5173524,0, 0,-0.00024180325,0.0015341907, 101.8253,260,-0.5391662,1.18574,-0.2632936,0,0,0)
#Yasso07 Experiment monthly: AWEN parameters are divided by 12, the precipitation parameter is multiplied by 12 (monthly
#precipitation input, yasso is not made for real monthly precipitation)
#a.skan.month<-c( -0.5172509/12.0,-3.551512/12.0,-0.3458914/12.0,-0.2660175/12.0,0.044852223,0.0029265443,0.9779027,0.6373951,0.3124745,0.018712098,	0.022490378,0.011738963,0.00099046889,0.3361765,0.041966144,0.089885026,0.089501545,-0.0022709155,0.17,-0.0015,0.17,-0.0015,	0.17,-0.0015,0,12.0*-2.935411,0,101.8253, 260,-0.080983594,-0.315179,-0.5173524,0, 0,-0.00024180325,0.0015341907, 101.8253,260,-0.5391662,1.18574,-0.2632936,0,0,0)
#Yasso15 parameter set, annual time step
a.y15 <- c(0.49,4.9,0.24,0.095,0.44,0.25,0.92,0.99,0.084,0.011,0.00061,0.00048,0.066,0.00077,0.1,0.65,-0.15,-0.02,-0.92,-0.0004,-0.00017,0.091,-0.00021,0.049,-0.000079,0.035,-0.00021,-1.8,-1.2,-13,0.0046,0.0013,-0.44,1.3,0.26)
#Yasso15 parameter set, one month time step experiment
#AWEN parameters are divided by 12, the precipitation parameters are multiplied by 12 (monthly
#precipitation input, yasso is not made for real monthly precipitation, c.f. yasso07 monthly experiment)
a.y15.month <- c(0.49/12.0,4.9/12.0,0.24/12.0,0.095/12.0,0.44,0.25,0.92,0.99,0.084,0.011,0.00061,0.00048,0.066,0.00077,0.1,0.65,-0.15,-0.02,-0.92,-0.0004,-0.00017,0.091,-0.00021,0.049,-0.000079,0.035,-0.00021,-1.8*12.0,-1.2*12.0,-13*12.0,0.0046,0.0013,-0.44,1.3,0.26)
#Mean monthly rainfall
monthly.rainfall.mean.2014<-57.98333333
monthly.rainfall.mean.2015<-57.925
#Annual rainfall sum
annual.rainfall.sum.2014<-695.8
annual.rainfall.sum.2015<-695.1
#Annual climate: Mean Temperature, Rainfall sum, Mean Temperature, Amplitude
climate.2015<-c(5.949810833,annual.rainfall.sum.2015,18.15)
#Read monthly data for litter infall and weather
awen.monthly.file<-"awen_monthly.txt"
weather.monthly.file<-"weather_monthly.txt"
awen.table<-read.table(awen.monthly.file,header=TRUE)
weather.table<-read.table(weather.monthly.file,header=TRUE)
#To compare monthly data with annual data
time.step<-1
#For the test start from 0
init.x0.2015<-c(0,0,0,0,0)
#Annual litter fall: A,W,E,N (H=0.0)
awen.annual<-c(1.82895834,0.409638388,0.236951784,0.90181951,0.0)
litter.size<-0
leaching<-0.0
result<-c(0,0,0,0,0)
#One year time step
#Uncomment this if you want the init.x0.2015 to be zero (i.e. c(0,0,0,0,0))
init.x0.2015<-yasso15(a.y15,1000,climate.2015,init.sf,awen.annual,litter.size,leaching,result)
one.year.2015<-yasso15(a.y15,time.step,climate.2015,init.x0.2015,awen.annual,litter.size,leaching,result)
#Year 2015 in 12 steps, annual data
one.month.2015.result<-matrix(0.0,12,5)
one.month.step<-yasso15(a.y15,1/12,climate.2015,init.x0.2015,awen.annual,litter.size,leaching,result)
one.month.2015.result[1,]<-one.month.step
for (i in 2:12){
    one.month.step<-yasso15(a.y15,1/12.0,climate.2015,one.month.step,awen.annual,litter.size,leaching,result)
    one.month.2015.result[i,]<-one.month.step
}
f<-"Yasso15OneMonthStep2015.txt"
write.table(as.table(one.month.2015.result),f,col.names=FALSE,row.names=FALSE)

#---------monthly try-------------------
#1. sine function, monthly 12 points from the sine function, AWEN monthly
sinusoid.temp<-read.table(file.temp.2015.sinusoid,header=FALSE)
one.month.sine.2015.matrix<-matrix(0,12,5)
#Mean Temperature, Rainfall sum, Amplitude
climate.month.sine.2015<-c(sinusoid.temp[1],weather.table[13,7],
                           weather.table[13,12])
awen.month.2015<-c(awen.table[1,2],awen.table[1,3],awen.table[1,4],awen.table[1,5],0.0)
one.month.2015<-yasso15month(a.y15.month,time.step,climate.month.sine.2015,init.x0.2015,awen.month.2015,
                             litter.size,leaching,result)
one.month.sine.2015.matrix[1,]<-one.month.2015
for (i in 2:12){
    climate.month.sine.2015<-c(sinusoid.temp[i],weather.table[13+(i-1),7],
                               weather.table[13+(i-1),12])   
    awen.month.2015<-c(awen.table[i,2],awen.table[i,3],awen.table[i,4],awen.table[i,5],0.0)
    one.month.2015<-yasso15month(a.y15.month,time.step,climate.month.sine.2015,
                                 one.month.2015,awen.month.2015,litter.size,leaching,result)
    one.month.sine.2015.matrix[i,]<-one.month.2015
}
file<-"Yasso15ResultsMonth.txt"
space<-matrix(' ',1,0)
cname<-c("A","W","E","N","H")
init.x0.2015.matrix<-matrix(0,1,5)
init.x0.2015.matrix[1,]<-init.x0.2015
colnames(init.x0.2015.matrix)<-cname
write("Year_2015",file)
write("Init",file,append=TRUE)
write.table(as.table(init.x0.2015.matrix),file,row.names=FALSE,append=TRUE)
write("Temperature Sinusoid One Month, AWEN One Month",file,append=TRUE)
colnames(one.month.sine.2015.matrix)<-cname
write.table(as.table(one.month.sine.2015.matrix),file,row.names=FALSE,append=TRUE)
write("One_Year",file,append=TRUE)
write.table(as.table(t(as.matrix(one.year.2015))),file,row.names=FALSE,col.names=FALSE,append=TRUE)
write("Difference(Year/Month)",file,append=TRUE)
write(space,file,append=TRUE)


#2. 2015: One month time step, monthly climate data, AWEN monthly
one.month.2015.matrix<-matrix(0.0,12,5)
#Mean Temperature, Rainfall sum, Amplitude
climate.month.2015<-c(weather.table[13,8],weather.table[13,7],
                      weather.table[13,12])
#Monthly litter fall: A,W,E,N (H=0.0)
awen.month.2015<-c(awen.table[1,2],awen.table[1,3],awen.table[1,4],awen.table[1,5],0.0)
one.month.2015<-yasso15month(a.y15.month,time.step,climate.month.2015,init.x0.2015,awen.month.2015,litter.size,leaching,result)
one.month.2015.matrix[1,]<-one.month.2015
for (i in 2:12){
    climate.month.2015<-c(weather.table[13+(i-1),8],weather.table[13+(i-1),7],
                          weather.table[13+(i-1),12])
    awen.month.2015<-c(awen.table[i,2],awen.table[i,3],awen.table[i,4],awen.table[i,5],0.0)
    print(awen.month.2015)
    #Previous month is the input to next month
    one.month.2015<-yasso15month(a.y15.month,time.step,climate.month.2015,one.month.2015,awen.month.2015,litter.size,leaching,result)
    one.month.2015.matrix[i,]<-one.month.2015
}
write("Climate One Month, AWEN One Month",file,append=TRUE)
colnames(one.month.2015.matrix)<-cname
write.table(as.table(one.month.2015.matrix),file,row.names=FALSE,append=TRUE)
write("One_Year",file,append=TRUE)
write.table(as.table(t(as.matrix(one.year.2015))),file,row.names=FALSE,col.names=FALSE,append=TRUE)
write("Difference(Year/Month)",file,append=TRUE)
write(space,file,append=TRUE)
print(awen.table)
