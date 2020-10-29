source('load_yasso15.r')
#Yasso init.sf from GAF/KHK to run Yasso to steady state
init.sf <- c(8.2359634,0.8931869,0.9298725,9.1277685,44.4027225)
#Yassom parameters
a.skan <- c( -0.5172509,-3.551512,-0.3458914,-0.2660175,0.044852223,0.0029265443,0.9779027,0.6373951,0.3124745,0.018712098,	0.022490378,0.011738963,0.00099046889,0.3361765,0.041966144,0.089885026,0.089501545,-0.0022709155,0.17,-0.0015,0.17,-0.0015,	0.17,-0.0015,0,-2.935411,0,101.8253, 260,-0.080983594,-0.315179,-0.5173524,0, 0,-0.00024180325,0.0015341907, 101.8253,260,-0.5391662,1.18574,-0.2632936,0,0,0)
a.y15 <- c(0.49,4.9,0.24,0.095,0.44,0.25,0.92,0.99,0.084,0.011,0.00061,0.00048,0.066,0.00077,0.1,0.65,-0.15,-0.02,-0.92,-0.0004,-0.00017,0.091,-0.00021,0.049,-0.000079,0.035,-0.00021,-1.8,-1.2,-13,0.0046,0.0013,-0.44,1.3,0.26)
print(length(a.skan))
print(length(a.y15))
#Amplitude or Amplitude/2 for yasso07?
amplitude<-2.0
#Try with amplitude (1) or no amplitude (0.0)
amplitude.factor<-1.0
#Mean monthly rainfall
monthly.rainfall.mean.2014<-57.98333333
monthly.rainfall.mean.2015<-57.925
#Annual rainfall sum
annual.rainfall.sum.2014<-695.8
annual.rainfall.sum.2015<-695.1
#Annual climate: Mean Temperature, Rainfall sum, Mean Temperature, Amplitude
climate.2014<-c(5.64228833,annual.rainfall.sum.2014,14.725/amplitude)
climate.2015<-c(5.949810833,annual.rainfall.sum.2015,10.325/amplitude)
#To compare monthly data with annual data
time.step<-1
#For the test start from 0
init.x0.2014<-c(0,0,0,0,0)
init.x0.2015<-c(0,0,0,0,0)
#Annual litter fall: A,W,E,N (H=0.0)
awen.annual<-c(1.82895834,0.409638388,0.236951784,0.90181951,0.0)
litter.size<-2
leaching<-0
result<-c(0,0,0,0,0)
#One year time step. First run the model into steady state, 1000 iterations should do
init.x0.2014<-yasso15(a.y15,1000,climate.2014,init.sf,
                      awen.annual,litter.size,leaching,result)
one.year.2014<-yasso15(a.y15,time.step,climate.2014,init.x0.2014,
                       awen.annual,litter.size,leaching,result)
