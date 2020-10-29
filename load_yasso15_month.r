yasso15path<-"./"
yasso15month <- function(theta, time, climate, init, b, d, leac, xt) {
    if(!is.loaded("yasso15month")){
        print("Loading")
        l <- dyn.load(paste(yasso15path,"y15_subroutine_month.so",sep=""))
    }
    pa = .Fortran("yasso15month", theta=as.double(theta), time=as.double(time),
        climate=as.double(climate), init=as.double(init),
        b=as.double(b), d=as.double(d),
        leac=as.double(leac),xt=as.double(xt))
    return(pa$xt)
} 
