yasso15path<-"./"
yasso15 <- function(theta, time, climate, init, b, d, leac, xt) {
    if(!is.loaded("yasso15")){
        print("Loading yasso15")
        l <- dyn.load(paste(yasso15path,"y15_subroutine.so",sep=""))
    }
    pa = .Fortran("yasso15", theta=as.double(theta),time=as.double(time),
        climate=as.double(climate),init=as.double(init),
        b=as.double(b), d=as.double(d),
        leac=as.double(leac), xt=as.double(xt))
    return(pa$xt)
} 
