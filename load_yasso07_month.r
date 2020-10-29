yasso07path<-"./"
yasso07month <- function(a, t, cl, init, inf, s, z) {
    if(!is.loaded("yasso07month")){
        print("Loading")
        l <- dyn.load(paste(yasso07path,"y07_subroutine_month.so",sep=""))
    }
    pa = .Fortran("yasso07month", a=as.double(a), t=as.double(t),cl=as.double(cl),
        init=as.double(init), inf=as.double(inf), s=as.double(s), z=as.double(z))
    return(pa$z)
} 
