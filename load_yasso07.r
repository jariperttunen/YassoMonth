yasso07path<-"./"
yasso07 <- function(a, t, cl, init, inf, s, z) {
    if(!is.loaded("yasso07")){
        print("Loading")
        l <- dyn.load(paste(yasso07path,"y07_subroutine.so",sep=""))
    }
    pa = .Fortran("yasso07", a=as.double(a), t=as.double(t),cl=as.double(cl),
        init=as.double(init), inf=as.double(inf), s=as.double(s), z=as.double(z))
    return(pa$z)
} 
