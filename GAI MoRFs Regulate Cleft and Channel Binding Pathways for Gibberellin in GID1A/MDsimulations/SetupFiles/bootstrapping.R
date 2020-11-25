
#------------------- READ RAMD DATA
getdata <- function(name) 
{
	s <- read.table(name,header = FALSE)
	x1 <- s[["V3"]]
     return (x1/500000.0)
}
#---------------------- GET RAMD effective dissociation TIME from RAMD data
gettime <- function(xtot) { #this is median calc, why not just return median????
      htot_norm <- sort(xtot, decreasing = TRUE)
      k<- 1
      for (i in seq_along(htot_norm)) #this just gets half the length conservatively
      {
           if (k > length (htot_norm)/2.0) break#so if uneven list....will avg actual median by position+1
           k <- k+1 #else
      }
      restime <- (htot_norm[k-1]+htot_norm[k+1])/2.0 #so just next of median divided by 2?
      return (restime)# i think this is just the freaking median after ranking?
}
#-------------------- BOOTSTRAPPING DATA
getbootstr  <- function(x) 
{
    data <- gettime(x)
    for (i in 1:100000) 
    {
            #y <-sample(x, length(x)*0.025,replace=T)
        y <-sample(x, length(x)*0.8,replace=T)
    		data <- c(data,gettime(y))
    }
    return (data)
}
# runtime
filelist <- list.files(path=".",pattern=paste("times","_",sep=""),full.names = TRUE)
out <- "FILE to be analyzed:"
out
filelist
list_out <- NULL
#pdfx=matrix()
#pdfy=matrix()
for (k in seq_along(filelist)) 
{
    data_inp <- getdata(filelist[k])
    x1=getbootstr(data_inp)
    out_distr=paste(MEAN=c(mean(x1),SD=sd(x1)))
    list_out <- c(list_out,mean(x1),sd(x1))
    h <- hist(x1, plot=FALSE)
    plot(h, col="grey") #plot hist
    xfit<-seq(min(x1),max(x1),length=400) 
    yfit<-dnorm(xfit,mean=mean(x1),sd=sd(x1)) 
    yfit <- yfit*diff(h$mids[1:2])*length(x1) 
    #pdfx[,k]=xfit
    #pdfy[,k]=yfit
    # ylim=c(0,1.5)
    # lines(xfit, yfit, col="blue", lwd=2)
    lines(xfit, yfit,  col="blue", lwd=2)
    # xlines <-seq(min(h$breaks),max(h$breaks),length=100) #seq of x for pdf
    # lines(x = xlines,y=dgamma(xlines,2,1)*length(x1)*diff(h$breaks)[1])
    #write(paste(h$breaks,h$density),paste("distribution.dat", toString(k),sep=""),ncolumns = 1, sep = "\t")
    write(paste(xfit,yfit),paste("distribution.dat", toString(k),sep=""),ncolumns = 1, sep = "\t")
    }
# plot(pdfx, pdfy, lwd=2)
out <- " MEAN /ns   SD "
out
list_out

#------------------
# data_inp <- times
# x1=getbootstr(data_inp)
# out_distr=paste(MEAN=c(mean(x1),SD=sd(x1)))
# list_out <- c(list_out,mean(x1),sd(x1))
# h <- hist(x1)
# write(paste(h$breaks,h$density),paste("distribution.dat", toString(k),sep=""),ncolumns = 1, sep = "\t")
# out <- " MEAN /ns   SD "
# out
# list_out
