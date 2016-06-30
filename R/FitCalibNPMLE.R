FitCalibNpmle <- function(w,w.res)
{
lr.for.fit <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
colnames(lr.for.fit) <- c("left","right")
lr.for.fit[lr.for.fit==Inf] <- 200
lr.for.fit[lr.for.fit==0] <- 0.0001
#df.cln <- lr.for.fit[apply(lr.for.fit,1,function(x)  sum(is.na(x)))<2,]
fit.npmple <- ic_np(cbind(left,right)~0,data = lr.for.fit)
if (inherits(fit.npmple, "error")) { 
  return(NA) 
} else {
return(fit.npmple)
}
}