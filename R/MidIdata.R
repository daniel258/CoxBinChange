MidIdata <- function(w, w.res, obs.tm, delta, Z)
{
  new.data <- matrix(nr = 2*length(obs.tm), nc = ncol(Z) + 5) # 5:1 for ID,  2 for (start,stop) time, one for X and one for delta
  colnames(new.data) <- c("ID", "start.time", "stop.time", "delta", "X", paste0("Z", 1:ncol(Z)))
  first.one <- apply(w.res, 1, function(x) Position(function(y) y==1, x))  # Find for each observation the first time w.res==1
  k <- 1
  # For each observation, set one row if X(t)==0 foal all observed t, and set two rows (one with X==0 and one with X==1) if a change was observed
  for (j in 1:n.sample)
  {
    if (is.na(first.one[j]))
    {
    new.data[k, ] <- c(j, 0, obs.tm[j], delta[j], 0, Z[j,]) 
    k <- k + 1
    } else 
      {
        if (first.one[j]==1) {change.point <- w[j, first.one[j]]/2} # if X(t)==1 was first observed for first questionire divide the point by two for the MidI mehtod
        else {change.point <- (w[j, first.one[j]] + w[j, first.one[j] - 1 ])/2} #  
        new.data[k, ] <- c(j, 0, change.point, 0, 0, Z[j, ])
        new.data[k + 1, ] <- c(j, change.point, obs.tm[j], delta[j], 1, Z[j, ])
        k <- k + 2
      }
  }
  new.data <- as.data.frame(new.data[1:(k - 1), ]) # k + 1 was the last row, but then we added +2 to k
  return(new.data)
}
