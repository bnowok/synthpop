###-----syn.ranger---------------------------------------------------------
# bagging when mtry = ncol(x) - using all predictors
syn.ranger <- function(y, x, xp, smoothing, proper = FALSE, ...) 
{ 
  #nodesize <- max(1, nodesize)  # safety
  #if (proper == TRUE) {
  #  s <- sample(length(y), replace = T); y <- y[s]
  #  x <- x[s, , drop = FALSE]
  #}  
  
  for (i in which(sapply(x,class) != sapply(xp,class))) xp[,i] <-
      eval(parse(text = paste0("as.", class(x[,i]), "(xp[,i])", sep = "")))
  
  if (is.factor(y)) {
    obslevels <- levels(y)
    y <- droplevels(y)
  }
  
  # fit a random forest
  # regression (mtry=p/3), classification (mtry=sqrt(p))
  rf.fit <- ranger(y ~ ., data = cbind.data.frame(y,x), ...)
  nodessyn <- predict(rf.fit, data = xp, type = "terminalNodes")$predictions
  nodesobs <- predict(rf.fit, data = x, type = "terminalNodes")$predictions
  
  ndonors <- vector("list", nrow(xp))
  n0      <- vector("list", ntree)
  for (j in 1:nrow(xp)) {
    for (i in 1:ntree) {
      n0[[i]] <- y[nodesobs[,i] == nodessyn[j,i]]
    }
    empty <- sapply(n0, length)
    ndonors[[j]] <- unlist(n0[empty != 0])
  }
  
  yhat <- sapply(ndonors, sample, size = 1)          
  
  if (is.factor(y)) yhat <- factor(yhat, levels = obslevels) 
  if (!is.factor(y) & smoothing == "density") yhat <- syn.smooth(yhat,y)
  
  return(list(res = yhat, fit = rf.fit))  # "rf"
}
