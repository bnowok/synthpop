# Functions for synthesising data, some of which are adapted
# from mice package by S. van Buuren and K. Groothuis-Oudshoorn,
# TNO Quality of Life


###-----.norm.fix.syn------------------------------------------------------

.norm.fix.syn <- function(y, x, ridge=0.00001, ...)
{
# Calculates regression coefficients + error estimate

  xtx <- t(x) %*% x
  pen <- ridge * diag(xtx)
  if (length(pen)==1) pen <- matrix(pen)
  v           <- solve(xtx + diag(pen))
  coef        <- t(y %*% x %*% v)
  residuals   <- y - x %*% coef
  sigma       <- sqrt((sum(residuals^2))/(length(y)-ncol(x)-1))
  parm        <- list(coef, sigma)
  names(parm) <- c("beta","sigma")
  return(parm)
}


###-----.norm.draw.syn-----------------------------------------------------

.norm.draw.syn <- function(y, x, ridge=0.00001, ...)
{
# Draws values of beta and sigma for Bayesian linear regression synthesis 
# of y given x according to Rubin p.167

  xtx <- t(x) %*% x
  pen <- ridge * diag(xtx)
  if (length(pen)==1) pen <- matrix(pen)
  v           <- solve(xtx + diag(pen))
  coef        <- t(y %*% x %*% v)
  residuals   <- y - x %*% coef
  sigma.star  <- sqrt(sum((residuals)^2)/rchisq(1, length(y) - ncol(x)))
  beta.star   <- coef + (t(chol((v + t(v))/2)) %*% rnorm(ncol(x))) * sigma.star
  parm        <- list(coef, beta.star, sigma.star)      
  names(parm) <- c("coef","beta","sigma")      
  return(parm)
}


###-----syn.norm-----------------------------------------------------------

syn.norm <- function(y, x, xp, proper = FALSE, ...)
{
  x   <- cbind(1, as.matrix(x))
  xp  <- cbind(1, as.matrix(xp))
  if (proper==FALSE){
    parm <- .norm.fix.syn(y, x, ...)
  } else {
    parm <- .norm.draw.syn(y, x, ...)
  }  
  res <- xp %*% parm$beta + rnorm(nrow(xp)) * parm$sigma
  res <- round(res,max(sapply(y,decimalplaces)))
  
  return(list(res = res, fit = parm))
}


###-----syn.lognorm--------------------------------------------------------

syn.lognorm <- function (y, x, xp, proper = FALSE, ...) 
{
  addbit <- FALSE
  if (any(y < 0)) stop("Log transformation not appropriate for negative values.\n")
  if (any(y == 0)) {y <- y + .5*min(y[y!=0]); y <- log(y); addbit <- TRUE}  ##  warning about this and above should be in check model
  else y <- log(y)
  x   <- cbind(1, as.matrix(x))
  xp  <- cbind(1, as.matrix(xp))
  if (proper==FALSE){
    parm <- .norm.fix.syn(y, x, ...)
  } else {
    parm <- .norm.draw.syn(y, x, ...)
  }
  res <- xp %*% parm$beta + rnorm(nrow(xp)) * parm$sigma
  if (addbit) {res <- res-.5*min(y[y!=0]); res[res<=0] <- 0}
  res <- exp(res)
  res <- round(res,max(sapply(y,decimalplaces)))
  
  return(list(res = res, fit = parm))
}


###-----syn.sqrtnorm-------------------------------------------------------

syn.sqrtnorm <- function (y, x, xp, proper = FALSE, ...) 
{
  addbit <- FALSE
  if (any(y < 0)) stop("Square root transformation not appropriate for negative values.\n")   ##  needs check in checkmodel
  else y <- sqrt(y)
  x   <- cbind(1, as.matrix(x))
  xp  <- cbind(1, as.matrix(xp))
  if (proper==FALSE){
    parm <- .norm.fix.syn(y, x, ...)
  } else {
    parm <- .norm.draw.syn(y, x, ...)
  }
  res <- xp %*% parm$beta + rnorm(nrow(xp)) * parm$sigma
  res <- res^2
  res <- round(res,max(sapply(y,decimalplaces)))
  
  return(list(res = res, fit = parm))
}


###-----syn.cubertnorm-----------------------------------------------------

syn.cubertnorm <- function (y, x, xp, proper = FALSE, ...) 
{
  addbit <- FALSE
  if (any(y < 0)) stop("Cube root transformation not appropriate for negative values.\n")   ##  needs check in checkmodel
  else y <- y^(1/3)
  x   <- cbind(1, as.matrix(x))
  xp  <- cbind(1, as.matrix(xp))
  if (proper==FALSE){
    parm <- .norm.fix.syn(y, x, ...)
  } else {
    parm <- .norm.draw.syn(y, x, ...)
  }
  res <- xp %*% parm$beta + rnorm(nrow(xp)) * parm$sigma
  res <- res^3
  res <- round(res,max(sapply(y,decimalplaces)))

  return(list(res = res, fit = parm))
}


###-----syn.normrank-------------------------------------------------------

syn.normrank <- function(y, x, xp, smoothing, proper = FALSE, ...)
{
  # Regression synthesis of y given x, with a fixed regression
  # line, and with random draws of the residuals around the line.
  # Adapted from norm by carrying out regression on Z scores from ranks
  # predicting new Z scores and then transforming back
  # similar to method by ? and ?
  #
  # First get approx rank position of vector in one of another length
  # so that result returned has correct length for xp
  # matters for sub-samples and missing data

  z    <- qnorm(rank(y)/(length(y)+1))

  x    <- cbind(1, as.matrix(x))
  xp   <- cbind(1, as.matrix(xp))

  if (proper==FALSE){
    parm <- .norm.fix.syn(z, x, ...)
  } else {
    parm <- .norm.draw.syn(z, x, ...)
  }
  
  pred <- (xp %*% parm$beta + rnorm(nrow(xp)) * parm$sigma)
  res  <- round(pnorm(pred)*(length(y)+1))
  res[res<1] <- 1
  res[res>length(y)] <- length(y)

  if (smoothing=="") res  <- sort(y)[res]

  if (smoothing=="density"){
    ydsamp <- y
    ys     <- 1:length(y)
    maxfreq <- which.max(table(y))
    maxcat  <- as.numeric(names(table(y))[maxfreq])
    if (table(y)[maxfreq]/sum(table(y))>.7) ys <- which(y!=maxcat)
    if (10*table(y)[length(table(y))-1] < 
      tail(table(y),n=1)-table(y)[length(table(y))-1]){
      ys <- ys[-which(y==max(y))]  
      maxy <- max(y)
    }   
    densbw <- density(y[ys],width="SJ")$bw
    ydsamp[ys] <- rnorm(length(ydsamp[ys]),
                    mean=sample(ydsamp[ys],length(ydsamp[ys]),replace=TRUE),
                    sd=densbw)
    if (!exists("maxy")) maxy <- max(y) + densbw
    ydsamp[ys] <- pmax(pmin(ydsamp[ys],maxy),min(y))
    res <- sort(ydsamp)[res]
  }

  return(list(res = res, fit = parm))
}


###-----syn.ranknorm-------------------------------------------------------
# check if this function is needed
syn.ranknorm <- function(y, x, xp, proper = FALSE, ...)
{
# Regression synthesis of y given x, with a fixed regression
# line, and with random draws of the residuals around the line.
# Adapted from norm to get ranks in real data
#
# First get approx rank position of vector in one of another length
# so that result returned has correct length for xp 
# matters for sub-samples and missing data

  newrank <- function(pred,oldn) {
    newn <- length(pred)
    res  <- round((1:newn)*(oldn/newn))
    res[res<1] <- 1
    res[res>oldn] <- oldn
    res[rank(pred)]
  }

  x  <- cbind(1, as.matrix(x))
  xp <- cbind(1, as.matrix(xp))

  if (proper==FALSE){
    parm <- .norm.fix.syn(y, x, ...)
  } else {
    parm <- .norm.draw.syn(y, x, ...)
  }
  
  pred <- (xp %*% parm$beta + rnorm(nrow(xp)) * parm$sigma)
  if (nrow(x)== nrow(xp)) rankpred <- rank(pred)
  else rankpred <- newrank(pred,length(y))
  res  <- sort(sample(y,replace=TRUE))[rankpred]  ##  note bootstrap added
  
  return(list(res = res, fit = parm))
}


###-----.pmm.match---------------------------------------------------------

.pmm.match <- function(z, yhat=yhat, y=y, donors=3, ...)
{
# Auxilary function for syn.pmm.
# z    = target predicted value (scalar)
# yhat = array of fitted values, to be matched against z
# y    = array of donor data values

# Finds the three cases for which abs(yhat-z) is minimal,
# and makes a random draw from these.

  d <- abs(yhat-z)
  m <- sample(y[rank(d, ties.method="random") <= donors], 1)
  return(m)
}


###-----syn.pmm------------------------------------------------------------

syn.pmm <- function (y, x, xp, proper = FALSE, ...)
{
# Synthesis of y by predictive mean matching
# Warning: can be slow for large data sets 
# for which syn.normrank may be a better choice
  x       <- cbind(1, as.matrix(x))
  xp      <- cbind(1, as.matrix(xp))
  if (proper==FALSE){
    parm <- .norm.fix.syn(y, x, ...)
  } else {
    parm <- .norm.draw.syn(y, x, ...)
  }
  yhatobs <- x  %*% parm$coef
  yhatmis <- xp %*% parm$beta
  res <- apply(as.array(yhatmis), 1, .pmm.match, yhat=yhatobs, y=y, ...)
  
  return(list(res = res, fit = parm))
}


###-----augment.syn--------------------------------------------------------

augment.syn <- function(y, x, ...)
{
  # define augmented data for stabilizing logreg and polyreg
  # by the ad hoc procedure of White, Daniel & Royston, CSDA, 2010
  # This function will prevent augmented data beyond the min and
  # the max of the data
  # Input:
  # x: numeric data.frame (n rows)
  # y: factor or numeric vector (length n)
  # Output:
  # return a list with elements y, x, and w with length n+2*(ncol(x))*length(levels(y))
  
  x    <- as.data.frame(x)
  icod <- sort(unique(unclass(y)))
  ki   <- length(icod)
  # if (ki>maxcat) stop(paste("Maximum number of categories (",maxcat,") exceeded", sep=""))
  p    <- ncol(x)

  # skip augmentation if there are no predictors
  if (p==0) return(list(y=y,x=x,w=rep(1,length(y))))
  
  # skip augmentation if there is only 1 missing value  
  if (length(y)==1) return(list(y=y,x=x,w=rep(1,length(y))))
    
  # calculate values to augment
  mean <- apply(x,2,mean)
  sd   <- sqrt(apply(x,2,var))
  minx <- apply(x,2,min)
  maxx <- apply(x,2,max)
  nr   <- 2 * p * ki
  a    <- matrix(mean, nrow=nr, ncol=p, byrow=TRUE)
  b    <- matrix(rep(c(rep(c(0.5,-0.5),ki),rep(0,nr)),length=nr*p), nrow=nr, ncol=p, byrow=FALSE)
  c    <- matrix(sd, nrow=nr, ncol=p, byrow=TRUE)
  d    <- a + b * c
  d    <- pmax(matrix(minx,nrow=nr,ncol=p, byrow=TRUE), d)
  d    <- pmin(matrix(maxx,nrow=nr,ncol=p, byrow=TRUE), d)
  e    <- rep(rep(icod, each=2), p)
  dimnames(d) <- list(paste("AUG",1:nrow(d),sep=""),dimnames(x)[[2]])
  xa   <- rbind.data.frame(x, d)
  # beware, concatenation of factors
  # this change needed to avoid reordering of factors                           
  # if (is.factor(y)) ya <- as.factor(levels(y)[c(y,e)]) else ya  <- c(y, e)
  if (is.factor(y)) ya <- addNA(factor(levels(y)[c(y,e)],levels=levels(y)),ifany=TRUE) else ya <- c(y,e)   
  wa <- c(rep(1,length(y)),rep((p+1)/nr,nr))

  return(list(y=ya, x=xa, w=wa))
}


###-----syn.logreg---------------------------------------------------------
                                                    
syn.logreg <- function(y, x, xp, denom = NULL, denomp = NULL, 
                       proper = FALSE, ...)            
{
  # Synthesis for binary or binomial response variables by
  # logistic regression model. See Rubin (1987, p. 169-170) for
  # a description of the method.
  

  # The method consists of the following steps:
  # 1. Fit a logit, and find (bhat, V(bhat))
  # 2. Draw BETA from N(bhat, V(bhat))
  # 3. Compute predicted scores for m.d., i.e. logit-1(X BETA)
  # 4. Compare the score to a random (0,1) deviate, and synthesise.

  #browser()
  if (is.null(denom)) {
    aug <- augment.syn(y, x, ...)
    # when no missing data must set xf to augmented version
    xf   <- aug$x
    y    <- aug$y
    w    <- aug$w
    xf   <- cbind(1, as.matrix(xf))
    xp   <- cbind(1, as.matrix(xp))
    expr <- expression(glm.fit(xf, y, family = binomial(link = logit), weights=w))
    fit  <- suppressWarnings(eval(expr))
    fit.sum <- summary.glm(fit)
    beta <- coef(fit)
    if (proper==TRUE){
      rv   <- t(chol(fit.sum$cov.unscaled))
      beta <- beta + rv %*% rnorm(ncol(rv))  
    }
    p   <- 1/(1 + exp(-(xp %*% beta)))  
    vec <- (runif(nrow(p))<= p)
    if (!is.logical(y)) vec <- as.numeric(vec)          
    if (is.factor(y)) vec <- factor(vec,c(0,1),labels=levels(y))
  
  } else {
    aug <- augment.syn(y, x, ...)
    # when no missing data must set xf to augmented version
    xf   <- aug$x
    y    <- aug$y
    w    <- aug$w
    xf   <- cbind(1, as.matrix(xf))
    xp   <- cbind(1, as.matrix(xp))
    den  <- w
    denind <- which(den==1)
    den[denind] <- denom
    yy   <- y/den        #GR denom give then average response
    yy[den<1]   <- mean(yy[denind]) 
    expr <- expression(glm.fit(xf, yy, family = binomial(link = logit), weights=den))
    fit  <- suppressWarnings(eval(expr))
    fit.sum <- summary.glm(fit)
    beta <- coef(fit.sum)[,1]
    if (proper==TRUE){
      rv   <- t(chol(fit.sum$cov.unscaled))
      beta <- beta + rv %*% rnorm(ncol(rv))  
    }
    p <- 1/(1 + exp(-(xp %*% beta)))  
    vec <- rbinom(nrow(p),denomp, p) 
  }
  return(list(res = vec, fit = fit$coefficients))  # fit = "logreg"
}


###-----syn.polyreg--------------------------------------------------------   
    
syn.polyreg <- function(y, x, xp, proper = FALSE, maxit = 1000, 
                        trace = FALSE, MaxNWts = 10000, ...)
{
# synthesis for categorical response variables by the Bayesian
# polytomous regression model. See J.P.L. Brand (1999), Chapter 4,
# Appendix B.
#
# The method consists of the following steps:
# 1. Fit categorical response as a multinomial model
# 2. Compute predicted categories
# 3. Add appropriate noise to predictions.
#
# This algorithm uses the function multinom from the libraries nnet and MASS
# (Venables and Ripley).
  
  x   <- as.matrix(x)
  xp  <- as.matrix(xp)
  
  if (proper==TRUE){ # bootstrap to make proper
    s   <- sample(length(y), replace=TRUE)
    x   <- x[s, , drop = FALSE]
    y   <- y[s]  
    y   <- factor(y)
  }
  
  aug <- augment.syn(y, x, ...)
  # yf and xf needed for augmented data to save x as non augmented  not now needed can tidy
  xf  <- aug$x
  yf  <- aug$y
  w   <- aug$w
  xfy <- cbind.data.frame(yf, xf)  
  fit <- multinom(formula(xfy), data = xfy, weights = w,
  maxit = maxit, trace = trace, MaxNWts = MaxNWts, ...)
  if(fit$convergence == 1) cat(" <- Reached max number of iterations for a multinomial model suggest rerunning with polyreg.maxit increased (default 1000)\n")             
  post <- predict(fit, xp, type = "probs")   
  if (length(y)==1) post <- matrix(post, nrow=1, ncol=length(post)) 
  if (!is.factor(y)) y <- as.factor(y)
  nc <- length(levels(yf))                    
  un <- rep(runif(nrow(xp)),each=nc)
  if (is.vector(post)) post <- matrix(c(1-post,post),ncol=2)
  draws <- un > apply(post,1,cumsum)
  idx   <- 1 + apply(draws,2,sum)
   
  res <- levels(yf)[idx]
  return(list(res = res, fit = fit$coefficients))   # fit = "polyreg"
}


###-----syn.polr-----------------------------------------------------------

syn.polr <- function(y, x, xp, proper = FALSE, maxit = 1000,
                     trace = FALSE, MaxNWts = 10000, ...)
{
  x   <- as.matrix(x)
  xp  <- as.matrix(xp)
  
  if (proper==TRUE){  # bootstrap to make proper
    s   <- sample(length(y), replace=TRUE)
    x   <- x[s,]
    y   <- y[s]
    y   <- factor(y)
  }
  
  aug <- augment.syn(y, x, ...)
  # yf, wf and xf needed for augmented data to save x as non augmented  GR
  xf  <- aug$x
  yf  <- aug$y
  wf  <- aug$w
  #xy  <- cbind.data.frame(y = y,  x = xp)
  xfy <- cbind.data.frame(yf, xf)

  ## polr may fail on sparse data. We revert to multinom in such cases. 
  fit <- try(suppressWarnings(polr(formula(xfy), data = xfy, weights=wf, ...)), silent = TRUE)
  if (inherits(fit, "try-error")) {
    fit <- multinom(formula(xfy), data = xfy, weights = wf,
                    maxit = maxit, trace = trace, MaxNWts = MaxNWts, ...)
    cat("\tMethod changed to multinomial")
    if(fit$convergence==1) cat("\nReached max number of iterations for a multinomial model\nRerun with polyreg.maxit increased (default 100)\n")
  }
  post  <- predict(fit, xp, type = "probs")
  if (length(y) == 1) post <- matrix(post, nrow = 1, ncol = length(post))
  y     <- as.factor(y)
  nc    <- length(levels(yf))                       
  un    <- rep(runif(nrow(xp)), each = nc)
  if (is.vector(post)) post <- matrix(c(1 - post, post), ncol = 2)
  draws <- un > apply(post, 1, cumsum)
  idx   <- 1 + apply(draws, 2, sum)
# this slightly clumsy code needed to ensure y retains its labels and levels
#  y[1:length(y)]<-(levels(y)[idx])

  res <- levels(yf)[idx]

  return(list(res = res, fit = fit$coefficients))  # fit = "polr" 
}


###-----syn.sample---------------------------------------------------

syn.sample <- function(y, xp, smoothing, cont.na, proper = FALSE, ...) 
{
  # Generates random sample from the observed y's
  # with bootstrap if proper == TRUE
  if (proper==TRUE) y <- sample(y,replace=TRUE)
  yp <- sample(y, size=xp, replace=TRUE)
  
  if (smoothing=="density") yp[!(yp %in% cont.na)] <- 
    syn.smooth(yp[!(yp %in% cont.na)],y[!(y %in% cont.na)])
  
  return(list(res = yp, fit = "sample"))
}


###-----syn.passive--------------------------------------------------------
syn.passive <- function(data, func)
{
  # Special elementary synthesis method for transformed data.
  res <- suppressWarnings(model.frame(as.formula(func), data, na.action=na.pass))	#BN 25/08 added suppressWarnings to avoid NAs by coersion for NAtemp
 
  return(list(res = res, fit = "passive"))
}


###-----syn.cart-----------------------------------------------------------

syn.cart <- function(y, x, xp, smoothing, proper = FALSE, 
                     minbucket = 5, cp = 1e-04, ...)
{
  if (proper==TRUE){
    s <- sample(length(y), replace=TRUE)
    x <- x[s,,drop=FALSE]
    y <- y[s]
  }
  
  #for (j in 1:ncol(x)){
  #  if(is.factor(x[,j])) { 
  #    attributes(x[,j])$contrasts <- NULL
  #    attributes(xp[,j])$contrasts <- NULL
  #  }
  #}
  minbucket <- max(1, minbucket)  # safety
  if (!is.factor(y) & !is.logical(y)) {
    fit <- rpart(y ~ ., data = as.data.frame(cbind(y, x)), method = "anova",
                 minbucket = minbucket, cp = cp, ...)
    # get leaf number for observed data
    leafnr  <- floor(as.numeric(row.names(fit$frame[fit$where,])))
    # replace yval with leaf number in order to predict later node number 
    # rather than yval (mean y for observations classified to a leaf) 
    fit$frame$yval <- as.numeric(row.names(fit$frame))
    # predict leaf number
    nodes       <- predict(object=fit, newdata=xp)
    uniquenodes <- unique(nodes)
    new  <- vector("numeric",nrow(xp))
    for(j in uniquenodes){
      donors        <- y[leafnr==j] # values of y in a leaf
      new[nodes==j] <- resample(donors,size=sum(nodes==j),replace=T)
    }
    if (smoothing=="density") new <- syn.smooth(new, y) 
    #donor <- lapply(nodes, function(s) y[leafnr == s])
    #new   <- sapply(1:length(donor),function(s) resample(donor[[s]], 1))
  } else {
    fit   <- rpart(y ~ ., data = as.data.frame(cbind(y, x)), method = "class",
                   minbucket = minbucket, cp = cp, ...)
    nodes <- predict(object = fit, newdata = xp)
    new   <- apply(nodes, MARGIN=1, FUN=function(s) resample(colnames(nodes),size=1,prob=s))
    if (is.logical(y)) {
      new   <- as.logical(new)
    } else {
      new   <- factor(new,levels=levels(y)) 
    }
  }
  
  return(list(res = new, fit = fit))
}


###-----syn.ctree----------------------------------------------------------

syn.ctree <- function(y, x, xp, smoothing, proper = FALSE, minbucket = 5, ... )
{ 
  if (proper==TRUE){
    s <- sample(length(y),replace=T)
    y <- y[s]
    x <- x[s,,drop=FALSE]
  }
  
  for (i in which(sapply(x,class)!=sapply(xp,class))) xp[,i] <-
  eval(parse(text=paste0("as.",class(x[,i]),"(xp[,i])",sep="")))
  # Fit a tree
  datact     <- ctree(y ~ ., data=as.data.frame(cbind(y,x)), 
                   controls=ctree_control(minbucket=minbucket, ...))
  fit.nodes  <- where(datact)
  nodes      <- unique(fit.nodes)
  no.nodes   <- length(nodes)
  pred.nodes <- where(datact,newdata=xp)
  # Get row numbers for predicted by sampling
  # with replacement from existing data
  rowno      <- 1:length(y)
  newrowno   <- vector("integer",nrow(xp))

  for (i in nodes){
    newrowno[pred.nodes==i] <- sample(rowno[fit.nodes==i],
                                      length(newrowno[pred.nodes==i]),
                                      replace=T)
  }
  new <- y[newrowno]
  if (!is.factor(y) & smoothing=="density") new <- syn.smooth(new,y)
  
  return(list(res = new, fit = datact))
}


###-----syn.survctree------------------------------------------------------

syn.survctree <- function(y, yevent, x, xp, proper = FALSE, minbucket = 5, ...)
# time, event - data column numbers
{      
  if (proper==TRUE){
    s <- sample(length(y),replace=T)
    y <- y[s]
    x <- x[s,,drop=FALSE]                        
    yevent <- yevent[s]
  }
  browser()
  # Fit a tree
  datact     <- ctree(Surv(y,yevent) ~ ., data=as.data.frame(cbind(y,yevent,x)),
                      controls=ctree_control(minbucket = minbucket, ...))
  fit.nodes  <- where(datact)
  nodes      <- unique(fit.nodes)
  no.nodes   <- length(nodes)
  pred.nodes <- where(datact,newdata=xp)
  # Get row numbers for predicted by sampling
  # with replacement from existing data
  rowno      <- 1:length(y)
  newrowno   <- rep(0,nrow(xp))
  for (i in nodes){
    newrowno[pred.nodes==i] <- sample(rowno[fit.nodes==i],
    length(newrowno[pred.nodes==i]),replace=T)
   }
  #Predicte node & sample time+event
  faketime  <- y[newrowno]
  fakeevent <- yevent[newrowno]

  return(list(syn.time = faketime, syn.event = fakeevent, fit = datact))
}


###-----syn.rf-------------------------------------------------------------
# bagging when mtry = ncol(x) - using all predictors
syn.rf <- function(y, x, xp, smoothing, proper = FALSE, ntree = 10, ...) 
{ 
  #nodesize <- max(1, nodesize)  # safety
  #if (proper == TRUE) {
  #  s <- sample(length(y), replace = T); y <- y[s]
  #  x <- x[s, , drop = FALSE]
  #}  

  for (i in which(sapply(x,class)!=sapply(xp,class))) xp[,i] <-
  eval(parse(text=paste0("as.",class(x[,i]),"(xp[,i])",sep="")))
  
  if (is.factor(y)){
    obslevels <- levels(y)
    y <- droplevels(y)
  }

  # fit a random forest
  # regression (mtry=p/3), classification (mtry=sqrt(p))
  rf.fit <- randomForest(y ~ ., data = cbind.data.frame(y,x), ntree = ntree, ...)
  nodessyn <- attr(predict(rf.fit, newdata = xp, nodes = T), "nodes")
  nodesobs <- attr(predict(rf.fit, newdata = x, nodes = T), "nodes")

  ndonors <- vector("list", nrow(xp))
  n0      <- vector("list", ntree)
  for (j in 1:nrow(xp)){
    for (i in 1:ntree) {
      n0[[i]] <- y[nodesobs[,i]==nodessyn[j,i]]
    }
    empty <- sapply(n0, length)
    ndonors[[j]] <- unlist(n0[empty!=0])
  }
  
  yhat <- sapply(ndonors, sample, size=1)          
  
  if (is.factor(y)) yhat <- factor(yhat, levels = obslevels) 
  if (!is.factor(y) & smoothing=="density") yhat <- syn.smooth(yhat,y)
    
  return(list(res = yhat, fit = rf.fit))  # "rf"
}


###-----syn.bag-------------------------------------------------------------
# bagging when mtry = ncol(x) - using all predictors
syn.bag <- function(y, x, xp, smoothing, proper = FALSE, ntree = 10, ...) 
{ 
  #nodesize <- max(1, nodesize)  # safety
  #if (proper == TRUE) {
  #  s <- sample(length(y), replace = T); y <- y[s]
  #  x <- x[s, , drop = FALSE]
  #}  

  for (i in which(sapply(x,class)!=sapply(xp,class))) xp[,i] <-
  eval(parse(text=paste0("as.",class(x[,i]),"(xp[,i])",sep="")))
  
  if (is.factor(y)){
    obslevels <- levels(y)
    y <- droplevels(y)
  }

  # fit a random forest
  # regression (mtry=p/3), classification (mtry=sqrt(p))
  rf.fit <- randomForest(y ~ ., data = cbind.data.frame(y,x), ntree = ntree, mtry = ncol(x), ...)
  nodessyn <- attr(predict(rf.fit, newdata = xp, nodes = T), "nodes")
  nodesobs <- attr(predict(rf.fit, newdata = x, nodes = T), "nodes")

  ndonors <- vector("list", nrow(xp))
  n0      <- vector("list", ntree)
  for (j in 1:nrow(xp)){
    for (i in 1:ntree) {
      n0[[i]] <- y[nodesobs[,i]==nodessyn[j,i]]
    }
    empty <- sapply(n0, length)
    ndonors[[j]] <- unlist(n0[empty!=0])
  }
  
  yhat <- sapply(ndonors, sample, size=1)
  
  if (is.factor(y)) yhat <- factor(yhat, levels = obslevels) 
  if (!is.factor(y) & smoothing=="density") yhat <- syn.smooth(yhat,y)
    
  return(list(res = yhat, fit = rf.fit))    # "bag"
}


###-----syn.cartbboot-----------------------------------------------------

syn.cartbboot <- function(y, x, xp, proper = FALSE, 
                           minbucket = 5, cp = 1e-04, ...) 
{
  if (proper==TRUE){
    s <- sample(length(y),replace=TRUE)
    y <- y[s]
    x <- x[s,,drop=FALSE]
  }
  
  minbucket <- max(1, minbucket) 
  if (!is.factor(y)) {
    fit <- rpart(y ~ ., data = as.data.frame(cbind(y, x)), method = "anova",
                 minbucket = minbucket, cp = cp, ...)
    leafnr   <- floor(as.numeric(row.names(fit$frame[fit$where,])))
    fit$frame$yval <- as.numeric(row.names(fit$frame))
    nodes    <- predict(object = fit, newdata = xp)
    donor    <- lapply(nodes, function(s) y[leafnr == s])
    bboot.p  <- lapply(unique(nodes),function(s) diff(c(0,sort(runif(length(y[leafnr == s])-1)),1))); names(bboot.p) <- unique(nodes)
    donor.p  <- lapply(nodes,function(s) bboot.p[[which(unique(nodes)==s)]])
    new      <- sapply(1:length(donor), function(s) resample(donor[[s]],1,p=donor.p[[s]]))
  } else {
    fit <- rpart(y ~ ., data = as.data.frame(cbind(y, x)), method = "class",
                 minbucket = minbucket, cp = cp, ...)
    fit$frame$yval <- as.numeric(row.names(fit$frame))
    pred.cc    <- predict(object=fit,newdata=xp,type="matrix")[,2:(nlevels(y)+1)]   # cc - class counts
    pred.nodes <- predict(object=fit,newdata=xp,type="vector")
    donor      <- sapply(rownames(pred.cc), function(s) rep(levels(y),pred.cc[s,]))
    bboot.p    <- lapply(unique(pred.nodes),function(s) diff(c(0,sort(runif(fit$frame$n[which(rownames(fit$frame)==s)]-1)),1))); names(bboot.p) <- unique(pred.nodes)
    donor.p    <- sapply(1:nrow(pred.cc),function(s) bboot.p[[as.character(pred.nodes[s])]]); names(donor.p) <- names(donor)
    new        <- sapply(1:length(donor),function(s) resample(donor[[s]],1,p=donor.p[[s]]))
    new        <- factor(new,levels=levels(y))
  }

  return(list(res = new, fit = fit))
}


###-----syn.nested---------------------------------------------------------
# function for allocating to subcategories (random sampling within groups)

syn.nested <- function (y, x, xp, ...)
{
  xr   <- apply(x, 1, function(x) paste(x, collapse = "-"))
  xpr  <- apply(xp, 1, function(x) paste(x, collapse = "-"))
  uxpr <- sort(unique(xpr))
  
  index  <- 1:length(y)
  indexp <- rep(0, nrow(xp))
  for (i in uxpr) {
    indexp[xpr == i] <- sample(index[xr == i], sum(xpr == i), TRUE)
  }
  yp <- y[indexp]
 
  return(list(res = yp, fit = "nested"))
}


###-----syn.satcat---------------------------------------------------------

syn.satcat <- function (y, x, xp, proper = FALSE, ...)
{
  # fits a saturated model to combinations of variables
  # method fails if the predictor variables generate 
  # a combination of variables not found in the original data
  
  xr  <- apply(x, 1, function(x) paste(x, collapse = "-"))
  syn.categories <- apply(xp, 1, function(x) paste(x, collapse = "-"))
  
  if (!all(names(table(syn.categories)) %in% names(table(xr)))) { 
    cat("\n\n")
    print(table(syn.categories)[!names(table(syn.categories)) %in% names(table(xr))])
    stop('The combined groups above for "satcat" have no records in original data.\n', call. = FALSE)
  }
  
  uxpr <- sort(unique(syn.categories))
  
  index  <- 1:length(y)
  indexp <- rep(0, nrow(xp))
  for (i in uxpr) {
    indexp[syn.categories == i] <- sample(index[xr == i], sum(syn.categories == i), TRUE)
    if (proper==TRUE) indexp[syn.categories == i] <- sample(indexp[syn.categories == i], sum(syn.categories == i), TRUE)
  }
  yp <- y[indexp]
  fit <- table(xr)
  return(list(res = yp, fit = fit))
}


###-----syn.constant-------------------------------------------------------

syn.constant <- function(y, xp, ...) 
{
  yp <- y
  length(yp) <- xp
  if (xp > length(y)){
    yp[(length(y) + 1):xp] <- names(which.max(table(y, exclude = NULL)))  # in case y is 'almost' constant 
    if (is.numeric(y)) yp <- as.numeric(yp)
    else if (is.logical(y)) yp <- as.logical(yp)
  }
  return(list(res = yp, fit = "constant"))
}


###-----syn.collinear------------------------------------------------------

syn.collinear <- function (y, x, xp, ...)
{
  x <- x[,1]                               #!BN to check
  xp <- xp[,1]                             #!BN to check
  indexp  <- match(xp,x)
  yp <- y[indexp]
  return(list(res = yp, fit = "collinear"))
}



# O T H E R   A U X I L I A R Y   F U N C T I O N S  
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| 

###-----is.passive---------------------------------------------------------
is.passive <- function(string) return("~"==substring(string,1,1))

  
###-----resample-----------------------------------------------------------
# used in syn.cart() and syn.cartbboot() instead of sample() 
# for safety reasons, i.e. to avoid sampling from 1:x when length(x)==1
resample   <- function(x, ...) x[sample.int(length(x),...)]


###-----decimalplaces------------------------------------------------------
# counts number of decimal places - used for rounding smoothed values
# approximate in some cases (as.character -> 15 significant digits; 
# scientific notation)
decimalplaces <- function (x) 
{
  x <- x - floor(x) # -> more digit numbers 
  if ((x%%1) != 0 & (round(x, 15)%%1 != 0)) {
    nchar(strsplit(sub("0+$","",as.character(x)),".",fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

###-----get.names----------------------------------------------------------
get.names <- function(formula,names){
  res <- match(all.vars(formula)[-1],names)
  return(res)
}


###-----addXfac------------------------------------------------------------
# function to add factors
addXfac <- function(x,...){
  df  <- cbind.data.frame(x,...)
  if (any(sapply(df,is.factor)+sapply(df,is.numeric)!=1)) stop("All arguments must be factors or numeric")
  fn  <- function(f){
    if (is.factor(f)) f <- as.numeric(levels(f))[f]
    else f <- f
  }
  df  <- as.data.frame(sapply(df,fn))
  add <- factor(rowSums(df))
  return(add)
  }
  
  
###-----syn.smooth--------------------------------------------------------- 
syn.smooth <- function(ysyn, yobs){
  ys <- 1:length(ysyn)
  # exclude from smoothing if freq for a single value higher than 70% 
  maxfreq <- which.max(table(ysyn))
  maxcat  <- as.numeric(names(table(ysyn))[maxfreq])
  if (table(ysyn)[maxfreq]/sum(table(ysyn))>.7) ys <- which(ysyn!=maxcat)
  # exclude from smoothing if data are top-coded - approximate check
  if (10*table(ysyn)[length(table(ysyn))-1] <
    tail(table(ysyn),n=1)-table(ysyn)[length(table(ysyn))-1]){
    ys   <- ys[-which(ysyn==max(yobs))]
    maxy <- max(yobs)
  }
  densbw  <- density(ysyn[ys],width="SJ")$bw
  ysyn[ys] <- rnorm(n=length(ysyn[ys]),mean=ysyn[ys],sd=densbw)
  if (!exists("maxy")) maxy <- max(yobs) + densbw
  ysyn[ys] <- pmax(pmin(ysyn[ys],maxy),min(yobs))
  ysyn[ys] <- round(ysyn[ys],max(sapply(yobs,decimalplaces)))      
  return(ysyn)
}

