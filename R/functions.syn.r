# Functions for synthesising data, some of which are adapted
# from mice package by S. van Buuren and K. Groothuis-Oudshoorn,
# TNO Quality of Life


###-----.norm.fix.syn------------------------------------------------------

.norm.fix.syn <- function(y, x, ridge = 0.00001, ...)
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

.norm.draw.syn <- function(y, x, ridge = 0.00001, ...)
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
  if (proper == FALSE){
    parm <- .norm.fix.syn(y, x, ...)
  } else {
    parm <- .norm.draw.syn(y, x, ...)
  }  
  res <- xp %*% parm$beta + rnorm(nrow(xp)) * parm$sigma
  res <- round(res, max(sapply(y, decimalplaces)))
  
  return(list(res = res, fit = parm))
}


###-----syn.lognorm--------------------------------------------------------

syn.lognorm <- function(y, x, xp, proper = FALSE, ...) 
{
  addbit <- FALSE
  if (any(y < 0)) stop("Log transformation not appropriate for negative values.\n", call. = FALSE)
  if (any(y == 0)) {y <- y + .5*min(y[y != 0]); y <- log(y); addbit <- TRUE}  ##  warning about this and above should be in check model
  else y <- log(y)
  x   <- cbind(1, as.matrix(x))
  xp  <- cbind(1, as.matrix(xp))
  if (proper == FALSE) {
    parm <- .norm.fix.syn(y, x, ...)
  } else {
    parm <- .norm.draw.syn(y, x, ...)
  }
  res <- xp %*% parm$beta + rnorm(nrow(xp)) * parm$sigma
  if (addbit) {res <- res - .5 * min(y[y != 0]); res[res <= 0] <- 0}
  res <- exp(res)
  res <- round(res, max(sapply(y, decimalplaces)))
  
  return(list(res = res, fit = parm))
}


###-----syn.sqrtnorm-------------------------------------------------------

syn.sqrtnorm <- function(y, x, xp, proper = FALSE, ...) 
{
  addbit <- FALSE
  if (any(y < 0)) stop("Square root transformation not appropriate for negative values.\n", call. = FALSE)   ##  needs check in checkmodel
  else y <- sqrt(y)
  x   <- cbind(1, as.matrix(x))
  xp  <- cbind(1, as.matrix(xp))
  if (proper == FALSE) {
    parm <- .norm.fix.syn(y, x, ...)
  } else {
    parm <- .norm.draw.syn(y, x, ...)
  }
  res <- xp %*% parm$beta + rnorm(nrow(xp)) * parm$sigma
  res <- res^2
  res <- round(res, max(sapply(y, decimalplaces)))
  
  return(list(res = res, fit = parm))
}


###-----syn.cubertnorm-----------------------------------------------------

syn.cubertnorm <- function(y, x, xp, proper = FALSE, ...) 
{
  addbit <- FALSE
  y <- sign(y)*abs(y)^(1/3)
  x   <- cbind(1, as.matrix(x))
  xp  <- cbind(1, as.matrix(xp))
  if (proper == FALSE) {
    parm <- .norm.fix.syn(y, x, ...)
  } else {
    parm <- .norm.draw.syn(y, x, ...)
  }
  res <- xp %*% parm$beta + rnorm(nrow(xp)) * parm$sigma
  res <- res^3
  res <- round(res, max(sapply(y, decimalplaces)))

  return(list(res = res, fit = parm))
}


###-----syn.normrank-------------------------------------------------------

syn.normrank <- function(y, x, xp, smoothing = "", proper = FALSE, ...)
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

  z  <- qnorm(rank(y)/(length(y) + 1))
  x  <- cbind(1, as.matrix(x))
  xp <- cbind(1, as.matrix(xp))

  if (proper == FALSE) {
    parm <- .norm.fix.syn(z, x, ...)
  } else {
    parm <- .norm.draw.syn(z, x, ...)
  }
  
  pred <- (xp %*% parm$beta + rnorm(nrow(xp)) * parm$sigma)
  res  <- round(pnorm(pred)*(length(y) + 1))
  res[res < 1] <- 1
  res[res > length(y)] <- length(y)
  res  <- sort(y)[res]

  if (smoothing != "") {
    res <- syn.smooth(res, y, smoothing = smoothing)
  }

  # if (smoothing == "") res  <- sort(y)[res]
  # 
  # if (smoothing == "density") {
  #   ydsamp <- y
  #   ys     <- 1:length(y)
  #   maxfreq <- which.max(table(y))
  #   maxcat  <- as.numeric(names(table(y))[maxfreq])
  #   if (table(y)[maxfreq]/sum(table(y)) > .7) ys <- which(y != maxcat)
  #   if (10 * table(y)[length(table(y)) - 1] < 
  #     tail(table(y), n = 1) - table(y)[length(table(y)) - 1]) {
  #     ys <- ys[-which(y == max(y))]  
  #     maxy <- max(y)
  #   }   
  #   densbw <- density(y[ys], width = "SJ")$bw
  #   ydsamp[ys] <- rnorm(length(ydsamp[ys]), 
  #     mean = sample(ydsamp[ys], length(ydsamp[ys]), replace = TRUE), sd = densbw)
  #   if (!exists("maxy")) maxy <- max(y) + densbw
  #   ydsamp[ys] <- pmax(pmin(ydsamp[ys],maxy),min(y))
  #   res <- sort(ydsamp)[res]
  # }

  return(list(res = res, fit = parm))
}


###-----.pmm.match---------------------------------------------------------

.pmm.match <- function(z, yhat = yhat, y = y, donors = 3, ...)
{
# Auxilary function for syn.pmm.
# z    = target predicted value (scalar)
# yhat = array of fitted values, to be matched against z
# y    = array of donor data values

# Finds the three cases for which abs(yhat-z) is minimal,
# and makes a random draw from these.

  d <- abs(yhat - z)
  m <- sample(y[rank(d, ties.method = "random") <= donors], 1)
  
  return(m)
}


###-----syn.pmm------------------------------------------------------------

syn.pmm <- function(y, x, xp, smoothing = "", proper = FALSE, ...)
{
# Synthesis of y by predictive mean matching
# Warning: can be slow for large data sets 
# for which syn.normrank may be a better choice
  x       <- cbind(1, as.matrix(x))
  xp      <- cbind(1, as.matrix(xp))
  if (proper == FALSE) {
    parm <- .norm.fix.syn(y, x, ...)
  } else {
    parm <- .norm.draw.syn(y, x, ...)
  }
  yhatobs <- x  %*% parm$coef
  yhatmis <- xp %*% parm$beta
  res <- apply(as.array(yhatmis), 1, .pmm.match, yhat = yhatobs, y = y, ...)
  
  if (smoothing != "") {
    res <- syn.smooth(res, y, smoothing = smoothing)
  }
  
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
  if (p == 0) return(list(y = y, x = x, w = rep(1, length(y))))
  
  # skip augmentation if there is only 1 missing value  
  if (length(y) == 1) return(list(y = y, x = x, w = rep(1, length(y))))
    
  # calculate values to augment
  mean <- apply(x,2,mean)
  sd   <- sqrt(apply(x,2,var))
  minx <- apply(x,2,min)
  maxx <- apply(x,2,max)
  nr   <- 2 * p * ki
  a    <- matrix(mean, nrow = nr, ncol = p, byrow = TRUE)
  b    <- matrix(rep(c(rep(c(0.5, -0.5), ki), rep(0,nr)), length = nr*p), 
                 nrow = nr, ncol = p, byrow = FALSE)
  c    <- matrix(sd, nrow = nr, ncol = p, byrow = TRUE)
  d    <- a + b * c
  d    <- pmax(matrix(minx, nrow = nr, ncol = p, byrow = TRUE), d)
  d    <- pmin(matrix(maxx, nrow = nr, ncol = p, byrow = TRUE), d)
  e    <- rep(rep(icod, each = 2), p)
  dimnames(d) <- list(paste("AUG", 1:nrow(d), sep = ""), dimnames(x)[[2]])
  xa   <- rbind.data.frame(x, d)
  # beware, concatenation of factors
  # this change needed to avoid reordering of factors                           
  # if (is.factor(y)) ya <- as.factor(levels(y)[c(y,e)]) else ya  <- c(y, e)
  if (is.factor(y)) ya <- addNA(factor(levels(y)[c(y, e)], 
    levels = levels(y)), ifany = TRUE) else ya <- c(y, e)   
  wa <- c(rep(1, length(y)),rep((p + 1)/nr, nr))

  return(list(y = ya, x = xa, w = wa))
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

  xmeans <- lapply(x, mean)                      ## x matrix centred
  x  <- mapply(function(x, y) x - y, x, xmeans)
  xp <- mapply(function(x, y) x - y, xp, xmeans) ## also xp to match
  
  if (is.null(denom)) {
    aug <- augment.syn(y, x, ...)
    # when no missing data must set xf to augmented version
    xf   <- aug$x
    y    <- aug$y
    w    <- aug$w
    xf   <- cbind(1, as.matrix(xf))
    xp   <- cbind(1, as.matrix(xp))
    expr <- expression(glm.fit(xf, y, family = binomial(link = logit), weights = w))
    fit  <- suppressWarnings(eval(expr))
    fit.sum <- summary.glm(fit)
    beta <- coef(fit)
    if (proper == TRUE) {
      rv   <- t(chol(fit.sum$cov.unscaled))
      beta <- beta + rv %*% rnorm(ncol(rv))  
    }
    p   <- 1/(1 + exp(-(xp %*% beta)))  
    vec <- (runif(nrow(p)) <= p)
    if (!is.logical(y)) vec <- as.numeric(vec)          
    if (is.factor(y)) vec <- factor(vec,c(0,1), labels = levels(y))
  } else {
    aug <- augment.syn(y, x, ...)
    # when no missing data must set xf to augmented version
    xf   <- aug$x
    y    <- aug$y
    w    <- aug$w
    xf   <- cbind(1, as.matrix(xf))
    xp   <- cbind(1, as.matrix(xp))
    den  <- w
    denind <- which(den == 1)
    den[denind] <- denom
    yy   <- y/den        #denom give then average response
    yy[den < 1]   <- mean(yy[denind]) 
    expr <- expression(glm.fit(xf, yy, family = binomial(link = logit), weights = den))
    fit  <- suppressWarnings(eval(expr))
    fit.sum <- summary.glm(fit)
    beta <- coef(fit.sum)[, 1]
    if (proper == TRUE) {
      rv   <- t(chol(fit.sum$cov.unscaled))
      beta <- beta + rv %*% rnorm(ncol(rv))  
    }
    p <- 1/(1 + exp(-(xp %*% beta)))  
    vec <- rbinom(nrow(p),denomp, p) 
  }
  return(list(res = vec, fit = fit.sum))
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
  
  if (proper == TRUE) { # bootstrap to make proper
    s   <- sample(length(y), replace = TRUE)
    x   <- x[s, , drop = FALSE]
    y   <- y[s]  
    y   <- factor(y)
  }
  aug <- augment.syn(y, x, ...)
  # yf and xf needed for augmented data to save x as non augmented  not now needed can tidy
  xf  <- aug$x
  yf  <- aug$y
  w   <- aug$w
    
  ### rescaling numeric to [0,1]
  toscale <- sapply(xf, function(z) (is.numeric(z) & (any(z < 0) | any(z > 1))))
  rsc <- sapply(xf[, toscale, drop = FALSE], range)
  xf_sc <- xf
  for (i in names(toscale[toscale == TRUE])) xf_sc[, i] <- (xf_sc[, i] - rsc[1,i])/(rsc[2,i] - rsc[1,i])
  for (i in names(toscale[toscale == TRUE])) xp[, i] <- (xp[, i] - rsc[1,i])/(rsc[2,i] - rsc[1,i])
  ###
  
  xfy <- cbind.data.frame(yf, xf_sc)  
  fit <- multinom(formula(xfy), data = xfy, weights = w,
    maxit = maxit, trace = trace, MaxNWts = MaxNWts, ...)
  if (fit$convergence == 1) cat("\nReached max number of iterations for a multinomial model\nsuggest rerunning with polyreg.maxit increased (default 1000)\n")             
  post <- predict(fit, xp, type = "probs") 
  if (length(y) == 1) post <- matrix(post, nrow = 1, ncol = length(post)) 
  if (!is.factor(y)) y <- as.factor(y)
  nc <- length(levels(yf))                    
  un <- rep(runif(nrow(xp)), each = nc)
  if (is.vector(post)) post <- matrix(c(1 - post, post), ncol = 2)
  draws <- un > apply(post, 1, cumsum)
  idx   <- 1 + apply(draws, 2, sum)
  res <- levels(yf)[idx]
  if (length(table(res)) == 1) {
    cat("\n***************************************************************************************")
    cat("\nWarning the polyreg fit produces only one category for the variable being synthesised." )
    cat("\nThis may indicate that the function multinom used in polyreg failed to iterate, possibly")
    cat("\nbecause the variable is sparse. Check results for this variable carefully.")
    cat("\n****************************************************************************************\n")
  } 
  fitted <- summary(fit)
  return(list(res = res, fit = fitted)) 
}


###-----syn.polr-----------------------------------------------------------

syn.polr <- function(y, x, xp, proper = FALSE, maxit = 1000,
                     trace = FALSE, MaxNWts = 10000, ...)
{
  x   <- as.matrix(x)
  xp  <- as.matrix(xp)
  
  if (proper == TRUE) {  # bootstrap to make proper
    s   <- sample(length(y), replace = TRUE)
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
  fit <- try(suppressWarnings(polr(formula(xfy), data = xfy, Hess = TRUE, weights = wf, ...)), silent = TRUE)

  if (inherits(fit, "try-error")) {
    fit <- multinom(formula(xfy), data = xfy, weights = wf,
                    maxit = maxit, trace = trace, Hess = TRUE, MaxNWts = MaxNWts, ...)
    cat("\tMethod changed to multinomial")
    if (fit$convergence == 1) cat("\nReached max number of iterations for a multinomial model\nRerun with polyreg.maxit increased (default 100)\n")
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
  fitted <- summary(fit)
  return(list(res = res, fit = fitted)) 
}


###-----syn.sample---------------------------------------------------

syn.sample <- function(y, xp, smoothing = "", cont.na = NA, proper = FALSE, ...) 
{
  # Generates random sample from the observed y's
  # with bootstrap if proper == TRUE
  if (proper == TRUE) y <- sample(y, replace = TRUE)
  yp <- sample(y, size = xp, replace = TRUE)
  
  if (smoothing != "") yp[!(yp %in% cont.na)] <- 
    syn.smooth(yp[!(yp %in% cont.na)], y[!(y %in% cont.na)], 
               smoothing = smoothing)
  
  return(list(res = yp, fit = "sample"))
}


###-----syn.passive--------------------------------------------------------
syn.passive <- function(data, func)
{
  # Special elementary synthesis method for transformed data.
  # SuppressWarnings to avoid message 'NAs by coercion for NAtemp  
  res <- suppressWarnings(model.frame(as.formula(func), data, 
                                      na.action = na.pass))	

  return(list(res = res, fit = "passive"))
}


###-----syn.cart-----------------------------------------------------------

syn.cart <- function(y, x, xp, smoothing = "", proper = FALSE, 
                     minbucket = 5, cp = 1e-08, ...)
{
  ylogical <- is.logical(y)
  
  if (proper == TRUE) {
    s <- sample(length(y), replace = TRUE)
    x <- x[s,,drop = FALSE]
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
    nodes       <- predict(object = fit, newdata = xp)
    # BN:16/06/20
    # node numbering: node * 2 + 0:1    
    notleaf <- setdiff(nodes, leafnr)
    # if (length(notleaf) > 0) {
    #   for (i in notleaf){
    #     nodes[which(nodes == i)] <- 2 * i + sample(0:1, 1)
    #   }
    # }
    if (length(notleaf) > 0) {
      for (i in notleaf){
        j <- i
        while(!(j %in% leafnr)){
          j <- 2 * j + sample(0:1, 1)
        }
        nodes[which(nodes == i)] <- j
      }
    }
    
    uniquenodes <- unique(nodes)
    new  <- vector("numeric",nrow(xp))
    for (j in uniquenodes) {
      donors <- y[leafnr == j] # values of y in a leaf
      new[nodes == j] <- resample(donors, size = sum(nodes == j), 
                                  replace = TRUE)
    }
 
    if (smoothing != "") new <- syn.smooth(new, y, smoothing = smoothing)
    
    #donor <- lapply(nodes, function(s) y[leafnr == s])
    #new   <- sapply(1:length(donor),function(s) resample(donor[[s]], 1))
  } else {
    y     <- factor(y)
    fit   <- rpart(y ~ ., data = as.data.frame(cbind(y, x)), method = "class",
                   minbucket = minbucket, cp = cp, ...)
    nodes <- predict(object = fit, newdata = xp)
    new   <- apply(nodes, MARGIN = 1, FUN = function(s) resample(colnames(nodes), 
                                                        size = 1, prob = s))
    if (ylogical) {
      new   <- as.logical(new)
    } else {
      new   <- factor(new, levels = levels(y)) 
    }
  }
  
  return(list(res = new, fit = fit))
}


###-----syn.ctree----------------------------------------------------------

syn.ctree <- function(y, x, xp, smoothing = "", proper = FALSE, minbucket = 5, 
                      mincriterion = 0.9, ...)
                      # teststat = "max", testtype = "Univariate", 
                      
{ 
  if (proper == TRUE) {
    s <- sample(length(y), replace = truehist())
    y <- y[s]
    x <- x[s, , drop = FALSE]
  }
  
  for (i in which(sapply(x, class) != sapply(xp,class))) xp[,i] <-
  eval(parse(text = paste0("as.", class(x[,i]), "(xp[,i])", sep = "")))
  # Fit a tree
  datact     <- ctree(y ~ ., data = as.data.frame(cbind(y,x)), 
    controls = ctree_control(minbucket = minbucket, mincriterion = mincriterion, 
                             # teststat = teststat, testtype = testtype, 
                             ...))
  fit.nodes  <- where(datact)
  nodes      <- unique(fit.nodes)
  no.nodes   <- length(nodes)
  pred.nodes <- where(datact, newdata = xp)
  # Get row numbers for predicted by sampling with replacement from existing data
  rowno      <- 1:length(y)
  newrowno   <- vector("integer", nrow(xp))

  for (i in nodes) {
    newrowno[pred.nodes == i] <- sample(rowno[fit.nodes == i],
                                      length(newrowno[pred.nodes == i]),
                                      replace = TRUE)
  }
  new <- y[newrowno]
  if (!is.factor(y) & smoothing != "") new <- 
    syn.smooth(new, y, smoothing = smoothing )
  
  return(list(res = new, fit = datact))
}


###-----syn.survctree------------------------------------------------------

syn.survctree <- function(y, yevent, x, xp, proper = FALSE, minbucket = 5, ...)
# time, event - data column numbers
{      
  if (proper == TRUE) {
    s <- sample(length(y), replace = TRUE)
    y <- y[s]
    x <- x[s, , drop = FALSE]                        
    yevent <- yevent[s]
  }

  for (i in which(sapply(x, class) != sapply(xp,class))) xp[,i] <-
      eval(parse(text = paste0("as.", class(x[,i]), "(xp[,i])", sep = "")))

  if (is.factor(yevent)) {
    yevent0 <- as.numeric(yevent) - 1
  } else {
    yevent0 <- yevent 
  }

  # Fit a tree  
  datact     <- ctree(Surv(y, yevent0) ~ ., 
                      data = as.data.frame(cbind(y, yevent0, x)),
                      controls = ctree_control(minbucket = minbucket, ...))
  fit.nodes  <- where(datact)
  nodes      <- unique(fit.nodes)
  no.nodes   <- length(nodes)
  pred.nodes <- where(datact, newdata = xp)
  # Get row numbers for predicted by sampling
  # with replacement from existing data
  rowno      <- 1:length(y)
  newrowno   <- rep(0,nrow(xp))
  for (i in nodes) {
    newrowno[pred.nodes == i] <- sample(rowno[fit.nodes == i],
    length(newrowno[pred.nodes == i]), replace = TRUE)
   }
  #Predicte node & sample time+event
  faketime  <- y[newrowno]
  fakeevent <- yevent[newrowno]

  return(list(syn.time = faketime, syn.event = fakeevent, fit = datact))
}


###-----syn.rf-------------------------------------------------------------
# bagging when mtry = ncol(x) - using all predictors
syn.rf <- function(y, x, xp, smoothing = "", proper = FALSE, ntree = 10, ...) 
{ 

  #nodesize <- max(1, nodesize)  # safety
  #if (proper == TRUE) {
  #  s <- sample(length(y), replace = T); y <- y[s]
  #  x <- x[s, , drop = FALSE]
  #}  

  for (i in which(sapply(x,class) != sapply(xp, class))) xp[,i] <-
    do.call(paste0("as.", class(x[,i])[1]), unname(xp[, i]))

  if (is.factor(y)) {
    obslevels <- levels(y)
    y <- droplevels(y)
  }

  # fit a random forest
  # regression (mtry = p/3), classification (mtry = sqrt(p))
  rf.fit <- randomForest(y ~ ., data = cbind.data.frame(y,x), ntree = ntree, ...)
  nodessyn <- attr(predict(rf.fit, newdata = xp, nodes = T), "nodes")
  nodesobs <- attr(predict(rf.fit, newdata = x, nodes = T), "nodes")

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
  if (!is.factor(y) & smoothing != "") yhat <- 
    syn.smooth(yhat, y, smoothing = smoothing)
    
  return(list(res = yhat, fit = rf.fit)) 
}


###-----syn.ranger---------------------------------------------------------
# bagging when mtry = ncol(x) - using all predictors
# contributed by Caspar J. van Lissa

syn.ranger <- function(y, x, xp, smoothing = "", proper = FALSE, ...) 
{ 
  dots <- list(...)
  dots[c("formula", "data")] <- NULL
  if("min.node.size" %in% names(dots)){
    dots[["min.node.size"]] <- max(1, dots[["min.node.size"]])  # safety
  }
  
  if (proper == TRUE) {
    s <- sample(length(y), replace = TRUE)
    y <- y[s]
    x <- x[s, , drop = FALSE]
  }
  
  for (i in which(sapply(x,class) != sapply(xp, class))) xp[,i] <-
    do.call(paste0("as.", class(x[,i])[1]), unname(xp[, i]))
  
  if (is.factor(y)) {
    obslevels <- levels(y)
    y <- droplevels(y)
  }
  
  # fit a random forest
  Args     <- c(list(formula = y ~ ., data = cbind.data.frame(y,x)), dots)
  rf.fit   <- do.call(ranger, Args)
  nodessyn <- predict(rf.fit, data = xp, type = "terminalNodes")$predictions
  nodesobs <- predict(rf.fit, data = x, type = "terminalNodes")$predictions
  ntree    <- rf.fit$num.trees
  ndonors  <- vector("list", nrow(xp))
  n0       <- vector("list", ntree)
  for (j in 1:nrow(xp)) {
    for (i in 1:ntree) {
      n0[[i]] <- y[nodesobs[,i] == nodessyn[j,i]]
    }
    empty <- sapply(n0, length)
    ndonors[[j]] <- unlist(n0[empty != 0])
  }
  
  yhat <- sapply(ndonors, sample, size = 1)          
  
  if (is.factor(y)) yhat <- factor(yhat, levels = obslevels) 
  if (!is.factor(y) & smoothing != "") yhat <- 
    syn.smooth(yhat, y, smoothing = "smoothing")
  
  return(list(res = yhat, fit = rf.fit))
}


###-----syn.bag-------------------------------------------------------------
# bagging when mtry = ncol(x) - using all predictors
syn.bag <- function(y, x, xp, smoothing = "", proper = FALSE, ntree = 10, ...) 
{ 
  #nodesize <- max(1, nodesize)  # safety
  #if (proper == TRUE) {
  #  s <- sample(length(y), replace = T); y <- y[s]
  #  x <- x[s, , drop = FALSE]
  #}  

  for (i in which(sapply(x,class) != sapply(xp, class))) xp[,i] <-
    do.call(paste0("as.", class(x[,i])[1]), unname(xp[, i]))
  
  if (is.factor(y)) {
    obslevels <- levels(y)
    y <- droplevels(y)
  }

  # fit a random forest
  # regression (mtry = p/3), classification (mtry = sqrt(p))
  rf.fit <- randomForest(y ~ ., data = cbind.data.frame(y,x), 
                         ntree = ntree, mtry = ncol(x), ...)
  nodessyn <- attr(predict(rf.fit, newdata = xp, nodes = T), "nodes")
  nodesobs <- attr(predict(rf.fit, newdata = x, nodes = T), "nodes")

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
  if (!is.factor(y) & smoothing != "") yhat <- 
    syn.smooth(yhat,y, smoothing = smoothing)
    
  return(list(res = yhat, fit = rf.fit))
}


###-----syn.nested---------------------------------------------------------
# function for allocating to subcategories (random sampling within groups)

syn.nested <- function(y, x, xp, smoothing = "", cont.na = NA, ...)
{
  xr   <- x[,1]
  xpr  <- xp[,1]
  uxpr <- sort(unique(xpr))
  
  index  <- 1:length(y)
  indexp <- rep(0, nrow(xp))
  for (i in uxpr) {
    indexp[xpr == i] <- sample(index[xr == i], sum(xpr == i), TRUE)
  }
  yp <- y[indexp]
  
  if (smoothing != "") yp[!(yp %in% cont.na)] <-
    syn.smooth(yp[!(yp %in% cont.na)], y[!(y %in% cont.na)], 
               smoothing = smoothing)
  
  return(list(res = yp, fit = "nested"))
}


###-----syn.satcat---------------------------------------------------------

syn.satcat <- function(y, x, xp, proper = FALSE, ...)
{
  # Fits a saturated model to combinations of variables.
  # Method fails if the predictor variables generate
  # a combination of variables not found in the original data.
  if (proper == TRUE) {
    s <- sample(length(y), replace = TRUE)
    y <- y[s]
    x <- x[s, , drop = FALSE]
  }

  xr  <- apply(x, 1, function(x) paste(x, collapse = "-"))
  syn.categories <- apply(xp, 1, function(x) paste(x, collapse = "-"))

  if (!all(names(table(syn.categories)) %in% names(table(xr)))) {
    cat("\n\n")
    print(table(syn.categories)[!names(table(syn.categories)) %in% names(table(xr))])
    stop('The combined groups above for "satcat" have no records in original data.\n
         Consider using grouped synthesis with "syn.catall" to overcome this', call. = FALSE)
  }

  uxpr <- sort(unique(syn.categories))

  index  <- 1:length(y)
  indexp <- rep(0, nrow(xp))
  for (i in uxpr) {
    indexp[syn.categories == i] <- sample(index[xr == i], sum(syn.categories == i), TRUE)
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
  if (xp > length(y)) {
    yp[(length(y) + 1):xp] <- names(which.max(table(y, exclude = NULL)))  # in case y is 'almost' constant 
    if (is.numeric(y)) yp <- as.numeric(yp)
    else if (is.logical(y)) yp <- as.logical(yp)
  }
  return(list(res = yp, fit = "constant"))
}


###-----syn.collinear------------------------------------------------------

syn.collinear <- function(y, x, xp, ...)
{
  x <- x[,1]                               #!BN to check
  xp <- xp[,1]                             #!BN to check
  indexp  <- match(xp, x)
  yp <- y[indexp]
  return(list(res = yp, fit = "collinear"))
}


###-----syn.catall---------------------------------------------------------

syn.catall <- function(x, k, proper = FALSE, priorn = 1, structzero = NULL, 
                       maxtable = 1e8, epsilon = 0, rand = TRUE, ...)
{
 # Fits a saturated model to combinations of variables
 # xp just holds number of synthetic records required

   levs <- sapply(x, function(x) {length(levels(x)) + any(is.na(x))})  # all NAtemp here already
 table.size <- prod(levs)   # exp(sum(log(levs)))
 if (table.size > maxtable) stop("Table has more than ", maxtable/1e6,
   " million cells (", round(table.size/1e6, 2),
   " millions),\nwhich may lead to memory problems.\nYou can rerun syn() with catall.maxtable increased (default: ", maxtable/1e6,
   " millions).\nAlternatively use a smaller group of variables for catall.", 
   sep = "", call. = FALSE)

 N <- dim(x)[1]
 if (proper == TRUE) x <- x[sample(1:N, replace = TRUE), ]
 tab <- table(x)
 n <- length(tab) 
 addon <- priorn/n
 # Set structural zero cells to zeros 
 if (!is.null(structzero)) {
   sz <- checksz(structzero, x) ## checks and converts var names to numbers
   if (sum(tab[sz]) > 0) cat("
\n************************************************************************
WARNING: Total of ", sum(tab[sz])," counts of original data in structural zero cells.
************************************************************************\n", sep = "")
 }
 # Add extra to prior
 tab <- (tab + addon)
 if (!is.null(structzero)) tab[sz] <- 0
 dt  <- dim(tab)
 dn  <- dimnames(tab)
  if (epsilon > 0) {
    if (rand == TRUE) {
      if (!is.null(structzero)) tab[sz] <- addlapn(tab[sz], epsilon) 
      else tab <- addlapn(tab, epsilon) 
      fit <- tab
      tab <- tab/sum(tab)   # get it as proportions
      tab <- rmultinom(1, k, tab)
    } else {
      if (!is.null(structzero)) tab[sz] <- addlapn(tab[sz], epsilon) 
      else tab <- addlapn(tab, epsilon ) #GR2022
      fit <- tab
      tab <- roundspec(tab*k/sum(tab))
    }
 } else {
    if (rand == TRUE) {
      fit <- tab
      tab <- tab/sum(tab)   # get it as proportions
      tab <- rmultinom(1, k, tab)
    }
    else stop("If you set rand = FALSE when epsilon = 0 synthpop will return the data unchanged", call. = FALSE) 
 }
 tab <- array(tab, dt)
 dimnames(tab) <- dn
 res <- array.to.frame(tab)
 for (i in 1:dim(res)[2]) res[, i] <- addNA(res[, i], ifany = TRUE)
 res <- res[sample(1:nrow(res)), ]
 return(list(res = res, fit = fit))
}


###-----syn.ipf------------------------------------------------------------
  
syn.ipf <- function(x, k, proper = FALSE, priorn = 1, structzero = NULL, 
                    gmargins = "twoway", othmargins = NULL, tol = 1e-3, max.its = 5000,
                    maxtable = 1e8, print.its = FALSE, epsilon= 0, rand = TRUE,...)
{
 # Fits log-linear model to combinations of variables
 # k just holds number of synthetic records required

 levs <- sapply(x, function(x) {length(levels(x)) + any(is.na(x))}) # all NAtemp here already
 table.size <- prod(levs)   # exp(sum(log(levs))) 
 if (table.size > maxtable) stop("Table has more than ", maxtable/1e6,
   " million cells (", round(table.size/1e6, 2),
   " millions),\nwhich may lead to memory problems.\nYou can rerun syn() with ipf.maxtable increased (default: ", maxtable/1e6,
   " millions).\nAlternatively use a smaller group of variables for ipf.", 
   sep = "", call. = FALSE)

 N  <- dim(x)[1]
 nv <- dim(x)[2] # number of variables
 if (proper == TRUE) x <- x[sample(1:N, replace = TRUE),]
 tab <- table(x, useNA = "ifany")
 n <- length(tab)
 # Add extra to prior
 addon <- priorn/n
 # Set structural zero cells to zeros 
 if (!is.null(structzero)) {
   sz <- checksz(structzero, x) # checks and converts var names to numbers
   if (sum(tab[sz]) > 0) cat("
\n************************************************************************
WARNING: Total of ", sum(tab[sz])," counts of original data in structural zero cells.
************************************************************************\n", sep = "")
 }
 if (!is.null(gmargins)) {
   if (gmargins == "twoway") {
     n_margins  <- nv*(nv - 1)/2
     mx_margins <- combn(1:nv, 2)
     margins <- split(mx_margins, col(mx_margins))
   } else if (gmargins == "oneway") {
     margins <- as.list(1:nv)
   } else stop("Only 'oneway' or 'twoway' are implemented for gmargins.\n", 
               call. = FALSE)
   if (!is.null(othmargins)) for (i in 1:length(othmargins)) {
     margins[[length(margins) + 1]] <- othmargins[[i]]
   }
 } else {
   if (!is.null(othmargins)) margins <- othmargins
   else stop("Need to specify some margins.\n", call. = FALSE)
 }
 
 umar <- unique(unlist(margins))
 missed <- (1:nv)[!(1:nv %in% umar)]
  
 if (length(missed) > 0) cat("\n
**************************************************
SEVERE WARNING: Margins ", missed, " not fitted.
This means they will be fitted as having
the same proportion in each level.
**************************************************\n", sep = "")
 if (epsilon > 0) {eps <- epsilon / length(margins)
   cat("Overall epsilon for DP is ",epsilon," divided equally between ", 
       length(margins)," margins to give ", eps, " each,\n")
 } 

# Get data for margins
 margins.data <- vector("list", length(margins))
 for (i in 1:length(margins)) {
   margins.data[[i]] <- table(x[, margins[[i]]], useNA = "ifany")
    margins.data[[i]] <- margins.data[[i]] + priorn/length(margins.data[[i]])
    if (epsilon > 0) {
      margins.data[[i]] <- addlapn(margins.data[[i]], eps)
    }
 }
 start <- array(1, dim(tab))
  
 results.1 <- suppressWarnings(Ipfp(start, margins, margins.data, 
                iter = max.its, print = print.its, tol = tol, ...))  # note larger than default to speed things up needs checking
 
 if (results.1$conv == TRUE) cat("\n['ipf' converged in ", 
   length(results.1$evol.stp.crit), " iterations]\n", sep = "")
 else cat("\n['ipf' failed to converge in ", length(results.1$evol.stp.crit),
   " iterations]\n\nYou can try to change parameters of Ipfp() function,\ne.g. syn(..., ipf.iter = 2000).", sep = "")
 
 exp1 <- array(N * results.1$p.hat , dim(tab))
 if (!is.null(structzero)) exp1[sz] <- 0
 
 exp1 <- exp1 / sum(exp1)   # get it as proportions
 if (epsilon == 0 & rand == FALSE) cat("WARNING: No DP noise or random noise added,\nData returned will be close to NULL expectation for model defined by margins.\n")
 if (rand == TRUE)  z    <- rmultinom(1, k, exp1)
 else z <- roundspec(exp1*k/sum(exp1))
 res  <- array(z, dim(exp1))
 dimnames(res) <- dimnames(tab)
 res  <- array.to.frame(res)
 for (i in 1:dim(res)[2]) res[,i] <- addNA(res[,i], ifany = TRUE)
 res <- res[sample(1:nrow(res)),]
 fit <- list(margins = margins, margins.data = margins.data)
 return(list(res = res, fit = fit))
}


# O T H E R   A U X I L I A R Y   F U N C T I O N S  
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| 

###-----is.passive---------------------------------------------------------
is.passive <- function(string) return("~" == substring(string,1,1))


###-----resample-----------------------------------------------------------
# used in syn.cart() and syn.cartbboot() instead of sample() 
# for safety reasons, i.e. to avoid sampling from 1:x when length(x)==1
resample <- function(x, ...) x[sample.int(length(x), ...)]


###-----decimalplaces------------------------------------------------------
# counts number of decimal places - used for rounding smoothed values
# approximate in some cases (as.character -> 15 significant digits; 
# scientific notation)
decimalplaces <- function(x) 
{
  x <- x - floor(x) # -> more digit numbers 
  if (!is.na(x) & (x %% 1) != 0 & (round(x, 15) %% 1 != 0)) {
    nchar(strsplit(sub("0+$", "", as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}


###-----get.names----------------------------------------------------------
get.names <- function(formula, names)
{
  res <- match(all.vars(formula)[-1], names)
  return(res)
}


###-----addXfac------------------------------------------------------------
# function to add factors
addXfac <- function(x,...)
{
  df  <- cbind.data.frame(x,...)
  if (any(sapply(df,is.factor) + sapply(df,is.numeric) != 1)) 
    stop("All arguments must be factors or numeric", call. = FALSE)
  fn  <- function(f){
    if (is.factor(f)) f <- as.numeric(levels(f))[f]
    else f <- f
  }
  df  <- as.data.frame(sapply(df, fn))
  add <- factor(rowSums(df))
  return(add)
}


###-----syn.smooth--------------------------------------------------------- 

syn.smooth <- function(ysyn, yobs = NULL, smoothing = "spline", window = 5, ...)
{
  if (!(smoothing %in% c("", "spline", "density", "rmean"))) 
    cat('Smoothing must be one of "spline", "density" or "rmean". No smoothing done.\n')
  if (any(is.na(ysyn))) stop("ysyn cannot contain missing values", call. = FALSE)
  
  else if (smoothing == "density") {
    ys <- 1:length(ysyn)
    # exclude from smoothing if freq for a single value higher than 70% 
    maxfreq <- which.max(table(ysyn))
    maxcat  <- as.numeric(names(table(ysyn))[maxfreq])
    if (table(ysyn)[maxfreq]/sum(table(ysyn)) > .7) ys <- which(ysyn != maxcat)
    # exclude from smoothing if data are top-coded - approximate check
    if (10 * table(ysyn)[length(table(ysyn)) - 1] <
        tail(table(ysyn), n = 1) - table(ysyn)[length(table(ysyn)) - 1]) {
          ys   <- ys[-which(ysyn == max(yobs))]
          maxy <- max(yobs)
    }
 
    densbw  <- density(ysyn[ys], width = "SJ")$bw
    ysyn[ys] <- rnorm(n = length(ysyn[ys]), mean = ysyn[ys], sd = densbw)
    if (!exists("maxy")) maxy <- max(yobs) + densbw
    ysyn[ys] <- pmax(pmin(ysyn[ys], maxy), min(yobs))
    ysyn[ys] <- round(ysyn[ys], max(sapply(yobs, decimalplaces))) 
  }    
  else if (smoothing == "rmean") {
    ord <- order(ysyn)
    ysyn <- runningmean(1:length(ysyn), ysyn[ord], window = window)
    ysyn[ord] <- round(ysyn, max(sapply(yobs, decimalplaces))) 
  }
  else if (smoothing == "spline") {
    ord <- order(ysyn)
    ysyn <- smooth.spline(sort(ysyn), all.knots = FALSE)$y
    ysyn[ord] <- round(ysyn, max(sapply(yobs, decimalplaces))) 
  }
  
  return(ysyn)
}


###-----checksz------------------------------------------------------------

checksz <- function(sz, x) 
{
  if (!is.list(sz)) stop("structzero needs to be a list.\n", call. = FALSE)
  if (!is.character(names(sz)) || !all(grepl("_", names(sz)))) stop("\nstructzero list elements must be named using variable names\nseperated by an underscore, e.g. sex_edu", call. = FALSE)
  
  list_sub_names <- names(do.call("c", sz))
  allvars <- unique(sub(".*\\.", "", list_sub_names))  
  
  if (!all(allvars %in% names(x))) stop("\nStructural zero variables must match names of variables in the data.", call. = FALSE)

  dd  <- as.data.frame(table(x, useNA = "ifany"))
  res <- rep(FALSE, nrow(dd))
  
  for (i in 1:length(sz)) {
    
    tempz <- rep(TRUE, nrow(dd))
    if (names(sz)[i] != paste0(names(sz[[i]]), collapse = "_")) stop("\nNames of structzero list elements must correspond to names of their sublists.", call. = FALSE)
    vars  <- names(sz[[i]])    
    vars  <- match(vars, names(x))
    nvars <- length(vars)
    if (!(is.list(sz[[i]]) & length(sz[[i]]) == nvars)) stop("Each element of structzero list must be a list of length\nequal to the number of variables used to define structural zeros.\n", call. = FALSE)
    
    for (j in 1:nvars) {
      if (!is.numeric(sz[[i]][[j]])) {
        if (!all(sz[[i]][[j]] %in% levels(x[,vars[j]]))) 
          stop("Structural zeros (element ", i, ", variable ", j, 
               "): level(s) of variable not in data.\n", call. = FALSE)
        sz[[i]][[j]] <- match(sz[[i]][[j]], levels(x[, vars[j]]))
      }
      if (!all(sz[[i]][[j]] %in% 1:nlevels(x[,vars[j]]))) 
        stop("Structural zeros (element ", i, ", variable ", j, 
             "): numeric level(s) of variable not in data.\n", sep = "", 
             call. = FALSE)
      tempz <- tempz & (as.numeric(dd[, vars[j]]) %in% sz[[i]][[j]])
    }
    res[tempz] <- TRUE
  }
  return(res)
}


###-----array.to.frame-----------------------------------------------------

array.to.frame <- function(x)
{
  df1 <- as.data.frame.table(x)
  # df1 <- df1[df1$Freq != 0,]
  res <- df1[rep(1:nrow(df1), df1$Freq), -ncol(df1)]
  dimnames(res)[[1]] <- 1:sum(x)
  return(res)
}


###-----roundspec----------------------------------------------------------

roundspec <- function(tab) {  ## special rounding to preserve total
  if (abs(round(sum(tab)) - sum(tab)) > 1e-6) 
    stop("This function assumes sum of tab is a whole number\n", .call= FALSE)
  diff <- round(sum(tab) - sum(round(tab)))
  
  if (diff == 0 & all(tab >= 0)) { 
    result <- round(tab)
  } else if (any(round(tab) <= 0)){
    if (diff < 0 ) {
      inds <- (1:length(tab))[round(tab) >0]
      vals <- tab[round(tab) > 0]
      ordinds <- inds[order(vals)]
      newtab <- round(tab)
      newtab[ordinds[1:(-diff)]]<- newtab[ordinds[1:(-diff)]]  - 1
    } else {
      inds <- (1:length(tab))
      newtab <- round(tab)
      newtab[1:(diff)]<- newtab[1:(diff)] + 1
      
    }
    result <- newtab
  } else {
    if (abs(diff) > length(tab)) cat("Unexpected problem with rounding.\n
Please report to maintainer of synthpop with example\n
including this output\n", tab, "\n")
    result <- round(tab)
    inds <- (1:length(result))[order(result)][1:abs(diff)]
    if (diff >0 ) result[inds] <- result[inds] + 1
    if (diff <0 ) result[inds] <- result[inds] - 1
  }
  return(result)
}


###-----addlapn------------------------------------------------------------

addlapn <- function(x, eps){
  # add Laplace noise with re-scaling to a total or sample size
  
  if (eps <= 0) stop("eps must be > 0\n")
  res <- x + rlaplace(length(x), 0, 1/eps)
  res <- makepos(res, sum(x))
  return(res)
}


###---------------makepos--------------------------------------------------

makepos <- function(lap, tot) {  ## makes positive summing to total
  olap <- order(-lap)
  lapo <- lap[olap]
  lap[olap[cumsum(abs(lapo)) >= tot]] <- 0
  abs(lap)
}





