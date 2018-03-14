###-----utility.gen-------------------------------------------------------- 

utility.gen <- function(object, data, method = "logit", maxorder = 1, 
                        tree.method = "rpart", resamp.method = NULL,
                        nperms = 50, cp = 1e-3, minbucket = 5, mincriterion = 0,  
                        vars = NULL, aggregate = FALSE, maxit = 200, ngroups = NULL, 
                        print.every = 10, digits = 2, print.zscores = FALSE, zthresh = 1.6, 
                        print.ind.results = TRUE, print.variable.importance = FALSE,  ...)
{
 
 m   <- object$m
 cna <- object$cont.na
  
 # Check input parameters 
 # ---
 if (is.null(method) || length(method) != 1 || is.na(match(method, c("cart", "logit")))) {
   stop("Invalid 'method' type - must be either 'logit' or 'cart'.", call. = FALSE) 
 }
 if (!is.null(resamp.method) && is.na(match(resamp.method, c("perm", "pairs")))) {
   stop("Invalid 'resamp.method' type - must be NULL, 'perm' or 'pairs'.", call. = FALSE)
 }
 if (aggregate == TRUE & method != "logit") stop("Aggregation only works for 'logit' method.", call. = FALSE) 
 if (is.null(data)) stop("Requires parameter 'data' to give name of the real data.\n",  call. = FALSE)
 if (!class(object) == "synds") stop("Object must have class 'synds'.\n", call. = FALSE)

 # Check selected variables and make observed and synthetic comparable
 if (!(is.null(vars))) {               
   if (!(all(vars %in% names(data))))	stop("Some 'vars' specified not in data.\n", call. = FALSE)
   data <- data[, vars]
   if (m == 1) {
     if (!all(vars %in% names(object$syn))) stop("Some 'vars' specified not in synthetic data.\n", call. = FALSE)
     else object$syn <- object$syn[, vars]
   } else {
     if (!all(vars %in% names(object$syn[[1]]))) stop("Some 'vars' specified not in synthetic data.\n", call. = FALSE)
     else object$syn <- lapply(object$syn, "[", vars) 
   }
 } else {
   if (m == 1) vars <- names(object$syn) else vars <- names(object$syn[[1]])
   if (!all(vars %in% names(data))) stop("Some variables in synthetic data not in data.\n", call. = FALSE)
   else data <- data[, vars]  # make data match synthetic
 } 
  
 # Check whether some variables are unsynthesised  
 incomplete <- FALSE
 nsynthd <- length(vars)
 unsyn.vars <- names(object$method)[object$method == ""]  # identify unsynthesised
 if (any(vars %in% unsyn.vars) & !is.null(unsyn.vars)) {
   notunsyn <- vars[!vars %in% unsyn.vars]  # synthesised vars
   incomplete <- TRUE
 }
  
 # Set default resampling according to completeness
 if (is.null(resamp.method)) {
   if (method == "cart") {
     if (incomplete == FALSE) resamp.method = "perm" 
     else if (m == 1) stop("Incomplete synthetic data with utility evaluated from \"cart\" require resamp.method = \"pairs\" that needs a synthesis with m > 1, m = 10 suggested.\n", call. = FALSE)
     else resamp.method = "pairs"
   }
 }
 else if (!is.null(resamp.method) && resamp.method == "pairs" & m == 1) stop("resamp.method = \"pairs\" needs a synthesis with m > 1, m = 10 suggested.\n", call. = FALSE)
 
 # Drop any single value columns 
 leneq1 <- function(x) length(table(as.numeric(x), useNA = "ifany")) %in% (0:1) 
 dout <- sapply(data,leneq1)
 if (m == 1) sout <- sapply(object$syn,leneq1)
 else  sout <- sapply(object$syn[[1]],leneq1)
 dout <- dout & sout
 if (any(dout == TRUE)) {
   cat("Some columns with single values or all missing values in original and synthetic\nexcluded from utility comparisons (excluded variables: ",
       paste(names(data)[dout], collapse = ", "), ").\n", sep = "")
   data <- data[,!dout]
   if (m == 1) object$syn <- object$syn[, !dout]
   else object$syn <- lapply(object$syn, "[", !dout)
 }

 # Categoricasl vars: make missing data part of factor
 catvars <- (1:dim(data)[2])[sapply(data, is.factor)]
 for (i in catvars) {
   data[,i] <- factor(data[,i])
   if (m == 1) object$syn[,i] <- factor(object$syn[,i])
   else for (j in 1:m) object$syn[[j]][,i] <- factor(object$syn[[j]][,i])  
   if (any(is.na(data[,i]))) {
     data[,i] <- addNA(data[,i])
     if (m == 1) object$syn[,i] <- addNA(object$syn[,i])
     else for (j in 1:m) object$syn[[j]][,i] <- addNA(object$syn[[j]][,i])
   }
 }

 # Numeric variables
 numvars <- (1:dim(data)[2])[sapply(data,is.numeric)]
 names(numvars) <- names(data)[numvars]
 # If ngroups != NULL divide numeric variables into ngroups
 data0 <- data  # to save if m > 1
 if (!is.null(ngroups)) {
   for (i in numvars) {
     if (m == 1) {
       groups <- group_num(data[,i], object$syn[,i], ngroups, cont.na = cna[[i]])
       data[,i] <- groups[[1]]
       object$syn[,i] <- groups[[2]]
     } else {
       for (j in 1:m) {
         groups <- group_num(data0[,i], object$syn[[j]][,i], ngroups, cont.na = cna[[i]])
         data[,i] <- groups[[1]]
         object$syn[[j]][,i] <- groups[[2]]
       }
     }
   }
 }   
 for (i in numvars) {
   if (anyNA(data[,i])) {
     newname <- paste(names(data)[i], "NA", sep = "_")
     data <- data.frame(data, 1*(is.na(data[,i])))
     names(data)[length(data)] <- newname
     data[is.na(data[,i]), i] <- 0
     if (m == 1) {
       object$syn <- data.frame(object$syn, 1*(is.na(object$syn[,i])))
       names(object$syn)[length(object$syn)] <- newname
       object$syn[is.na(object$syn[,i]), i] <- 0
     } else {
       for (j in 1:m) {
         object$syn[[j]] <- data.frame(object$syn[[j]], 1*(is.na(object$syn[[j]][,i])))
         names(object$syn[[j]])[length(object$syn[[j]])] <- newname
         object$syn[[j]][is.na(object$syn[[j]][,i]),i] <- 0
       }
     }
   }
 }

 nnosplits <- NA  # for 'logit'
 
 # Function for getting propensity scores 
 # --------------------------------------
 propcalcs <- function(syndata, data) {
    
   df.prop <- rbind(syndata, data)  # make data frame for calculating propensity score
   n1 <- dim(data)[1]
   n2 <- dim(syndata)[1]
   df.prop <- data.frame(df.prop, t = c(rep(1,n2), rep(0,n1)))
    
   if (aggregate == TRUE) {
     aggdat <- aggregate(df.prop[,1], by = df.prop, FUN = length)
     wt <- aggdat$x
     aggdat <- aggdat[, -dim(aggdat)[2]]
   }
    
   if (method == "logit") {
     if (maxorder > 4) cat("Maximum order of interactions over 4 may cause computational problems.\n")
     if (maxorder >= 1) logit.int <- as.formula(paste("t ~ .^", maxorder + 1))
     else logit.int <- as.formula(paste("t ~ ."))
     
     if (aggregate == TRUE) fit <- glm(logit.int, data = aggdat, family = "binomial", 
                                       control = list(maxit = maxit), weights = wt)
     else fit <- glm(logit.int, data = df.prop, family = "binomial",
                     control = list(maxit = maxit))
     if (fit$converged == FALSE) cat("Warning: Logistic model did not converge in ",
                                    maxit, " iterations.\nYou should increase parameter 'maxit'.\n", sep = "")

     # Get number of parameters that involve synthesised variables
     score <- predict(fit, type = "response")
     if (incomplete == FALSE) km1 <- length(fit$coefficients) - 1
     else {
       namescoef <- names(fit$coefficients)
       coefOK <- rep(FALSE, length(namescoef))
       for (nn in notunsyn) coefOK[grepl(nn, namescoef)] <- TRUE
       km1 <- sum(coefOK)
       if (m == 1 || (m > 1 & j == 1)) cat("Expectation of utility uses only coefficients involving synthesised variables: ",
                                           km1, " from ", length(fit$coefficients) - 1, "\n", sep = "")
     }

     if (aggregate == TRUE) { 
       utilVal <- (sum(wt*(score - n2/(n1 + n2))^2, na.rm = T))*(n1 + n2)^3/n1^2/n2
     } else {
       utilVal <- (sum(   (score - n2/(n1 + n2))^2, na.rm = T))*(n1 + n2)^3/n1^2/n2
     } 
      
     utilR     <- utilVal / km1
     utilExp   <- km1
     utilStd   <- (utilVal - km1) / sqrt(km1*2)
     fit$data  <- NULL  # to save space
     fit$model <- NULL  # to save space
     # res.ind <- list(utilVal = utilVal, utilExp = km1, utilR = utilR, utilStd = utilStd, fit = fit)
    
   } else if (method == "cart") {
      
     if (tree.method == "rpart") {
       fit <- rpart(t ~ ., data = df.prop, method = 'class', 
         control = rpart.control(cp = cp, minbucket = minbucket))
       score <- predict(fit)[, 2]
     } else if (tree.method == "ctree") {
       fit <- ctree(t ~ ., data = df.prop, 
         controls = ctree_control(mincriterion = mincriterion, minbucket = minbucket))
       score <- predict(fit)
     }
     utilVal <- (sum((score - n2/(n1 + n2))^2, na.rm = T))*(n1 + n2)^3/n1^2/n2
   }
   
   # Permutation test 
   if (!is.null(resamp.method) && resamp.method == "perm") { # to allow resamp for logit models
     simutil <- rep(0, nperms)
     if (m == 1) j <- 1
     if (j == 1 ) {
       if (print.every == 0) cat("Running ", nperms, " permutations to get NULL pMSE.\n", sep = "")
       else  cat("Running ", nperms," permutations to get NULL pMSE and printing every ", print.every, "th.", sep = "")
     }
     if (print.every > 0) cat("\nsynthesis ", j, ": ", sep = "")
     
     for (i in 1:nperms) {
       if (print.every > 0 & floor(i/print.every) == i/print.every)  cat(i, " ", sep = "")
       pdata <- df.prop
       pdata$t <- sample(pdata$t)
       
       if (method == "cart") {
         if (tree.method == "rpart") {
           sfit <- rpart(t ~ ., data = pdata, method = 'class', control = rpart.control(cp = cp, minbucket = minbucket))
           score <- predict(sfit)[,2]
         } else if (tree.method == "ctree") {
           sfit <- ctree(t ~ ., data = pdata, 
                         controls = ctree_control(mincriterion = mincriterion, minbucket = minbucket))
           score <- predict(sfit)
         }
         simutil[i] <- (sum((score - n2/(n1 + n2))^2, na.rm = T))*(n1 + n2)^3/n1^2/n2
       
       } else if (method == "logit") {
         if (maxorder >= 1) logit.int <- as.formula(paste("t ~ .^", maxorder + 1))
         else logit.int <- as.formula(paste("t ~ ."))
         
         if (aggregate == TRUE) {
           aggdat1 <- aggregate(pdata[,1], by = pdata, FUN = length)
           wt <- aggdat1$x
           aggdat1 <- aggdat1[, -dim(aggdat1)[2]]
           sfit <- glm(logit.int, data = aggdat1, family = "binomial", 
                       control = list(maxit = maxit), weights = wt)
         } else sfit <- glm(logit.int, data = pdata, family = "binomial", control = list(maxit = maxit))
         
         if (sfit$converged == FALSE) cat("Warning: Logistic model did not converge in ",
                                         maxit, " iterations.\nYou should increase parameter 'maxit'.\n", sep = "")
         score <- predict(sfit, type = "response")
         if (aggregate == TRUE) { 
           simutil[i] <- (sum(wt*(score - n2/(n1 + n2))^2, na.rm = T))*(n1 + n2)^3/n1^2/n2/2  # reduced by factor of 2 for logit model
         } else {
           simutil[i] <- (sum(   (score - n2/(n1 + n2))^2, na.rm = T))*(n1 + n2)^3/n1^2/n2/2  # reduced by factor of 2 for logit model
         } 
       }
     }
     nnosplits <- c(sum(simutil < 1e-8), length(simutil))
     utilExp   <- mean(simutil)
     utilR     <- utilVal / utilExp
     utilStd   <- (utilVal - utilExp) / sqrt(var(simutil))
   }
    
   if (!is.null(resamp.method) && resamp.method == "pairs") res.ind <- list(utilVal = utilVal, utilExp = NA, utilR = NA, 
                                                                            utilStd = NA, fit = fit, nnosplits = NA)
   else res.ind <- list(utilVal = utilVal, utilExp = utilExp, utilR = utilR,
                        utilStd = utilStd, fit = fit, nnosplits = nnosplits)
   return(res.ind)
 }
 # -------------------------------------- 
 
  n1 <- nrow(data)

  if (m == 1) {
    n2 <- nrow(object$syn)
    res.ind <- propcalcs(object$syn, data)
    pMSE <- res.ind$utilVal * n1^2 * n2/(n1 + n2)^4
    res <- list(call = match.call(), m = m, method = method, tree.method = tree.method, 
                resamp.method = resamp.method, maxorder = maxorder, vars = vars, 
                aggregate = aggregate, maxit = maxit, ngroups = ngroups,
                mincriterion = mincriterion, nperms = nperms, 
                pMSE = pMSE, utilVal = res.ind$utilVal, utilExp = res.ind$utilExp, utilR = res.ind$utilR, 
                utilStd = res.ind$utilStd, fit = res.ind$fit, nnosplits = res.ind$nnosplits, 
                digits = digits, print.ind.results = print.ind.results, print.zscores = print.zscores,
                zthresh = zthresh, print.variable.importance = print.variable.importance)
  } else {
    n2 <- nrow(object$syn[[1]])  
    pMSE <- utilVal <- utilExp <- utilR <- utilStd <- rep(NA, m)
    fit <- nnosplits <- as.list(1:m)
    if (!is.null(resamp.method) && resamp.method == "pairs") {
      kk <- 0
      simvals <- rep(NA, m*(m - 1)/2)
    }
    for (j in 1:m) {
      res.ind <- propcalcs(object$syn[[j]], data)
      if (method == "logit" & is.null(resamp.method)) {
        if (j == 1) cat("Fitting syntheses: ")
        cat(j, " ", sep = "")
      }
      utilVal[j] <- res.ind$utilVal; fit[[j]] <- res.ind$fit; nnosplits[[j]] <- res.ind$nnosplits
      if (!is.null(resamp.method) && resamp.method == "pairs") {
        if (j == 1) {
          if (print.every == 0) cat("Simulating NULL pMSE from ", m*(m - 1)/2, " pairs.", sep = "")
          else cat("Simulating NULL pMSE from ", m*(m - 1)/2, " pairs, printing every ", print.every, "th:\n", sep = "")
          if (m*(m - 1)/2 < 6 ) cat("Number of pairs too low, suggest increasing number of syntheses (m).\n")
        }
        if (j < m) {
          for (jj in (j + 1):(m)) {
            kk <- kk + 1
            if (print.every > 0 & floor(kk/print.every) == kk/print.every) cat(kk," ",sep = "")
            simvals[kk] <- propcalcs(object$syn[[j]], object$syn[[jj]])$utilVal
          }
        }
      } else { 
        utilExp[j] <- res.ind$utilExp
        utilR[j]   <- res.ind$utilR 
        utilStd[j] <- res.ind$utilStd
      }
    }
    
    if (!is.null(resamp.method) && resamp.method == "pairs") {
      nnosplits <- c(sum(simvals < 1e-8), length(simvals))      
      if (method == "logit") simvals <- simvals / 2
      exp <- mean(simvals)
      expsd <- sqrt(var(simvals))
      for (j in 1:m) {
        utilExp[j] <- exp
        utilR[j]   <- utilVal[j] / exp
        utilStd[j] <- (utilVal[j] - exp) / expsd
      }
    }
    for (j in 1:m) pMSE[j] <- utilVal[j] * n1^2 * n2/(n1 + n2)^4
    
    res <- list(call = match.call(), m = m, method = method, tree.method = tree.method,
                resamp.method = resamp.method, maxorder = maxorder, vars = vars,
                aggregate = aggregate, maxit = maxit, ngroups = ngroups,
                mincriterion = mincriterion, nperms = nperms,
                pMSE = pMSE, utilVal = utilVal, utilExp = utilExp, utilR = utilR,
                utilStd = utilStd, fit = fit, nnosplits = nnosplits, 
                digits = digits, print.ind.results = print.ind.results, print.zscores = print.zscores, 
                zthresh = zthresh, print.variable.importance = print.variable.importance)
  }

  class(res) <- "utility.gen"
  return(res)
}


###-----utility.tab--------------------------------------------------------

utility.tab <- function(object, data, vars = NULL, ngroups = 5, 
                        print.tables = TRUE, print.zdiff = FALSE, 
                        digits = 2, ...) 
{
  # CHECKS 
  #---------
  if (is.null(data)) 
    stop("Requires parameter 'data' to give name of the real data.\n", 
         call. = FALSE)
  if (!class(data) == "data.frame") 
    stop("Data must have class 'data.frame'.\n", call. = FALSE)
  if (!class(object) == "synds") 
    stop("Object must have class 'synds'.\n", call. = FALSE)
  if (is.null(vars)) stop("Need to set variables with vars parameter.\n", call. = FALSE) else if 
    (!(all(vars %in% names(data)))) stop("Unrecognized variable(s) in vars parameter: ", 
     paste(vars[!(vars %in% names(data))], collapse = ", "), call. = FALSE)
  #---------

  my_cont.na <- object$cont.na[match(vars,names(data))]
  data  <- data[, vars, drop = FALSE]
  nvars <- ncol(data)  
  data.orig <- data

  m <- object$m
  if (m == 1) syndata <- list(object$syn) else syndata <- object$syn   
  syndata <- lapply(syndata,'[',vars)
  
  #if (nrow(data) != nrow(syndata[[1]])) 
  #  stop("\nThis function is not for case when sizes of original and synthetic data differ.\n", call.=FALSE)

  df <- ratioFT <- ratioVW <- pvalVW <- pvalFT <- 
    UtabFT <- UtabVW <- nempty <- vector("numeric", m)
  tab.syn <- tab.zdiff <- tabd <- vector("list", m) 
  
  # nobs.missings <- sum(apply(is.na(data),1,sum)>0) 
  # nsyn.missings <- vector("numeric", m)

  for (i in 1:m) {
    data <- data.orig
    # nsyn.missings[i] <- sum(apply(is.na(syndata[[i]]),1,sum)>0) 
    # make all variables into factors
    for (j in 1:nvars) {
      if (is.numeric(data[,j])) {
        grpd <- group_num(data[,j], syndata[[i]][,j], 
          n = ngroups, cont.na = my_cont.na[[j]], ...)
        data[,j] <- grpd[[1]]; syndata[[i]][,j] <- grpd[[2]]
      } else if (is.character(data[,j])) {
        data[,j] <- factor(data[,j])
        syndata[[i]][,j] <- factor(syndata[[i]][,j], 
                                   levels = levels(data[,j]))
      }
      if (any(is.na(data[,j]))) {
        # makes missings into part of factors if present 
        data[,j] <- addNA(data[,j]) 
        syndata[[i]][,j] <- addNA(syndata[[i]][,j]) 
      }
    }
    
    # Adjustment for synthetic data that have different size than original data
    tabd[[i]]       <- table(data) * (nrow(syndata[[i]])/nrow(data))
    
    totcells        <- length(tabd[[i]])
    tab.syn[[i]]    <- table(syndata[[i]])
    nempty[i]       <- sum(tabd[[i]] + tab.syn[[i]] == 0)  
    df[i]           <- totcells - nempty[i] - 1
    diff            <- (tab.syn[[i]] - tabd[[i]])
    expect          <- (tab.syn[[i]] + tabd[[i]])/2
    UtabFT[i]       <- 4*sum((tab.syn[[i]]^(0.5) - tabd[[i]]^(0.5))^2)
    tabsq           <- diff^2/expect
    tabsq[expect == 0] <- 0
    UtabVW[i]       <- sum(tabsq)
    ratioFT[i]      <- UtabFT[i] / df[i]
    pvalFT[i]       <- 1 - pchisq(UtabFT[i], df[i])
    # stdFT[i] <-(UtabFT[i] - df[i]) / sqrt(2*df[i])
    ratioVW[i]      <- UtabVW[i] / df[i]
    pvalVW[i]       <- 1 - pchisq(UtabVW[i], df[i])
    # stdVW[i] <-(UtabVW[i] - df[i]) / sqrt(2*df[i])
    tab.zdiff[[i]]  <- diff/sqrt(expect)

  }

  if (m == 1) {
    tab.syn   <- tab.syn[[1]]
    tab.zdiff <- tab.zdiff[[1]]
    tabd      <- tabd[[1]]  
  } 
  
  # If all frequency tables for original data are the same, keep only one 
  if (m > 1 && (all(sapply(object$syn, nrow) == sum(object$n)) || all(sapply(object$syn, nrow) == object$n))) tabd <- tabd[[1]]

  res <- list(m = m, 
              UtabFT  = UtabFT, 
              UtabVW  = UtabVW, 
              df      = df, 
              ratioFT = ratioFT, 
              pvalFT  = pvalFT,
              ratioVW = ratioVW, 
              pvalVW  = pvalVW,
              nempty  = unlist(nempty), 
              # nobs.missings = nobs.missings,
              # nsyn.missings = nsyn.missings,
              tab.obs = tabd, 
              tab.syn = tab.syn, 
              tab.zdiff = tab.zdiff,
              digits    = digits, 
              print.zdiff  = print.zdiff, 
              print.tables = print.tables, 
              n = object$n)
  
  class(res) <- "utility.tab"
  return(res)
}


###-----group_num----------------------------------------------------------
# function to categorise continuous variables

group_num <- function(x1, x2, n = 5, style = "fisher", cont.na = NA, ...) {
  
  # Categorise 2 continous variables into factors of n groups 
  # with same groupings determined by the first one
  
  if (!is.numeric(x1) | !is.numeric(x2)) stop("x1 and x2 must be numeric.\n", call. = FALSE)
  
  # Select non-missing(nm) values 
  x1nm <- x1[!(x1 %in% cont.na)]
  x2nm <- x2[!(x2 %in% cont.na)]
  # Derive breaks
  my_breaks <- classIntervals(rbind(x1nm, x2nm), n = n, style = style, ...)$brks
  my_levels <- c(levels(cut(x1nm, breaks = my_breaks, 
                 dig.lab = 8, right = FALSE, include.lowest = TRUE)),
                 cont.na[!is.na(cont.na)])
  # Apply groupings to non-missing data
  x1[!(x1 %in% cont.na)] <- as.character(cut(x1nm, breaks = my_breaks, 
                            dig.lab = 8, right = FALSE, include.lowest = TRUE))
  x2[!(x2 %in% cont.na)] <- as.character(cut(x2nm, breaks = my_breaks, 
                            dig.lab = 8, right = FALSE, include.lowest = TRUE))
  x1 <- factor(x1, levels = my_levels)
  x2 <- factor(x2, levels = my_levels)
  return(list(x1,x2))  
  
}  




