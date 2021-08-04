###-----utility.gen--------------------------------------------------------
utility.gen <- function(object, data, ...) UseMethod("utility.gen")


###-----utility.gen.default------------------------------------------------
utility.gen.default <- function(object, ...)
  stop("No compare method associated with class ", class(object), call. = FALSE)


###-----utility.gen.data.frame---utility.gen.list--------------------------
utility.gen.data.frame <- utility.gen.list <- 
    function(object, data, 
             not.synthesised = NULL, cont.na = NULL,
             method = "logit", maxorder = 1, k.syn = FALSE, 
             tree.method = "rpart", resamp.method = NULL, nperms = 50,   
             cp = 1e-3, minbucket = 5, mincriterion = 0, vars = NULL,  
             aggregate = FALSE, maxit = 200, ngroups = NULL, print.flag = TRUE, 
             print.every = 10, digits = 6, print.zscores = FALSE, zthresh = 1.6,    
             print.ind.results = TRUE, print.variable.importance = FALSE,  ...)
  {
  
  if (is.null(data)) stop("Requires parameter 'data' to give name of the real data.\n\n",  call. = FALSE)
  if (is.null(object)) stop("Requires parameter 'object' to give name of the synthetic data.\n\n",  call. = FALSE)   
  
  if (is.list(object) & !is.data.frame(object)) m <- length(object)
  else if (is.data.frame(object)) m <- 1
  else stop("object must be a data frame or a list of data frames.\n", call. = FALSE)
  
  # sort out cont.na to make it into a complete named list
  cna <- cont.na
  cont.na <- as.list(rep(NA, length(data)))
  names(cont.na) <- names(data)
  if (!is.null(cna)) {
    if (!is.list(cna) | any(names(cna) == "") | is.null(names(cna))) 
      stop("Argument 'cont.na' must be a named list with names of selected variables.", call. = FALSE)  
    if (any(!names(cna) %in% names(data))) stop("Names of the list cont.na must be variables in data.\n", call. = FALSE)
    for (i in 1:length(cna)) {
      j <- (1:length(data))[names(cna)[i] == names(data)]
      cont.na[[j]] <- unique(c(NA,cna[[i]]))
    }
  }
  
  syn.method = rep("ok", length(data))
  if (!is.null(not.synthesised)) {
    if (!is.null(not.synthesised) && !all(not.synthesised %in% names(data))) stop("not.synthesised must be names of variables in data.\n", call. = FALSE)
    syn.method[names(data) %in% not.synthesised] <- ""
  }
  
  object <- list(syn = object, m = m, strata.syn = NULL, method = syn.method, cont.na = cont.na)
  class(object ) <- "synds"
  
  res <- utility.gen.synds(object = object, data = data, method = method, maxorder = maxorder, print.flag = print.flag, 
                     tree.method = tree.method, resamp.method = resamp.method, k.syn = k.syn,    
                     nperms = nperms, cp = cp, minbucket = minbucket, mincriterion = mincriterion,  
                     vars = vars, aggregate = aggregate, maxit = maxit, ngroups = ngroups, 
                     print.every = print.every, digits = digits, print.zscores = print.zscores, zthresh = zthresh, 
                     print.ind.results = print.ind.results, print.variable.importance = print.variable.importance )
  res$call <- match.call()
  return(res)
}


###-----utility.gen-------------------------------------------------------- 
utility.gen.synds <- function(object, data, method = "logit", maxorder = 1,  
                        k.syn = FALSE, tree.method = "rpart", 
                        resamp.method = NULL,nperms = 50, cp = 1e-3, 
                        minbucket = 5, mincriterion = 0, vars = NULL,   
                        aggregate = FALSE, maxit = 200, ngroups = NULL, 
                        print.flag = TRUE, print.every = 10, 
                        digits = 6, print.zscores = FALSE, zthresh = 1.6, 
                        print.ind.results = TRUE, 
                        print.variable.importance = FALSE,  ...)
{

 m   <- object$m

 # Check input parameters 
 if (is.null(method) || length(method) != 1 || is.na(match(method, c("cart", "logit")))) {
   stop("Invalid 'method' type - must be either 'logit' or 'cart'.\n", call. = FALSE) 
 }
 if (!is.null(resamp.method) && is.na(match(resamp.method, c("perm", "pairs")))) {
   stop("Invalid 'resamp.method' type - must be NULL, 'perm' or 'pairs'.\n", call. = FALSE)
 }
 if (aggregate == TRUE & method != "logit") stop("Aggregation only works for 'logit' method.\n", call. = FALSE) 
 if (is.null(data)) stop("Requires parameter 'data' to give name of the real data.\n",  call. = FALSE)
 if (!class(object) == "synds") stop("Object must have class 'synds'.\n", call. = FALSE)
 if (k.syn & !is.null(resamp.method) && resamp.method == "pairs") stop('\nresamp.method = "pairs" will give the wrong answer when k.syn is TRUE.\n', call. = FALSE)

 # Check selected variables and make observed and synthetic comparable
 if (!(is.null(vars))) {               
   if (!(all(vars %in% names(data))))	stop("Some 'vars' specified not in data.\n", call. = FALSE)
   data <- data[, vars, drop = FALSE]
   if (m == 1) {
     if (!all(vars %in% names(object$syn))) stop("Some 'vars' specified not in synthetic data.\n", call. = FALSE)
     else object$syn <- object$syn[, vars, drop = FALSE ]
   } else {
     if (!all(vars %in% names(object$syn[[1]]))) stop("Some 'vars' specified not in synthetic data.\n", call. = FALSE)
     else object$syn <- lapply(object$syn, "[", vars)
   }
 } else {
   if (m == 1) vars <- names(object$syn) else vars <- names(object$syn[[1]])
   if (!all(vars %in% names(data))) stop("Some variables in synthetic data not in original data.\n", call. = FALSE)
   else data <- data[, vars]  # make data match synthetic
 } 

 # get cont.na and method parameters for stratified synthesis    
 if (!is.null(object$strata.syn)) {
   cna <- object$cont.na[1,]
   syn.method <- object$method[1,] 
 } else {
   cna <- object$cont.na 
   syn.method <- object$method
 } 

 cna <- cna[names(cna) %in% vars]

 for ( i in 1:length(cna)) { 
   nm <- names(cna)[i]
   vals <- unique(cna[[i]][!is.na(cna[[i]])])  # get variables with cont.na other than missing
   if (length(vals) > 0){
     for (j in 1:length(vals))
       n_cna <- sum(vals[j] == data[,nm] & !is.na(data[,nm]))
     if (n_cna == 0) stop("\nValue ", vals[j], " identified as denoting a special or missing in cont.na for ",nm, " is not in data.\n",sep = "", call. = FALSE)
     else if (n_cna < 10 & print.flag) cat ("\nWarning: Only ",n_cna ," record(s) in data with value ",vals[j]," identified as denoting a missing value in cont.na for ",nm, "\n\n", sep = "")
   }
 }
 # Check whether some variables are unsynthesised  
 incomplete <- FALSE
 nsynthd <- length(vars)
 unsyn.vars <- names(syn.method)[syn.method == ""]  # identify unsynthesised
 if (any(vars %in% unsyn.vars) & !is.null(unsyn.vars)) {
   notunsyn <- vars[!vars %in% unsyn.vars]  # synthesised vars
   if (!all(unsyn.vars %in% vars)) stop("Unsynthesised variables must be a subset of variables contributing to the utility measure.\n", call. = FALSE)
   if ( all(vars %in% unsyn.vars)) stop("Utility measure impossible if all in vars are unsynthesised.\n", call. = FALSE)
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
 else if (!is.null(resamp.method) && resamp.method == "pairs" & m == 1) stop('resamp.method = "pairs" needs a synthesis with m > 1, m = 10 suggested.\n', call. = FALSE)
 
 # Drop any single value columns 
 leneq1 <- function(x) length(table(as.numeric(x[!is.na(x)]), useNA = "ifany")) %in% (0:1) 
 
 dchar <- sapply(data,is.character)   ##GR
 if (any(dchar == TRUE)) for ( i in 1:dim(data)[2]) if (dchar[i] == TRUE) data[,i] <- factor(data[,i])
 dout <- sapply(data,leneq1)
 if (m == 1) sout <- sapply(object$syn,leneq1)
 else  sout <- sapply(object$syn[[1]],leneq1)
 dout <- dout & sout
 if (any(dout == TRUE) & print.flag) {
   cat("Some columns with single values or all missing values in original and synthetic\nexcluded from utility comparisons (excluded variables: ",
       paste(names(data)[dout], collapse = ", "), ").\n", sep = "")
   data <- data[,!dout]
   if (m == 1) object$syn <- object$syn[, !dout, drop = FALSE]
   else object$syn <- lapply(object$syn, "[", !dout)
 }

 # Numeric variables
 numvars <- (1:dim(data)[2])[sapply(data,is.numeric)]
 names(numvars) <- names(data)[numvars]
 # If ngroups != NULL divide numeric variables into ngroups
 data0 <- data  # to save if m > 1
 if (!is.null(ngroups)) {
   for (i in numvars) {
     if (m == 1) {
       groups <- group_num(data[,i], object$syn[,i], ngroups, cont.na = cna)
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
 
 for (i in numvars) {
   if (anyNA(data[,i]) & is.null(ngroups)) {
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
   if (any(!is.na(cna[[i]]))  & is.null(ngroups)) {               
     cna[[i]] <- cna[[i]][!is.na(cna[[i]])]
     for (j in 1:length(cna[[i]])) {
       newname <- paste(names(data)[i], "cna",j, sep = "_")
       data <- data.frame(data, 1*(data[,i] == cna[[i]][j]))
       data[data[,i] == cna[[i]][j], i] <- 0
       names(data)[length(data)] <- newname
     }
     if (m == 1) {
       for (j in 1:length(cna[[i]])) {
         newname <- paste(names(object$syn)[i], "cna",j, sep = "_")
         object$syn <- data.frame(object$syn, 1*(object$syn[,i] == cna[[i]][j]))
         object$syn[object$syn[,i] == cna[[i]][j], i] <- 0
         names(object$syn)[length(object$syn)] <- newname
       }
     } else {
       for (k in 1:m) {
         for (j in 1:length(cna[[i]])) {
           newname <- paste(names(object$syn[[k]])[i], "cna",j, sep = "_")
           object$syn[[k]] <- data.frame(object$syn[[k]], 1*(object$syn[[k]][,i] == cna[[i]][j]))
           object$syn[[k]][object$syn[[k]][,i] == cna[[i]][j], i] <- 0
           names(object$syn[[k]])[length(object$syn[[k]])] <- newname
         }
       }
     }
   }  
 }

 nnosplits <- NA  # for 'logit' model with no perm or pairs
 
 # Function for getting propensity scores 
 # --------------------------------------
 propcalcs <- function(syndata, data) {
 
   n1 <- dim(data)[1]
   n2 <- dim(syndata)[1]
   N <- n1 + n2                         
   cc <- n2 / N                         
   if (k.syn) cc <- 0.5                 
    
   df.prop <- rbind(syndata, data)  # make data frame for calculating propensity score
   df.prop <- data.frame(df.prop, t = c(rep(1,n2), rep(0,n1)))
   
   # remove any levels of factors that don't exist in data or syndata
   catvars <- ( 1:(dim(df.prop)[2]) )[ sapply(df.prop,is.factor)]  
   for (i in catvars) {
     if (any(table(df.prop[,i]) == 0)) {
       df.prop[,i] <- as.factor(as.character(df.prop[,i]))
       if (print.flag) cat("Empty levels of factor(s) for variable ", names(df.prop)[i]," removed \n" )
     }
   }
    
   if (aggregate == TRUE) {
     aggdat <- aggregate(df.prop[,1], by = df.prop, FUN = length)
     wt <- aggdat$x
     aggdat <- aggdat[, -dim(aggdat)[2]]
   }
    
   if (method == "logit" ) {
     if (maxorder > 4 & print.flag) cat("Maximum order of interactions over 4 may cause computational problems.\n")
     if (maxorder >= 1) logit.int <- as.formula(paste("t ~ .^", maxorder + 1))
     else logit.int <- as.formula(paste("t ~ ."))
     
     if (aggregate == TRUE) fit <- glm(logit.int, data = aggdat, family = "binomial", 
                                       control = list(maxit = maxit), weights = wt)
     else fit <- suppressWarnings(glm(logit.int, data = df.prop, family = "binomial",
                     control = list(maxit = maxit)))                            ## Error changed to warning
     if (fit$converged == FALSE) cat("\nConvergence failed.\n")
     
     # Get number of parameters that involve synthesised variables
     score <- predict(fit, type = "response")
     if (incomplete == FALSE) km1 <- length(fit$coefficients[!is.na(fit$coefficients)]) - 1  # To allow for non-identified coefficients
     else {
       namescoef <- names(fit$coefficients)
       coefOK <- rep(FALSE, length(namescoef))
       for (nn in notunsyn) coefOK[grepl(nn, namescoef)] <- TRUE
       km1 <- sum(coefOK & print.flag)
       if (m == 1 || (m > 1 & j == 1)) cat("Expectation of utility uses only coefficients involving synthesised variables: ",
                                           km1, " from ", length(fit$coefficients) - 1, "\n", sep = "")
     }
     # one more cofficient (intercept needed if k.syn TRUE)
     if (k.syn) km1 <- km1 + 1                                              
     if (aggregate == TRUE) { 
       pMSE <- (sum(wt*(score - cc)^2, na.rm = T)) / N  
     } else {
       pMSE <- (sum(   (score - cc)^2, na.rm = T)) / N 
     } 
     if (aggregate == TRUE) pct.correct <- sum(wt[( score > cc & df.prop$t ==1) | ( score <= cc & df.prop$t ==0) ])/N
     else  pct.correct <- sum(( score > cc & df.prop$t ==1) | ( score <= cc & df.prop$t ==0) )/N
     

     pval <- 1 - pchisq(pMSE * N  / cc / (1-cc)^2, km1) 
     pMSEExp   <- km1 * (1 - cc)^2 * cc / N
     utilR     <- pMSE / pMSEExp                          
     fit$data  <- NULL  # to save space
     fit$model <- NULL  # to save space 
     #fit$residuals <- fit$weights <- fit$R <- fit$linear.predictors <- fit$qr <-  NULL   # to save space

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
     pMSE <- (sum((score - cc )^2, na.rm = T)) / N  
     pct.correct <- sum(( score > cc & df.prop$t ==1) | ( score <= cc & df.prop$t ==0) )/N
   }
   
   # Permutation test 
   if (!is.null(resamp.method) && resamp.method == "perm") { # to allow resamp for logit models
     simutil <- rep(0, nperms)
     if (m == 1) j <- 1
     if (j == 1 ) {
       if (print.every == 0 & print.flag) cat("Running ", nperms, " permutations to get NULL pMSE.\n", sep = "")
       else  cat("Running ", nperms," permutations to get NULL pMSE and printing every ", print.every, "th.", sep = "")
     }
     if (print.every > 0 & print.flag) cat("\nsynthesis ", j, ": ", sep = "")
     
     for (i in 1:nperms) {
       if (print.every > 0 & floor(i/print.every) == i/print.every & print.flag)  cat(i, " ", sep = "")
       pdata <- df.prop
       if (!k.syn) pdata$t <- sample(pdata$t)
       else pdata$t <- rbinom(N, 1, 0.5)    ##GR0221#
       
       if (method == "cart") {
         if (tree.method == "rpart") {
           sfit <- rpart(t ~ ., data = pdata, method = 'class', control = rpart.control(cp = cp, minbucket = minbucket))
           score <- predict(sfit)[,2]
         } else if (tree.method == "ctree") {
           sfit <- ctree(t ~ ., data = pdata, 
                         controls = ctree_control(mincriterion = mincriterion, minbucket = minbucket))
           score <- predict(sfit)
         }
         simutil[i] <- (sum((score -cc)^2, na.rm = T)) / N / 2  
       
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
         
         if (sfit$converged == FALSE & print.flag) cat("Warning: Logistic model did not converge in ",
                                         maxit, " iterations.\nYou could try increasing parameter 'maxit'.\n", sep = "")
         score <- predict(sfit, type = "response")
         if (aggregate == TRUE) { 
           simutil[i] <- (sum(wt*(score - cc)^2, na.rm = T)) / N / 2  # reduced by factor of 2 
         } else {
           simutil[i] <- (sum(   (score - cc)^2, na.rm = T)) / N / 2  # reduced by factor of 2 
         } 
       }
     }
     nnosplits <- c(sum(simutil < 1e-8), length(simutil))
     pMSEExp   <- mean(simutil) 
     utilR     <- pMSE / pMSEExp
     pval <- sum(simutil > pMSE) / length(simutil)
     #utilStd <- (pMSE - pMSEExp) / sqrt(var(simutil))  
   }
   if (!is.null(resamp.method) && resamp.method == "pairs") res.ind <- list(pMSE = pMSE, pct.correct = pct.correct * 100, 
            pMSEExp = NA, utilR = NA, fit = fit, nnosplits = nnosplits, pval = NA, df = NA) 
   else if (!is.null(resamp.method) && resamp.method == "perm") res.ind <- list(pMSE = pMSE, pct.correct = pct.correct * 100,
            pMSEExp = pMSEExp, utilR = utilR, fit = fit, nnosplits = nnosplits, pval = pval, df = NA) 
   else res.ind <- list(pMSE = pMSE, pct.correct = pct.correct *100, 
            pMSEExp = pMSEExp, utilR = utilR, fit = fit, nnosplits = nnosplits, pval = pval, df = km1) 

   return(res.ind)
 }
 # -------------------------------------- 
 
  n1 <- nrow(data)

  if (m == 1) {
    n2 <- nrow(object$syn)
    res.ind <- propcalcs(object$syn, data)
    res <- list(call = match.call(), m = m, method = method, tree.method = tree.method, 
                resamp.method = resamp.method, maxorder = maxorder, vars = vars,
                k.syn = k.syn, aggregate = aggregate, maxit = maxit, 
                ngroups = ngroups, mincriterion = mincriterion, 
                nperms = nperms, df = res.ind$df, incomplete = incomplete,       
                pMSE = res.ind$pMSE, pMSEExp = res.ind$pMSEExp, utilR = res.ind$utilR, 
                pval = res.ind$pval, pct.correct = res.ind$pct.correct,
                fit = res.ind$fit, nnosplits = res.ind$nnosplits,   
                digits = digits, print.ind.results = print.ind.results, 
                print.zscores = print.zscores, zthresh = zthresh, 
                print.variable.importance = print.variable.importance)
  } else {
    n2 <- nrow(object$syn[[1]]) 
    pMSE <- pMSEExp <- utilR  <- pval <- pct.correct <- rep(NA, m)
    fit <- nnosplits <- as.list(1:m)
    if (!is.null(resamp.method) && resamp.method == "pairs") {
      kk <- 0
      simvals <- rep(NA, m*(m - 1)/2)
    }
    for (j in 1:m) {
      res.ind <- propcalcs(object$syn[[j]], data)
      if (method == "logit" & is.null(resamp.method)) {
        if (j == 1 & print.flag) cat("Fitting syntheses: ")
        if (print.flag) cat(j, " ", sep = "")
      }

      nnosplits[[j]] <- res.ind$nnosplits   # to stop failure for m > 1 logit models
      pMSE[j] <- res.ind$pMSE;  pval[j] <- res.ind$pval; pMSEExp[j] <- res.ind$pMSEExp    
      utilR[j] <- res.ind$utilR ; pct.correct[j] <- res.ind$pct.correct
      
      if (!is.null(resamp.method) && resamp.method == "pairs") {
        if (j == 1) {
          if (print.every == 0 & print.flag) cat("Simulating NULL pMSE from ", m*(m - 1)/2, " pairs.", sep = "")
          else cat("Simulating NULL pMSE from ", m*(m - 1)/2, " pairs, printing every ", print.every, "th:\n", sep = "")
          if (m*(m - 1)/2 < 6 ) cat("Number of pairs too low, suggest increasing number of syntheses (m).\n")
        }
        if (j < m) {
          for (jj in (j + 1):(m)) {
            kk <- kk + 1
            if (print.every > 0 & floor(kk/print.every) == kk/print.every & print.flag) cat(kk," ",sep = "")
            simvals[kk] <- propcalcs(object$syn[[j]], object$syn[[jj]])$pMSE
          }
        }
      } 
      fit[[j]] <- res.ind$fit
    }
    
    if (!is.null(resamp.method) && resamp.method == "pairs") {
      nnosplits<- c(sum(simvals < 1e-8), length(simvals))
      #if (method == "logit")
      simvals <- simvals / 2
      pval[j] <- sum(simvals > pMSE[j])/ length(simvals) 
      exp <- mean(simvals)
      for (j in 1:m) {
        pMSEExp[j] <- exp 
        utilR[j]   <- pMSE[j] / exp
      }
    }
    res <- list(call = match.call(), m = m, method = method, tree.method = tree.method,
                resamp.method = resamp.method, maxorder = maxorder, vars = vars, 
                k.syn = k.syn, aggregate = aggregate, maxit = maxit, 
                ngroups = ngroups, mincriterion = mincriterion, 
                nperms = nperms, df = res.ind$df, incomplete = incomplete,
                pMSE = pMSE, pMSEExp = pMSEExp, utilR = utilR, pct.correct = pct.correct ,
                pval = pval, fit = fit, nnosplits = nnosplits,  
                digits = digits, print.ind.results = print.ind.results, print.zscores = print.zscores, 
                zthresh = zthresh, print.variable.importance = print.variable.importance)
  }

  class(res) <- "utility.gen"
  res$call <- match.call()
  return(res)
}


###-----utility.tab--------------------------------------------------------
utility.tab <- function(object, data, ...) UseMethod("utility.tab")


###-----utility.tab.default------------------------------------------------
utility.tab.default <- function(object, ...)
  stop("No compare method associated with class ", class(object), call. = FALSE)


###-----utility.tab.data.frame---utility.tab.list--------------------------
utility.tab.data.frame <- utility.tab.list <- 
        function(object, data, vars = NULL, cont.na = NULL, 
                 ngroups = 5, useNA = TRUE, 
                 print.tables = length(vars) < 4, 
                 print.stats = 'VW', print.zdiff = FALSE, 
                 digits = 3, k.syn = FALSE, ...) 
{
  
  if (is.null(data)) stop("Requires parameter 'data' to give name of the real data.\n\n",  call. = FALSE)
  if (is.null(object)) stop("Requires parameter 'object' to give name of the synthetic data.\n\n",  call. = FALSE)   
  
  if (is.list(object) & !is.data.frame(object)) m <- length(object)
  else if (is.data.frame(object)) m <- 1
  else stop("object must be a data frame or a list of data frames.\n", call. = FALSE)
  
  # sort out cont.na to make it into a complete named list
  cna <- cont.na
  cont.na <- as.list(rep(NA, length(data)))
  names(cont.na) <- names(data)
  if (!is.null(cna)) {
    if (!is.list(cna) | any(names(cna) == "") | is.null(names(cna))) 
      stop("Argument 'cont.na' must be a named list with names of selected variables.", call. = FALSE)  
    if (any(!names(cna) %in% names(data))) stop("Names of the list cont.na must be variables in data.\n", call. = FALSE)
    for (i in 1:length(cna)) {
      j <- (1:length(data))[names(cna)[i] == names(data)]
      cont.na[[j]] <- unique(c(NA,cna[[i]]))
    }
  }
  
  object <- list(syn = object, m = m, cont.na = cont.na)
  class(object ) <- "synds"
  
  res <- utility.tab.synds(object = object, data = data, vars = vars, ngroups = ngroups,
                     useNA = useNA, print.tables = print.tables,
                     print.stats = print.stats, print.zdiff = print.zdiff, 
                     digits = digits, k.syn = k.syn, ...) 
  return(res)
}


###-----utility.tab--------------------------------------------------------
utility.tab.synds <- function(object, data, vars = NULL, ngroups = 5, 
                              useNA = TRUE, print.tables = length(vars) < 4, 
                              print.stats = "VW", print.zdiff = FALSE,
                              digits = 3, k.syn = FALSE, ...) 
{
  vars <- unique(vars)
    
  # CHECKS 
  #---------
  if (is.null(data)) 
    stop("Requires parameter 'data' to give name of the real data.\n", 
         call. = FALSE)
  if (!is.data.frame(data)) 
    stop("Data must have class 'data.frame'.\n", call. = FALSE)
  if (!class(object) == "synds") 
    stop("Object must have class 'synds'.\n", call. = FALSE)
  if (is.null(vars)) stop("Need to set variables with vars parameter.\n", call. = FALSE) else if 
    (!(all(vars %in% names(data)))) stop("Unrecognized variable(s) in vars parameter: ", 
     paste(vars[!(vars %in% names(data))], collapse = ", "), call. = FALSE)
  #---------

  data  <- data[, vars, drop = FALSE]
  nvars <- ncol(data)  
  data.orig <- data

  # get cont.na parameters for stratified synthesis    
  # --------
  if (!is.null(object$strata.syn)) {
    # cna <- apply(object$cont.na, 2, function(y) {unlist(unique(y))})
    cna <- object$cont.na[1,]
  } else {
    cna <- object$cont.na
  } 
  cna <- cna[vars]  

  m <- object$m
  if (m == 1) syndata <- list(object$syn) else syndata <- object$syn   
  syndata <- lapply(syndata,'[',vars)
  
  df <- ratioFT <- pvalFT <- ratioVW <- pvalVW <- ratioG <- pvalG <- 
    JSD <- pct.correct <- UtabG <- UtabFT <- UtabVW <- 
    nempty <- vector("numeric", m)
  tab.syn <- tab.zdiff <- tabd <- vector("list", m) 
  
  # nobs.missings <- sum(apply(is.na(data), 1, sum)>0) 
  # nsyn.missings <- vector("numeric", m)

  for (i in 1:m) {
    data <- data.orig
    # nsyn.missings[i] <- sum(apply(is.na(syndata[[i]]), 1, sum)>0) 
    # make all variables into factors
    for (j in 1:nvars) {
      if (is.numeric(data[,j])) {
        grpd <- group_num(data[,j], syndata[[i]][,j], 
          n = ngroups, cont.na = cna[[j]], ...)
        data[,j] <- grpd[[1]]; syndata[[i]][,j] <- grpd[[2]]
      } else if (is.character(data[,j])) {
        data[,j] <- factor(data[,j])
        syndata[[i]][,j] <- factor(syndata[[i]][,j], 
                                   levels = levels(data[,j]))
      }
      if (any(is.na(data[,j])) & useNA) {
        # makes missings into part of factors if present 
        data[,j] <- addNA(data[,j]) 
        syndata[[i]][,j] <- addNA(syndata[[i]][,j]) 
      }
    }
    
    # Adjustment for synthetic data that have different size than original data
    if (!k.syn) tabd[[i]] <- table(data) * (nrow(syndata[[i]])/nrow(data))
           else tabd[[i]] <- table(data)
    totcells        <- length(tabd[[i]])
    tab.syn[[i]]    <- table(syndata[[i]])
    nempty[i]       <- sum(tabd[[i]] + tab.syn[[i]] == 0)  
    if (!k.syn) df[i] <- totcells - nempty[i] - 1
           else df[i] <- totcells - nempty[i]
    diff            <- (tab.syn[[i]] - tabd[[i]])
    expect          <- (tab.syn[[i]] + tabd[[i]])/2
    UtabFT[i]       <- 4*sum((tab.syn[[i]]^(0.5) - tabd[[i]]^(0.5))^2)
    tabsq           <- diff^2/expect
    tabsq[expect == 0] <- 0
    UtabVW[i]       <- sum(tabsq)
    ratioFT[i]      <- UtabFT[i] / df[i]
    pvalFT[i]       <- 1 - pchisq(UtabFT[i], df[i])
    # stdFT[i] <- (UtabFT[i] - df[i]) / sqrt(2*df[i])
    ratioVW[i]      <- UtabVW[i] / df[i]
    pvalVW[i]       <- 1 - pchisq(UtabVW[i], df[i])
    # stdVW[i] <- (UtabVW[i] - df[i]) / sqrt(2*df[i])
    tab.zdiff[[i]]  <- diff/sqrt(expect)
    ptabd <- tabd[[i]] / sum(tabd[[i]])
    ptabs <- tab.syn[[i]] / sum(tab.syn[[i]])
    phalf <- (ptabd + ptabs) *0.5
    ## Jensen-shannon divergence
    JSD[i] <- (sum((ptabd * log2(ptabd/phalf))[ptabd > 0])/2 +
               sum((ptabs * log2(ptabs/phalf))[ptabs > 0])/2)^0.5
    ## calculate pct.correct
    predsyn <- (ptabs > ptabd)
    # if (!k.syn) pct.correct[i] <- (sum(ptabs[predsyn]) + sum(ptabd[!predsyn])) * 50 else      
    pct.correct[i] <- (sum(tab.syn[[i]][predsyn]) + sum(tabd[[i]][!predsyn])) / (sum(tab.syn[[i]]) + sum(tabd[[i]])) * 100

    ## symmetric likelihood ratio chisq
    UtabG[i]  <-  sum((tabd[[i]] * log(tabd[[i]]/tab.syn[[i]]))[tabd[[i]] > 0 &  tab.syn[[i]] > 0] ) +   
      sum((tab.syn[[i]] * log(tab.syn[[i]]/tabd[[i]]))[tabd[[i]] > 0 &  tab.syn[[i]] > 0])
    ratioG[i] <- UtabG[i] / df[i]
    pvalG[i]  <- 1 - pchisq(UtabG[i], df[i])
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
              UtabG   = UtabG,
              df      = df, 
              ratioFT = ratioFT, 
              pvalFT  = pvalFT,
              ratioVW = ratioVW, 
              pvalVW  = pvalVW,
              ratioG  = ratioG,
              pvalG   = pvalG,
              JSD     = JSD,
              pct.correct = pct.correct,
              nempty   = unlist(nempty), 
              # nobs.missings = nobs.missings,
              # nsyn.missings = nsyn.missings,
              tab.obs = tabd, 
              tab.syn = tab.syn, 
              tab.zdiff = tab.zdiff,
              digits    = digits, 
              print.stats  = print.stats, 
              print.zdiff  = print.zdiff, 
              print.tables = print.tables, 
              n = sum(object$n),
              k.syn = k.syn)
  
  class(res) <- "utility.tab"
  return(res)
}


###-----group_num----------------------------------------------------------
# function to categorise continuous variables

group_num <- function(x1, x2, n = 5, style = "quantile", cont.na = NA, ...) {
  
  # Categorise 2 continuous variables into factors of n groups 
  # with same groupings determined by the first one
  
  if (!is.numeric(x1) | !is.numeric(x2)) stop("x1 and x2 must be numeric.\n", 
                                              call. = FALSE)
  
  # Select non-missing(nm) values 
  x1nm <- x1[!(x1 %in% cont.na) & !is.na(x1)]
  x2nm <- x2[!(x2 %in% cont.na) & !is.na(x2)]

  # Derive breaks
  my_breaks <- unique(suppressWarnings(classIntervals(c(x1nm, x2nm), 
                                          n = n, style = style, ...))$brks)

  my_levels <- c(levels(cut(x1nm, breaks = my_breaks, 
                 dig.lab = 8, right = FALSE, include.lowest = TRUE)),
                 cont.na[!is.na(cont.na)])
  
  # Apply groupings to non-missing data
  x1[!(x1 %in% cont.na) & !is.na(x1)] <- as.character(cut(x1nm, 
    breaks = my_breaks, dig.lab = 8, right = FALSE, include.lowest = TRUE))
  x2[!(x2 %in% cont.na) & !is.na(x2)] <- as.character(cut(x2nm, 
    breaks = my_breaks, dig.lab = 8, right = FALSE, include.lowest = TRUE))
  x1 <- factor(x1, levels = my_levels)
  x2 <- factor(x2, levels = my_levels)
  
  return(list(x1,x2))  
}  




