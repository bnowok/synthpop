###-----print.synds--------------------------------------------------------

print.synds <- function(x, ...){
  cat("Call:\n($call) ")
  print(x$call)
  cat("\nNumber of synthesised data sets: \n($m) ",x$m,"\n")  
  if (x$m == 1) {
    cat("\nFirst rows of synthesised data set: \n($syn)\n")
    print(head(x$syn))
  } else {
    cat("\nFirst rows of first synthesised data set: \n($syn)\n")
    print(head(x$syn[[1]]))
  }    
  cat("...\n")
  cat("\nSynthesising methods: \n($method)\n")
  print(x$method)
  cat("\nOrder of synthesis: \n($visit.sequence)\n")
  print(x$visit.sequence)
  cat("\nMatrix of predictors: \n($predictor.matrix)\n")
  print(x$predictor.matrix)     
  invisible(x)
}


###-----summary.synds------------------------------------------------------

summary.synds <- function(object, msel = NULL, 
  maxsum = 7, digits = max(3, getOption("digits") - 3), ...){
  if (!is.null(msel) & !all(msel %in% (1:object$m))) 
    stop("Invalid synthesis number(s)", call. = FALSE)

  sy <- list(m = object$m, msel = msel, method = object$method)

  if (object$m == 1) {
    sy$result <- summary(object$syn,...)
  } else if (is.null(msel)) {
    zall <- vector("list",object$m) 
    for (i in 1:object$m) zall[[i]] <- lapply(object$syn[[i]], summary,
      maxsum = maxsum, digits = digits, ...)
    zall.df <- Reduce(function(x,y) mapply("rbind",x,y),zall)
    meanres <- lapply(zall.df, function(x) apply(x,2,mean))
    sy$result <- summary.out(meanres)
  } else if (length(msel) == 1) {
    sy$result <- summary(object$syn[[msel]],...)
  } else {
    for (i in (1:length(msel))) {
      sy$result[[i]] <- summary(object$syn[[msel[i]]],...)
    }
  }
  class(sy) <- "summary.synds"
  return(sy)
}


###-----print.summary.synds------------------------------------------------

print.summary.synds <- function(x, ...){

 if (x$m == 1) {
   cat("Synthetic object with one synethesis using methods:\n")
   print(x$method)
   cat("\n")
   print(x$result)
 } else if (is.null(x$msel)) {
   cat("Synthetic object with ",x$m," syntheses using methods:\n",sep = "")
   print(x$method)
   cat("\nSummary (average) for all synthetic data sets:\n",sep = "")
   print(x$result)  
 } else if (length(x$msel) == 1) {
   cat("Synthetic object with ",x$m," syntheses using methods:\n",sep = "")
   print(x$method)
   cat("\nSummary for synthetic data set ",x$msel,":\n",sep = "")
   print(x$result)
 } else {
   cat("Synthetic object with ",x$m," syntheses using methods:\n",sep = "")
   print(x$method)
   for (i in (1:length(x$msel))) {
     cat("\nSummary for synthetic data set ",x$msel[i],":\n",sep = "")
     print(x$result[[i]])
   }
 }
 invisible(x)
}


###-----mcoefvar--------------------------------------------------
# Arrange coefficients from all m syntheses in a matrix
# (same with their variances). 
# [used in lm.synds and glm.synds function]

mcoefvar <- function(analyses, ...) {
  m <- length(analyses)
  if (m == 1) {
    matcoef <- mcoefavg <- analyses[[1]]$coefficients[,1]
    matvar  <- mvaravg  <- analyses[[1]]$coefficients[,2]^2
  } else {
    namesbyfit <- lapply(lapply(analyses,coefficients),rownames)
    allnames <- Reduce(union,namesbyfit)
    matcoef <- matvar <- matrix(NA, m, length(allnames))
    dimnames(matcoef)[[2]] <- dimnames(matvar)[[2]] <- allnames
    for (i in 1:m) {
      pos <- match(namesbyfit[[i]],allnames)
      matcoef[i,pos] <- analyses[[i]]$coefficients[,1]
      matvar[i,pos] <- analyses[[i]]$coefficients[,2]^2
    }
    mcoefavg <- apply(matcoef, 2, mean, na.rm = TRUE)
    mvaravg  <- apply(matvar,  2, mean, na.rm = TRUE)
    #bm <- apply(matcoef,2,var) not needed xpt for partial synthesis
  }
  if (m > 1) rownames(matcoef) <- rownames(matvar) <- paste0("syn=", 1:m)
  return(list(mcoef    = matcoef,  mvar    = matvar, 
              mcoefavg = mcoefavg, mvaravg = mvaravg))
}


###-----lm.synds-----------------------------------------------------------

lm.synds <- function(formula, data, ...)
{
 if (!class(data) == "synds") stop("Data must have class synds\n", call. = FALSE)
 if (is.matrix(data$method)) data$method <- data$method[1,]
 if (is.matrix(data$visit.sequence)) data$visit.sequence <- data$visit.sequence[1,]
 if (data$m > 1) vars <- names(data$syn[[1]])  else  vars <- names(data$syn)  
 n <- sum(data$n)
 if (is.list(data$k)) k <- sum(data$k[[1]]) else k <- sum(data$k)  
 
 call <- match.call()
 fitting.function <- "lm"
 analyses <- as.list(1:data$m)

 # Do the repeated analysis, store the result without data
 if (data$m == 1) {
   analyses[[1]] <- summary(lm(formula, data = data$syn, ...))
 } else {
   for (i in 1:data$m) {
     analyses[[i]] <- summary(lm(formula, data = data$syn[[i]], ...))
   }
 }
 
 # Check validity of inference from vars not in visit sequence or with method ""
 incomplete <- checkcomplete(vars, formula, data$visit.sequence, data$method) 
 
 # Get matrices from coefficients
 allcoefvar <- mcoefvar(analyses = analyses)
      
 # Return the complete data analyses as a list of length m
 object <- list(call = call, mcoefavg = allcoefvar$mcoefavg, 
                mvaravg = allcoefvar$mvaravg, analyses = analyses,  
                fitting.function = fitting.function,
                n = n, k = k, proper = data$proper, 
                m = data$m, method = data$method, incomplete = incomplete,
                mcoef = allcoefvar$mcoef, mvar = allcoefvar$mvar)
 class(object) <- "fit.synds"
 return(object)
}


###-----glm.synds----------------------------------------------------------

glm.synds <- function(formula, family = "binomial", data, ...)
{
 if (!class(data) == "synds") stop("Data must have class synds\n", call. = FALSE)
 if (is.matrix(data$method)) data$method <- data$method[1,]
 if (is.matrix(data$visit.sequence)) data$visit.sequence <- data$visit.sequence[1,]
 if (data$m > 1) vars <- names(data$syn[[1]])  else  vars <- names(data$syn)  
 n <- sum(data$n)
 if (is.list(data$k)) k <- sum(data$k[[1]]) else k <- sum(data$k)  
 
 call <- match.call()
 fitting.function <- "glm"
 analyses <- as.list(1:data$m)
 
 # Do the repeated analysis, store the result without data
 if (data$m == 1) {
   analyses[[1]] <- summary(glm(formula,data = data$syn, family = family, ...))
 } else {
   for (i in 1:data$m) {
     analyses[[i]] <- summary(glm(formula,data = data$syn[[i]], family = family, ...))
   }
 }
 
 # Check completeness for inference from vars not in visit sequence or with method ""
 incomplete <- checkcomplete(vars, formula, data$visit.sequence, data$method) 
 
 # Get matrices from coefficients
 allcoefvar <- mcoefvar(analyses = analyses)

 # Return the complete data analyses as a list of length m
 object <- list(call = call, mcoefavg = allcoefvar$mcoefavg, 
                mvaravg = allcoefvar$mvaravg, analyses = analyses,  
                fitting.function = fitting.function,
                n = n, k = k, proper = data$proper, 
                m = data$m, method = data$method, incomplete = incomplete,
                mcoef = allcoefvar$mcoef, mvar = allcoefvar$mvar)
 class(object) <- "fit.synds"
 return(object)
}


###-----polr.synds-----------------------------------------------------

polr.synds <- function(formula, data, ...)
{
  if (!class(data) == "synds") stop("Data must have class 'synds'.\n", call. = FALSE)
  if (is.matrix(data$method)) data$method <- data$method[1,]
  if (is.matrix(data$visit.sequence)) data$visit.sequence <- data$visit.sequence[1,]
  if (data$m > 1) vars <- names(data$syn[[1]]) else  vars <- names(data$syn)  
  n <- sum(data$n)
  if (is.list(data$k)) k <- sum(data$k[[1]]) else k <- sum(data$k)  
  
  call <- match.call()
  fitting.function <- "polr"
  analyses <- as.list(1:data$m)
  
  # Do the repeated analysis, store the result without data
  for (i in 1:data$m) {
    if (data$m == 1) fit <- polr(formula, data = data$syn, Hess = TRUE, ...)
    else fit <- polr(formula, data = data$syn[[i]], Hess = TRUE, ...)
    ss <- summary(fit)
    analyses[[i]] <- ss
  }

  # Check validity of inference from vars not in visit sequence or with method ""
  incomplete <- checkcomplete(vars, formula, data$visit.sequence, data$method) 
  
  # Get matrices from coefficients
  allcoefvar <- mcoefvar(analyses = analyses)
 
  # Return the complete data analyses as a list of length m
  object <- list(call = call, mcoefavg = allcoefvar$mcoefavg, 
                 mvaravg = allcoefvar$mvaravg, analyses = analyses,  
                 fitting.function = fitting.function,
                 n = n, k = k, proper = data$proper, 
                 m = data$m, method = data$method, incomplete = incomplete,
                 mcoef = allcoefvar$mcoef, mvar = allcoefvar$mvar)
  class(object) <- "fit.synds"
  return(object)
}


###-----multinom.synds-----------------------------------------------------

multinom.synds <- function(formula, data, ...)
{
  if (!class(data) == "synds") stop("Data must have class 'synds'.\n", call. = FALSE)
  if (is.matrix(data$method)) data$method <- data$method[1,]
  if (is.matrix(data$visit.sequence)) data$visit.sequence <- data$visit.sequence[1,]
  if (data$m > 1) vars <- names(data$syn[[1]]) else  vars <- names(data$syn)  
  n <- sum(data$n)
  if (is.list(data$k)) k <- sum(data$k[[1]]) else k <- sum(data$k)  
  
  call <- match.call()
  fitting.function <- "multinom"
  analyses <- as.list(1:data$m)
  
  # Do the repated analysis, store the result without data
  for (i in 1:data$m) {
    if (data$m == 1) fit <- multinom(formula, data = data$syn, Hess = TRUE, ...)
    else fit <- multinom(formula, data = data$syn[[i]], Hess = TRUE, ...)
    ss <- summary(fit)
    analyses[[i]] <- list(coefficients = cbind(as.vector(t(ss$coefficients)),
                                               as.vector(t(ss$standard.errors)),
                                               as.vector(t(ss$coefficients)/t(ss$standard.errors))))
    dd <- dimnames(t(ss$coefficients))
    dimnames(analyses[[i]]$coefficients) <- list(paste(rep(dd[[2]], each = length(dd[[1]])),
                                                       rep(dd[[1]], length(dd[[2]])), sep = ":"),
                                                       c("Estimate", "se", "z value"))
  }

  # Check validity of inference from vars not in visit sequence or with method ""
  incomplete <- checkcomplete(vars, formula, data$visit.sequence, data$method) 
  
  # Get matrices from coefficients
  allcoefvar <- mcoefvar(analyses = analyses)
  
  # Return the complete data analyses as a list of length m
  object <- list(call = call, mcoefavg = allcoefvar$mcoefavg, 
                 mvaravg = allcoefvar$mvaravg, analyses = analyses,  
                 fitting.function = fitting.function,
		 n = n, k = k, proper = data$proper, 
                 m = data$m, method = data$method, incomplete = incomplete,
                 mcoef = allcoefvar$mcoef, mvar = allcoefvar$mvar)
  class(object) <- "fit.synds"
  return(object)
}


###-----checkcomplete------------------------------------------------------
# Used in lm.synds, glm.synds, multinom.synds, and polr.synds

checkcomplete <- function(vars, formula, vs, method)
{  
  inform <- all.vars(formula) # get all variables in formula
  if ("." %in% inform) inform <- vars
  if (any(!inform %in% vars)) stop("Variable(s) in formula (model to be fitted) are not in synthetic data: ",
                                   paste(inform[!inform %in% vars], collapse = ", "), call. = FALSE)
  methin <- method[names(method) %in% inform]
  if (all(methin == "")) cat("No variables in your formula (model to be fitted) have been synthesised.\nIf the data contain exactly the same observations as the original ones,\nresults will be identical.\n")

  order_vs  <- match(names(vs), names(methin))
  order_vs  <- order_vs[!is.na(order_vs)]
  order_oth <- setdiff(1:length(methin), order_vs)
  methin_order <- methin[c(order_vs, order_oth)] 

  blankmeths <- (1:length(methin_order))[methin_order == ""]

  if (!all(blankmeths == (1:length(blankmeths)))){ 
cat(
"**********************************************************",
"\nWARNING: Some variable(s) in formula (model to be fitted)  
are not synthesised and not used in synthesising models
for all other variables:", 
paste(names(methin_order)[blankmeths][!(blankmeths == (1:length(blankmeths)))], collapse = ", "), 
"\nMethods in synthesis order are:\n")
print(methin_order)
cat("Results may not be correct.
**********************************************************\n")
}
  
  incomplete <- !all(blankmeths == (1:length(blankmeths)))
  # incomplete <- length(blankmeths) > 0
  return(incomplete)
}


###-----print.fit.synds----------------------------------------------------

print.fit.synds <- function(x, msel = NULL, ...)
{
  if (!is.null(msel) & !all(msel %in% (1:x$m))) stop("Invalid synthesis number(s): `msel` must be selected from 1:", x$m, call. = FALSE, sep = "")

  #if (x$n != x$k | x$m > 1) cat("Note: To get a summary of results you would expect from the original data\nor for population inference use the summary() function on your fit.\nSee vignette on inference to get more details.\n") 
  if (x$n != x$k | x$m > 1) cat("Note: To get more details of the fit see vignette on inference.\n") 

  cat("\nCall:\n")
  print(x$call)
  if (is.null(msel) & x$m > 1) {
    cat("\nAverage coefficient estimates from", x$m, "syntheses:\n")
    print(x$mcoefavg)
  } else if (x$m == 1) {
    cat("\nCoefficient estimates from a single synthesis:\n")
    print(x$mcoefavg)
  } else {
    cat("\nCoefficient estimates for selected synthetic data set(s):\n")
    print(x$mcoef[msel, , drop = FALSE])
  }
  invisible(x)
}


###-----summary.fit.synds--------------------------------------------------

summary.fit.synds <- function(object, population.inference = FALSE, msel = NULL, 
                              real.varcov = NULL, ...)
{ # df.residual changed to df[2] because didn't work for lm 
  if (!class(object) == "fit.synds") stop("Object must have class fit.synds\n", call. = FALSE)
  m <- object$m
  n <- object$n
  k <- object$k
  incomplete <- object$incomplete

  coefficients <- object$mcoefavg  # mean of coefficients (over m syntheses)
  if (!is.null(real.varcov))  vars <- diag(real.varcov)
  else  vars <- object$mvaravg * k/n  # mean of variances (over m syntheses) * adjustment

## Checks, messages and warnings for population inference
#---
  if (population.inference == TRUE) {
    if (incomplete == TRUE & m == 1) {
      cat("Warning: You have selected population inference when some variables in your model are",
          "\nnot synthesised and when only a single synthetic data set has been created (m = 1).",
          "\nThe correct method for this case requires m > 1, ideally m > 5.",
          "\nTo provide some results calculations proceed as if all variables had been synthesised.\n\n")
      incomplete <- FALSE}
    else if (incomplete == TRUE & m < 5) {
      cat("Note: You have selected population inference for incompletely synthesised data with m = ", m,
          ",\nwhich is smaller than the minimum of 5 recommended. The estimated standard errors of your\ncoefficients may be inaccurate.\n\n", sep = "")
    }
  }
#--- 

## Inference to Q hat
#---
 if (population.inference == FALSE) { 
   result <- cbind(coefficients,
                   sqrt(vars),
                   coefficients/sqrt(vars),
                   2*pnorm(-abs(coefficients/sqrt(vars))))
   colnames(result) <- c("xpct(Beta)", "xpct(se.Beta)", "xpct(z)", "Pr(>|xpct(z)|)")
#--- 

## Population inference to Q
#---   
 } else { 
  
  ## incomplete method  
   if (incomplete == TRUE) {
     bm <- apply(object$mcoef, 2, var)
     result <- cbind(coefficients,
                     sqrt(bm/m + vars),
                     coefficients/sqrt(bm/m + vars),
                     2*pnorm(-abs(coefficients/sqrt(bm/m + vars))))

  ## simple synthesis   
    } else {
      if (object$proper == FALSE) Tf <- vars*(1 + n/k/m) else Tf <- vars*(1 + (n/k + 1)/m)
      result <- cbind(coefficients, 
                      sqrt(Tf), 
                      coefficients/sqrt(Tf),
                      2*pnorm(-abs(coefficients/sqrt(Tf)))) 
    }
    colnames(result) <- c("Beta.syn","se.Beta.syn","z.syn","Pr(>|z.syn|)")
  }
#---
  
 res <- list(call = object$call, proper = object$proper,
             population.inference = population.inference,
             incomplete = incomplete, 
             fitting.function = object$fitting.function,
             m = m, coefficients = result, n = n, k = k, 
             analyses = object$analyses, msel = msel)
 class(res) <- "summary.fit.synds"
 return(res)
}


###-----print.summary.fit.synds--------------------------------------------

print.summary.fit.synds <- function(x, ...) {
 
 if (!is.null(x$msel) & !all(x$msel %in% (1:x$m))) stop("Invalid synthesis number(s)", call. = FALSE)
 #cat("\nNote that all these results depend on the synthesis model being correct.\n")  
 #cat("\nFor details, see package vignette on inference.\n")

 if (x$m == 1) {
   cat("Fit to synthetic data set with a single synthesis. ")
 } else {
   cat("Fit to synthetic data set with ", x$m, " syntheses. ", sep = "")
 }

 if (x$population.inference) {
   if (x$incomplete == TRUE) cat("Inference to population coefficients when\nsome variables in the model are not synthesised. Methods for incomplete/partial\nsynthesis are used.\n")
   else cat("Inference to population\ncoefficients when all variables in the model are synthesised.\n")
 } else {
   cat("Inference to coefficients\nand standard errors that would be obtained from the original data.\n")
   if (x$k != x$n) 
     cat("\nThe synthetic data have a different size (", x$k, ") from the original data (", x$n, "),",
         "\nso the standard errors of the coefficients have been adjusted to estimate",
         "\nthe standard errors from the original data.\n", sep = "")
 }

 cat("\nCall:\n")
 print(x$call)
 cat("\nCombined estimates:\n")
 printCoefmat(x$coefficients)

 if (!is.null(x$msel)) {
   allcoef <- lapply(lapply(x$analyses[x$msel], "[[", "coefficients"), as.data.frame)
   
   estimates <- lapply(allcoef, "[", "Estimate")
   allestimates <- do.call(cbind, estimates)
   
   zvalues <- lapply(allcoef, "[", "z value")
   allzvalues <- do.call(cbind, zvalues)
   
   colnames(allestimates) <- colnames(allzvalues) <- paste0("syn=",x$msel)
   
   cat("\nEstimates for selected syntheses contributing to the combined estimates:\n")
   
   cat("\nCoefficients:\n")
   print(allestimates)
   cat("\nz values:\n")
   print(allzvalues)
   
   # for(i in x$msel) {          
   #   cat("\nsyn=",i,"\n",sep = "")
   #   print(x$analyses[[i]]$coefficients)
   # }
 }      
 invisible(x)
}


###-----print.compare.fit.synds--------------------------------------------

print.compare.fit.synds <- function(x, print.coef = x$print.coef, ...){

  cat("\nCall used to fit models to the data:\n")
  print(x$call)
  if (print.coef == TRUE) {
    cat("\nEstimates for the observed data set:\n")
    print(x$coef.obs)
    cat("\nCombined estimates for the synthesised data set(s):\n")
    print(x$coef.syn)
  }  

  cat("\nDifferences between results based on synthetic and observed data:\n")
    print(cbind.data.frame(x$coef.diff,x$ci.overlap))                                 
  if (x$m == 1) {
    cat("\nMeasures for one synthesis and ", x$ncoef, " coefficients", sep = "") 
  } else {
    cat("\nMeasures for ", x$m, " syntheses and ", x$ncoef, " coefficients", sep = "") 
  }   
  cat("\nMean confidence interval overlap: ", x$mean.ci.overlap)
  cat("\nMean absolute std. coef diff: ", x$mean.abs.std.diff)
  if (!is.null(x$lack.of.fit)){
    cat("\nLack-of-fit: ", x$lack.of.fit,"; p-value ", round(x$lof.pval, 4), " for test that synthesis model is compatible ", sep = "")
    if (x$incomplete == FALSE) cat("\nwith a chi-squared test with ", x$ncoef, " degrees of freedom\n", sep = "")
    else cat("\nwith an F distribution with ",x$ncoef," and ",x$m - x$ncoef," degrees of freedom\n", sep = "") 
  }
  if (!is.null(x$ci.plot)) {
    cat("\nConfidence interval plot:\n")
    print(x$ci.plot)
  }
  invisible(x)
}


###-----print.compare.synds------------------------------------------------

print.compare.synds <- function(x, ...) {
  if (x$stat == "counts") cat("\nComparing counts observed with synthetic\n\n") 
  else cat("\nComparing percentages observed with synthetic\n\n")
  if (class(x$plots)[1] == "gg") {
    if (x$table) print(x$tables) 
    if (x$plot) print(x$plots)
  } else {
    for (i in 1:length(x$tables)) {
      if (x$table) print(x$tables[[i]]) 
      if (x$plot) {
        print(x$plots[[i]])
        if (i < length(x$tables)) {
          cat("Press return for next variable(s): ")
          ans <- readline()
        }
      }
    }
  }
 invisible(x)
}


###-----print.utility.gen--------------------------------------------------

print.utility.gen <- function(x, digits = NULL, zthresh = NULL,    
                              print.zscores = NULL,  
                              print.ind.results = NULL,                               
                              print.variable.importance = NULL, ...){

  if (x$method == "logit"){  
    converged <- ifelse(x$m == 1, x$fit$converged, x$fit[[1]]$converged)
    length.y <- ifelse(x$m == 1, length(x$fit$y)/2, length(x$fit[[1]]$y)/2)
    length.coef <- ifelse(x$m == 1, length(x$fit$coefficients), length(x$fit[[1]]$coefficients)) 
    if (converged == FALSE)  cat("\nWarning: Logistic model did not converge in ", x$maxit,
    " iterations. Check results.\nThis can be due to cells with small numbers or too complicated a model.
    In first case utility measures are probably fine but coefficients may be unstable. 
    In the latter case you may be able to overcome this by increasing the parameter 'maxit'.
    Your combined data has ", length.y, " observations and your model has ",
    length.coef, " parameters.\n", sep = "")
  }
  
  if (is.null(digits))  digits <- x$digits
  if (is.null(zthresh)) zthresh <- x$zthresh
  if (is.null(print.zscores))  print.zscores <- x$print.zscores
  if (is.null(print.ind.results)) print.ind.results <- x$print.ind.results   
  if (is.null(print.variable.importance)) print.variable.importance <- x$print.variable.importance
  
  cat("\nUtility score calculated by method: ", x$method, "\n", sep = "")
  cat("\nCall:\n")
  print(x$call)
  
  if (!is.null(x$resamp.method)) {
    if (x$resamp.method == "perm") cat("\nNull pMSE simulated from a permutation test with ", x$nperm," replications.\n", sep = "")
    else if (x$resamp.method == "pairs") cat("\nNull pMSE simulated from ", x$m*(x$m - 1)/2," pairs of syntheses.\n", sep = "")
    
    if (!is.list(x$nnosplits)) { 
      if (x$nnosplits[1] > 0) cat(
"\n***************************************************************
Warning: null pMSE resamples failed to split ", x$nnosplits[1], " times from ", x$nnosplits[2],
"\n***************************************************************\n", sep = "")
    } else {
      for (ss in 1:x$m) {
        if (x$nnosplits[[ss]][1] > 0) cat("\nSynthesis ", ss, 
          " null pMSE resamples failed to split ", x$nnosplits[[ss]][1],
          " times from ", x$nnosplits[[ss]][2], sep = "")
      }
      cat("\n")
    }
  }
  
  if (x$m > 1) {
    cat("\nMean utility results from ", x$m, " syntheses:\n", 
        "pMSE (propensity score mean square error): ", round(mean(x$pMSE), digits),
        "\nExpected value: ", round(mean(x$pMSEExp), digits), 
        "\nRatio to expected: ", round(mean(x$utilR), digits),
        "\nMean of % correct (percentage correctly predicted): ", round(mean(x$pct.correct), 2),
        "\nMedian of p-values: ", round(median(x$pval), 4), "\n\n", sep = "") 

    if (print.ind.results == TRUE) {
      cat("Individual utility results from ", x$m, " syntheses\n", sep = "")
      tabres <- data.frame(round(x$pMSE, digits), round(x$pMSEExp, digits),
                           round(x$utilR, digits), round(x$pval, 4), round(x$pct.correct, 2))  
      names(tabres) <- c("pMSE", "Expected", "Ratio", "p-value", "% correct")
      print(tabres)
      cat("\n")
    }
    
  } else {
    cat("\nUtility results:",
        "\npMSE (propensity score mean square error): ", round(x$pMSE, digits),              
        "\nExpected value: ", round(x$pMSEExp, digits),  
        "\nRatio to expected: ", round(x$utilR, digits),
        "\n% correct (percentage correctly predicted): ", round(x$pct.correct, 2),
        "\np-value: ", round(x$pval, 4),"\n\n", sep = "")                 
  }
  
  if (print.zscores == TRUE) {
    if (x$method == "cart") {
      cat("\nz-scores not available for CART models\n")
    } else {
      if (x$m > 1) {
        allzscores <- vector("list", x$m) 
        for (i in 1:x$m) allzscores[[i]] <- summary(x$fit[[i]])$coefficients[ ,3] 
        allnames <- unique(unlist(lapply(allzscores, names)))
        allzscores.NA <- lapply(allzscores, "[", allnames) 
        allzscores.NA.df <- do.call(cbind, allzscores.NA)
        zscores <- apply(allzscores.NA.df, 1, mean, na.rm = TRUE)  
        names(zscores) <- allnames
      } else {
        zscores <- summary(x$fit)$coefficients[ ,3]
      }
      
      if (!is.na(zthresh)) { 
        zscores <- zscores[abs(zscores) > zthresh]
        if (length(zscores) == 0) {
          cat("\nNo z-scores (or mean z-scores if m > 1) above threshold of +/-", zthresh,"\n", sep = "")
        } else {
          cat("\nz-scores (or mean z-scores if m > 1) greater than the threshold of +/- ", zthresh, "\n", sep = "")
          print(zscores)
        }  
      } else {
        cat("\nAll z-scores (or mean z-scores if m > 1)\n")
        print(zscores)
      }
    }
  }
  
  if (print.variable.importance == TRUE) {
    if (x$method != "cart" | x$tree.method != "rpart") {
      cat("\nVariable importance only available for CART models using function 'rpart'\n")
    } else {
      cat("\nRelative importance of each variable scaled to add to 100\n" )
      if (x$m == 1) {
        variable.importance <- x$fit$variable.importance
        variable.importance <- round(variable.importance/sum(variable.importance)*100, digits)
        print(variable.importance)
      } else {
        cat("(results for ", x$m, " syntheses)\n", sep = "")
        variable.importance <- vector("list", x$m)
        for (i in 1:x$m) {
          if (is.null(x$fit[[i]]$variable.importance)) x$fit[[i]]$variable.importance <- NA
          variable.importance[[i]] <- x$fit[[i]]$variable.importance
          variable.importance[[i]] <- round(variable.importance[[i]]/sum(variable.importance[[i]])*100, digits)
        }
        allnames <- unique(unlist(lapply(variable.importance, names)))
        all.vars <- lapply(variable.importance, "[", allnames) 
        all.vars.importance <- do.call(rbind, all.vars)
        colnames(all.vars.importance) <- allnames
        rownames(all.vars.importance) <- 1:x$m
        print(all.vars.importance)
      }
    }
  }
  invisible(x)
  
}


###-----print.utility.tab--------------------------------------------------

print.utility.tab <- function(x, print.tables = NULL,
                              print.zdiff = NULL, 
                              print.stats = NULL, 
                              digits = NULL, ...){
  
  if (is.null(print.stats)) print.stats <- x$print.stats
  if (is.null(print.zdiff)) print.zdiff <- x$print.zdiff
  if (is.null(print.tables)) print.tables <- x$print.tables
  if (is.null(digits)) digits <- x$digits
  
  if (print.tables == TRUE) {
    if (is.table(x$tab.obs)) {
      if (x$k.syn) cat("\nObserved not adjusted to match the size of the synthetic data since 'k.syn' set to TRUE: \n($tab.obs)\n")
      else if (sum(x$tab.obs) != x$n & !x$k.syn) { 
        cat("\nObserved adjusted to match the size of the synthetic data: \n($tab.obs)\n")
        print(round(x$tab.obs, digits))
      } else {
        cat("\nObserved: \n($tab.obs)\n")
        print(x$tab.obs)
      }  
    } else {
      #if (sum(x$tabd)/length(x$tabd)!= x$n) 
      cat("\nMean of ",x$m," observed tables ($tab.obs) adjusted to match the size of synthetic data:\n", sep = "")
      meantabd <- apply(simplify2array(x$tab.obs), c(1,2), mean)
      print(round(meantabd, digits))
    } 
    
    if (x$m == 1) {
      cat("\nSynthesised: \n($tab.syn)\n")
      print(x$tab.syn) 
    } else {
      meantab <- apply(simplify2array(x$tab.syn), c(1,2), mean)
      cat("\nMean of ", x$m, " synthetic tables ($tab.syn):\n", sep = "")
      print(round(meantab, digits))
    }
  }
  
  if (print.zdiff == TRUE) {
    cat("\nTable of z-scores for differences: \n($tab.zdiff)\n")
    if (x$m == 1) {
      print(round(x$tab.zdiff, digits)) 
    } else {
      meanzdiff <- apply(simplify2array(x$tab.zdiff), c(1,2), mean)
      cat("\nMean of ",x$m," z-score tables:\n", sep = "")
      print(round(as.table(meanzdiff), digits))
    }
  }
  if (print.stats[1] == "all") print.stats <- c("VW", "FT", "G", "JSD", "pct.correct")
  if (x$k.syn) celldiff <- 0 else celldiff <- 1
  
  if (x$m == 1) {
    cat("\nNumber of cells in each table: ", 
        x$df[1] + x$nempty[1] + celldiff,
        "; Number of cells contributing to utility measures: ", 
        x$df[1] + celldiff,"\n", sep = "")
    cat("\nUtility score results\n")
    if ("VW" %in% print.stats) cat("Voas Williamson (VW): ", round(x$UtabVW, digits), ";",
      " Ratio to degrees of freedom (df): ", round(x$ratioVW, digits), ";",
      " p-value: ", round(x$pvalVW, 4), "\n", sep = "")
    if ("FT" %in% print.stats) cat("Freeman Tukey (FT): ", round(x$UtabFT, digits), ";",
      " Ratio to degrees of freedom (df): ", round(x$ratioFT, digits), ";",
      " p-value: ", round(x$pvalFT, 4), "\n", sep = "")
    if ("G" %in% print.stats) cat("Likelihood ratio (G): ", round(x$UtabG, digits), ";",                                    
                                  " Ratio to degrees of freedom (df): ", round(x$ratioG, digits), ";",
                                  " p-value: ", round(x$pvalG, 4), "\n", sep = "")
    if ("JSD" %in% print.stats) cat("Jensen-Shannon Distance: ", round(x$JSD, digits),"\n", sep = "")  
    if ("pct.correct" %in% print.stats) cat("Percent correctly predicted (pct.correct): ",  round(x$pct.correct, 1), "\n", sep = "")     

  } else if (x$m > 1) {
    cat("\nAverage results for ", x$m, " syntheses\n", sep = "")
    cat("\nNumber of cells in each table: ", 
        round(mean(x$df[1] + x$nempty[1] + celldiff), digits),
        "; Number of cells contributing to utility measures: ", 
        round(mean(x$df + celldiff), digits),"\n", sep = "")
    cat("\nAverage utility score results\n")
    if ("VW" %in% print.stats) cat("Voas Williamson (VW): ", round(mean(x$UtabVW), digits), ";",  ## round all pval to 4  
                                   " Ratio to degrees of freedom (df): ", round(mean(x$ratioVW), digits), ";",
                                   " Median p-value: ",round(median(x$pvalVW), 4), "\n", sep = "")    
    if ("FT" %in% print.stats) cat("Freeman Tukey (FT): ", round(mean(x$UtabFT), digits), ";",
                                   " Ratio to degrees of freedom (df): ", round(mean(x$ratioFT), digits), ";",
                                   " Median p-value: ", round(median(x$pvalFT), 4), "\n", sep = "")
    if ("G" %in% print.stats) cat("Likelihood ratio (G): ", round(mean(x$UtabG), digits), ";",                                    
                                  " Ratio to degrees of freedom (df): ", round(mean(x$ratioG), digits), ";",
                                  " Median p-value: ", round(median(x$pvalG), 4),"\n", sep = "")
    if ("JSD" %in% print.stats) cat("Jensen-Shannon Distance: ", round(mean(x$JSD), digits),"\n", sep = "")  
    if ("pct.correct" %in% print.stats) cat("Mean percent correctly predicted (pct.correct): ",  round(mean(x$pct.correct), 1), "\n", sep = "")  
    
    cat("\nResults from individual syntheses\n")
    tab.res <- cbind.data.frame(
      round(x$ratioVW, digits), round(x$pvalVW, 4),
      round(x$ratioFT, digits), round(x$pvalFT, 4),
      round(x$ratioG, digits),  round(x$pvalG, 4), 
      round(x$JSD, digits), round(x$pct.correct, 1)) 
    colnames(tab.res) <- c("VW Ratio", "VW p-value",
                           "FT Ratio", "FT p-value",
                           "G Ratio", "G p-value",
                           "JSD", "% correct")
    out <- NULL
    if (!("VW" %in% print.stats)) out <- 1:2               
    if (!("FT" %in% print.stats)) out <- c(out, 3:4)
    if (!("G" %in% print.stats)) out <- c(out, 5:6)
    if (!("JSD" %in% print.stats)) out <- c(out, 7)
    if (!("pct.correct" %in% print.stats)) out <- c(out, 8)
    if (!is.null(out)) tab.res <- tab.res[, -out]
    print(tab.res)
  }
  
  invisible(x)
}


###-----summary.out--------------------------------------------------------
summary.out <- function(z, digits = max(3L, getOption("digits") - 3L), ...)
{
    ncw <- function(x) {
        zz <- nchar(x, type = "w")
        if (any(na <- is.na(zz))) {
            zz[na] <- nchar(encodeString(zz[na]), "b")
        }
        zz
    }
    nv <- length(z)
    nm <- names(z)
    lw <- numeric(nv)
    nr <- if (nv)
        max(unlist(lapply(z, NROW)))
    else 0
    for (i in seq_len(nv)) {
        sms <- z[[i]]
        if (is.matrix(sms)) {
            cn <- paste(nm[i], gsub("^ +", "", colnames(sms),
                useBytes = TRUE), sep = ".")
            tmp <- format(sms)
            if (nrow(sms) < nr)
                tmp <- rbind(tmp, matrix("", nr - nrow(sms),
                  ncol(sms)))
            sms <- apply(tmp, 1L, function(x) paste(x, collapse = "  "))
            wid <- sapply(tmp[1L, ], nchar, type = "w")
            blanks <- paste(character(max(wid)), collapse = " ")
            wcn <- ncw(cn)
            pad0 <- floor((wid - wcn)/2)
            pad1 <- wid - wcn - pad0
            cn <- paste0(substring(blanks, 1L, pad0), cn, substring(blanks,
                1L, pad1))
            nm[i] <- paste(cn, collapse = "  ")
            z[[i]] <- sms
        }
        else {
            sms <- format(sms, digits = digits)
            lbs <- format(names(sms))
            sms <- paste0(lbs, ":", sms, "  ")
            lw[i] <- ncw(lbs[1L])
            length(sms) <- nr
            z[[i]] <- sms
        }
    }
    if (nv) {
        z <- unlist(z, use.names = TRUE)
        dim(z) <- c(nr, nv)
        if (anyNA(lw))
            warning("probably wrong encoding in names(.) of column ",
                paste(which(is.na(lw)), collapse = ", "))
        blanks <- paste(character(max(lw, na.rm = TRUE) + 2L),
            collapse = " ")
        pad <- floor(lw - ncw(nm)/2)
        nm <- paste0(substring(blanks, 1, pad), nm)
        dimnames(z) <- list(rep.int("", nr), nm)
    }
    else {
        z <- character()
        dim(z) <- c(nr, nv)
    }
    attr(z, "class") <- c("table")
    z
}
