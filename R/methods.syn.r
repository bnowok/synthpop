###-----print.synds--------------------------------------------------------

print.synds <- function(x, ...){
  cat("Call:\n($call) ")
  print(x$call)
  cat("\nNumber of synthesised data sets: \n($m) ",x$m,"\n")  
  if (x$m==1){
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
  maxsum = 7, digits = max(3, getOption("digits")-3), ...){
  if (!is.null(msel) & !all(msel %in% (1:object$m))) stop("Invalid synthesis number(s)", call. = FALSE)

  sy <- list(m = object$m, msel = msel, method = object$method)

  if (object$m == 1){
    sy$result <- summary(object$syn,...)
  } else if (is.null(msel)){
    zall <- vector("list",object$m) 
    for (i in 1:object$m) zall[[i]] <- lapply(object$syn[[i]], summary,
      maxsum = maxsum, digits = digits, ...)
    zall.df <- Reduce(function(x,y) mapply("rbind",x,y),zall)
    meanres <- lapply(zall.df, function(x) apply(x,2,mean))
    sy$result <- summary.out(meanres)
  } else if (length(msel)==1){
    sy$result <- summary(object$syn[[msel]],...)
  } else {
    for (i in (1:length(msel))){
      sy$result[[i]] <- summary(object$syn[[msel[i]]],...)
    }
  }
  class(sy) <- "summary.synds"
  return(sy)
}


###-----print.summary.synds------------------------------------------------

print.summary.synds <- function(x, ...){

 if (x$m==1){
   cat("Synthetic object with one synethesis using methods:\n")
   print(x$method)
   cat("\n")
   print(x$result)
 } else if (is.null(x$msel)){
   cat("Synthetic object with ",x$m," syntheses using methods:\n",sep="")
   print(x$method)
   cat("\nSummary (average) for all synthetic data sets:\n",sep="")
   print(x$result)  
 } else if (length(x$msel)==1){
   cat("Synthetic object with ",x$m," syntheses using methods:\n",sep="")
   print(x$method)
   cat("\nSummary for synthetic data set ",x$msel,":\n",sep="")
   print(x$result)
 } else {
   cat("Synthetic object with ",x$m," syntheses using methods:\n",sep="")
   print(x$method)
   for (i in (1:length(x$msel))){
     cat("\nSummary for synthetic data set ",x$msel[i],":\n",sep="")
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
    for (i in 1:m){
      pos <- match(namesbyfit[[i]],allnames)
      matcoef[i,pos] <- analyses[[i]]$coefficients[,1]
      matvar [i,pos] <- analyses[[i]]$coefficients[,2]^2
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
 if (!class(data)=="synds") stop("Data must have class synds\n")
 call <- match.call()
 fitting.function <- "lm"
 analyses <- as.list(1:data$m)

 # do the repated analysis, store the result without data
 if (data$m==1) {
   analyses[[1]] <- summary(lm(formula,data=data$syn,...))
 } else {
   for(i in 1:data$m) {
     analyses[[i]] <- summary(lm(formula,data=data$syn[[i]],...))
   }
 }
 allcoefvar <- mcoefvar(analyses = analyses)
      
 # return the complete data analyses as a list of length m
 object <- list(call=call, mcoefavg = allcoefvar$mcoefavg, 
             mvaravg = allcoefvar$mvaravg, proper = data$proper, m = data$m, 
             analyses = analyses, fitting.function = fitting.function,
             n = data$n, k=data$k, mcoef = allcoefvar$mcoef,
             mvar = allcoefvar$mvar)
 class(object) <- "fit.synds"
 return(object)
}


###-----glm.synds----------------------------------------------------------

glm.synds <- function(formula, family="binomial", data, ...)
{
 if (!class(data)=="synds") stop("Data must have class synds\n")
 call <- match.call()
 fitting.function <- "glm"
 analyses <- as.list(1:data$m)
 
 # do the repated analysis, store the result without data
 if (data$m==1) {
   analyses[[1]] <- summary(glm(formula,data=data$syn,family=family,...))
 } else {
   for(i in 1:data$m) {
     analyses[[i]] <- summary(glm(formula,data=data$syn[[i]],family=family,...))
   }
 }
 allcoefvar <- mcoefvar(analyses = analyses)
 
 # return the complete data analyses as a list of length m
 object <- list(call=call, mcoefavg = allcoefvar$mcoefavg, 
             mvaravg = allcoefvar$mvaravg, proper = data$proper, m = data$m,
             analyses = analyses, fitting.function = fitting.function,
             n = data$n, k = data$k, mcoef = allcoefvar$mcoef,
             mvar = allcoefvar$mvar)
 class(object) <- "fit.synds"
 return(object)
}


###-----print.fit.synds----------------------------------------------------

print.fit.synds <- function(x, msel = NULL, ...)
{
  if (!is.null(msel) & !all(msel %in% (1:x$m))) stop("Invalid synthesis number(s)", call. = FALSE)
  cat("\nCall:\n")
  print(x$call)
  if (is.null(msel)){
    cat("\nCombined coefficient estimates:\n")
    print(x$mcoefavg)
  } else {
    cat("\nCoefficient estimates for selected syntheses:\n")
    print(x$mcoef[msel,,drop=FALSE])
  }
  invisible(x)
}


###-----summary.fit.synds--------------------------------------------------

summary.fit.synds <- function(object, population.inference = FALSE, msel = NULL, partly = FALSE, ...)
{ # df.residual changed to df[2] because didn't work for lm - check if that's ok
  if (!class(object) == "fit.synds") stop("Object must have class fit.synds\n")
  m <- object$m
  n <- sum(object$n)                                                     #!BN-15/08/2016, strata
  if (is.list(object$k)) k <- sum(object$k[[1]]) else k <- sum(object$k) #!BN-15/08/2016, strata 
    
  coefficients <- object$mcoefavg     # mean of coefficients (over m syntheses)
  vars <- object$mvaravg              # mean of variances (over m syntheses)

  if (population.inference == F){ ## inf to Q hat

    if(partly == TRUE){
      bm = c(0)
      for(i in 1:m){
        bm = bm + ((object$mcoef[i,] - object$mcoefavg)^2 / (m - 1))
      }
      result <- cbind(coefficients,
                      sqrt(bm/m),
                      sqrt(vars*k/n),
                      coefficients/sqrt(vars*k/n),
                      #sqrt((1 + coefficients^2/vars/2/object$analyses[[1]]$df[2])*n/k/m)
                      sqrt(1 + coefficients^2/vars/4/object$analyses[[1]]$df[2] * n/k/m))  #!?????
      dimnames(result)[[2]] <- c("B.syn","se(B.syn)","se(Beta).syn","Z.syn","se(Z.syn)")
    }
    else if (object$proper == F){
      ## simple synthesis
      result <- cbind(coefficients,
                      sqrt(vars/m),
                      sqrt(vars*k/n),
                      coefficients/sqrt(vars*k/n),
                      #sqrt((1 + coefficients^2/vars/2/object$analyses[[1]]$df[2])*n/k/m)
                      sqrt(1 + coefficients^2/vars/4/object$analyses[[1]]$df[2] * n/k/m))  #!?????
      dimnames(result)[[2]] <- c("B.syn","se(B.syn)","se(Beta).syn","Z.syn","se(Z.syn)")
    } else {
      ## proper synthesis
      result <- cbind(coefficients,
                      sqrt(vars*(1+k/n)/m), 
                      sqrt(vars*k/n),
                      coefficients/sqrt(vars*k/n),
                      #sqrt((1 + k/n + coefficients^2/vars/2/object$analyses[[1]]$df[2])*n/k/m)
                      sqrt(1 + coefficients^2/vars/4/object$analyses[[1]]$df[2]*n/k/m))    #!?????
      dimnames(result)[[2]] <- c("B.syn","se(B.syn)","se(Beta).syn","Z.syn","se(Z.syn)")
    }

  } else { ## pop inference to Q

  	if(partly == TRUE){
      bm = c(0)
      for(i in 1:m){
        bm = bm + ((object$mcoef[i,] - object$mcoefavg)^2 / (m - 1))
      }
      result <- cbind(coefficients,
                      sqrt(bm/m + vars),
                      coefficients/sqrt(bm/m + vars))
                      #sqrt((1 + coefficients^2/vars/2/object$analyses[[1]]$df[2])*n/k/m)
      dimnames(result)[[2]] <- c("B.syn","se(B.syn)","Z.syn")
    }
    else if (object$proper == F){
      ## simple synthesis
      Tf <- vars*(k/n+1/m)
      result <- cbind(coefficients,sqrt(Tf),coefficients/sqrt(Tf))
      dimnames(result)[[2]] <- c("B.syn","se(B.syn)","Z.syn")
    }
    else {
      ## proper synthesis
      Tf <- vars*(k/n+(1+k/n)/m)
      result <- cbind(coefficients,sqrt(Tf),coefficients/sqrt(Tf))
      dimnames(result)[[2]] <- c("B.syn","se(B.syn)","Z.syn")
    }
  }
  res <- list(call = object$call, proper = object$proper,
              population.inference = population.inference,
              fitting.function = object$fitting.function,
              m = m, coefficients = result, n = n, k = k, 
              analyses = object$analyses, msel = msel)
  class(res) <- "summary.fit.synds"
  return(res)
}


###-----print.summary.fit.synds--------------------------------------------

print.summary.fit.synds <- function(x, ...) {
  if (!is.null(x$msel) & !all(x$msel %in% (1:x$m))) stop("Invalid synthesis number(s)", call. = FALSE)
  
  if (is.null(x$msel)){
    if (x$m==1) {
      cat("\nFit to synthetic data set with a single synthesis.\n")
    } else {
      cat("\nFit to synthetic data set with ",x$m," syntheses.\n",sep="")
    }
    if (x$population.inference) {
      cat("Inference to population coefficients.\n")
    } else {
      cat("Inference to coefficients and standard errors\nthat would be obtained from the observed data.\n")
    }
    cat("\nCall:\n")
    print(x$call)
    cat("\nCombined estimates:\n")
    if (x$population.inference){
      print(x$coefficients[,c("B.syn","se(B.syn)","Z.syn")])
    } else {
      print(x$coefficients[,c("B.syn","se(Beta).syn","Z.syn")])
    }
  } else {
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficient estimates for selected syntheses:\n")
    for(i in x$msel) {          
      cat("\nsyn=",i,"\n",sep="")
      print(x$analyses[[i]]$coefficients)
    }
  }      
  invisible(x)
}


###-----print.compare.fit.synds--------------------------------------------

print.compare.fit.synds <- function(x, print.coef = x$print.coef, ...){

  cat("\nCall used to fit models to the data:\n")
  print(x$call)
  if (print.coef == TRUE){
    cat("\nEstimates for the observed data set:\n")
    print(x$coef.obs)
    cat("\nCombined estimates for the synthetised data set(s):\n")
    print(x$coef.syn[,c("B.syn","se(Beta).syn","se(B.syn)","Z.syn")])
  }  
    
  cat("\nDifferences between results based on synthetic and observed data:\n")
    print(cbind.data.frame(x$coef.diff,x$ci.overlap))                                 
  if (x$m==1) {
    cat("\nMeasures for one synthesis and ", x$ncoef, " coefficients", sep="") 
  } else {
    cat("\nMeasures for ", x$m, " syntheses and ", x$ncoef, " coefficients", sep="") 
  }   
  cat("\nMean confidence interval overlap: ", x$mean.ci.overlap)
  cat("\nMean absolute std. coef diff: ", x$mean.abs.std.diff)
  cat("\nMean lack-of-fit: ", x$mean.lof )
  cat("\nMean lack-of-fit ratio to expected: ", x$mean.lof.exp)
  cat("\nStandardised lack-of-fit: ", x$std.lof ,"\n")

  if (!is.null(x$ci.plot)){
    cat("\nConfidence interval plot:\n")
    print(x$ci.plot)
  }
  invisible(x)
}


###-----print.compare.synds------------------------------------------------

print.compare.synds <- function(x, ...) {
  cat("\nComparing percentages observed with synthetic\n\n")
  if (class(x$plots)[1]=="gg"){
    print(x$tables) 
    print(x$plots)
  } else {
    for (i in 1:length(x$tables)) {
      print(x$tables[[i]]) 
      print(x$plots[[i]])
      if (i < length(x$tables)) {
        cat("Press return for next plot: ")
        ans <- readline()
      }
    }
  }
 invisible(x)
}


###-----print.utility.gen------------------------------------------------

print.utility.gen <- function(x,print.zscores =x$print.zscores,digits= x$digits,
  usethresh = x$usethresh,zthresh=x$zthresh, print.variable.importance=x$print.variable.importance, ...){
  #
	cat("\nUtility score calculated by method: ",x$method,"\n")
  cat("\nCall was: ","\n")
  print(x$call)
	  if (x$method=="cart") cat("\nNull utility sumulated from a permutation test with ",x$nperm,"replications\n")
  if (x$m>1) {
	mnU<-mean(x$utilVal); mnE<-mean(x$utilExp); mnR<-mean(x$utilR); mnS<-mean(x$utilStd)
      cat("\nMean utility score from ",x$m," syntheses \nUtility ",round(mnU,x$digits)," Expected value",round(mnE,x$digits),
      " Ratio to expected ",round(mnR,x$digits)," Standardised ",round(mnS,x$digits),"\n")
      if (x$m<11) {
        cat("\nIndividual utility score results from ",x$m,"syntheses\n\n")
        tabres<-rbind(round(x$utilVal,x$digits),round(x$utilExp,x$digits),round(x$utilR,x$digits), round(x$utilStd,x$digits))
        dimnames(tabres)<-list(c("Utility","Expected","Ratio","Standardised"),paste("syn",1:length(x$utilVal)))
        print(tabres)
      }
  }
  else {
     cat("\nUtility score results\nUtility ",round(x$utilVal,x$digits)," Expected value",round(x$utilExp,x$digits),
      " Ratio to expected ",round(x$utilR,x$digits)," Standardised ",round(x$utilStd,x$digits),"\n")}
  if (x$print.zscores==TRUE) {
    if(x$method=="cart") cat("\nZscores not available for CART models\n")
    else {
      if (x$m==1) {
        zscores=summary(x$fit)$coefficients[,3]
        if (x$usethresh==TRUE) { 
          cat("\nShowing Zscores greater than the threshold of +/- ",x$zthresh,"\n")
          zscores<-zscores[abs(zscores)>x$zthresh]
          if (length(zscores)==0) cat("No z scores above threshold of ",x$zthresh,"\n")
          else print(zscores)
        }
        else print(zscores)
      }
      else{
        if (x$usethresh==TRUE) { 
          cat("\nShowing Zscores greater than the threshold of +/- ",x$zthresh,"\n")
          for (i in 1:x$m){
            zscores=summary(x$fit[[i]])$coefficients[,3]
            cat("Synthesis ",i,"\n")
            zscores<-zscores[abs(zscores)>x$zthresh]
            if (length(zscores)==0) cat("No z scores above threshold of ",x$zthresh,"\n")
            else print(zscores)
          }
        }
        else  { 
          for (i in 1:x$m){
            zscores=summary(x$fit[[i]])$coefficients[,3]
            cat("Synthesis ",i,"\n")
            print(zscores)
          }
        }
      }
    }
  }
  
 if (x$print.variable.importance==TRUE) {
  if(x$method!="cart" & x$tree.method=="rpart") cat("Variable importance only available for CART models using function rpart\n")
  else {
  cat("\nRelative importance of each variable scaled to add to 100\n" )
    if (x$m==1) {
      variable.importance<-x$fit$variable.importance
      variable.importance<-round(variable.importance/sum(variable.importance)*100,x$digits)
      print(variable.importance)
    }
    else{
       for (i in 1:x$m){
       cat("Synthesis ",i,"\n")
         variable.importance<-x$fit[[i]]$variable.importance
         variable.importance<-round(variable.importance/sum(variable.importance)*100,x$digits)
         print(variable.importance)
    }
  }
 }
}
	invisible(x)
	
}


###-----print.utility.tab--------------------------------------------------

print.utility.tab <- function(x, print.tables = x$print.tables,  
  print.zdiff = x$print.zdiff, digits = x$digits,...){

  if(print.tables == TRUE) {
    if (is.table(x$tab.obs)) {
      if (sum(x$tab.obs)!= x$n) {
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

    if (x$m==1) {
      cat("\nSynthesised: \n($tab.syn)\n")
	    print(x$tab.syn) 
    } else {
      meantab <- apply(simplify2array(x$tab.syn), c(1,2), mean)
      cat("\nMean of ",x$m," synthetic tables ($tab.syn):\n", sep="")
      print(round(meantab, digits))
    }
  }
  
  if(print.zdiff == TRUE) {
    cat("\nTable of Z scores for differences: \n($tab.zdiff)\n")
    if (x$m==1) {
      print(round(x$tab.zdiff, digits)) 
    } else {
      meanzdiff <- apply(simplify2array(x$tab.zdiff), c(1,2), mean)
      cat("\nMean of ",x$m," Z score tables:\n", sep="")
      print(round(as.table(meanzdiff), digits))
    }
  }
  
  if (x$m==1){
    cat("\nNumber of cells in each table: ", 
        x$df[1] + x$nempty[1] + 1,
        "; Number of cells contributing to utility measures: ", 
        x$df + 1,"\n", sep="")
    cat("\nUtility score results\n")
    cat("Freeman Tukey (FT): ", round(x$UtabFT,digits), ";",
        " Ratio to degrees of freedom (df): ", round(x$ratioFT,digits), ";",
        " p-value: ", x$pvalFT, "\n", sep="")
    cat("Voas Williamson (VW): ", round(x$UtabVW,digits), ";",
        " Ratio to degrees of freedom (df): ", round(x$ratioVW,digits), ";",
        " p-value: ", x$pvalVW, "\n", sep="")
  } else if (x$m>1){
    cat("\nAverage results for ", x$m, " syntheses\n", sep="")
    cat("\nNumber of cells in each table: ", 
        round(mean(x$df[1] + x$nempty[1] + 1),digits),
        "; Number of cells contributing to utility measures: ", 
        round(mean(x$df + 1),digits),"\n", sep="")
    cat("\nUtility score results\n")
    cat("Freeman Tukey (FT): ", round(mean(x$UtabFT),digits), ";",
        " Ratio to degrees of freedom (df): ", round(mean(x$ratioFT),digits),"\n", sep="")
    cat("Voas Williamson (VW): ", round(mean(x$UtabVW), digits), ";",
        " Ratio to degrees of freedom (df): ",  round(mean(x$ratioVW), digits),"\n", sep="")
    
    cat("\nResults from individual syntheses\n")
    tab.res <- cbind.data.frame(x$df,
    round(x$UtabFT,digits), round(x$pvalFT,digits),
    round(x$UtabVW,digits), round(x$pvalVW,digits))
    colnames(tab.res) <- c("df", 
                           "FT Utility","FT p-value",
                           "VW Utility","VW p-value")
    print(tab.res)
  }

 	invisible(x)
}


###-----summary.out--------------------------------------------------------
summary.out <- function (z, digits = max(3L, getOption("digits") - 3L), ...)
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

