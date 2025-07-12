synorig.compare <- function(syn,orig, print.flag = TRUE){

 needsfix <- FALSE
 unchanged <- TRUE
  ## to convert any tibbles or matrices
  if (!is.data.frame(syn)) { syn <- data.frame(syn) ; unchanged = FALSE}
  if (!is.data.frame(orig)) { orig <- data.frame(orig); unchanged = FALSE}


 
  if (any(!(names(syn)) %in% names(orig))) {
  cat("Variables",names(syn)[!(names(syn)%in% names(orig))],"in synthetic but not in original\n")
  out <- names(syn)[!(names(syn) %in% names(orig))]
  syn <- syn[,(names(syn) %in% names(orig))]
  cat(out,"dropped from syn\n\n")
  }
## reduce syn and orig to common vars in same order
common <- names(orig)[names(orig) %in% names(syn)]
len.common <- length(common)
if (print.flag) cat(len.common, "variables in common out of" , length(names(syn)), "in syn out of", length(names(orig)),"in orig\n")
## reorder to match up
orig <- orig[,names(orig) %in% common] 
syn <- syn[,names(syn) %in% common] 
syn <- syn[, names(orig)]

##-------------- change common  variables  that are numeric in syn and factors in orig to factors in syn---------------------------
 nch_syn <- 0 ; nch_orig <- 0
  for( i in 1:len.common ) {
    if ( is.numeric( syn[,i] ) & is.factor( orig[,i]) ) {
      syn[,i] <- factor( syn[,i] ) 
      nch_syn <- nch_syn +1
      unchanged = FALSE
    }
    if ( is.factor( syn[,i] ) & is.numeric( orig[,i]) ) {
      orig[,i] <- factor( orig[,i] ) 
      nch_orig <- nch_orig +1
      unchanged = FALSE
    }
  }
 if (!unchanged) cat("\nVariables changed from numeric to factor",nch_syn,"in syn",nch_orig,"in original\n\n")

##-------------- change common character variables to factors---------------------------
nch_syn <- 0 ; nch_orig <- 0; unchanged2 <- TRUE
for( i in 1:len.common ) {
  if ( is.character( syn[,i] ) ) {
    syn[,i] <- factor( syn[,i] ) 
    nch_syn <- nch_syn +1
    unchanged2 = FALSE ; unchanged = FALSE
  }
  if ( is.character( orig[,i] )) {
    orig[,i] <- factor( orig[,i] ) 
    nch_orig <- nch_orig +1
    unchanged2 = FALSE ; unchanged = FALSE
  }
}
if (!unchanged2) cat("\nVariables changed from character to factor x",nch_syn,"in syn","and",nch_orig,"in orig\n\n")

  ##--------------- check data types match in common variables----------------------------------

  for (i in 1:len.common){
    if ( is.integer(syn[,i])  &  (is.numeric(orig[,i]) & !is.integer(orig[,i])) ) { 
      syn[,i] <-  as.numeric(syn[,i])
      cat(names(syn)[i],"changed from integer to numeric in synthetic to match original\n")
      unchanged <- FALSE
    }
    else if ( is.integer(orig[,i])  & (is.numeric(syn[,i]) & !is.integer(syn[,i])) ) { 
      orig[,i] <-  as.numeric(orig[,i])
      cat(names(orig)[i],"changed from integer to numeric in original to match synthetic\n")
      unchanged <- FALSE
    }
    else if ( length(class(syn[,i])) !=  length(class(orig[,i])) ||
         !all(class(syn[,i])  ==  class(orig[,i]) )	) {
       cat("\nDifferent classes for",names(syn)[i],"in syn:",class(syn[,i]),"in orig:",class(orig[,i]),"\n")
      needsfix <- TRUE
    }

  } 
  ##--------------------------- compare missingness and levels for factors------------------------------

  for (i in 1:length(common)){
    #cat(names(syn)[i], names(orig)[i],"NAMES\n\n")
    if (!any(is.na(orig[,i])) & any(is.na(syn[,i])) ) { 
      cat("\n\nMissing data for common variable", names(syn)[i], "in syn but not in orig\nThis looks wrong check carefully\n")
      orig[,i] <- addNA(orig[,i])
      print(levels(syn[,i]))
      cat("NA added to factor ",names(orig)[i],"in orig\n\n\n")
      unchanged <- FALSE
    }
    if ( is.factor(syn[,i])  & is.factor(orig[,i]) ) {  
    if ( any(is.na(orig[,i])) & !any(is.na(syn[,i])) ) {
        syn[,i] <- addNA(syn[,i])
        cat("NA added to factor ",names(syn)[i],"in syn\n")
        unchanged <- FALSE
    }
      lev1 <- levels(syn[,i]) ;lev2 <- levels(orig[,i])
      
      if ( length(lev1) != length(lev2)  ||  !all(lev1 == lev2) ) {
        cat("\nFactor levels don't match for", names(syn)[i],"levels combined\n",
               "syn levels ", lev1,"\n","orig levels", lev2,"\n")

        syn[,i] <- factor(as.character(syn[,i]), exclude = NULL, levels = unique(c(lev1,lev2)))
        orig[,i] <- factor(as.character(orig[,i]), exclude = NULL, levels = unique(c(lev1,lev2)))
         unchanged = FALSE
      }
    }
  }

  if (needsfix)  cat("\n***********************************************************************************\n",
                      "STOP: you may need to change the original or synthetic data to make them match:\n")

           
  if (!unchanged) { cat("\n*****************************************************************\n",
                    "Differences detected and corrections attempted check output above.\n")
  }
  else {cat("Synthetic and original data checked with synorig.compare,\n looks like no adjustment needed\n\n")}

     res <- list(syn=syn, orig = orig, needsfix = needsfix, unchanged = unchanged)
}
