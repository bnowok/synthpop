# Source functions for synthpop library to create synthetic data
# for the SYLLS project.

# Structure and some functions based on code from MICE package
# by S. van Buuren and K. Groothuis-Oudshoorn

syn <- function(data, method = "cart",
                visit.sequence = (1:ncol(data)),
                predictor.matrix = NULL,  
                m = 1, k = nrow(data), proper = FALSE, 
                minnumlevels = -1, maxfaclevels = 60,
                rules = NULL, rvalues = NULL,                                  
                cont.na = NULL, semicont = NULL,
                smoothing = NULL, event = NULL, denom = NULL,
                drop.not.used = FALSE, drop.pred.only = FALSE,                                            
                default.method = c("normrank", "logreg", "polyreg", "polr"),
                numtocat = NULL, catgroups = rep(5, length(numtocat)),
                models = FALSE,
                print.flag = TRUE,
                seed = "sample",
                ...)
{
#----------------------------------------------------------------------
# the code for the following checking functions is included within syn
# so as to allow them to access the local objects
#----------------------------------------------------------------------
obs.vars <- names(data)

#if (k == 0) m <- 0 #!BN-12/08/2016

# set default method for everything to cart and to blank (which will get defaults) if method is "parametric"
# if (all(method == "")) method = "cart"  # change to "ctree"?
# else if (length(method)==1 && method=="parametric") method=rep("",dim(data)[2])

if (!is.null(attr(data,"filetype")) && attr(data,"filetype") == "sav") {
  var.lab <- attr(data,"variable.labels")
  val.lab <- attr(data,"label.table")
} else {
  var.lab <- val.lab <- NULL
}

# if visit sequence includes variable names change them into column indecies 
if (is.character(visit.sequence)) {
  nametoind <- match(visit.sequence, colnames(data))
    if (any(is.na(nametoind))) stop("Unrecognized variable(s) in visit.sequence: ", 
      paste(visit.sequence[is.na(nametoind)], collapse = ", "), call. = FALSE)
    else visit.sequence <- nametoind
} else if (!all(visit.sequence %in% 1:length(data))) stop("Column indices in visit.sequence must be between 1 and ", 
                                                          length(data), sep = "", call. = FALSE)  


##-----------------------code for numtocat------------------------------------------
if (!is.null(numtocat)) {
  
  # if numtocat is numeric change to column names and check names if not
  if (is.numeric(numtocat)) {
    if (!all(numtocat %in% 1:length(data))) stop("Column indices in visit.sequence must be between 1 and ", 
      length(data), sep = "", call. = FALSE)  
    numtocat <- names(data)[numtocat]
  } else if (!all(numtocat %in% names(data))) stop("numtocat variable(s): ",
      paste(numtocat[!numtocat %in% names(data)], collapse = ", ")," are not in data.", sep = "", call. = FALSE)
 
 # check if numtocat in visit sequence and remove from numtocat if not
 if (!all(numtocat %in% names(data)[visit.sequence])) {
   cat("\nVariable(s) in numtocat (", paste(numtocat[!numtocat %in% names(data)[visit.sequence]], collapse = ", "),
     ") not in visit.sequence and have been removed from numtocat.\n" , sep = "")
   if (length(catgroups) > 1) catgroups <- catgroups[numtocat %in% names(data)[visit.sequence]]
   numtocat <- numtocat[numtocat %in% names(data)[visit.sequence]]
   if (length(numtocat) == 0) stop("\nAll variables in numtocat removed, perhaps try without this parameter.\n", call. = FALSE)
 }
 
 # get cont.na for numtocat vars
 cna <- cont.na
 if (!any(names(cna) %in% numtocat)) cna <- NULL
 if (!is.null(cna) & any(names(cna) %in% numtocat)) {
   cna <- cna[names(cna) %in%  numtocat]
   if (length(cna) == 0 ) cna <- NULL
 }

 # make numtocat object incorporating some checks
 numtocat.obj <- numtocat.syn(data, numtocat = numtocat, cont.na = cna, 
                              catgroups = catgroups, print.flag = FALSE)
 
 # adjust visit.sequence
 visit.sequence <- c(visit.sequence, length(data) + 1:length(numtocat))
  
 # replace data with categorised with orig vars at end
 data <- cbind(numtocat.obj$data, numtocat.obj$orig)
 
 if (length(method) == 1) {
   if (method == "parametric") stop("'parametric' not available, when numtocat is used, methods must be specified in detail.\n", call. = FALSE )
   else method <- rep(method, length(data) - length(numtocat))
 }
  
 # assign method for final columns as nested and give changed parameters correct names
 method <- c(method, paste("nested", numtocat, sep = "."))
#browser()
 # checks on methods for numtocat variables
 if (any(method[numtocat.obj$ind] %in% 
     c("norm", "normrank", "lognorm", "sqrtnorm", "cubertnorm"))) {
   nummth <- method[numtocat.obj$ind] %in% 
     c("norm", "normrank", "lognorm", "sqrtnorm", "cubertnorm")
   stop("Method(s) (", paste(method[numtocat.obj$ind][nummth], collapse = ", "),
        ") assigned to variable(s) in numcat (", 
        paste(numtocat[nummth], collapse = ", "),
        ")\nunsuitable for categorical data.", sep = "", call. = FALSE)
 }  

 # modify visit sequence and predictor matrix to synthesis numtocat 
 # original variables after the others
 # and adjust other parameters to match if not null
 if (!is.null(predictor.matrix)) {
   predictor.matrix <- cbind(predictor.matrix, matrix(0, dim(data)[2], length(numtocat)))
   predictor.matrix <- rbind(predictor.matrix, matrix(0, length(numtocat), dim(data)[2] + length(numtocat)))
   for (i in 1:length(numtocat)) { 
     predictor.matrix[dim(data)[2] + i, numtocat.obj$ind[i]] <- 1
   }
 }

 # move these parameters to numeric versions
 newnames <- function(x) { 
   # mini function to change names of lists
   xn <- names(x)
   ind <- match(numtocat, xn)
   ind <- ind[!is.na(ind)]
   names(x)[ind] <- paste("orig", names(x)[ind], sep = ".")
   return(x)
 }
 if (!is.null(rules))     rules     <- newnames(rules)
 if (!is.null(rvalues))   rvalues   <- newnames(rvalues)
 if (!is.null(cont.na))   cont.na   <- newnames(cont.na)
 if (!is.null(smoothing)) smoothing <- newnames(smoothing) 

 if (print.flag == TRUE) cat(
"**************************************************************
 The numeric variable(s): ", 
 paste(names(data)[numtocat.obj$ind], collapse = ", "), 
"\n will been synthesised as grouped variables and their numeric
 values generated from boostrap samples within categories.
**************************************************************\n", sep = "") 
}
##-----------------end of--code for numtocat----------------------------


##-----------------------check.visit.sequence.syn-----------------------

check.visit.sequence.syn <- function(setup){

 vis      <- setup$visit.sequence
 nvar     <- setup$nvar
 varnames <- setup$varnames
 method   <- setup$method

 # visit.sequence can include column indices only
 # not all columns have to be given - variables
 # not in the visit sequence won't be synthesised
 
 # stop if variable in visit.sequence more than once
 if (any(duplicated(vis))) stop("Visit sequence includes repeated variable names/numbers.\n", call. = FALSE)

 # remove any visitSequnce members outside 1:nvar
 if (any(!(vis %in% 1:nvar))) {
   cat("Element(s): {",paste(vis[!(vis %in% 1:nvar)],
       collapse = ", "),"} of the 'visit.sequence' removed as not valid. No such column.\n\n", sep = "")
   vis <- as.numeric(vis[vis %in% 1:nvar])
 }

 # remove event indicator(s) from visit.sequence, if present
 event.in.vis <- !is.na(match(vis,event))
 if (!is.null(event) & any(event.in.vis)) {
   cat("Element(s) {",paste(vis[event.in.vis],collapse = ", "),
       "} of the 'visit.sequence' with method(s) {",
       paste(method[vis[event.in.vis]],collapse = ", "),
       "} removed because used as event indicator(s).\nAny event indicators will be synthesised along with the corresponding survival time(s). \n\n")
   vis <- vis[!event.in.vis]
   if (length(vis) < 2) stop("Visit sequence now of length ",
       length(vis),". No synthesis done.")
 } 
                                                            #GRdenom new code
 #!BN adjusted to allow visit sequence with selected vars only 
 #! have to add a condition when denominator is not in visit seq at all;
 #! sampler has to be changed still
 #!---                                                                   
 #  check that denominator comes before the count for a binomial with denom
 #if (any(denom>0)) {
 #   denom.in.vis<-(1:nvar)[denom>0]
 #       for (j in denom.in.vis){
 #          posj<-(1:length(vis))[match(j,vis)]
 #          kj <-denom[j]
 #          posk<-(1:length(vis))[match(kj,vis)]
 #      if (posj<=posk) 
 #         stop("\n Denominator ",varnames[j]," for ",varnames[kj]," must be synthesisied before it\n")
 #   }
 #
 # }                                                               
 
 # check that denominator comes before the count for a binomial with denom 
 for (j in which(denom[vis] > 0)) { 
   denom.pos <- match(denom[vis][j],vis)
   if (j < denom.pos) stop("Denominator ",varnames[denom[vis][j]]," for ",
                           varnames[vis[j]]," must be synthesisied before it\n",
                           call. = FALSE)
 }
 #!---
                                                                        
 # stop if visit.sequence is empty
 if (length(vis) == 0) stop(paste("Seems no variables being synthesised.\nCheck parameter 'visit.sequence'."))

 # All variables used in passive synthesis have to be synthesised BEFORE
 # the variables they apply to
 for (j in which(is.passive(method[vis]))) {  #  
   var.present <- match(all.vars(as.formula(method[vis][j])),varnames) 
   var.in.vis  <- match(var.present,vis)
   if (j < max(var.in.vis) | any(is.na(var.in.vis))) stop("Variable(s) used in passive synthesis for ",
     varnames[vis][j]," has/have to be synthesised BEFORE the variables they apply to.")
 }

 setup$visit.sequence <- vis
 return(setup)
}
##-----------------end of--check.visit.sequence.syn---------------------


##-----------------------check.predictor.matrix.syn---------------------

check.predictor.matrix.syn <- function(setup){
 ## checks the predictor.matrix
 ## makes consistency edits of the predictor.matrix

 pred     <- setup$predictor.matrix
 nvar     <- setup$nvar
 varnames <- setup$varnames
 vis      <- setup$visit.sequence
 method   <- setup$method
 denom    <- setup$denom                     #GRdenom new
    
 # set up default predictor matrix (if not provided by a user)
 # to lower diagonal in order of visitSequnce but with
 # elements for variables not to be synthesised set to 0
 
 pred.dt           <- matrix(0, nvar, nvar)
 pred.dt[vis, vis] <- outer(1:length(vis), 1:length(vis), ">")
 if (is.null(pred)) pred <- pred.dt

 # basic corrections for a default matrix or the one provided by a user
 dimnames(pred)   <- list(varnames, varnames)
 diag(pred)       <- 0

 # select from visit.sequence variables that are synthesised
 # (=method different than "")
 vis.syn <- vis
 if (!all(method == "") & length(method) > 1) vis.syn <- intersect(vis, which(method != ""))
 # removing predictors for variables with "" method
 if (length(vis.syn) < length(vis)) {
   vis.blank        <- setdiff(vis,vis.syn)
   pred[vis.blank,] <- 0
 }
 # removing predictors for variables not in visit.sequence
 pred[setdiff(1:nvar, vis),] <- 0
 
 # removing predictors for variables with "sample" method
 for (j in which(method == "sample")) pred[j,] <- 0
 
 # removing survival time from predictors
 for (j in which(method == "survctree")) pred[,j] <- 0

 # removing event indicator from predictors
 for (j in which(method == "survctree" & event > 0)) pred[,event[j]] <- 0
                                                                     #GRdenom new lines
 #  remove denom from prediction of its numerator
    for  (j in which(method == "logreg")) {
      if (denom[j] > 0) pred[j, denom[j]] <- 0
    }
                                                                    # to here
 # checking consistency between visit.sequence and predictor matrix
 # provided by a user: dropping from predictors variables that are
 # not already synthesised
 preddel <- which((pred[, vis.syn, drop = FALSE] == 1 & 
                   pred.dt[, vis.syn, drop = FALSE] == 0), arr.ind = TRUE)
 if (length(vis) > 1) {
 	 pred[,vis.syn] <- ifelse((pred[,vis.syn] == 1 & pred.dt[, vis.syn] == 0),
                             0, pred[, vis.syn])
 	 if (nrow(preddel) > 0) cat(paste("Not synthesised predictor ",
                         varnames[vis.syn][preddel[, 2]],
                         " removed from predictor.matrix for variable ",
                         varnames[preddel[, 1]], ".\n", sep = ""))
 }
 setup$predictor.matrix <- pred
 setup$visit.sequence   <- vis
 return(setup)
}
##-----------------end of--check.predictor.matrix.syn-------------------


##-----------------------check.method.syn------------------------------

check.method.syn <- function(setup, data, proper) {
 ## check method, set defaults if appropriate

 method         <- setup$method
 default.method <- setup$default.method
 vis            <- setup$visit.sequence
 nvar           <- setup$nvar
 varnames	      <- setup$varnames
 pred		        <- setup$predictor.matrix
 event          <- setup$event
 denom          <- setup$denom                    

 # check that all ipf and allcat are at start of visit sequence   
 mcatall <- (method %in% "catall")[vis]
 mipf    <- (method %in% "ipf")[vis]
 if (any(mipf) & any(mcatall)) stop("Methods 'ipf' and 'catall' cannot both be used.\nIf you want all margins fitted for a set of variables,\nthen you could use 'ipf' and specify othmargins appropriately.\n", call. = FALSE)
 
 if (any(mcatall)) {
   if (any(mcatall != mcatall[order(!mcatall)])) stop("All variables with method 'catall' must be together at start of visit sequence.\n", call. = FALSE)
   if (sum(mcatall) == 1) {
     method[1] <- "sample"
     cat("First method changed to 'sample' from 'catall' as set for a single variable only.\n", call. = FALSE)
   }
 }
 if (any(mipf)) {
   if (any(mipf != mipf[order(!mipf)])) stop("All variables with method 'ipf' must be together at start of visit sequence.\n", call. = FALSE)
   if (sum(mipf) == 1) {
     method[1] <- "sample"
     cat("First method changed to 'sample' from 'ipf' as set for a single variable only.\n", call. = FALSE)
   }
 }
 
 # change method for constant variables but leave passive variables untouched
 # factors and character variables with missing data won't count,
 # as NA is made into an additional level
 for (j in 1:nvar) {
   if (!is.passive(method[j]) & method[j] != "ipf" & method[j] != "catall") {
     if (is.numeric(data[,j])) {
       v <- var(data[,j], na.rm = TRUE)
       if (!is.na(v)) constant <- (v < 1000 * .Machine$double.eps) else
       constant <- is.na(v) | v < 1000 * .Machine$double.eps
     } else {
       constant <- all(duplicated(data[,j])[-1L])
     }

     if (constant) {
       if (any(vis == j)) {
	       method[j] <- "constant" 
         if (print.flag == T) cat('Variable ', varnames[j], 
           ' has only one value so its method has been changed to "constant".\n', sep = "")
	       pred[j, ] <- 0
       }
       if (any(pred[, j] != 0)) { 
         if (print.flag == T) cat("Variable ", varnames[j], 
           " removed as predictor because only one value.\n", sep = "")
         pred[, j] <- 0
       }
     }
   }
 } 
 
 # check nested variables
 #---
 nestmth.idx <- grep("nested", method)
 gr.vars <- vector("character", length(method))
 gr.vars[nestmth.idx] <- substring(method[nestmth.idx], 8)

 if (length(nestmth.idx) > 0) { 
   for (i in nestmth.idx) {
     # check if provided grouped var exists
     if (gr.vars[i] == "") stop("No variable name provided for 'nested' method for ", 
       varnames[i] ,".\nSet method as 'nested.varname' instead of 'nested'.\n", call. = FALSE)
     if (!(gr.vars[i] %in% varnames)) stop("Unrecognized variable ", gr.vars[i], 
       " provided for 'nested' method for ", varnames[i] ,"\n", call. = FALSE)
     if (gr.vars[i] == varnames[i]) stop("Variable ", varnames[i], 
       " can not be predicted by itself.\n", call. = FALSE) 
       
     # check if var nested in gr.var
     #? tabvars   <- table(data[,i], data[,gr.vars[i]]) 
     tabvars <- table(data[,i], data[,gr.vars[i]], useNA = "ifany") 
     tabvars01 <- ifelse(tabvars > 0, 1, 0)
     ymulti <- rowSums(tabvars01) > 1
     if ("NAtemp" %in% names(ymulti)) ymulti["NAtemp"] <- FALSE  
     ymulti[names(ymulti) %in% cont.na[[i]]] <- FALSE  # missing values and cont.na are excluded 
     if (any(ymulti)) cat("\nNOTE: Variable ", varnames[i], 
       " is not nested within its predictor ", gr.vars[i], ".\nCheck values of ", 
       varnames[i], ": ", paste0(rownames(tabvars01)[ymulti], collapse = ", "), 
       "\n\n", sep = "")
   
   # adjust predictor matrix
   pred[i, -match(gr.vars[i], varnames)] <- 0  # remove all predictors except the group var
   pred[-match(varnames[i], gr.vars), i] <- 0  # exclude detailed var from predictors except when used for nested method
   }
   if (m > 0) method[nestmth.idx] <- "nested"
 }
 #---

 # check if var has predictors
 if (sum(pred) > 0) has.pred <- apply(pred != 0, 1, any)   # GR condition added
 else has.pred <- rep(0, nvar) 
 
 if (any(method == "parametric")) { # the default argument
   # set method for first in visit.sequence to "sample"
   # method[vis[1]] <- "sample"  
   # change to default methods for variables with predictors

   if (length(vis) > 1) {
	   for (j in vis[-1]) {
	     if (has.pred[j]) {
	       y <- data[,j]
	       if (is.numeric(y))        method[j] <- default.method[1]
	       else if (nlevels(y) == 2) method[j] <- default.method[2]
	       else if (is.ordered(y) & nlevels(y) > 2) method[j] <- default.method[4]
	       else if (nlevels(y) > 2)  method[j] <- default.method[3]
	       else if (is.logical(y))   method[j] <- default.method[2]
	       else if (nlevels(y) != 1) stop("Variable ",j," ",varnames[j],
		     " type not numeric or factor.") # to prevent a constant values failing
	     } else if (method[j] != "constant") method[j] <- "sample" 
	   }
	 }
 }
  
 # check whether the elementary syhthesising methods are available
 # on the search path
 active    <- !is.passive(method) & !(method == "") & !(method == "constant")   #BN-20/10/2016 "constant" added
 if (sum(active) > 0) {
   # fullNames <- paste("syn", method[active], sep=".") #!GR-29/8/16
   fullNames <- method[active]                                #!GR-29/8/16
   if (m == 0) fullNames[grep("nested",fullNames)] <- "nested"  #!GR-29/8/16
   fullNames <- paste("syn", fullNames, sep = ".")               #!GR-29/8/16
   notFound  <- !(fullNames %in% c('syn.bag', 'syn.cart', 'syn.cartbboot', 'syn.collinear', 
      'syn.ctree', 'syn.cubertnorm', 'syn.lognorm', 'syn.logreg', 'syn.nested', 'syn.norm', 
      'syn.normrank', 'syn.pmm', 'syn.polr', 'syn.polyreg', 'syn.ranknorm', 'syn.rf', 
      'syn.sample', 'syn.satcat', 'syn.smooth', 'syn.sqrtnorm', 'syn.survctree') | 
      sapply(fullNames, exists, mode = "function", inherit = TRUE)) 
   if (any(notFound)) stop(paste("The following functions were not found:",
                           paste(unique(fullNames[notFound]), collapse = ", ")), call. = FALSE)
 }

 # type checks on built-in  methods 

 for (j in vis) {
 	 y     <- data[,j]
   vname <- colnames(data)[j]
   mj    <- method[j]
   mlist <- list(m1 = c("logreg","polyreg","polr","ipf","catall"), #GRBN
                 m2 = c("norm","normrank","survctree"),
                 m3 = c("norm","normrank","survctree","logreg"))
   # In case of type mismatch stop execution

#                                                       #GRdenom lines changed
# check for logistic with denominator
# 
   if (denom[j] > 0) {
     if (!(mj %in% c("logreg"))) {
       method[j] <- "logreg"
       cat("Variable ", vname," has denominator (", colnames(data[denom[j]]), 
       ") and method ", mj, " has been changed to logreg\n", sep = "")
     }
   #if (!(mj %in% c("logreg"))) stop("Variable ", vname," has denominator (",
   #  colnames(data[denom[j]]), ") and method should be set to logreg and not ",mj,"\n", 
   #  call. = FALSE)
#
#  check all integers
#
  #browser()
  if (!((is.integer(y) | all((y - round(y)) == 0, na.rm = TRUE)) &  #!!!!!! address missing data issue
     (is.integer(data[denom[j]]) | all((data[denom[j]] - round(data[denom[j]]) == 0), na.rm = TRUE)))) #!!!!!! address missing data issue   
     stop("Variable ", vname," and denominator ", colnames(data[denom[j]]),
     " must be integers\n", call. = FALSE)
   if (any((data[denom[j]] - y) < 0, na.rm = TRUE)) stop("Variable ", vname,    #!!!!!! address missing data issue
     " must be less than or equal denominator ",
     colnames(data[denom[j]]),"\n", call. = FALSE)
   } else {
     if (is.numeric(y) & (mj %in% mlist$m1) & !(j %in% numtocat)) {                                             #!GRipf    numtocat added
       stop('Type mismatch for variable ', vname,
            '.\nSynthesis method "', mj, 
            '" is for categorical data unless grouped with numtocat.',
            sep = "", call. = FALSE)
     } else if (is.factor(y) & nlevels(y) == 2 & (mj %in% mlist$m2)) {
       stop('Type mismatch for variable ', vname,
       	    '.\nSyhthesis method "', mj, '" is not for factors.',
            sep = "", call. = FALSE)
     } else if (is.factor(y) & nlevels(y) > 2 & (mj %in% mlist$m3)) {
       stop('Type mismatch for variable ', vname,
    	      '.\nSyhthesis method "', mj,
            '" is not for factors with three or more levels.',
            sep = "", call. = FALSE)
     }                                               
   }
 }

 # check method for variables without predictors
 # set to "sample" if method is not valid
 
 # check if var has predictors (re-compute it)
 if (sum(pred) > 0) has.pred <- apply(pred != 0, 1, any)  #  GR condition added
 else has.pred <- rep(0, sqrt(length(pred)))            # this needed in case pred now has dimension 1

 for (j in vis) {
   if (!has.pred[j] & substr(method[j],1,6) != "nested" & is.na(any(match(method[j],
      c("","constant","sample","sample.proper","catall","ipf", 
      "norm","norm.proper","logreg","logreg.proper")))))
   {
     if (print.flag == TRUE) cat('\nMethod "',method[j],
     '" is not valid for a variable without predictors (',
     names(data)[j],')\nMethod has been changed to "sample"\n\n', sep = "")
     method[j] <- "sample"
   }
 }

 
 # chck survival method	and events are consistent
 if (any(method == "survctree")) {
   for (j in vis) {   # checks for survival variables
     y     <- data[,j]
     vname <- colnames(data)[j]
     mj    <- method[j]
     if (mj == "survctree") {
   	   if (!is.numeric(y)) stop("Variable ",vname,
         " should be a numeric survival time.")
    	 if (any(!is.na(y) & y < 0)) stop("Variable ",vname,          #!!!!!!!!!!!! BN added !is.na 
         " should be a non-negative survival time.")

       if (is.na((match(event[j],1:nvar)))) {
         cat("Variable ",vname," is a survival time. Corresponding event not in data, assuming no censoring.\n\n", sep = "")
         event[j] <- -1      # used to indicate no censoring
       }
    	 else {
         tabEI <- table(data[,event[j]])
         if (length(tabEI) != 2) {
           if (length(tabEI) == 1 & all(tabEI == 1)) cat("Variable ", vname,
             " is a survival time with all cases having events.\n", sep = "")
           else if (length(tabEI) == 1 & all(tabEI == 0)) stop("Variable ",
             vname," is a survival time with no cases having events.\n",
             "Estimation not possible.", sep = "")
           else stop("Event must be binary 0/1 but has values {", paste(
                names(tabEI), collapse = ", "),"}.\nNo data synthesised.")
         }
         if (!all(as.character(names(tabEI)) == c("0","1")) &&
             !is.logical(data[,event[j]])) {
           stop("Event must be binary 0/1 but it is a ", class(data[,event[j]]),
                " variable with values {", paste(names(tabEI)[1:2],
                collapse = ", "),"}.", sep = "")
         }
       }
     } else {
       # checks for non-survival variables
       if (event[j] != 0) {
         cat("Event for variable ", vname, " set to ", event[j],
 	           ' although method is "', mj, '". Event reset to 0.\n', sep = "")
         event[j] <- 0
       }
     }# end checking events for non  survival
   }# end of j loop
 }# end of checking when some surv variables
 else if (!all(event == 0)) {
   cat("No variables have a survival method so event vector which was \n{",
   paste(event, collapse = ","),"} set to 0s.\n\n", sep = "")
   event <- rep(0, nvar)
 }

 # change names for proper imputations and check
 #for(j in unique(vis)){
 #  if(proper==T & method[j]!="") method[j] <- paste(method[j],
 #                                                   ".proper",sep="")
 #}

 # check collinearity of variables
 if (sum(pred > 0)  & m > 0) {                                     
   inpred <- apply(pred != 0, 1, any) | apply(pred != 0, 2, any)
   if (any(inpred)) {
     collout <- collinear.out(data[, inpred, drop = FALSE])
     if (length(collout) > 0) {
       for (i in 1:length(collout)) {
         if (print.flag) cat("\nVariables ", paste(collout[[i]], collapse = ", ")," are collinear.\n", sep = "")
         vars <- match(collout[[i]], varnames[vis])
			   vfirst <- collout[[i]][vars == min(vars)]
			   nfirst <- match(vfirst,varnames)
         nall   <- match(collout[[i]],varnames)
         if (print.flag) cat("\nVariables later in visit.sequence are derived from ",vfirst,".\n", sep = "")
         for (ii in nall) {
           if (ii != nfirst) {
             method[ii] <- "collinear"
             pred[ii,]  <- 0
             pred[,ii]  <- 0
             pred[ii, nfirst] <- 1
           }
			   }
       }
     }                 
   } 
 } 

 setup$event            <- event
 setup$method           <- method
 setup$predictor.matrix <- pred
 setup$visit.sequence   <- vis
 setup$denom            <- denom                  #GRdenom new
 
 return(setup)
}
##--------------------end-of--check.method.syn-------------------------
 

##------------------check.rules.syn------------------------------------

check.rules.syn <- function(setup, data) {

 rules      <- setup$rules
 rvalues    <- setup$rvalues
 pred       <- setup$predictor.matrix
 nvar       <- setup$nvar                                                               
 varnames   <- setup$varnames
 method     <- setup$method
 vis        <- setup$visit.sequence
  
 # Check the syntax
 #------------------
 # check the length of the rules and corresponding values
 if (any(sapply(rules,length) != sapply(rvalues,length)))
   stop("The number of data rules for each variable should equal the number of corresponding values.\n  Check variable(s): ",
     paste(varnames[sapply(rules,length) != sapply(rvalues,length)], collapse = ", "), ".")

 # special characters 
 char.allowed <- c("","|","||","&","&&","==",">=","<=","<",">",
   "!=","==-",">=-","<=-","<-",">-","!=-","=='",".",")","(",";","-",
   "'","\"","\"(",")\"","'(",")'") #### . ( and ) added
 char.present <- paste(gsub("\\w"," ",unlist(rules)),collapse = " ") # remove word characters and concatenate
 char.present <- strsplit(char.present,"[[:space:]]+")[[1]]    # split into seperate characters
 char.wrong   <- !(char.present %in% char.allowed)             # identify unxepected characters
 #if (any(char.wrong)) stop("Unexpected character(s) in rules: ",paste(char.present[char.wrong],collapse=" "),".")

 # variables names (=a string before a special character) must be in varnames 
 rule.sep <- lapply(sapply(rules,strsplit,"[|&]"),unlist)       # split into seperate conditions
 get.vars <- lapply(rule.sep,function(x) gsub("[<>=!].*","",x)) # remove evrything after a special character
 get.vars <- lapply(get.vars,function(x) gsub(" ","",x))        # remove spaces
 get.vars <- lapply(get.vars,function(x) gsub("[\\(\\)]","",x)) # remove brackets  
 get.vars <- lapply(get.vars,function(x) gsub("is.na","",x))    # remove function name 
 get.vars <- lapply(get.vars,function(x) x[x != ""])            # remove empty strings  ?? why this
 
 vars.in.rules <- unique(unlist(get.vars))
 vars.wrong <- !(vars.in.rules %in% varnames)                   # identify unxepected variables
 
 if (any(vars.wrong)) stop("Unexpected variable(s) in rules: ",
   paste(vars.in.rules[vars.wrong], collapse = ", "), ".", call. = FALSE)
 # remove rules with warning for ipf and catall

 vars.with.rules <- varnames[rules != ""]
 if (any(method[varnames %in% vars.with.rules] %in% c("catall","ipf"))){
   cat("\nRules cannot be used for variables synthesised by ipf or catall")
   cat("\nbut values can be restricted by defining structural zero cells\nwith ipf.structzero or catall.structzero parameter.\n")
   rules[method %in% c("catall","ipf")] <-  rvalues[method %in% c("catall","ipf")] <- ""
   cat("\nRules defined for variable(s) ",
       paste0(varnames[method %in% c("catall","ipf") & varnames %in% vars.with.rules], collapse = ", "),
       " have been deleted.\n\n", sep = "")
   setup$rules <- rules
   setup$rvalues <- rvalues
   if (all(rules == "")) {
     return(setup)
   }
 }
 
 if (any(char.wrong)) {
   cat("One of rules may not be correct. If this is the case compare your rules and Error below.\nOtherwise rules have been applied.\n") 
   rs <- unlist(rules); names(rs) <- varnames
   rs <- cbind(rs[rs != ""]); colnames(rs) <- ""
   cat("\nYour rules are:")
   print(rs); cat("\n")
 }

 # Check that missingness in the data obeys the rules in rules
 nonmissing <- vector("list", nvar)
 isfactor   <- sapply(data, is.factor)
 yes.rules <- sapply(rules, function(x) any(x != ""))
 lth.rules <- sapply(rules, length)
 for (i in 1:nvar) {
   if (yes.rules[i]) {
     for (r in 1:lth.rules[i]) {
       if (is.na(rvalues[[i]][r]) & !isfactor[i]) {
         nonmissing[[i]][r] <- with(data,sum(!is.na(data[eval(parse(text = rules[[i]][r])), i])))
       } else if (is.na(rvalues[[i]][r]) & isfactor[i]) {    # different for factors because <NA> is treated as a level
         #nonmissing[[i]][r] <- with(data,sum(!is.na(as.character(data[eval(parse(text=rules[[i]][r])),i]))))
         nonmissing[[i]][r] <- with(data,sum(as.character(data[eval(parse(text = rules[[i]][r])),i]) != "NAtemp" &
                               as.character(data[eval(parse(text = rules[[i]][r])),i]) != "NAlogical"))       
       } else {
         nonmissing[[i]][r] <- with(data,sum(data[eval(parse(text = rules[[i]][r])),i] != rvalues[[i]][r] |
                               is.na(data[eval(parse(text = rules[[i]][r])),i])))
       }
     }
   }
 }
 any.nonmissing <- sapply(nonmissing, function(x) any(x > 0))
 if (any(any.nonmissing) > 0) cat("\nUnexpected values (not obeying the rules) found for variable(s): ",
     paste(varnames[any.nonmissing > 0], collapse = ", "),
     ".\nRules have been applied but make sure they are correct.\n", sep = "")

 # Check visit sequence 
 # all variables used in missing data rules have to be synthesised BEFORE 
 # the variables they apply to
 var.position <- lapply(get.vars, function(x) match(unique(x),varnames))
 var.in.vis   <- lapply(var.position, function(x) if (length(x) == 0) {
                                        x <- 0
                                        } else if (any(is.na(match(x,vis)))) {
                                        x[!is.na(match(x, vis))] <- match(x, vis)
                                        x[is.na(match(x, vis))]  <- nvar
                                        } else {
                                        x <- match(x,vis)})
 max.seq      <- sapply(var.in.vis, max, na.rm = T)
 not.synth    <- match(1:nvar,vis)[!is.na( match(1:nvar,vis))] <= max.seq[!is.na( match(1:nvar,vis))]
 if (any(not.synth,na.rm = TRUE)) stop("Variable(s) used in missing data rules for ",
       paste(varnames[!is.na( match(1:nvar,vis))][not.synth & !is.na(not.synth)], collapse = " "),
       " have to be synthesised BEFORE the variables they apply to.")

 # Check if a variable with missing values predicts other variables only if its
 # missing values are a subset of the missing values of the predicted variables
 # and remove from a prediction matrix if not. 
 # It refers to missing values coded as NA, otherwise variable can be used as 
 # a predictor without restrictions.
  
 #for (i in 1:nvar){
 #  if (!is.na(rvalues[i])) data[with(data,eval(parse(text=rules[i]))),i] <- NA
 #}
 patternRules <- matrix(0, nrow = nrow(data), ncol = ncol(data))
 for (i in 1:nvar) {
   if (yes.rules[i]) {
     for (r in 1:lth.rules[i]) {
       if (is.na(rvalues[[i]][r])) patternRules[with(data,eval(parse(text = rules[[i]][r]))), i] <- 1
     }
   }
 }
 patternNA <- is.na(data) + 0
 patternNA <- ifelse(patternRules == patternNA, patternNA, 0)
 diffNAij  <- function(i, j, dataNA) sum(dataNA[, i] - dataNA[, j] < 0)
 diffNA    <- Vectorize(diffNAij, vectorize.args = list("i", "j"))
 predNA    <- outer(1:nvar, 1:nvar, diffNA, dataNA = patternNA)
   
 # predNAwrong <- which ((pred==1 & predNA>0),arr.ind=TRUE)
 # pred        <- ifelse((pred==1 & predNA>0),0,pred)
 # if(nrow(predNAwrong)>0) cat(paste("Missing values of variable ",
 # varnames[predNAwrong[,2]]," are not a subset of missing values of variable ",
 # varnames[predNAwrong[,1]]," and cannot be used as its predictor (removed).\n",sep=""),
 # "\n",sep="")

 setup$predictor.matrix <- pred

 return(setup)
}
##-----------------end of--check.rules.syn----------------------------


##------------------namedlist------------------------------------
# check args that should be provided as a named list 
# and create list with elements for each variable 
 namedlist <- function(x, varnames = colnames(data), 
                       nvars = length(varnames), 
                       missval = NA, argname, argdescription = "", 
                       asvector = FALSE){
   if (is.null(x)) {
     x <- as.list(rep(missval,nvars))
   } else if (!is.list(x) | any(names(x) == "") | is.null(names(x))) {
     stop("Argument '", argname,"' must be a named list with names of selected ", 
     argdescription, " variables.", call. = FALSE)  
   } else {
     x.missval <- as.list(rep(missval,nvars))
     x.ind <- match(names(x), varnames)
     if (any(is.na(x.ind))) stop("Unrecognized variable names in '", 
       argname,"': ",paste(names(x)[is.na(x.ind)], collapse = ", "), call. = FALSE)
     # For 'event' and 'denom' check if denominators' name exist and 
     # change them to column indecies   
     if (argname %in% c("denom", "event") & is.character(argname)) {
       denom.ind <- lapply(x,match,varnames)
       if (any(is.na(denom.ind))) stop("Unrecognized variable(s) provided as ", argname, "(s): ", 
       paste(unlist(x)[is.na(denom.ind)], collapse = ", "), call. = FALSE) 
       x <- denom.ind
     }
     x.missval[x.ind] <- x
     x <- x.missval 
   }
   names(x) <- varnames
   if (asvector) x <- unlist(x)
   return(x)  
 }
##-----------------end of--namedlist-----------------------------


#----------------------- now syn continues here ----------------------
# Basic checks of provided parameters:
# dimensions, consistency, replication, ...

 call <- match.call()
 nvar <- ncol(data)
 if (!is.na(seed) & seed == "sample") {
   seed <- sample.int(1e9, 1)
   # cat("No seed has been provided and it has been set to ", seed,".\n\n", sep="")
 }
 if (!is.na(seed)) set.seed(seed)

 if (!(is.matrix(data) | is.data.frame(data)))
    stop("Data should be a matrix or data frame.")
 if (nvar < 2) stop("Data should contain at least two columns.")
                              
 # S U B S A M P L E   S I Z E
 if (k != nrow(data) & print.flag == TRUE & m > 0) {
  # if (k > nrow(data)) {
  #   cat("Warning: Subpopulation size (k=",k,") cannot be greater than the population size (",
  #       nrow(data),").\n","Synthetic data sets of same size as data will be produced.\n\n",sep="")
  #       k <- nrow(data)
  # } else
   cat("Sample(s) of size ", k, " will be generated from original data of size ",
         nrow(data),".\n\n", sep = "")
 }

 # M E T H O D S
 method <- gsub(" ","",method) # remove any spaces in or around method
 # expand user's syhthesising method (single string) to all variables
 if (length(method) == 1) {
   if (is.passive(method)) stop("Cannot have a passive syhthesising method for every column.")
   method <- rep(method, nvar)
   if (!(method[1] %in% c("catall","ipf"))) method[visit.sequence[1]] <- "sample" #GRBN
   # set method to "" for vars not in visit.sequence
   method[setdiff(1:nvar,visit.sequence)] <- ""
 }
 
 # if user specifies multiple methods, check the length of the argument
 # methods must be given for all columns in the data
 if (length(method) != nvar) stop(paste("The length of method (", length(method),
    ") does not match the number of columns in the data (",nvar,").", sep = ""), 
    call. = FALSE)

 # P R E D I C T O R   M A T R I X
 if (!is.null(predictor.matrix)) {
   if (!is.matrix(predictor.matrix)) {
     stop("Argument 'predictor.matrix' is not a matrix.")
   } else if (nvar != nrow(predictor.matrix) | nvar != ncol(predictor.matrix))
     stop(paste("The 'predictor.matrix' has ",nrow(predictor.matrix),
          " row(s) and ", ncol(predictor.matrix),
          " column(s). \nBoth should match the number of columns in the data (",
          nvar, ").", sep = ""))
 }

 data     <- as.data.frame(data)
 varnames <- dimnames(data)[[2]]
 
 # Named lists: check args and create list with elements for each variables 
 # C O N T I N O U S  V A R S  W I T H  M I S S I N G  D A T A  C O D E S
 # S E M I - C O N T I N O U S  V A R S 
 semicont <- namedlist(semicont, missval = NA, argname = "semicont", 
   argdescription = "semi-continuous")
 cont.na  <- namedlist(cont.na, missval = NA, argname = "cont.na", 
   argdescription = "")
 # combine cont.na and semicont lists  
 cont.na.ini <- cont.na
 cont.na <- mapply(c, cont.na, semicont, SIMPLIFY = FALSE)
 cont.na <- lapply(cont.na, unique)
 # R U L E S  and  R V A L U E S 
 rules   <- namedlist(rules, missval = "", argname = "rules", 
                      argdescription = "")
 rvalues <- namedlist(rvalues, missval = NA, argname = "rvalues", 
                      argdescription = "")
 # S M O O T H I N G
 smoothing <- namedlist(smoothing, missval = "", argname = "smoothing", 
   argdescription = "", asvector = TRUE)
 if (any(smoothing != "")) {
   varsmoothind <- which(smoothing != "")
   varnumind    <- which(sapply(data,is.numeric))
   smoothnumind <- match(varsmoothind,varnumind)
   if (any(is.na(smoothnumind)) & print.flag == TRUE)
   cat("\nSmoothing can only be applied to numeric variables.\nNo smoothing will be applied to variable(s): ",
     paste(varnames[varsmoothind[is.na(smoothnumind)]], collapse = ", "), "\n", sep = "")   
   smoothing[varsmoothind[is.na(smoothnumind)]] <- "" 
 }
 
 # D E N O M   
 denom <- namedlist(denom, missval = 0, argname = "denom", asvector = TRUE)
 # E V E N T
 event <- namedlist(event, missval = 0, argname = "event", asvector = TRUE)
 
 # Perform various validity checks on the specified arguments
 setup <- list(visit.sequence = visit.sequence,
               method = method,
               default.method = default.method,
               predictor.matrix = predictor.matrix,
               nvar = nvar,
               varnames = varnames, 
               rules = rules,
               rvalues = rvalues,
               cont.na = cont.na,                                         
               event = event,
               denom = denom)                   #GRdenom new
              
 setup <- check.visit.sequence.syn(setup)
 setup <- check.predictor.matrix.syn(setup)


# C H A N G E  D A T A  T Y P E  &  M O D I F Y  F A C T O R  L E V E L S
#---
 # apply only if in predictor matrix 
                                                                               # GR added condition and else
 if (!is.null(setup$predictor.matrix) & sum(setup$predictor.matrix > 0)) {
   inpred <- apply(setup$predictor.matrix != 0, 1, any)*(!(method %in% c("","sample"))) |                # GR added to allow null methods not affected
             apply(setup$predictor.matrix != 0, 2, any)  # if anywhere in predictor.matrix
 } else {
   inpred <- rep(FALSE, sqrt(length(setup$predictor.matrix)))
 }
 notevent <- is.na(match(1:nvar,setup$event))       # if not in event list

 # Convert any character variables into factors for variables in pred
 ischar    <- sapply(data,is.character)
 chartofac <- (ischar * inpred) > 0
 if (sum(chartofac) > 0) {
   cat("\nVariable(s):",paste0(varnames[chartofac], collapse = ", "),
       "have been changed from character to factor.\n", sep = " ")
   for (j in (1:nvar)[chartofac]) data[,j] <- as.factor(data[,j]) 
 }

 # Changing numeric variables with fewer than 'minnumlevels' into factors
 #  Default for this now set to -1 (no changes are made) 
 #  Warning if numeric vars with < 5 levels (20 too many as picks up months)
 #  Also only need to do this if variable in predictionMatrix
 #  and any inappropriate methods are changed to the default for factors
 nlevel      <- sapply(data, function(x) length(table(x)))
 ifnum       <- sapply(data, is.numeric)
 innumtocat  <- rep(FALSE,length(data)) #GRBN
 
 if (minnumlevels < 5 & any(nlevel > minnumlevels & nlevel <= 5 & ifnum & inpred & notevent)) {
   cat("Warning: In your synthesis there are numeric variables with 5 or less levels:\n",
       paste0(varnames[nlevel > minnumlevels & nlevel <= 5 & ifnum & inpred & notevent], collapse = ", "), ".",
       " Consider changing them to factors.\nYou can do it using parameter `minnumlevels`.\n\n", sep = "")  
 }

 vartofactor <- which(nlevel <= minnumlevels & ifnum & inpred & notevent)
 for (j in vartofactor) data[,j] <- as.factor(data[,j])
 if (length(vartofactor) > 0) {
   cat("\nVariable(s): ", paste0(varnames[vartofactor], collapse = ", "),
       "numeric but with fewer than", minnumlevels, 
       "levels turned into factor(s).\n\n", sep = " ")
   for (j in vartofactor) {
     if (setup$method[j] %in% c("norm","norm.proper",
                                "normrank","normrank.proper")) {
       if (nlevel[j] == 2) setup$method[j] <- default.method[2]
       else setup$method[j] <- default.method[3]
  	   cat("Method for ",varnames[j]," changed to ",setup$method[j],"\n\n")
     }
   }
 }

 # Modifies a factor by turning NA into an extra level
 isfactor  <- sapply(data,is.factor)
 for (j in (1:nvar)[isfactor & inpred & notevent]) {
   data[,j] <- addNA(data[,j], ifany = TRUE)
   levels(data[,j])[is.na(levels(data[,j]))] <- "NAtemp"          #!BN 10/08/15 
 } 

 islogicalNA <- sapply(data, function(x) (is.logical(x) & any(is.na(x))))   #!BN 20/04/16 
 for (j in (1:nvar)[islogicalNA & inpred & notevent]) {
   data[,j] <- addNA(data[,j], ifany = TRUE)
   levels(data[,j])[is.na(levels(data[,j]))] <- "NAlogical"          
 }

 
#---
 # browser()
 setup  <- check.method.syn(setup, data, proper)
 if (any(rules != "")) setup <- check.rules.syn(setup, data)

 method           <- setup$method
 predictor.matrix <- setup$predictor.matrix
 visit.sequence   <- setup$visit.sequence
 event            <- setup$event
 rules            <- setup$rules
 rvalues          <- setup$rvalues
 cont.na          <- setup$cont.na
 default.method   <- setup$default.method
 denom            <- setup$denom                       #GRdenom new
 ############################################################!GRipf added this as seemed right
 #
 method[!(1:length(method) %in% visit.sequence)] <- ""  #GR?BN?

 # Identify any factors with > maxfaclevels levels that are in visit.sequence
 no.fac.levels <- sapply(data, function(x) length(levels(x)))
 too.many.levels <- no.fac.levels > maxfaclevels
 notsampling <- !(grepl("nested", method) | grepl("sample", method) | grepl("satcat", method) | grepl("constant", method))
 if (any(inpred & too.many.levels & notsampling)) {
   stop("Factor(s) with more than ",maxfaclevels," levels: ",
     paste0(varnames[inpred & too.many.levels]," (",
     no.fac.levels[inpred & too.many.levels],")", collapse = ", "), call. = FALSE,
     "\nYou can increase 'maxfaclevels' to continue but it may cause computational problems. \nConsider removing factor(s) from prediction matrix, combining categories or using 'nested' method.\n\n")
 }

 # Not used variables are identified and dropped if drop.not.used==T
 # reclculate inpred & notevent in case they have changed after
 # check.method and check.data
 if (sum(predictor.matrix) > 0) {                            # GR condition added
   inpred      <- apply(predictor.matrix != 0, 1, any) | apply(predictor.matrix != 0, 2, any) # if anywhere in predictor.matrix
	 ispredictor <- apply(predictor.matrix != 0, 2, any)    # if used as predictor
 }
 else inpred <- ispredictor <- rep(0, sqrt(length(predictor.matrix))) 

 notinvs     <- is.na(match(1:nvar,visit.sequence)) # if not in visit.sequence
 notsynth    <- notinvs | (!notinvs & method == "")  # if not synthesised
 notevent    <- is.na(match(1:nvar,event))         # if not in event list

 # identify columns not used as events or predictors or in visitSequnce
 out <- !inpred & notevent & notsynth

 if (any(out) & print.flag == TRUE) {
   cat("\nVariable(s):", paste0(varnames[out], collapse = ", "),
       "not synthesised or used in prediction.\n", sep = " ")
   if (drop.not.used == T) cat("The variable(s) will be removed from data and not saved in synthesised data.\n\n")
   else cat("CAUTION: The synthesised data will contain the variable(s) unchanged.\n\n")
 }

#browser()   
# remove columns not used from data and replace predictor matrix, visit sequence, nvar and others
 if (any(out) & drop.not.used == T) {
   if (sum(!out) == 0) stop("No variables left to be synthesised") ######to stop if all data excluded 
   newnumbers        <- rep(0,nvar)
   newnumbers[!out]  <- 1:sum(!out)
   visit.sequence    <- newnumbers[visit.sequence]
   visit.sequence    <- visit.sequence[!visit.sequence == 0]
   predictor.matrix  <- predictor.matrix[!out,!out]
   
   event[event != 0] <- newnumbers[event[event != 0]]
   event             <- event[!out]
   denom[denom != 0] <- newnumbers[denom[denom != 0]]
   denom             <- denom[!out]
   
   data              <- data[,!out]
   nvar              <- sum(!out)
   method            <- method[!out]
   varnames          <- varnames[!out]

   if (nvar == 1) {                             #  GR added  note having to reassign character vector
   	 cl <- class(data)
  	 data <- data.frame(data)
  	 if (cl == "character") data[,1] <- as.character(data[,1])
  	 names(data) <- varnames
   }
   cont.na          <- cont.na[!out]
   cont.na.ini      <- cont.na.ini[!out]      #BN13/11  
   semicont         <- semicont[!out]         #BN13/11  
   smoothing        <- smoothing[!out]
   rules            <- rules[!out]
   rvalues          <- rvalues[!out]
   var.lab          <- var.lab[!out]
   val.lab          <- val.lab[!out]
   
   # recalculate these
   if (sum(predictor.matrix > 0)) {                                 # GR condition added
	   inpred <- apply(predictor.matrix != 0, 1, any) |
	             apply(predictor.matrix != 0, 2, any)      # if anywhere in predictor.matrix
	   ispredictor <- apply(predictor.matrix != 0, 2, any) # if used as predictor
   }
   else inpred <- ispredictor <- rep(0,sqrt(length(predictor.matrix))) 

   notinvs  <- is.na(match(1:nvar,visit.sequence)) # if not in visit.sequence
   notsynth <- notinvs | (!notinvs & method == "") # if not synthesised
   notevent <- is.na(match(1:nvar,event))          # if not in event list
 }

 # Print out info on variables not synthesised but used in prediction
 pred.not.syn <- (ispredictor & notsynth)
 if (sum(pred.not.syn ) > 0 & drop.pred.only == FALSE) pred.not.syn[pred.not.syn == TRUE] <- FALSE
 
 if (sum(pred.not.syn ) > 0 & print.flag == TRUE) {
   cat("\nVariable(s):", paste0(varnames[ispredictor & notsynth], collapse = ", "),
       "used as predictors but not synthesised.\n", sep = " ")
   if (drop.pred.only == T) {
     cat("The variable(s) will not be saved with the synthesised data.\n")
   } else {
     cat("CAUTION: The synthesised data will contain the variable(s) unchanged.\n")
   }
 } 
 #browser()                                                     #  GR condition added
 if (sum(predictor.matrix) > 0) {
	 
   pm <- padMis.syn(data, method, predictor.matrix, visit.sequence,
			   nvar, rules, rvalues, default.method, cont.na, smoothing, event, denom)

	 # Pad the Syhthesis model with dummy variables for the factors
	 # p <- padModel.syn(data, method, predictor.matrix, visit.sequence,
	 #                   nvar, rules, rvalues)
	 p  <- padModel.syn(pm$data, pm$method, pm$predictor.matrix, pm$visit.sequence,
			   pm$nvar, pm$rules, pm$rvalues, pm$factorNA, pm$smoothing, pm$event, pm$denom)
   if (k != dim(data)[1]) {
     # create a non-empty data frame in case some variables are kept unsynthesised
     p$syn <- p$syn[sample(1:nrow(data), k, replace = TRUE),]
     dimnames(p$syn)[[1]] <- 1:k                             #!BN-15/06/2016
   }
   if (sum(duplicated(names(p$data))) > 0)
     stop("Column names of padded data should be unique.")
 
	 p$cont.na <- pm$cont.na  #!bn
	 
 } else {
 
   p <- list(data = data,                             #  GR this added
             syn = data,
             predictor.matrix = predictor.matrix, 
             method = method, 
             visit.sequence = visit.sequence, 
             rules = rules,
             rvalues = rvalues,
             cont.na = cont.na, 
             event = event,                  #GRdenom new
             denom = denom,                  #GRdenom new
             categories = NULL,
             smoothing = smoothing)         
   
   if (k != dim(data)[1]) {                                  #!BN-27/06/2018
     p$syn <- p$syn[sample(1:nrow(data), k, replace = TRUE),]   
     dimnames(p$syn)[[1]] <- 1:k
   }
 }

 if (m > 0) {
#  syn <- list(m)
#	 for(i in 1:m){
#     syn[[i]] <- data
#     if (k != dim(data)[1]) syn[[i]] <- syn[[i]][sample(1:dim(data)[1], k, replace = TRUE), ]
#   }    
   syn <- vector("list",m)                        #!BN-15/06/2016
   for (i in 1:m) {                               #!BN-15/06/2016
     syn[[i]] <- data[0,]                         #!BN-15/06/2016
     #if (k > 0){                                 #!BN-12/08/2016 - for syn.strata
       syn[[i]] <- syn[[i]][1:k,]                 #!BN-15/06/2016
       dimnames(syn[[i]])[[1]] <- 1:k             #!BN-15/06/2016
     #}                                           #!BN-12/08/2016
   } 
 }
 else syn <- NULL

 synall <- sampler.syn(p, data, m, syn, visit.sequence, rules, rvalues, 
    event, proper, print.flag, k, pred.not.syn, models, numtocat, ...)
 
 syn <- synall$syn
 fits <- synall$fits
 
 if (m == 1) {
   syn <- syn[[1]]
   fits <- fits[[1]]
 } 
 

#-----------------------------------------------------------------------
# restore the original NA's in the data
# for(j in p$visit.sequence) p$data[(!r[,j]),j] <- NA
 
 names(method)         <- varnames
 names(visit.sequence) <- varnames[visit.sequence]

# Put grouped numeric variables and their synthesising parameters 
# back into their correct positions
#---
 if (!is.null(numtocat)) { 
   out <- (length(method) - length(numtocat) + 1:length(numtocat))
   if (m == 1) {
     tocat <- match(numtocat, names(syn))
     syn[, tocat] <- syn[, out]
     syn <- syn[, -out]
   } else {
     for (i in 1:m) {
       tocat <- match(numtocat, names(syn[[1]]))
       syn[[i]][, tocat] <- syn[[i]][, out]
       syn[[i]] <- syn[[i]][,-out]
     }
   }
   
   # move their parameters to numeric versions
   method             <- method[-out]
   visit.sequence     <- visit.sequence[-out]
   smoothing[tocat]   <- smoothing[out]
   smoothing          <- smoothing[-out]
   cont.na.ini[tocat] <- cont.na.ini[out]
   cont.na.ini        <- cont.na.ini[-out]
   rules[tocat]       <- rules[out]
   rules              <- rules[-out]
   rvalues[tocat]     <- rvalues[out]
   rvalues            <- rvalues[-out]
   predictor.matrix   <- predictor.matrix[-out,-out]
 }

# Tidy models
#---
 if (models) {
   if (m == 1) { 
     # drop fits from dummies and nulls
     fitout <- sapply(fits, function(x) is.null(x) || (is.character(x) && x =="dummy")) | 
               grepl("orig.", names(fits))
     if (any(fitout)) fits <- fits[!fitout]
     # move the models for non-missing values to original position
     n.0 <- grepl(".0", names(fits))      
     if (any(n.0)) {
       fits[gsub(".0", "", names(fits[n.0]))] <- fits[n.0]
       fits <- fits[!n.0]
     }
   }
   if (m > 1) { 
     for (j in 1:m) {
       fitout <- sapply(fits[[j]], function(x) is.null(x) || (is.character(x) && x =="dummy")) |
                 grepl("orig.", names(fits))
       if (any(fitout)) fits[[j]] <- fits[[j]][!fitout]
       n.0 <- grepl(".0",names(fits[[j]]))
       if ( any(n.0)) {
         fits[[j]][gsub(".0","",names(fits[[j]][n.0]))] <- fits[[j]][n.0]
         fits[[j]] <- fits[[j]][!n.0]
       }
     }
   }
 }
 
# Save, and return, but don't return data
#---
 syndsobj <- list(call = call,
                  m = m,
                  syn = syn,
                  method = method,
                  visit.sequence = visit.sequence,
                  predictor.matrix = predictor.matrix,
                  smoothing = smoothing,
                  event = event,
                  denom = denom,                  
#                 minbucket = minbucket,
                  proper = proper,
                  n = nrow(data),
                  k = k,
                  rules = rules,
                  rvalues = rvalues,
                  cont.na = cont.na.ini,          
                  semicont = semicont,            
                  drop.not.used = drop.not.used,
                  drop.pred.only = drop.pred.only,
                  models = fits,
                  seed = seed,
                  var.lab = var.lab,
                  val.lab = val.lab,
                  obs.vars = obs.vars,
                  numtocat = numtocat,                      #!GRipf 2 new parameters
                  catgroups = catgroups)
 # if (diagnostics) syndsobj <- c(syndsobj, list(pad = p))
 class(syndsobj) <- "synds"
 return(syndsobj)
}


###-----collinear.out------------------------------------------------------

collinear.out <- function(x, threshold = 0.99999) {
  nvar     <- ncol(x)
  x        <- data.matrix(x)
  varnames <- dimnames(x)[[2]]
  z        <- suppressWarnings(cor(x, use = "pairwise.complete.obs"))
  z[is.na(z)] <- 0    #GR-09/2016 to handle all NAs
  hit      <- abs(z) >= threshold
  mysets   <- !duplicated(hit) & rowSums(hit) > 1
  if (sum(mysets) == 0) setlist <- NULL
  else {
    setlist  <- vector("list",sum(mysets))
    for (i in 1:sum(mysets)) setlist[[i]] <- colnames(hit)[hit[mysets, , drop = FALSE][i,]]  
  } 
  return(setlist)
}
