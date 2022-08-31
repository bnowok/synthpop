###-----utility.gen--------------------------------------------------------
utility.gen <- function(object, data, ...) UseMethod("utility.gen")


###-----utility.gen.default------------------------------------------------
utility.gen.default <- function(object, ...)
  stop("No compare method associated with class ", class(object), call. = FALSE)


###-----utility.gen.data.frame---utility.gen.list--------------------------
utility.gen.data.frame <- utility.gen.list <-
    function(object, data, 
             not.synthesised = NULL, cont.na = NULL,
             method = "cart", maxorder = 1,
             k.syn = FALSE, tree.method = "rpart",
             max.params = 400, print.stats = c("pMSE", "S_pMSE"),
             resamp.method = NULL, nperms = 50, cp = 1e-3, 
             minbucket = 5, mincriterion = 0, vars = NULL,
             aggregate = FALSE, maxit = 200, ngroups = NULL, 
             print.flag = TRUE, print.every = 10, 
             digits = 6, print.zscores = FALSE, zthresh = 1.6,
             print.ind.results = FALSE,
             print.variable.importance = FALSE, ...)
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

  res <- utility.gen.synds(object = object, data = data, 
                           method = method, maxorder = maxorder, 
                           k.syn = k.syn, tree.method = tree.method,
                           max.params = max.params, print.stats = print.stats,
                           resamp.method = resamp.method, nperms = nperms, cp = cp, 
                           minbucket = minbucket, mincriterion = mincriterion, 
                           vars = vars, aggregate = aggregate, maxit = maxit, 
                           ngroups = ngroups, print.flag = print.flag, 
                           print.every = print.every, digits = digits, 
                           print.zscores = print.zscores, zthresh = zthresh, 
                           print.ind.results = print.ind.results, 
                           print.variable.importance = print.variable.importance)
  res$call <- match.call()
  return(res)
}


###-----utility.gen--------------------------------------------------------
utility.gen.synds <- function(object, data, 
                              method = "cart", maxorder = 1,
                              k.syn = FALSE, tree.method = "rpart", 
                              max.params = 400, print.stats = c("pMSE", "S_pMSE"),
                              resamp.method = NULL, nperms = 50, cp = 1e-3,
                              minbucket = 5, mincriterion = 0, vars = NULL,
                              aggregate = FALSE, maxit = 200, ngroups = NULL,
                              print.flag = TRUE, print.every = 10,
                              digits = 6, print.zscores = FALSE, 
                              zthresh = 1.6, print.ind.results = FALSE,
                              print.variable.importance = FALSE, ...)
{
 m  <- object$m

 # Check input parameters
 if (is.null(method) || length(method) != 1 || is.na(match(method, c("cart", "logit"))))
   stop("Invalid 'method' type - must be either 'logit' or 'cart'.\n", call. = FALSE)
 if (is.null(print.stats) || any(is.na(match(print.stats, c("pMSE", "SPECKS", "PO50", "U", "S_pMSE", "S_SPECKS", "S_PO50", "S_U", "all")))))
   stop("Invalid 'print.stats'. Can only include 'pMSE', 'SPECKS', 'PO50', 'U', 'S_pMSE', 'S_SPECKS', 'S_PO50', 'S_U'.\nAternatively it can be set to 'all'.\n", call. = FALSE)
 if (!is.null(resamp.method) && is.na(match(resamp.method, c("perm", "pairs", "none"))))
   stop("Invalid 'resamp.method' type - must be NULL, 'perm', 'pairs' or 'none'.\n", call. = FALSE)
 if (aggregate == TRUE & method != "logit") stop("Aggregation only works for 'logit' method.\n", call. = FALSE)
 if (is.null(data)) stop("Requires parameter 'data' to give name of the real data.\n",  call. = FALSE)
 if (!inherits(object, "synds")) stop("Object must have class 'synds'.\n", call. = FALSE)
 if (k.syn & !is.null(resamp.method) && resamp.method == "pairs") stop('\nresamp.method = "pairs" will give the wrong answer when k.syn is TRUE.\n', call. = FALSE)
 if (is.null(tree.method) || length(tree.method) != 1 || is.na(match(tree.method, c("rpart", "ctree"))))
   stop("Invalid 'tree.method' - must be either 'rpart' or 'ctree'.\n", call. = FALSE)
 
 
 # Check selected variables and make observed and synthetic comparable
 if (!(is.null(vars))) {
   if (is.numeric(vars)){
     if (!(all(vars %in% 1:length(data)))) stop("Column indices of 'vars' must be in 1 to length(data).\n", call. = FALSE)
   } else if (!(all(vars %in% names(data)))) stop("Some 'vars' specified not in data.\n", call. = FALSE)
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

 # Set default resampling according to completeness and print.stats (incl. S_SPECKS or S_PO50 or S_U)
 if (is.null(resamp.method)) {
   if ("S_SPECKS" %in% print.stats || "S_PO50" %in% print.stats || "S_U" %in% print.stats || incomplete) {
     resamp.method <- "pairs"
     cat('Resampling method set to "pairs" because S_SPECKS or S_PO50 or S_U in print.stats or incomplete = TRUE.\n') 
   } else if (method == "cart") resamp.method <- "perm"
 } else {
   if (incomplete & resamp.method == "perm")
     stop('Incomplete synthesis requires resamp.method = "pairs".\n', call. = FALSE)
   if (any(c("S_SPECKS", "S_PO50", "S_U") %in% print.stats) & resamp.method == "perm")
     stop('Stat SPECKS, PO50, and U requires resamp.method = "pairs" to get S_SPECKS, S_PO50, and S_U respectively.\n', call. = FALSE)
   if (resamp.method == "pairs" & m == 1) 
     stop('resamp.method = "pairs" needs a synthesis with m > 1, m = 10 suggested.\n', call. = FALSE)
 }

 # Drop any single value columns
 leneq1 <- function(x) length(table(as.numeric(x[!is.na(x)]), useNA = "ifany")) %in% (0:1)

 dchar <- sapply(data,is.character)
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
 numvars <- (1:dim(data)[2])[sapply(data, is.numeric)]
 names(numvars) <- names(data)[numvars]
 # If ngroups != NULL divide numeric variables into ngroups
 data0 <- data  # to save if m > 1

 if (!is.null(ngroups)) {
   for (i in numvars) {
     if (m == 1) {
       groups <- group_num(data[,i], object$syn[,i], object$syn[,i], 
                           ngroups, cont.na = cna, ...)
       data[,i] <- groups[[1]]
       object$syn[,i] <- groups[[2]]
     } else {
       syn0 <- c(sapply(object$syn, '[[', i)) 
       for (j in 1:m) {
         groups <- group_num(data0[,i], object$syn[[j]][,i], syn0,
                             ngroups, cont.na = cna[[i]], ...)
         data[,i] <- groups[[1]]
         object$syn[[j]][,i] <- groups[[2]]
       }
     }
   }
 }
 
 # Categorical vars: make missing data part of factor
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
   catvars <- (1:(dim(df.prop)[2]))[sapply(df.prop,is.factor)]
   for (i in catvars) {
     if (any(table(df.prop[,i]) == 0)) {
       df.prop[,i] <- as.factor(as.character(df.prop[,i]))
       if (print.flag) cat("Empty levels of factor(s) for variable ", names(df.prop)[i]," removed.\n" )
     }
   }

   if (aggregate == TRUE) {
     aggdat <- aggregate(df.prop[,1], by = df.prop, FUN = length)
     wt <- aggdat$x
     aggdat <- aggdat[, -dim(aggdat)[2]]
   }

   if (method == "logit" ) {
     
     if (maxorder >= dim(data)[2])
       stop("maxorder cannot be greater or equal to the number of variables.\n", call. = FALSE)

     # cheking for large models
     levs <- sapply(data, function(x) length(levels(x)))
     levs[levs == 0] <- 2
     tt1 <- apply(combn(length(levs), 1), 2, function(x) {prod(levs[x] - 1)})
     if (maxorder == 0) nparams <- 1 + sum(tt1)
     else {
       tt2 <- apply(combn(length(levs), 2), 2, function(x) {prod(levs[x] - 1)})
       if (maxorder == 1) nparams <- 1 + sum(tt1) + sum(tt2)
       else {
         tt3 <- apply(combn(length(levs), 3), 2, function(x) {prod(levs[x] - 1)})
         if (maxorder == 2) nparams <- 1 + sum(tt1) + sum(tt2) + sum(tt3)
         else {
           tt4 <- apply(combn(length(levs), 4), 2, function(x) {prod(levs[x] - 1)})
           if (maxorder == 3) nparams <- 1 +  sum(tt1) + sum(tt2) + sum(tt3) + sum(tt4)
           else {
             tt5 <- apply(combn(length(levs), 5), 2, function(x) {prod(levs[x] - 1)})
             if (maxorder == 4) nparams <- 1 +  sum(tt1) + sum(tt2) + sum(tt3) + sum(tt4) + sum(tt5)
           }  
         }  
       }  
     }
     if (nparams > max.params) stop("You will be fitting a large model with ", nparams, 
       " parameters that may take a long time and fail to converge.
Have you selected variables with vars?
You can try again, if you really want to, by increasing max.params.\n", sep = "", call. = FALSE)
     else if (nparams > dim(data)[[1]]/5) cat("You will be fitting a large model with ", nparams,
       " parameters and only ", dim(data)[[1]], " records
that may take a long time and fail to converge.
Have you selected variables with vars?\n")

     if (maxorder >= 1) logit.int <- as.formula(paste("t ~ .^", maxorder + 1))
     else logit.int <- as.formula(paste("t ~ ."))

     if (aggregate == TRUE) fit <- glm(logit.int, data = aggdat, family = "binomial",
                                       control = list(maxit = maxit), weights = wt)
     else fit <- suppressWarnings(glm(logit.int, data = df.prop, family = "binomial",
                                      control = list(maxit = maxit)))
     #if (fit$converged == FALSE) cat("\nConvergence failed.\n")

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
     # one more coefficient (intercept needed if k.syn TRUE)
     if (k.syn) km1 <- km1 + 1
     if (aggregate == TRUE) {
       pMSE <- (sum(wt*(score - cc)^2, na.rm = T)) / N
       KSt <- suppressWarnings(ks.test(rep(score[aggdat$t == 1], wt[aggdat$t == 1]),
                                       rep(score[aggdat$t == 0], wt[aggdat$t == 0])))
       SPECKS <- KSt$statistic
       PO50 <- sum(wt[(score > 0.5 & df.prop$t == 1) | ( score <= 0.5 & df.prop$t == 0)])/N*100 - 50
       U      <- suppressWarnings(wilcox.test(rep(score[aggdat$t == 1], wt[aggdat$t == 1]),
                                               rep(score[aggdat$t == 0], wt[aggdat$t == 0]))$statistic) 
     } else {
       pMSE <- (sum((score - cc)^2, na.rm = T)) / N
       KSt <- suppressWarnings(ks.test(score[df.prop$t == 1], score[df.prop$t == 0]))
       SPECKS <- KSt$statistic
       PO50 <- sum((score > 0.5 & df.prop$t == 1) | ( score <= 0.5 & df.prop$t == 0))/N*100 - 50
       U      <- suppressWarnings(wilcox.test(score[df.prop$t == 1], score[df.prop$t == 0])$statistic) 
     }
     pMSEExp <- km1 * (1 - cc)^2 * cc / N
     S_pMSE  <- pMSE / pMSEExp

     # to save space
     fit$data <- NULL
     # fit$model <- fit$residuals <- fit$y <- NULL ?

   } else if (method == "cart") {
     km1 <- NA
     if (tree.method == "rpart") {
       fit <- rpart(t ~ ., data = df.prop, method = 'class',
                    control = rpart.control(cp = cp, minbucket = minbucket))
       score <- predict(fit)[, 2]
     } else if (tree.method == "ctree") {
       fit <- ctree(t ~ ., data = df.prop,
                    controls = ctree_control(mincriterion = mincriterion, minbucket = minbucket))
       score <- predict(fit)
     }
     pMSE <- sum((score - cc)^2, na.rm = T) / N
     KSt <- suppressWarnings(ks.test(score[df.prop$t == 1], score[df.prop$t == 0]))
     SPECKS <- KSt$statistic
     PO50 <- sum((score > 0.5 & df.prop$t == 1) | ( score <= 0.5 & df.prop$t == 0))/N*100 - 50
     U <- suppressWarnings(wilcox.test(score[df.prop$t == 1], score[df.prop$t == 0])$statistic)
   }

   # Permutations
   if (!is.null(resamp.method) && resamp.method == "none") S_pMSE <- NA
   else if (!is.null(resamp.method) && resamp.method == "perm") { # to allow resamp for logit models
     S_pMSE <- rep(NA, m)
     simpMSE <- rep(0, nperms)
     if (m == 1) j <- 1
     if (j == 1 & print.flag) {
       if (print.every == 0 | print.every >= nperms) cat("Running ", nperms, " permutations to get NULL utilities.", sep = "")
       else cat("Running ", nperms, " permutations to get NULL utilities and printing every ", print.every, "th.", sep = "")
     }
     #if (print.flag) cat("\nsynthesis ", j, "   ", sep = "")
     if (print.flag) cat("\nsynthesis ")
     
     for (i in 1:nperms) {
       if (print.every > 0 & nperms > print.every & floor(i/print.every) == i/print.every & print.flag)  cat(i, " ", sep = "")
       pdata <- df.prop
       if (!k.syn) pdata$t <- sample(pdata$t)
       else pdata$t <- rbinom(N, 1, 0.5)

       if (method == "cart") {
         if (tree.method == "rpart") {
           sfit <- rpart(t ~ ., data = pdata, method = 'class', control = rpart.control(cp = cp, minbucket = minbucket))
           score <- predict(sfit)[,2]
         } else if (tree.method == "ctree") {
           sfit <- ctree(t ~ ., data = pdata,
                         controls = ctree_control(mincriterion = mincriterion, minbucket = minbucket))
           score <- predict(sfit)
         }
         simpMSE[i] <- (sum((score - cc)^2, na.rm = T)) / N / 2
       
       } else if (method == "logit") {
         if (maxorder >= 1) logit.int <- as.formula(paste("t ~ .^", maxorder + 1))
         else logit.int <- as.formula(paste("t ~ ."))

         if (aggregate == TRUE) {
           aggdat1 <- aggregate(pdata[,1], by = pdata, FUN = length)
           wt <- aggdat1$x
           aggdat1 <- aggdat1[, -dim(aggdat1)[2]]
           sfit <- glm(logit.int, data = aggdat1, family = "binomial",
                       control = list(maxit = maxit), weights = wt)
         } else sfit <- glm(logit.int, data = pdata, family = "binomial", 
                            control = list(maxit = maxit))

         if (sfit$converged == FALSE & print.flag) cat("Warning: Logistic model did not converge in ",
           maxit, " iterations.\nYou could try increasing parameter 'maxit'.\n", sep = "")
         score <- predict(sfit, type = "response")
         if (aggregate == TRUE) {
           simpMSE[i] <- sum(wt*(score - cc)^2, na.rm = T) / N / 2  # reduced by factor of 2
         } else {
           simpMSE[i] <- sum((score - cc)^2, na.rm = T) / N / 2  # reduced by factor of 2
         }
       }
     }
     nnosplits <- c(sum(simpMSE < 1e-8), length(simpMSE))
     S_pMSE <- pMSE/mean(simpMSE)
   }
   if (!is.null(resamp.method) && resamp.method == "pairs") 
     res.ind <- list(pMSE = pMSE, SPECKS = SPECKS, PO50 = PO50, U = U, 
                     S_pMSE= NA, S_SPECKS = NA,  S_PO50 = NA, S_U = NA, 
                     fit = fit, nnosplits = NA, df = NA)
   else if (!is.null(resamp.method) && resamp.method == "perm") 
     res.ind <- list(pMSE = pMSE, SPECKS = SPECKS, PO50 = PO50,U = U,
                     S_pMSE= S_pMSE, S_SPECKS = NA, S_PO50 = NA, S_U = NA,
                     fit = fit, nnosplits = nnosplits, df = NA)
   else res.ind <- list(pMSE = pMSE, SPECKS = SPECKS, PO50 = PO50, U =U,
                        S_pMSE = S_pMSE, S_SPECKS = NA, S_PO50 = NA, S_U = NA,
                        fit = fit, nnosplits = NA, df = km1) ## changed to NA
   return(res.ind)
 }
 # --------------------------------------
 # end propcalcs

  n1 <- nrow(data)

  if (m == 1) {
    n2 <- nrow(object$syn)
    res.ind <- propcalcs(object$syn, data)
    res <- list(call = match.call(), m = m, method = method, tree.method = tree.method,
                resamp.method = resamp.method, maxorder = maxorder, vars = vars,
                k.syn = k.syn, aggregate = aggregate, maxit = maxit,
                ngroups = ngroups, mincriterion = mincriterion,
                nperms = nperms, df = res.ind$df, incomplete = incomplete,
                pMSE = res.ind$pMSE, S_pMSE = res.ind$S_pMSE,
                S_SPECKS = res.ind$S_SPECKS, S_PO50 = res.ind$S_PO50,S_U = res.ind$S_U,
                SPECKS = res.ind$SPECKS, PO50 = res.ind$PO50, U = res.ind$U,
                print.stats = print.stats,
                fit = res.ind$fit, nnosplits = res.ind$nnosplits,
                digits = digits, print.ind.results = print.ind.results,
                print.zscores = print.zscores, zthresh = zthresh,
                print.variable.importance = print.variable.importance)
  } else {
    n2 <- nrow(object$syn[[1]])
    pMSE <- SPECKS <- PO50 <- U <- S_pMSE <- S_SPECKS <- S_PO50 <- S_U <- rep(NA, m)
    fit <- nnosplits <- as.list(1:m)
    if (!is.null(resamp.method) && !(resamp.method == "none") && resamp.method == "pairs") {
      kk <- 0
      simpMSE <- simKS <- simPO50 <- simU <- rep(NA, m*(m - 1)/2)
    }
    for (j in 1:m) {
      res.ind <- propcalcs(object$syn[[j]], data)
      pMSE[j] <- res.ind$pMSE
      SPECKS[j] <- res.ind$SPECKS
      PO50[j] <- res.ind$PO50
      U[j] <- res.ind$U
      fit[[j]] <- res.ind$fit

      if (resamp.method == "none" || (method == "logit" & (is.null(resamp.method)))) {
        if (j == 1 & print.flag) cat("Fitting syntheses: ")
        if (print.flag) {
          cat(j, " ", sep = "")
          if (res.ind$fit$converged == FALSE) cat("Convergence failed.\n")
        }
        if (j == m ) cat("\n")
        S_pMSE[j] <- res.ind$S_pMSE
      }
      
      if (!is.null(resamp.method) && resamp.method == "pairs") {
        if (j == 1 & print.flag) {
          if (print.every == 0 | m*(m - 1)/2 <= print.every) cat("Simulating NULL pMSE from ", m*(m - 1)/2, " pair(s).", sep = "")
          else cat("Simulating NULL pMSE from ", m*(m - 1)/2, " pairs, printing every ", print.every, "th:\n", sep = "")
          if (m*(m - 1)/2 < 6 ) cat("\nNumber of pairs too low, we suggest increasing number of syntheses (m).\n")
        }
        if (j < m) {
          for (jj in (j + 1):(m)) {
            kk <- kk + 1
            if (print.every > 0 & print.every < m*(m - 1)/2 & floor(kk/print.every) == kk/print.every & print.flag) cat(kk," ",sep = "")
            simvals <- propcalcs(object$syn[[j]], object$syn[[jj]])
            simpMSE[kk] <- simvals$pMSE
            simKS[kk] <- simvals$SPECKS
            simPO50[kk] <- simvals$SPECKS
            simU[kk] <- simvals$U
          }
        }
        nnosplits<- c(sum(simpMSE < 1e-8), length(simpMSE))
        for (j in 1:m) {
          S_pMSE[j] <- pMSE[j] *2 /mean(simpMSE)
          S_SPECKS[j] <- SPECKS[j] *2 /mean(simKS)
          S_PO50[j] <- PO50[j] *2 /mean(simPO50)
          S_U[j] <- U[j] *2 /mean(simU)
        }
        
      } else {
        nnosplits[[j]] <- res.ind$nnosplits  
        S_pMSE[j] <- res.ind$S_pMSE
        S_SPECKS[j] <- res.ind$S_SPECKS
        S_PO50[j] <- res.ind$S_PO50
        S_U[j] <- res.ind$S_U
      }
   }   
   res <- list(call = match.call(), m = m, method = method, tree.method = tree.method,
               resamp.method = resamp.method, maxorder = maxorder, vars = vars,
               k.syn = k.syn, aggregate = aggregate, maxit = maxit, 
               ngroups = ngroups, mincriterion = mincriterion,
               nperms = nperms, df = res.ind$df, incomplete = incomplete,
               pMSE = pMSE,  S_pMSE = S_pMSE,
               S_SPECKS = S_SPECKS, S_PO50 = S_PO50, S_U = S_U,
               SPECKS = SPECKS, PO50 = PO50, U = U,
               print.stats = print.stats,
               fit = fit, nnosplits = nnosplits,
               digits = digits, print.ind.results = print.ind.results,
               print.zscores = print.zscores, zthresh = zthresh,
               print.variable.importance = print.variable.importance)
 
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
                 ngroups = 5, useNA = TRUE, max.table = 1e6,
                 print.tables = length(vars) < 4,
                 print.stats = c("pMSE", "S_pMSE", "df"), 
                 print.zdiff = FALSE, print.flag = TRUE,
                 digits = 4, k.syn = FALSE, ...)
{

  if (is.null(data)) stop("Requires parameter 'data' to give name of the real data.\n",  call. = FALSE)
  if (is.null(object)) stop("Requires parameter 'object' to give name of the synthetic data.\n",  call. = FALSE)

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

  res <- utility.tab.synds(object = object, data = data, vars = vars, 
                           ngroups = ngroups, useNA = useNA, 
                           print.tables = print.tables,
                           print.stats = print.stats, 
                           print.zdiff = print.zdiff,
                           print.flag = print.flag,
                           digits = digits, k.syn = k.syn, ...)
  return(res)
}


###-----utility.tab--------------------------------------------------------
utility.tab.synds <- function(object, data, vars = NULL, ngroups = 5,
                              useNA = TRUE, max.table = 1e6, 
                              print.tables = length(vars) < 4,
                              print.stats = c("pMSE", "S_pMSE", "df"), 
                              print.zdiff = FALSE, print.flag = TRUE,
                              digits = 4, k.syn = FALSE, ...)
{
  vars <- unique(vars)

  # CHECKS
  #---------
  if (is.null(data))
    stop("Requires parameter 'data' to give name of the real data.\n", call. = FALSE)
  if (!is.data.frame(data))
    stop("Data must have class 'data.frame'.\n", call. = FALSE)
  if (!inherits(object, "synds"))
    stop("Object must have class 'synds'.\n", call. = FALSE)
  if (is.null(vars)) stop("You need to set variables with vars parameter.\n", call. = FALSE) else if
    (!(all(vars %in% names(data)))) stop("Unrecognized variable(s) in vars parameter: ",
     paste(vars[!(vars %in% names(data))], collapse = ", "), call. = FALSE)
  if (!all(print.stats %in% c("VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", "MabsDD", "dBhatt","S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE", "df", "dfG", "all")))
    stop('print.stats must be set to "all" or selected from "VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", "MabsDD", "dBhatt", "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE", "df" or "dfG".\n', call. = FALSE)
  #---------

  data <- data[, vars, drop = FALSE]
  nvars <- ncol(data)
  data.orig <- data

  # get cont.na parameters for stratified synthesis
  # --------
  if (!is.null(object$strata.syn)) {
    # cna <- apply(object$cont.na, 2, function(y) {unlist(unique(y))})
    cna <- object$cont.na[1, ]
  } else {
    cna <- object$cont.na
  }
  cna <- cna[vars]

  m <- object$m
  if (m == 1) syndata <- list(object$syn) else syndata <- object$syn
  syndata <- lapply(syndata, '[', vars)

  pMSE <- S_pMSE  <- df <- dfG <- VW <- S_VW  <- FT <- S_FT  <- G  <-  S_G  <- 
    JSD  <- U <- S_JSD <- MabsDD <- WMabsDD <- S_WMabsDD <- SPECKS <-  
    dBhatt <- PO50 <- nempty <- vector("numeric", m)
  tab.syn <- tab.obs <- tab.zdiff <- vector("list", m)
  
  syn.mvar <- vector("list", nvars)
  for (j in 1:nvars) {
    if (is.numeric(syndata[[1]][, j])) syn.mvar[[j]] <- c(sapply(syndata, '[[', j))
  }  
  
  for (i in 1:m) {
    data <- data.orig
    # make all variables into factors
    for (j in 1:nvars) {
      if (is.numeric(data[, j])) {
        grpd <- group_num(data[, j], syndata[[i]][, j], syn.mvar[[j]],
                          n = ngroups, cont.na = cna[[j]], ...)
        data[, j] <- grpd[[1]]; syndata[[i]][, j] <- grpd[[2]]
      } else if (is.character(data[, j])) {
        data[, j] <- factor(data[, j])
        syndata[[i]][, j] <- factor(syndata[[i]][, j],
                                   levels = levels(data[, j]))
      }
      if (any(is.na(data[, j])) & useNA) {
        # makes missings into part of factors if present
        data[, j] <- addNA(data[, j])
        syndata[[i]][, j] <- addNA(syndata[[i]][, j])
      }
    }

    ## check table size
    table.size <- prod(sapply(data, function(x) length(levels(x))))
    if (table.size > max.table)
    stop("Table size ", round(table.size), " exceeds max.table limit of ", round(max.table),".",
         "\nYou could try increasing max.table but memory problems are likely.\n", call. = FALSE)
    else if (i == 1 & table.size > dim(data)[1]/2 & print.flag) cat("Warning: You are creating tables with ", table.size, 
         " cells from ", dim(data)[1], " observations.\nResults from sparse tables may be unreliable.\n", sep = "")
    ## make tables
    if (useNA){
      tab.obs[[i]] <- table(data, useNA = "ifany", deparse.level = 0)
      tab.syn[[i]] <- table(syndata[[i]], useNA = "ifany", deparse.level = 0)
    } else {
      tab.obs[[i]] <- table(data, useNA = "no", deparse.level = 0)
      tab.syn[[i]] <- table(syndata[[i]], useNA = "no", deparse.level = 0)
    }

   ## remove cells all zeros
    nempty[i] <-   sum(tab.obs[[i]] + tab.syn[[i]] == 0)
    td <- tab.obs[[i]][tab.obs[[i]] + tab.syn[[i]]  > 0]
    ts <- tab.syn[[i]][tab.obs[[i]] + tab.syn[[i]]  > 0]
    totcells <- length(td)

   ## calculate utility measures
    if (!k.syn) df[i] <- totcells - 1 else df[i] <- totcells
    cc      <- sum(ts) / sum(ts + td)
    N       <- sum(ts + td)
    sumos   <- ts + td
    expect  <- sumos * cc
    diff    <- ts - td * cc / (1 - cc)
    VW[i]   <- sum(diff^2 / expect)
    FT[i]   <- 4*sum((ts^(0.5) - (cc / (1 - cc) * td)^(0.5))^2)
    S_FT[i] <- FT[i] / df[i]
    S_VW[i] <- S_pMSE[i] <- VW[i] / df[i]
    pMSE[i] <- VW[i] * cc * (1 - cc)^2 / N
   ## standardized difference (diff/sqrt(expect))
    tab.zdiff[[i]] <- suppressWarnings((tab.syn[[i]] - tab.obs[[i]] * cc/(1-cc)) / 
                                        sqrt((tab.syn[[i]] + tab.obs[[i]]) * cc))
   ## Jensen-Shannon divergence 
    ptabd    <- td / sum(td)
    ptabs    <- ts / sum(ts)
    phalf    <- (ptabd + ptabs) *0.5
    JSD[i]   <- sum((ptabd * log2(ptabd/phalf))[ptabd > 0])/2 +        
                sum((ptabs * log2(ptabs/phalf))[ptabs > 0])/2
    S_JSD[i] <- JSD[i]*2*N/df[i]/log(2)
   ## Symmetric likelihood ratio chisq
    sok     <- ts[ts > 1e-8 & td > 1e-8]                              
    dok     <- td[ts > 1e-8 & td > 1e-8]
    if (!k.syn) dfG[i] <- length(dok) - 1 else dfG[i] <- length(dok)
    G[i]    <-  2 *sum(sok*log(sok/sum(sok)/dok*sum(dok)))
    S_G[i]  <- G[i] / dfG[i]
   ## Kolmogorov-Smirnov
    score   <- ts / (ts + td)                                                 
    kst     <- suppressWarnings(ks.test(rep(score, ts), rep(score, td)))
    SPECKS[i] <- kst$statistic
   ## Wilcoxon statistic 
    Ut      <- suppressWarnings(wilcox.test(rep(score, ts), rep(score, td))) 
    U[i]    <- Ut$statistic
   ## Calculate PO50
    predsyn <- (ptabs > ptabd)                                        
    PO50[i] <- (sum(ts[predsyn]) + sum(td[!predsyn])) / (sum(ts) + sum(td)) * 100 - 50
    
    MabsDD[i]  <- sum(abs(diff))/sum(ts)
    WMabsDD[i] <- sum(abs(diff)/sqrt(expect))/sqrt(2/pi)
    S_WMabsDD[i] <- WMabsDD[i]/df[i]
    
    dBhatt[i] <- sqrt(1 - sum(sqrt(ptabd*ptabs)))
  }

  tab.obs <- tab.obs[[1]]  

  if (m == 1) {
    tab.syn <- tab.syn[[1]]
    tab.zdiff <- tab.zdiff[[1]]
  }
  
  res <- list(m = m,
              VW = VW, 
              FT = FT,
              JSD = JSD, 
              SPECKS = SPECKS,
              WMabsDD = WMabsDD,
              U = U,
              G = G,
              pMSE = pMSE, 
              PO50 = PO50,
              MabsDD = MabsDD,
              dBhatt = dBhatt, 
              S_VW = S_VW,
              S_FT = S_FT,
              S_JSD = S_JSD,
              S_WMabsDD = S_WMabsDD,
              S_G = S_G,
              S_pMSE = S_pMSE,
              df = df,
              dfG = dfG,
              nempty = unlist(nempty),
              tab.obs = tab.obs,
              tab.syn = tab.syn,
              tab.zdiff = tab.zdiff,
              digits = digits,
              print.stats = print.stats,
              print.zdiff = print.zdiff,
              print.tables = print.tables,
              n = sum(object$n),
              k.syn = k.syn)

  class(res) <- "utility.tab"
  return(res)
}


###-----group_num----------------------------------------------------------
# function to categorise continuous variables

group_num <- function(x1, x2, xsyn, n = 5, style = "quantile", cont.na = NA, ...) {

  # Categorise 2 continuous variables into factors of n groups
  # with same groupings determined by the first one
  # xsyn - all synthetic values (for m syntheses)

  if (!is.numeric(x1) | !is.numeric(x2) | !is.numeric(xsyn)) 
    stop("x1, x2, and xsyn must be numeric.\n", call. = FALSE)

  # Select non-missing(nm) values
  x1nm <- x1[!(x1 %in% cont.na) & !is.na(x1)]
  x2nm <- x2[!(x2 %in% cont.na) & !is.na(x2)]
  xsynnm <- xsyn[!(xsyn %in% cont.na) & !is.na(xsyn)]
  
  # Derive breaks
  my_breaks <- unique(suppressWarnings(classIntervals(c(x1nm, xsynnm),
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

