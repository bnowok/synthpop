# Strata can be provided as a vector (=flag) or as variable names 
# If all stratifying values have the same values within strata they will be
# automatically removed; if not they will be taken into account

# ????? 
#! probably NOT possible and syn.strata needs to be kept as a separate FUN
# creating generic syn function
# changing syn to syn.default
# and if strata != NULL calling syn.strata
# ?????


# To check:
# 1. Defaults agree with syn() ?

# To do:
# 1. NAs in a stratifying variable/indicator: combine with another stratum? 

# Other issues:
# - multiple printout (per strata, per m)
# - print strata names?
# - compare() by strata?

# 12/08/2016:error when zero number of observation for strata in synthetic data

syn.strata <- function(data, strata = NULL, 
                minstratumsize = 10 + 10 * length(visit.sequence), 
                tab.strataobs = TRUE, tab.stratasyn = FALSE,
                method = vector("character", length = ncol(data)),
                visit.sequence = (1:ncol(data)),
                predictor.matrix = NULL,
                m = 1, k = nrow(data), proper = FALSE,
                minnumlevels = 5, maxfaclevels = 60,
                rules = NULL, rvalues = NULL,
                cont.na = NULL, semicont = NULL,
                smoothing = NULL, event = NULL, denom = NULL,
                drop.not.used = FALSE, drop.pred.only = FALSE,
                default.method = c("normrank","logreg","polyreg","polr"),
                models = FALSE,
                print.flag = TRUE,
                seed = "sample",
                ...){
 
 m0 <- max(1, m)

 # CHECKS
 #--------
 if (is.null(strata)) stop("Argument strata is missing.", call. = FALSE) 
 # If strata given as variable names (check if they exist) 
 # -> change into one factor with strata names 
 if (is.character(strata) & any(!duplicated(strata))) {
   varindex <- match(strata, colnames(data))
   if (any(is.na(varindex))) stop("Unrecognized variable(s) in strata: ", 
     paste(strata[is.na(varindex)],collapse=", "), call. = FALSE)
   else {
     dstrataNA  <- lapply(data[, strata, drop = FALSE], addNA, ifany = TRUE) # change NA to a factor level
     strata.lab <- interaction(dstrataNA, drop = TRUE, sep = "_")
     #strata.varnames <- paste0(strata, collapse="_") 
   }
 } else {
 # check length of strata vector if given as vector; check for missing
   if (length(strata) != nrow(data)) stop(paste("The length of strata index (",
     length(strata), ") does not match the number of rows in the data (",
     nrow(data),").",sep=""), call. = FALSE)
   if (any(is.na(strata))) stop("Strata indicator cannot have missing values.", 
     call. = FALSE)
   strata.lab <- factor(strata)
 }
 #--------
 
 # make sure stratification variables are include in visit.sequence
 # important when drop.not.used==T
 if (is.character(strata) & any(!duplicated(strata))){   #GR-20/10/2016 drop.not.used == TRUE removed from the condition
   strata <- match(strata, colnames(data))
   if (is.character(visit.sequence)) visit.sequence <- match(visit.sequence, colnames(data))
   if (any(is.na(visit.sequence))) stop("Unrecognized variable(s) in visit.sequence.", call. = FALSE)
   visit.sequence <- c(visit.sequence, strata[!(strata %in% visit.sequence)])
 }
 
 stratalev.lab <- levels(strata.lab) 
 strata.n.obs  <- table(strata.lab)   
 nstrata       <- dim(strata.n.obs)

 if (tab.strataobs == TRUE) {
   cat("Number of observations in strata (original data):")
   print(table(strata.lab, deparse.level = 0))
 }

 # check min number of observations in strata
 stratasize.stop <- minstratumsize         
 stratasize.warn <- 100 + 10*length(visit.sequence)

 smallstrata <- sum(strata.n.obs < stratasize.stop)
 if (smallstrata > 5) stop("In the original data multiple strata do not have enough observations.\nEach should have at least ",
   stratasize.stop, " observations (minstratumsize which by default is 10 + 10 * no. of variables used in prediction).\n", 
   sep="", call. = FALSE)
 if (smallstrata > 0) stop("In the original data some strata (", 
   paste(stratalev.lab[strata.n.obs < stratasize.stop], collapse=", "), 
   ") do not have enough observations.\nEach should have at least ",
   stratasize.stop, " observations (minstratumsize which by default is 10 + 10 * no. of variables used in prediction).\n", 
   sep="", call. = FALSE)
 if (any(strata.n.obs < stratasize.warn) & print.flag == TRUE) {
   cat("CAUTION: In the original data some strata (", 
   paste(stratalev.lab[strata.n.obs < stratasize.warn], collapse=", "), 
   ") have limited number of observations.\nThere should be at least ", 
   stratasize.warn, 
   " observations (100 + 10 * no. of variables used in prediction).\n", sep="")
 }
 
 synds.names <- c("call", "m", "syn", "method", "visit.sequence", 
   "predictor.matrix", "smoothing", "event", "denom", "proper", "n", "k", 
   "rules", "rvalues", "cont.na", "semicont", "drop.not.used", "drop.pred.only", 
   "models", "seed", "var.lab", "val.lab", "obs.vars", "strata.syn", "strata.lab")
 synds <- list(setNames(vector("list",length(synds.names)),synds.names)) 
 synds <- rep(synds, m0)
 sel.names <- match(c("call", "m", "predictor.matrix", "proper", "strata.syn",
   "strata.lab", "models"), synds.names)
 same.by.m <- c("call", "m", "method", "visit.sequence", "predictor.matrix", 
   "smoothing", "event", "denom", "proper", "n", "rules", "rvalues", "cont.na", 
   "semicont", "drop.not.used", "drop.pred.only", "models", "seed", "var.lab", 
   "val.lab", "obs.vars", "strata.lab")
 same.by.m.idx <- match(same.by.m, synds.names) 
  
 syn.args <- as.list(match.call()[-1])
   
 strata.n.syn <- vector("list", m0) 
 for (j in 1:m0){ 
   synds[[j]]$strata.syn <- sort(factor(sample(stratalev.lab, k, replace = TRUE, 
     prob = strata.n.obs), levels = stratalev.lab))   
   synds[[j]]$strata.lab <- stratalev.lab
   strata.n.syn[[j]] <- table(synds[[j]]$strata.syn, deparse.level = 0)
   if (tab.stratasyn == TRUE) {  
     cat("\nNumber of observations in strata (synthetic data, m = ", j, "):", sep="")
     print(strata.n.syn[[j]])
   }
 }
 #Different way of printing (all syn in one table)  
 #cat("\nNumber of observations in strata (synthetic data):\n")
 #starta.n.syn.df <- do.call("rbind", strata.n.syn) 
 #rownames(starta.n.syn.df) <- paste0("m = ", 1:m)
 #print(starta.n.syn.df)

 for (j in 1:m0){  
   synds.ind <- vector("list", nstrata) # results by stratum  
   # exclude args that are not in syn()
   syn.args$strata <- syn.args$tab.stratasyn <- syn.args$tab.strataobs <- 
     syn.args$minstratumsize <- NULL
   syn.args$m <- 1
   syn.args$visit.sequence <- visit.sequence
   
   for (i in 1:nstrata) {
     if (print.flag) cat("\nm = ",j,", strata = ", stratalev.lab[i],
       "\n-----------------------------------------------------\n", sep="")
     syn.args$data   <- data[strata.lab == stratalev.lab[i],]
     syn.args$k      <- strata.n.syn[[j]][i]; names(syn.args$k) <- "strata"
     if (syn.args$k == 0 | m==0) syn.args$m <- 0 else syn.args$m <- 1
     synds.ind[[i]]  <- do.call("syn", syn.args)
   }
   synds[[j]]$call <- match.call()
   synds[[j]]$m <- m
   synds[[j]]$proper <- proper
   synds[[j]]$predictor.matrix <- eval(parse(text = paste0("list(", 
     paste0("synds.ind[[", 1:nstrata, "]]$predictor.matrix", collapse = ", "),")")))
   synds[[j]]$models <- eval(parse(text = paste0("list(", 
     paste0("synds.ind[[", 1:nstrata, "]]$models", collapse = ", "),")")))
   synds[[j]][-sel.names] <- eval(parse(text = paste0("Map(rbind, ", 
     paste0("synds.ind[[", 1:nstrata, "]][-sel.names]", collapse = ", "),")")))
   if (k > 0 & m > 0) rownames(synds[[j]]$syn) <- 1:k
 }
 
 if (m==1 | m==0) {
   synds <- synds[[1]] 
 } else {
   synds <- eval(parse(text = paste0("Map(list, ", paste0("synds[[", 1:m, "]]", 
     collapse = ", "),")")))
   synds[same.by.m.idx] <- lapply(same.by.m.idx, function(x) synds[[x]][[1]])
 }

 if (m==0) cat("\nCAUTION: method, visit.sequence and predictor.matrix are lists that may vary by stratum and\nshould not be used to initialise syn. For initialising rerun with m = 0 without stratification.\n\n") 
 
 class(synds) <- "synds"
 return(synds) 
}
