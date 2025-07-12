###-----compare------------------------------------------------------------
compare <- function(object, data, ...) UseMethod("compare")

###-----compare.default----------------------------------------------------
compare.default <- function(object, ...)
  stop("No compare method associated with class ", class(object), call. = FALSE)


###-----compare.synds------------------------------------------------------
compare.synds <- function(object, data, vars = NULL, msel = NULL,  
                          stat = "percents", breaks = 20, ngroups =5, 
                          nrow = 2, ncol = 2, rel.size.x = 1,
                          utility.stats = c("pMSE", "S_pMSE", "df"),
                          utility.for.plot = "S_pMSE",  
                          cols = c("#1A3C5A","#4187BF"),
                          plot = TRUE, table = FALSE, print.flag = TRUE, ...){   

 if (is.null(data)) stop("Requires parameter data to give name of the real data.\n", call. = FALSE)
 if (!is.data.frame(data)) stop("Argument data must be a data frame.\n", call. = FALSE) 
  
 if (any(class(data) %in% c("tbl", "tbl_df"))) data <- as.data.frame(data)  
  
 if (!inherits(object, "synds")) stop("Object must have class synds.\n", call. = FALSE )                                                                    
 if (!is.null(msel) & !all(msel %in% (1:object$m))) stop("Invalid synthesis number(s).", call. = FALSE)
 if (!all(utility.stats %in% c("VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", "MabsDD", "dBhatt",
                              "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE", "df", "all"))) 
  stop('utility.stats must be set to "all" or selected from "VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", "MabsDD", "dBhatt", "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE" or "df".\n', call. = FALSE)  
 
  if (!is.null(utility.for.plot) &&
    !(length(utility.for.plot) == 1 & utility.for.plot %in% c("VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", "MabsDD", "dBhatt",
                                "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE"))) 
    stop('utility.for.plot must be one of "VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", "MabsDD", "dBhatt", "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G" or "S_pMSE" .\n', call. = FALSE)  
  if (!is.null(utility.for.plot) && !utility.for.plot %in% utility.stats ) cat(utility.for.plot, "used on plots added to table of results.\n")   
  utility.stats <- unique(c(utility.stats, utility.for.plot))
   
 if (!(length(stat) == 1 & stat %in% c("percents", "counts"))) { 
   cat('Parameter stat must be "percents" or "counts".\n' )
   stat <- "percents"
   cat('Changed to default value of "percents".\n')
 }

 # single / pooled synthetic data sets                                  
 if (object$m == 1) {
   syndsall <- object$syn 
 } else if (length(msel) == 1) {
   syndsall <- object$syn[[msel]]
 } else if (length(msel) > 1 | is.null(msel)) {                               
   syndsall <- do.call(rbind,object$syn)
 }
 # list of synthetic data sets for non-pooled results 
 if (length(msel) > 1) {
   synds <- vector("list",length(msel))
   for (i in 1:length(msel)) synds[[i]] <- object$syn[[msel[i]]]
 }
 synnames    <- names(syndsall)
 realnames   <- names(data)
 commonnames <- synnames[match(realnames,synnames)]
 commonnames <- commonnames[!is.na(commonnames)]

 if (!is.null(vars)) {
   if (!(all(vars %in% synnames))) stop("Variable(s) ", 
     paste0(vars[is.na(match(vars,synnames))], collapse = ", "),
     " not in synthetic data \n", call. = FALSE)
   if (!(all(vars %in% realnames))) stop("Variable(s) ", 
     paste0(vars[is.na(match(vars,realnames))], collapse = ", "),
     " not in observed data \n", call. = FALSE)
   commonnames <- commonnames[match(vars,commonnames)]
 }
 
 if (!(all(synnames %in% realnames))) cat("Warning: Variable(s)", 
   paste0(synnames[is.na(match(synnames,realnames))], collapse = ", "),
   "in synthetic object but not in observed data\n",
   " Looks like you might not have the correct data for comparison\n")
 
 if ((length(commonnames) == 0) && (typeof(commonnames) == "character"))        #! when would it apply?
   stop("None of variables selected for comparison in data", call. = FALSE)

 df.obs    <- data[, commonnames, drop = FALSE]
 df.synall <- syndsall[, commonnames, drop = FALSE]
 if (length(msel) > 1) {
   df.syn <- vector("list",length(msel))
   for (i in 1:length(msel)) df.syn[[i]] <- synds[[i]][, commonnames, drop = FALSE]
 }

 # change any numeric variables with < 6 distinct values to factors
 for (i in 1:length(commonnames)) {
   if (is.numeric(df.obs[,i]) && length(table(df.obs[,i])) < 6){
     df.obs[,i] <- as.factor(df.obs[,i])
     df.synall[,i] <- as.factor(df.synall[,i])
   }
 }

 num <- sapply(df.synall, is.numeric) | sapply(df.synall, is.integer)  
 fac <- sapply(df.synall, function(x) is.factor(x) | is.logical(x))  

 if (is.null(vars)) {
   if (object$m == 1) vars <- names(object$syn)
   else vars <- names(object$syn[[1]])
 }

 if (!is.null(utility.stats) | !is.null(utility.for.plot)) {
   utility.list <- as.list(1:length(vars))
   names(utility.list) <- vars
   if (utility.stats[1] == "all") utility.stats <- 
     c("VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", 
       "MabsDD", "dBhatt","S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", 
       "S_pMSE", "df")
   tab.utility <- matrix(NA, length(vars), length(unique(c(utility.stats, utility.for.plot))))
   dimnames(tab.utility) <- list(vars, unique(c(utility.stats, utility.for.plot)))

   for (i in 1:length(vars)) {
     utility.list[[i]] <- utility.tab(object, data, vars = vars[i], ngroups = ngroups)
     if (i == 1) tab.ind <- match(utility.stats, names(utility.list[[i]]))
     tab.utility[i, ] <- sapply(utility.list[[i]][tab.ind], mean)
     if (print.flag) cat("Calculations done for", vars[i],"\n")
   }

   if (!is.null(utility.for.plot)) utilvals.for.plot <- tab.utility[, match(utility.for.plot, dimnames(tab.utility)[[2]])]
 } else {
   if (is.null(utility.stats)) tab.utility <- NULL
   if (is.null(utility.for.plot)) utilvals.for.plot <- NULL
 }
 # frequency tables for factors
 if (sum(fac) > 0) {
   any.fac.na <- unlist(apply(df.obs[,fac,drop = FALSE],2,function(x) any(is.na(x))))    
   per.obs.fac <- ggfac(df.obs[, fac, drop = FALSE], anyna = any.fac.na, stat = stat)  
   if (length(msel) <= 1) {
     per.syn.facall <- ggfac(df.synall[, fac, drop = FALSE], 
     name = "synthetic", anyna = any.fac.na, stat = stat)
     if (stat == "counts") {
       per.syn.facall$perdf$Count <- per.syn.facall$perdf$Count/object$m
       per.syn.facall$perlist <- lapply(per.syn.facall$perlist,"/",object$m)
     }
   }
   if (length(msel) > 1) {
     per.syn.fac <- vector("list",length(msel))
     for (i in 1:length(msel)) per.syn.fac[[i]] <- ggfac(df.syn[[i]][, fac, drop = FALSE], 
       name = paste0("syn=", msel[i]), anyna = any.fac.na, stat = stat)  
   }
 } else {                                       
   per.obs.fac    <- NULL
   per.syn.facall <- NULL
 }

 # frequency tables for numeric variables
 if (sum(num) > 0) {
   cont.index <- match(colnames(df.obs[, num, drop = FALSE]), colnames(syndsall))    
   na <- object$cont.na[cont.index] 
  # to exclude from summaries if no missing in data
   any.na <- unlist(apply(df.obs[, num, drop = FALSE], 2, function(x) any(is.na(x))))  

   lbreaks    <- as.list(rep(breaks, length(na)))  

  ## get limits(=breaks) from both observed and synthetic
  #--
   df.both  <- rbind.data.frame(df.obs, df.synall)                
   per.both <- ggnum(df.both[, num, drop = FALSE], na = na, 
     breaks = lbreaks, anyna = any.na, stat = stat) ##GR stat added 
  #--  
   per.obs.num <- ggnum(df.obs[, num, drop = FALSE], na = na, 
     breaks = per.both$hbreaks, anyna = any.na, stat = stat) ##GR stat added
   
   if (length(msel) <= 1) {
     per.syn.numall <- ggnum(df.synall[, num, drop = FALSE], 
       name = "synthetic", na = na, breaks = per.both$hbreaks, 
       anyna = any.na, stat = stat) 
       if (stat == "counts") {
         per.syn.numall$perdf$Count <- per.syn.numall$perdf$Count/object$m
         per.syn.numall$perlist <- lapply(per.syn.numall$perlist,"/",object$m)
       }
   } 
   if (length(msel) > 1) {
     per.syn.num <- vector("list",length(msel))
     for (i in 1:length(msel)) per.syn.num[[i]] <- ggnum(df.syn[[i]][, num, drop = FALSE], 
       name =  paste0("syn=",msel[i]), na = na, breaks = per.both$hbreaks, 
       anyna = any.na, stat = stat) 
   } 
 } else {
   per.obs.num <- NULL
   per.syn.numall <- NULL
 }
 # data frame for plotting 
 if (length(msel) <= 1) per.fac <- rbind.data.frame(per.obs.fac$perdf, 
   per.obs.num$perdf, per.syn.facall$perdf, per.syn.numall$perdf )

 if (length(msel) > 1) {
   per.fac <- rbind.data.frame(per.obs.fac$perdf, per.obs.num$perdf)
   for (i in 1:length(msel)) {
     if (sum(fac) > 0) temp.fac <- per.syn.fac[[i]]$perdf else temp.fac <- NULL
     if (sum(num) > 0) temp.num <- per.syn.num[[i]]$perdf else temp.num <- NULL  
     per.fac <- rbind.data.frame(per.fac, temp.fac, temp.num)
   }
 }
  per.fac$Variable <- factor(per.fac$Variable, levels = commonnames, 
   ordered = T, exclude = NULL)
  per.fac$Value <- factor(per.fac$Value, levels = unique(per.fac$Value), 
   ordered = T, exclude = NULL)
 
 # list of result tables
 if (length(msel) <= 1) {
   os.table.fac <- mapply(rbind, obs = per.obs.fac$perlist, 
     syn = per.syn.facall$perlist, SIMPLIFY = FALSE)  
   os.table.num <- mapply(rbind, obs = per.obs.num$perlist, 
     syn = per.syn.numall$perlist, SIMPLIFY = FALSE)  
 }
 if (length(msel) > 1) {
   os.table.fac <- per.obs.fac$perlist 
   os.table.num <- per.obs.num$perlist 
   for (i in 1:length(msel)) {
     if (sum(fac) > 0) temp.fac <- per.syn.fac[[i]]$perlist else temp.fac <- NULL
     if (sum(num) > 0) temp.num <- per.syn.num[[i]]$perlist else temp.num <- NULL
     os.table.fac <- mapply(rbind, os.table.fac, temp.fac, SIMPLIFY = FALSE)    
     os.table.num <- mapply(rbind, os.table.num, temp.num, SIMPLIFY = FALSE)     
   }
  }
 os.table <- c(os.table.fac, os.table.num)
 if (is.null(msel)) {
   for (i in 1:length(os.table)) dimnames(os.table[[i]])[[1]] <- c("observed","synthetic")
 } else {
   for (i in 1:length(os.table)) dimnames(os.table[[i]])[[1]] <- c("observed", paste0("syn=", msel))
 }
 Value <- Percent <- Count <- Data <- NULL    ## otherwise 'no visible binding for global variables'
 # sorts the factor labels in the right order for numeric vars
 per.fac$Value <- as.character(per.fac$Value)
 vals <- unique(per.fac$Value)
 valsnum <- unique(per.fac$Value[per.fac$Variable %in% names(num[num == TRUE])])
 valsnum.nonmiss <- sort(as.numeric(vals[vals %in% valsnum & substr(vals, 1, 4) != "miss"]))
 valsnum.nonmiss <- format(valsnum.nonmiss, scientific = FALSE,
                           trim = TRUE, drop0trailing = TRUE)
 valsnum.miss <- sort(vals[vals %in% valsnum & substr(vals, 1, 4) == "miss"])
 vals[vals %in% valsnum] <- c(valsnum.nonmiss, valsnum.miss)
 per.fac$Value <- factor(as.character(per.fac$Value), levels = vals)
 
 if (!is.null(utility.for.plot)){
   levels(per.fac$Variable) <- paste0(levels(per.fac$Variable), ": ", 
     utility.for.plot, " = ", round(utilvals.for.plot, 2))
   commonnames_lab <- paste0(commonnames, ": ",
     utility.for.plot, " = ", round(utilvals.for.plot, 2))
 } else {
   commonnames_lab <- commonnames
 }

 # get different plots in order of data
 nperplot <- nrow*ncol
 nplots   <- ceiling(length(commonnames)/nperplot)
 plots    <- vector("list", nplots)
 tables   <- vector("list", nplots)
 
 for (i in 1:nplots) {
   min <- (i - 1)*nperplot + 1
   max <- min(length(commonnames), (i - 1)*nperplot + nperplot)
   # tables
   ttables <- vector("list", (max - min + 1))
   names(ttables) <- commonnames[min:max]
   for (j in commonnames[min:max]) {
     ttables[[j]] <- os.table[[j]]
   } 
   tables[[i]] <- ttables

  # plots
   per.fact <- per.fac[per.fac$Variable %in% commonnames_lab[min:max],]
   if (stat == "percents") p <- ggplot(data = per.fact, aes(x = Value, y = Percent, fill = Data))
   else p <- ggplot(data = per.fact, aes(x = Value, y = Count, fill = Data))
   p <- p + geom_bar(position = "dodge", colour = cols[1], stat = "identity") + 
      facet_wrap(~ Variable, scales = "free", ncol = ncol)  
   p <- p + guides(fill = guide_legend(override.aes = list(colour = NULL))) + 
        theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1, size = rel(rel.size.x)), 
              legend.position = "top", 
              legend.key = element_rect(colour = cols[1])) 
   p <- p + theme(legend.title = element_blank())
   if (length(msel) > 1) p <- p + scale_fill_manual(values = c(cols[1], rep(cols[2], length(msel))))
   if (length(msel) <= 1) p <- p + scale_fill_manual(values = cols)
   plots[[i]] <- p
 }
 
 if (length(tables) == 1) {
   tables <- tables[[1]]
   plots  <- plots[[1]]
 }

 res <- list(tables = tables, plots = plots, stat = stat, vars = vars,
             tab.utility = tab.utility, table = table, plot = plot)

 class(res) <- "compare.synds"
 return(res)
}


###-----compare.data.frame---compare.list----------------------------------    
compare.data.frame <- compare.list <- function(object, data, vars = NULL, cont.na = NULL,         
                                   msel = NULL, stat = "percents", breaks = 20, ngroups = 5, 
                                   nrow = 2, ncol = 2, rel.size.x = 1,
                                   utility.stats = c("pMSE", "S_pMSE", "df"),
                                   utility.for.plot = "S_pMSE",
                                   cols = c("#1A3C5A","#4187BF"),   
                                   plot = TRUE, table = FALSE ,print.flag = TRUE,
                                   compare.synorig = TRUE, ...){
  
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
  
  if (compare.synorig){
    if (m ==1) adjust.data <- synorig.compare(object,data, print.flag = FALSE) else
      if (m > 1) adjust.data <- synorig.compare(object[[1]],data, print.flag = FALSE)
    
    if (adjust.data$needsfix) stop("Synthetic data and/or original data needs more fixing before you can
      run the disclosure functions - see output. Use function synorig,compare() to check.", call. = FALSE)
    else if (!adjust.data$unchanged) {
      object <- adjust.data$syn
      data <- adjust.data$orig
      cat("Synthetic data or original or both adjusted with synorig.compare to try to make them comparable.\n")
      if (m > 1) cat("only first element of the list has been adjusted and will be used here\n")
      m <- 1 }
  }
  
  object <- list(syn = object, m = m, cont.na = cont.na)
  class(object ) <- "synds"
  
  res <- compare.synds(object = object, data = data, vars = vars, 
                       msel = msel, stat = stat, breaks = breaks, ngroups = ngroups,
                       nrow = nrow, ncol = ncol, rel.size.x = rel.size.x,
                       utility.stats = utility.stats,
                       utility.for.plot = utility.for.plot,
                       cols = cols, plot = plot, table = table,print.flag = print.flag, ...) 
  res$call <- match.call()
  return(res)
}


###-----pertable-----------------------------------------------------------
# calculate counts/percentages
pertable <- function(x, stat, ...) { ##GR stat added
 res <- table(x, useNA = "always")   
 if (stat == "percents") res <- res*100/sum(res) 
 return(res)
}


###-----ggfac--------------------------------------------------------------
# calculate percentages for factors and store in a data frame (a long format)
ggfac <- function(data, name = "observed", anyna, stat){ 
  data <- as.data.frame(data)
  perlist  <- lapply(data, pertable, stat = stat)
  for (i in 1:length(perlist)) {
    if (anyna[i] == FALSE) perlist[[i]] <- perlist[[i]][-length(perlist[[i]])]  
  }
  if (stat == "percents") Percent  <- unlist(perlist, use.names = FALSE) 
  else Count <- unlist(perlist, use.names = FALSE)
  Value <- unlist(lapply(perlist,names), use.names = FALSE)
  Variable <- rep(names(perlist), sapply(perlist,length))
  if (stat == "percents") perdf <- cbind.data.frame(Percent, Value, Variable) 
  else  perdf <- cbind.data.frame(Count, Value, Variable)
  
  perdf$Data <- name
  return(list(perdf = perdf, perlist = perlist))
} 

###-----ggnum--------------------------------------------------------------
# calculate percentages for numeric variables (non-missing values and 
# missing data categories seperately) and store in a data frame (a long format)
ggnum <- function(data, name = "observed", na = as.list(rep(NA,ncol(data))), 
                  breaks = as.list(rep(30,ncol(data))), anyna, 
                  stat = stat){ 
  data <- as.data.frame(data)

# non-missing values  
  nvar <- ncol(data)
  perlist <- vector("list", nvar)
  hbreaks <- vector("list", nvar)
  
  for (i in 1:nvar) {
  ## counts for non-missing values  
    vardata <- data[!(data[,i] %in% na[[i]]),i]                  
    hh      <- hist(vardata, breaks = breaks[[i]], plot = FALSE)
    counts  <- hh$counts
    names(counts) <- format(hh$breaks[-length(hh$breaks)],
                            scientific = FALSE, trim = TRUE,
                            drop0trailing = TRUE)
    hbreaks[[i]]  <- hh$breaks
  ## counts for missing values  
    dataNA <- data[data[,i] %in% c(NA,na[[i]]),i]
    
    if (length(dataNA) == 0) {
      if (anyna[i] == TRUE) {
        NAcounts <- 0   
        names(NAcounts) <- NA 
      } else {
        NAcounts <- NULL
      }
    } else {
      if (anyna[i] == TRUE) {
        NAcounts <- table(dataNA, useNA = "always")
      } else {
        NAcounts <- table(dataNA)  
      }
    } 
    if (!is.null(NAcounts)) names(NAcounts) <- paste("miss", names(NAcounts), sep = ".") 
    counts <- c(counts, NAcounts)
    if (stat == "percents") perlist[[i]] <- counts*100/length(data[,i])
    else perlist[[i]] <- counts
  }
  names(perlist) <- colnames(data)

# create data frame in a long format  
  if (stat == "percents") Percent  <- unlist(perlist, use.names = F) 
  else Count  <- unlist(perlist, use.names = F)
  Value    <- unlist(lapply(perlist,names), use.names = F)
  Variable <- rep(names(perlist),sapply(perlist,length))
  if (stat == "percents") perdf    <- cbind.data.frame(Percent, Value, Variable)
  else perdf    <- cbind.data.frame(Count, Value, Variable)
  perdf$Data <- name
  
return(list(perdf = perdf, perlist = perlist, hbreaks = hbreaks))
}


###-----dfNA---------------------------------------------------------------
dfNA <- function(data, na){
  all  <- length(data)
  data <- data[data %in% unlist(na)]
  NA.counts <- table(data, exclude = NULL)*100/all
  return(NA.counts)
}


###-----compare.fit.synds--------------------------------------------------
compare.fit.synds <- function(object, data, plot = "Z",
  print.coef = FALSE, return.plot = TRUE, plot.intercept = FALSE, 
  lwd = 1, lty = 1, lcol = c("#1A3C5A","#4187BF"),
  dodge.height = .5, point.size = 2.5, 
  population.inference = FALSE, ci.level = 0.95, ...) {   # c("#132B43", "#56B1F7")

 # Compares and plots fits to synthetic and original data
 # First parameter must be a fit to synthetic data from glm.synds(),
 # lm.synds(), polr.synds() or multinom.synds()

 value <- "Value"
 coefficient <- c("Coefficient", "Model")  
  
 if (!inherits(object, "fit.synds")) stop("Object must have class fit.synds.\n")
 if (!is.data.frame(data)) stop("Data must be a data frame.\n")  # theoretically can be a matrix (?)
 if (ci.level <= 0 | ci.level > 1) stop("ci.level must be beteen 0 and 1.\n")

 m <- object$m
 n <- object$n
 k <- object$k
 fitting.function <- object$fitting.function

##?? syn.coef    <- object$mcoefavg
 
 call <- match.call()

 # get fit to real data
 if (fitting.function %in% c("multinom", "polr")) {
   real.fit.0 <- do.call(object$fitting.function,
                         args = list(formula = formula(object),
                         Hess = TRUE, data = call$data))
   real.fit <- summary(real.fit.0)
 } else if (fitting.function %in% c("lm")) {
   real.fit <- summary(do.call(object$fitting.function,
                       args = list(formula = formula(object), 
                       data = call$data)))
 } else { ##  for glm
   real.fit <- summary(do.call(object$fitting.function,
                       args = list(formula = formula(object),
                       family = object$call$family, 
                       data = call$data)))
 }
 
 if (object$fitting.function == "multinom") { 
   real.varcov <- vcov(real.fit.0) 
   dd <- dimnames(t(real.fit$coefficients))
   real.fit$coefficients <- cbind(as.vector(t(real.fit$coefficients)), 
                                  as.vector(t(real.fit$standard.errors)),
                                  as.vector(t(real.fit$coefficients))/as.vector(t(real.fit$standard.errors))) 
   dimnames(real.fit$coefficients) <- list(
                                   paste(rep(dd[[2]], each = length(dd[[1]])), 
                                 rep(dd[[1]], length(dd[[2]])), sep = " : "), c("Estimate","Std error","t stat"))
 } else if (object$fitting.function == "polr") { 
   real.varcov <- vcov(real.fit.0) 
 } else if (object$fitting.function == "lm") {
   real.varcov <- real.fit$cov.unscaled * real.fit$sigma^2  ##??
 } else { # for "glm"
   real.varcov <- real.fit$cov.scaled   
 }

 syn.fit  <- summary.fit.synds(object, real.varcov = real.varcov,
                               population.inference = population.inference)
 incomplete <- syn.fit$incomplete

 # detailed results
 res.obs <- real.fit$coefficients[,1:3]
 colnames(res.obs) <- c("Beta","se(Beta)","Z")
 res.syn  <- syn.fit$coefficients[,1:3] 
 res.syn  <- res.syn[order(match(rownames(res.syn), rownames(res.obs))), ]
 res.overlap <- overlap.CI(res.syn, res.obs, ci.level = ci.level, intercept = TRUE)
 ncoef <- nrow(res.obs) 

 res.diff <-  cbind(res.syn[,1], res.obs[,1], 
		    res.syn[,1] - res.obs[,1], 
 		   (res.syn[,1] - res.obs[,1]) / res.obs[,2])
 dimnames(res.diff)[[2]] <- c("Synthetic", "Observed", "Diff", "Std. coef diff")
#? pval <- round(2 * (1 - pnorm(abs(res.syn[,1] - res.obs[,1])/sqrt(diag(lof.varcov)))), 3)  # "p value"
 
 if (incomplete == TRUE) {   
   if (object$m < ncoef) {
     cat("\n\nWnen some variables are not synthesised m  (= ",m,") must exceed number of",
         "\ncoefficients (= ",ncoef,") for lack of fit test. No test can be reported.\n\n")  
     lack.of.fit <- NULL
     lof.pvalue  <- NULL
   } else {
   QB <- matrix(NA, m, length(object$mcoefavg))
     for (i in 1:m) {
       QB[i,] <- object$mcoef[i,] - object$mcoefavg
     }
     lof.varcov <- t(QB) %*% QB/(m - 1)/m
     lack.of.fit <- t(res.diff[,3]) %*% ginv(lof.varcov) %*% res.diff[,3]
     lack.of.fit <- lack.of.fit * (object$m - ncoef)/ncoef/(object$m - 1) # Hotellings T square
     lof.pvalue  <- 1 - pf(lack.of.fit, ncoef, object$m - ncoef)          
   }
 } else {
   lof.varcov <- real.varcov*n/k/m
   lack.of.fit <- t(res.diff[,3]) %*% ginv(lof.varcov) %*% res.diff[,3]  #! multiply by m for combined test
   lof.pvalue  <- 1 - pchisq(lack.of.fit, ncoef)
 }

##?? if (object$proper == TRUE) lof.varcov <- real.varcov * (1 + n/k)/m

# Calculate summary measures
 mean.ci.overlap   <-  mean(res.overlap[,1])
 mean.abs.std.diff <-  mean(abs(res.diff[,4]))
 
 if (return.plot == TRUE) {
   yvar <- as.character(formula(object)[[2]])

 if (plot == "Z") {
     BetaCI <- dfCI(real.fit, Z = TRUE, ci.level = ci.level)
     
     if (population.inference == FALSE) {  ## get interval from real var
       BsynCI <- BetaCI
       for (i in c(1,3,4)) BsynCI[,i] <- BsynCI[,i] + (res.syn[,1] - res.obs[,1])/res.obs[,2]
       BsynCI[,5] <- "synthetic"
     } else {
       BsynCI <- dfCI(syn.fit, Z = TRUE, name.Z = "Z.syn", 
                      model.name = "synthetic", ci.level = ci.level)
     }
     xlab = "Z value"
     title = paste0("Z values for fit to ",yvar)
   
   } else {
     BetaCI <- dfCI(real.fit, Z = FALSE, ci.level = ci.level)  
     
     if (population.inference == FALSE) {  ## get interval from real var
       BsynCI <- BetaCI
       for (i in c(1,3,4)) BsynCI[,i] <- BsynCI[,i] + (res.syn[,1] - res.obs[,1])
       BsynCI[,5] <- "synthetic"
     } else {
       BsynCI <- dfCI(syn.fit, Z = FALSE, name.Z = "syn.coef", 
                      model.name = "synthetic", ci.level = ci.level)
     }

     xlab = "Value"
     title = paste0("Coefficients for fit to ",yvar)
   }
   
   modelCI <- rbind.data.frame(BetaCI, BsynCI)
   rownames(modelCI) <- 1:nrow(modelCI)
   
   if (!plot.intercept) modelCI <- modelCI[modelCI$Coefficient != "(Intercept)",]

   CI.geom <- geom_errorbar(aes_string(ymin = "LowCI", ymax = "HighCI",
     color = "Model", linetype = "Model"), data = modelCI, width = 0, 
     lwd = lwd, lty = lty, position = position_dodge(width = dodge.height))

   point.geom <- geom_point(aes_string(#ymin = value, ymax = value,  #BN-03/02/2017 commented
     color = "Model", shape = "Model"), data = modelCI, 
     size = point.size, position = position_dodge(width = dodge.height))

   p <- ggplot(data = modelCI, aes_string(x = "Coefficient", y = "Value"))
   p <- p + geom_hline(yintercept = 0, colour = "grey", linetype = 2, lwd = 1)
   p <- p + CI.geom + point.geom + labs(title = title, y = xlab)
   p <- p + scale_shape_manual(values = c(17:16), breaks = c("synthetic","observed")) +
            scale_colour_manual(values = lcol[2:1], breaks = c("synthetic", "observed"))
   p <- p + coord_flip()
   # p <- p + theme_bw()
   # scale_colour_manual(values = rev(brewer.pal(3,"Blues")))
   # scale_colour_grey(start = 0, end = .6)
   p
 } else p <- NULL

 res <- list(call = object$call, coef.obs = res.obs, coef.syn = res.syn, 
   coef.diff = res.diff, mean.abs.std.diff = mean.abs.std.diff,
   ci.overlap = res.overlap, mean.ci.overlap =  mean.ci.overlap, 
   lack.of.fit = lack.of.fit, lof.pvalue = lof.pvalue, 
   ci.plot = p, print.coef = print.coef,       
   m = object$m, ncoef = ncoef,
   incomplete = incomplete,
   population.inference = population.inference)  

 class(res) <- "compare.fit.synds"
 return(res)
}


###-----dfCI---------------------------------------------------------------
# extract info for plotting confidence intervals
dfCI <- function(modelsummary, names.est.se = c("Estimate","Std. Error"),
          model.name = "observed",  ci.level = 0.95, Z = FALSE, 
          name.Z = colnames(modelsummary$coefficients)[3]){
  
  CI <- qnorm(1 - (1 - ci.level)/2)
  if (!Z) {
    #msCI <- as.data.frame(modelsummary$coefficients[,names.est.se])
    msCI <- as.data.frame(modelsummary$coefficients[,1:2]) 
    names(msCI) <- c("Value", "SE")
  } else {
    # msCI <- as.data.frame(modelsummary$coefficients[,name.Z])
    msCI <- as.data.frame(modelsummary$coefficients[,3]) 
    names(msCI) <- c("Value")
    msCI$SE <- 1
  }  
  msCI$Coefficient <- rownames(msCI)

  msCI$HighCI <- msCI$Value + CI*msCI$SE
  msCI$LowCI  <- msCI$Value - CI*msCI$SE
  msCI$SE     <- NULL
  msCI$Model  <- model.name
  msCI$Coefficient <- factor(msCI$Coefficient, levels = rev(msCI$Coefficient))  #!BN290416, rev added 

  return(msCI)
}


###-----overlap.CI---------------------------------------------------------
overlap.CI <- function(synthetic, observed, ci.level, intercept, ...)
{
 CI <- qnorm(1- (1 - ci.level)/2)
 ##Initiate
 if(nrow(observed) > nrow(synthetic)){
   numVar <- nrow(synthetic); rNames <- rownames(synthetic)
 } else{
   numVar <- nrow(observed); rNames <- rownames(observed)
 } 
 CIoverlap <- matrix(NA, nrow = numVar, ncol = 1, 
   dimnames = list(rNames, "CI overlap"))
	
##Calculate CIoverlap
 for(i in 1:numVar) {

   ##Store CIs
   syn.upper <- synthetic[i, 1] + (CI * synthetic[i, 2])
   syn.lower <- synthetic[i, 1] - (CI * synthetic[i, 2])
   obs.upper <- observed[i, 1]  + (CI * observed[i, 2])
   obs.lower <- observed[i, 1]  - (CI * observed[i, 2])
	  
   ## CI overlap
   overlap.lower <- max(obs.lower, syn.lower)       
   overlap.upper <- min(obs.upper, syn.upper)       

   CIoverlap[i, 1] <- 0.5 * 
     (((overlap.upper - overlap.lower) / (obs.upper - obs.lower)) + 
	    ((overlap.upper - overlap.lower) / (syn.upper - syn.lower)))
 }

 if(intercept == FALSE) {
   CIoverlap = CIoverlap[-1, , drop = FALSE]
 }
  
 return(as.data.frame(CIoverlap))
}

 
