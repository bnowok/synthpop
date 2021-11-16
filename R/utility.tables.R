###-----utility.tables-----------------------------------------------------
utility.tables <- function(object, data, ...) UseMethod("utility.tables")


###-----utility.tables.default---------------------------------------------
utility.tables.default <- function(object, ...)
  stop("No compare method associated with class ", class(object), call. = FALSE)


###-----utility.tables.data.frame---utility.tables.list--------------------
utility.tables.data.frame <- utility.tables.list <- 
                  function(object, data, 
                           cont.na = NULL, not.synthesised = NULL, 
                           tables = "twoway", maxtables = 5e4, 
                           vars = NULL, third.var = NULL,
                           useNA = TRUE, ngroups = 5,
                           tab.stats = c("pMSE", "S_pMSE", "df"), 
                           plot.stat = "S_pMSE", plot = TRUE,
                           print.tabs = FALSE, digits.tabs = 4,
                           max.scale = NULL, min.scale = 0, plot.title = NULL,
                           nworst = 5, ntabstoprint = 0, k.syn = FALSE, ...) {

                     if (is.null(data)) stop("Requires parameter 'data' to give name of the real data.\n", call. = FALSE)
 if (is.null(object)) stop("Requires parameter 'object' to give name of the synthetic data.\n", call. = FALSE)   
  
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
 if  (!is.null(not.synthesised)) {
   if ( !is.null(not.synthesised) && !all(not.synthesised %in% names(data))) stop("not.synthesised must be names of variables in data\n\n", call. = FALSE)
   syn.method[names(data) %in% not.synthesised] <- ""
 }
  
 object <- list(syn = object, m = m, strata.syn = NULL, 
                method = syn.method, cont.na = cont.na)
 class(object ) <- "synds"
  
 res <- utility.tables.synds(object = object, data = data, 
                       tables = tables, maxtables = maxtables, 
                       vars = vars, third.var = third.var, 
                       useNA = useNA, ngroups = ngroups, 
                       tab.stats = tab.stats, plot.stat = plot.stat, 
                       plot = plot, print.tabs = print.tabs, 
                       digits.tabs = digits.tabs, max.scale = max.scale, 
                       min.scale = min.scale, plot.title = plot.title, 
                       nworst = nworst, ntabstoprint = ntabstoprint, 
                       k.syn = k.syn, ...)
 
 res$call <- match.call()
 return(res)
}


###-----utility.tables-----------------------------------------------------
utility.tables.synds <- function(object, data, 
                                 tables = "twoway", maxtables = 5e4, 
                                 vars = NULL, third.var = NULL,
                                 useNA = TRUE, ngroups = 5,  
                                 tab.stats = c("pMSE", "S_pMSE", "df"), 
                                 plot.stat = "S_pMSE", plot = TRUE,
                                 print.tabs = FALSE, digits.tabs = 4,
                                 max.scale = NULL, min.scale = 0, plot.title = NULL,
                                 nworst = 5, ntabstoprint = 0, k.syn = FALSE, ...){

 if (is.null(object)) stop("Requires parameter 'object' to give name of the synthetic data object.\n", call. = FALSE)   
 if (is.null(data)) stop("Requires parameter 'data' to give name of the original data.\n", call. = FALSE)
 if (class(object) != "synds") stop("'object' must be of class 'synds', a synthetic data object", call. = FALSE)
 if (!is.data.frame(data)) stop("'data' must be a data.frame.\n", call. = FALSE)
 if (!(tables %in% c("twoway", "oneway", "threeway"))) stop("Argument tables must be 'oneway', 'twoway' or 'threeway.'\n", call. = FALSE)

 if (is.null(vars) ) {
   if (object$m == 1) vars <- names(object$syn)
   else vars <- names(object$syn[[1]])
   vno <- 1:length(vars)
 } else if (is.character(vars)) {
   if (!all(vars %in% names(data))) stop("vars must be in names of original data.\n", call. = FALSE)
   if (object$m == 1){
     if (!all(vars %in% names(object$syn))) stop("vars must be in names of synthetic data.\n", call. = FALSE)
     vno <- match(vars, names(object$syn))
   } else if (object$m > 1){
     if (!all(vars %in% names(object$syn[[1]]))) stop("vars must be in names of synthetic data.\n", call. = FALSE)
     vno <- match(vars, names(object$syn[[1]]))
   }
 } else if (is.numeric(vars)) {
   vno <- vars
   if (object$m == 1)  {
     if(!all(vars %in% 1:length(object$syn))) stop("vars must be in 1:length(object$syn).\n", call. = FALSE)
     vars <- names(object$syn)[vno]
   } else if (object$m > 1) {  
     if(!all(vars %in% 1:length(object$syn[[1]]))) stop("vars must be in 1:length(object$syn[[1]]).\n", call. = FALSE)
     vars <- names(object$syn[[1]])[vno]
   }
 }
 ord <- order(vno)
 vars <- vars[ord]; vno <- vno[ord]  ## get in right order in syn
 # now vars is character, vno is numeric, both refer to position in synthetic data

# check tab.stats (for tabs) and plot.stat (for plots)
   if (is.null(tab.stats)) tab.stats <- c("pMSE", "S_pMSE", "df")
   else if(!all(tab.stats %in%  c("VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", 
     "MabsDD", "dBhatt", "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE", "df", "dfG", "all"))) stop('Parameter tab.stats can only include:
"VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", 
"MabsDD", "dBhatt", "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE", "df", "dfG", "all".\n', call. = FALSE, sep = "")
   if (any(tab.stats == "all")) tab.stats <- c("VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", 
     "MabsDD", "dBhatt", "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE", "df", "dfG")
 
   if (is.null(plot.stat)) plot.stat <- "S_pMSE"
   else if (!(length(plot.stat) == 1 && plot.stat %in%  c("VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", 
     "MabsDD", "dBhatt", "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE"))) 
   stop('Parameter plot.stat must be just one of:\n"VW", "FT", "JSD", "SPECKS", "WMabsDD", "U", "G", "pMSE", "PO50", 
"MabsDD", "dBhatt", "S_VW", "S_FT", "S_JSD", "S_WMabsDD", "S_G", "S_pMSE".\n', call. = FALSE)
# End of check of input parameters 

# Add labels to names.vars to get plots in right order
 nv <- length(vars)
 names.vars <- paste(ntoc(vno), vars, sep = ".")

# Code to make list of table utility values table of results and values to plot
 if (tables %in% c("twoway")) {
   npairs    <- nv*(nv - 1)/2
   pairs.num <- combn(1:nv, 2)
   pairs <- combn(vars, 2)
 } else if (tables == "threeway") {
   npairs <- nv*(nv - 1)*(nv - 2)/6
   pairs  <- combn(vars, 3)
   pairs.num <- combn(1:nv, 3)
   pair.names <- apply(pairs, 2, function(x) paste(as.vector(x), collapse=":"))
 } else if (tables == "oneway") {
   utility.list <- as.list(1:nv)
   npairs <- nv
   pairs <- matrix(vars, 1, nv)
   pairs.num <- matrix(1:nv, 1, nv)
 }

 if (npairs > maxtables) {
    sel <- sample(1:npairs, maxtables)
    pairs <- pairs[, sel, drop = FALSE]
    pairs.num <- pairs.num[, sel, drop = FALSE]
    npairs <- maxtables
    cat("Total tables requested exceeds limit of maxtables = ", maxtables,
        ", so utility is measured for a sample of ", maxtables, ".\n", sep = "")
    if (tables == "threeway") {
      plot <- FALSE    
      cat("You cannot produce plots of threeway tables from sampled tables.\n", sep = "")
    }
 }
 
 X1 <- X2 <- X3 <- rep("", npairs)
 for (i in 1:npairs){
   X1[i] <- names.vars[pairs.num[1, i]]
   if (!tables == "oneway")  X2[i] <- names.vars[pairs.num[2, i]]
   if (tables == "threeway") X3[i] <- names.vars[pairs.num[3, i]]
 }
  
 utility.list <- as.list(1:npairs)
   # now make tabs of results
   tabs <- matrix(NA, npairs, length(tab.stats)) 
   colnames(tabs) <- tab.stats
   if (tables == "oneway") dimnames(tabs)[[1]] <- paste(X1)
   else if (tables == "twoway") dimnames(tabs)[[1]] <- paste(X1, X2, sep = ":")
   else if (tables == "threeway") dimnames(tabs)[[1]] <- paste(X1, X2, X3, sep = ":")
   
   val <- rep(NA, npairs)
   
   for (i in 1:npairs) {
     utility.list[[i]] <- utility.tab(object, data, vars = pairs[, i], 
                                      ngroups = ngroups, k.syn = k.syn, 
                                      useNA = useNA, print.flag = FALSE, ...)
     if (i == 1) {
       tab.ind <- match(tab.stats, names(utility.list[[i]]))
       val.ind <- match(plot.stat, names(utility.list[[i]]) )
     }
     tabs[i,] <- sapply(utility.list[[i]][tab.ind], mean)
     val[i]      <- sapply(utility.list[[i]][val.ind], mean) # value for plotting
   }  
 
 # find worst set of variables
 nworst <- min(npairs, nworst)
 worstn <- val[order(-val)][1:nworst]
 names(worstn) <- dimnames(tabs)[[1]][order(-val)][1:nworst]
 worstnind <- (1:npairs)[order(-val)][1:nworst]
 worsttabs <- as.list(1:nworst)
 for (i in 1:nworst) {
   worsttabs[[i]] <- list(tab.obs = utility.list[[tab.obs = worstnind[i]]]$tab.obs, 
                          tab.syn = utility.list[[tab.obs = worstnind[i]]]$tab.syn)
 }
 names(worsttabs) <- names(worstn)

 # calculate vars with highest scores 
 var.scores <- NULL
 nn <- dimnames(tabs)[[1]]
 num <- score <- rep(NA, nv)
 names(num) <- names(score) <- names.vars
 for (i in 1:nv) {  
   num[i]   <- sum(grepl(names.vars[i], nn))
   score[i] <- sum(val[grepl(names.vars[i], nn)])
 }
 var.scores <- sort(score/num, decreasing = TRUE)
 if (tables == "threeway"){
   if (is.null(third.var)) third.var <- names(var.scores)[1]  # select worst as third.var if not specified
   else third.var <- names.vars[vars == third.var]
 } 
 
 if (tables == "oneway") toplot <- data.frame(X1 = X1, X2 = X2, val = val)
 else if (tables == "twoway") toplot <- rbind(data.frame(X1 = X1, X2 = X2, val = val), 
                                              data.frame(X1 = X2, X2 = X1, val = val))
 else if (tables == "threeway") {
   toplot <- data.frame(X1 = X1, X2 = X2, X3 = X3, val = val)
   toplot <- toplot[(third.var == toplot$X1 | third.var == toplot$X2 | third.var == toplot$X3), ]
   toplot[third.var == toplot$X1, 1:2] <- toplot[third.var == toplot$X1, 2:3]
   toplot[third.var == toplot$X2,   2] <- toplot[third.var == toplot$X2,   3]
   v2 <- toplot[, -3]
   v2$X1 <- toplot$X2
   v2$X2 <- toplot$X1
   toplot <- rbind(toplot[, c(1, 2, 4)], v2)
 }
    
 if (is.null(plot.title)){  
   if (tables == "twoway") plot.title <- bquote("Two-way utility:"~bold(.(plot.stat))~"for pairs of variables")
   else if (tables == "oneway") plot.title <- bquote("One-way utility:"~bold(.(plot.stat))~"for each variable")
   else if (tables == "threeway") plot.title <- bquote("Three-way utility:"~bold(.(plot.stat))~"for pairs along with"~bold(.(third.var)))
 }

 if (!is.null(max.scale)) {
   if (max(toplot$val, na.rm = TRUE) >  max.scale) {
     cat("Maximum of plot scale set to ", max.scale, 
         " (lower than maximum in results ",
         max(toplot$val, na.rm = TRUE), ").\n", sep = "")
     toplot$val[toplot$val > max.scale] <- max.scale
   }
 } else max.scale <- max(toplot$val, na.rm = TRUE) 
 
 if (!is.null(min.scale) & min.scale != 0 & (min(toplot$val, na.rm = TRUE) < min.scale)) {
   cat("Minimum of plot scale set to ", min.scale, " (greater than minimum in results ",
       min(toplot$val, na.rm = TRUE), ").\n", sep = "")
   toplot$val[toplot$val < min.scale] <- min.scale
 }

 utility.plot <- ggplot(toplot, aes(x = X2, y = X1)) + 
   geom_raster(aes(fill = val)) + 
   scale_fill_gradient(low = "grey92", high = "#E41A1C", limits = c(0, max.scale)) + 
   labs(x = "", y = "", title = plot.title) +
   theme_minimal() + theme(axis.text.x = element_text(size = 10, angle = 90), # vjust = 0.7),
                      axis.text.y = element_text(size = 10),
                      title = element_text(size = 11),
                      legend.title = element_blank(),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
 
 res <- list(tabs = tabs,
             plot.stat = plot.stat, 
             tables = tables, 
             third.var = third.var, 
             utility.plot = utility.plot,  
             var.scores = var.scores, plot = plot,
             print.tabs = print.tabs, 
             digits.tabs = digits.tabs, 
             plot.title = plot.title, 
             max.scale = max.scale, min.scale = min.scale,
             ntabstoprint = ntabstoprint, nworst = nworst, 
             worstn = worstn, worsttabs = worsttabs)

 class(res) <- "utility.tables"
 return(res)
}


###-----ntoc---------------------------------------------------------------
# to make labels for variables of constant length from an integer
ntoc <- function(x)
{
  nch <- nchar(max(x))
  res <- as.character(x)
  res <- str_pad(res, nch, "left", "0")
  return(res)
}


###-----print.utility.tables-----------------------------------------------
print.utility.tables <- function(x, print.tabs = NULL, digits.tabs = NULL,  
                                 plot = NULL, plot.title = NULL, 
                                 max.scale = NULL, min.scale = NULL,
                                 nworst = NULL, ntabstoprint = NULL, ...) {

  if (is.null(print.tabs)) print.tabs <- x$print.tabs
  if (is.null(digits.tabs)) digits.tabs <- x$digits.tabs
  if (is.null(plot)) plot <- x$plot
  if (is.null(plot.title)) plot.title <- x$plot.title 
  if (is.null(ntabstoprint)) ntabstoprint <- x$ntabstoprint
  if (is.null(nworst)) nworst <- x$nworst 
  if (is.null(max.scale)) max.scale <- x$max.scale
  if (is.null(min.scale)) min.scale <- x$min.scale
  
  if (x$tables == "twoway") cat("\nTwo-way utility: ", x$plot.stat, " value plotted for ", 
                                dim(x$tabs)[[1]], " pairs of variables.\n", sep = "")
  if (x$tables == "oneway") cat("\nUnivariate utility: ", x$plot.stat, " value plotted for ", 
                                dim(x$tabs)[[1]], " variables.\n", sep = "")
  if (x$tables == "threeway") {
    cat("\nThree-way utility (total of ", dim(x$tabs)[[1]]," variable combinations):\n", sep = "")
    cat("\nAverage of 3-way scores ", x$plot.stat,
        " (ordered) for 3-way tables including each variable.\n", sep = "")
    print(x$var.scores)
    cat("\nVariable with highest average score, ", x$third.var,
        ", selected to make plots.\nTo see others, set parameter 'third.var'.\n", sep = "")
  }

  cat("\nVariable combinations with worst ", nworst ,
      " utility scores (", x$plot.stat,"):\n", sep = "")
  print(round(x$worstn, digits.tabs))
  
  if (ntabstoprint > nworst) {
    ntabstoprint <- nworst
    cat("Only ", nworst, 
    " tables can be printed. For more rerun with 'nworst' set to a larger number (up to ", 
    dim(x$tabs)[[1]], ").\n", sep = "")
  }  
  if (ntabstoprint > 0) {
    cat("\nPrinting table of observed and synthetic data for the ", 
        ntabstoprint, " table with the worst utility\n", sep = "")
    for (i in 1:ntabstoprint) {
      cat("Tables of ", names(x$worsttabs[i]), "\nOriginal data:\n", sep = "")
      print(x$worsttabs[[i]]$tab.obs)
      cat("Synthetic data:\n")
      print(x$worsttabs[[i]]$tab.syn)
    }
  }
  
  if (plot) {
    if (!is.null(max.scale)) x$utility.plot$scales$scales[[1]]$limits[2] <- max.scale
    if (!is.null(min.scale)) x$utility.plot$scales$scales[[1]]$limits[1] <- min.scale
    print(x$utility.plot)
  }

  if (!print.tabs) {
    cat("\nMedians and maxima of selected utility measures for all tables compared\n")
    Medians <- apply(x$tabs, 2, function(x) median(x, na.rm = TRUE))
    Maxima  <- apply(x$tabs, 2, function(x) max(x, na.rm = TRUE))
    print(round(data.frame(Medians, Maxima), digits.tabs))
    cat("\nFor more details of all scores use print.tabs = TRUE.\n")
  } else {
    cat("\nTable of selected utility measures\n")
    print(round(x$tabs, digits.tabs))
  }

  invisible(x)
}
