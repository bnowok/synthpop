###-----utility.tables-----------------------------------------------------
utility.tables <- function(object, data, ...) UseMethod("utility.tables")


###-----utility.tables.default---------------------------------------------
utility.tables.default <- function(object, ...)
  stop("No compare method associated with class ", class(object), call. = FALSE)


###-----utility.tables.data.frame---utility.tables.list--------------------
utility.tables.data.frame <- utility.tables.list <- 
                  function(object, data, 
                           cont.na = NULL, not.synthesised = NULL, 
                           tables = "twoway", method = "tab", 
                           vars = NULL, third.var = NULL,
                           useNA = TRUE, ngroups = 6,  
                           utility.stats = NULL, plot.stat = NULL,   
                           convergence.action = "setNA", plot = TRUE, 
                           print.tab.res = FALSE,  digits.tab.res = 2,
                           maxit = 100, max.scale = NULL, plottitle = NULL,  
                           nworst = 5, ntabstoprint  = 0, k.syn = FALSE, ...) {

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
                       tables = tables, method = method, 
                       vars = vars, third.var = third.var, 
                       useNA = useNA, ngroups = ngroups,
                       utility.stats = utility.stats, plot.stat = plot.stat,  
                       convergence.action = convergence.action, plot = plot,
                       print.tab.res = print.tab.res, digits.tab.res = digits.tab.res,
                       maxit = maxit, max.scale = max.scale, plottitle = plottitle, 
                       nworst = nworst, ntabstoprint = ntabstoprint, k.syn = k.syn, ...)
 
 res$call <- match.call()
 return(res)
}


###-----utility.tables-----------------------------------------------------
utility.tables.synds <- function(object, data, 
                           tables = "twoway", method = "tab", 
                           vars = NULL, third.var = NULL, 
                           useNA = TRUE, ngroups = 6,  
                           utility.stats = NULL, plot.stat = NULL, 
                           convergence.action = "setNA", plot = TRUE,
                           print.tab.res = FALSE, digits.tab.res = 2,
                           maxit = 100, max.scale = NULL, plottitle = NULL,
                           nworst = 5, ntabstoprint = 0, k.syn = FALSE, ...){

 if (is.null(object)) stop("Requires parameter 'object' to give name of the synthetic data object.\n", call. = FALSE)   
 if (is.null(data)) stop("Requires parameter 'data' to give name of the real data.\n", call. = FALSE)
 if (class(object) != "synds") stop("'object' must be of class 'synds', a synthetic data object", call. = FALSE)
 if (!is.data.frame(data)) stop("'data' must be a data frame.\n", call. = FALSE)
 if (!(method %in% c("gen","tab"))) stop("'method' must be 'gen' or 'tab'.\n", call. = FALSE)
 if (method == "tab" & !(tables %in% c("twoway", "oneway", "threeway"))) stop("tables for method 'tab' must be 'oneway', 'twoway' or 'threeway.'\n", call. = FALSE)
 if (method == "gen" & !(tables %in% c("twoway", "oneway", "threeway"))) stop("tables for method 'gen' must be 'oneway' or 'twoway.'\n", call. = FALSE)

 any.not.converged <- FALSE
 if (!(convergence.action %in% c("setNA", "none"))) cat('convergence.action must be "setNA" or "none".\n', call. = FALSE)
 
 if (is.null(vars) ) {
   if (object$m == 1) vars <- names(object$syn)
   else vars <- names(object$syn[[1]])
   vno <- 1:length(vars)
 }
 else if (is.character(vars)) {
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
 vars <- vars[ord]; vno <- vno[ord] ## get in right order in syn
 nv <- length(vars)
# now vars is character, vno is numeric, both refer to position in synthetic data

# check utility.stats (for tab.res) and plot.stat (for plots)
 if (method == "tab"){
   if (is.null(utility.stats)) utility.stats <- c("VW", "df")
   else if(!all(utility.stats %in%  c("VW", "FT", "G", "JSD", "pct.correct", "df", "all"))) stop('Parameter utility.stats for method tab can only include:\n"VW", "FT", "G", "JSD", "pct.correct", "df", and "all".\n', call. = FALSE)
   if (is.null(plot.stat)) plot.stat <- "VW"
   else if (!(length(plot.stat) == 1 && plot.stat %in%  c("VW", "FT", "G", "JSD", "pct.correct", "df"))) stop('Parameter utility.stats for method tab can only include one of:\n"VW", "FT", "G", "JSD", "pct.correct").\n', call. = FALSE)
   if (any(utility.stats == "all")) utility.stats <- c("VW", "FT", "G", "JSD", "pct.correct", "df")
 }
 
 if (method == "gen"){
   if (is.null(utility.stats)) utility.stats <- c("utilR", "df", "pval", "pct.correct")
   else if (!all(utility.stats %in%  c("pMSE", "pMSEExp", "utilR", "pval", "df", "pct.correct", "all"))) stop('Parameter utility.stats for method gen can only include:\n"pMSE", "pMSEExp", "utilR", "pval", "df", "pct.correct", "all".\n', call. = FALSE)
   if (is.null(plot.stat)) plot.stat <- "utilR"
   else if (!(length(plot.stat) == 1 && plot.stat %in% c("utilR", "pct.correct"))) stop('Parameter plot.statfor method tab can only include one of: "utilR", "pct.correct".\n', call. = FALSE)
   if (any(utility.stats == "all")) utility.stats <- c("pMSE", "pMSEExp", "utilR", "pval", "df", "pct.correct")
 }  

# end of check input parameters now add labels to names.vars to get plots in right order
 nv <- length(vars)
 names.vars <- paste(ntoc(vno), vars)

# code to make list of table utility values table of results and values to plot
 if (tables %in% c("twoway")) {
   npairs  <- nv*(nv - 1)/2
   pairs.num <- combn(1:nv, 2)
   pairs <- combn(vars, 2)
 } else if (tables %in% c("threeway")) {
   npairs  <- nv*(nv - 1)*(nv-2)/6
   pairs <- combn(vars, 3)
   pairs.num <- combn(1:nv, 3)
   pair.names <- apply(pairs, 2, function(x) paste(as.vector(x), collapse="_"))
 } else if (tables == "oneway") {
   utility.list <- as.list(1:nv)
   npairs <- nv
   pairs <- matrix(vars, 1, nv)
   pairs.num <- matrix(1:nv, 1, nv)
 }
  
 X1 <- X2 <- X3 <- rep("", length(npairs))
 for (i in 1:npairs){
   X1[i] <- names.vars[pairs.num[1, i]]
   if (!tables == "oneway")  X2[i] <- names.vars[pairs.num[2, i]]
   if (tables == "threeway") X3[i] <- names.vars[pairs.num[3, i]]
 }
  
 utility.list <- as.list(1:npairs)
 if (method == "tab"){   
   # now make tab.res of results
   ratioVW <- pvalVW <- ratioFT <- pvalFT <- ratioG <- pvalG <- 
     JSD <- pct.correct <- df <- rep(NA, length(npairs))
     
   for (i in 1:npairs) {
     utility.list[[i]] <- utility.tab(object, data, vars = pairs[,i], 
                                      ngroups = ngroups, k.syn = k.syn, 
                                      useNA = useNA )
     ratioVW[i]  <- mean(utility.list[[i]]$ratioVW)
     pvalVW[i]   <- median(utility.list[[i]]$pvalVW)
     ratioFT[i]  <- mean(utility.list[[i]]$ratioFT)
     pvalFT[i]   <- median(utility.list[[i]]$pvalFT)
     ratioG[i]   <- mean(utility.list[[i]]$ratioG)
     pvalG[i]    <- median(utility.list[[i]]$pvalG)
     JSD[i]      <- mean(utility.list[[i]]$JSD)
     pct.correct[i] <- mean(utility.list[[i]]$pct.correct) - 50
     df[i]       <- median(utility.list[[i]]$df)
   }   

   # results tab.res
   tab.res <- cbind.data.frame(ratioVW, pvalVW, ratioFT, pvalFT,
                               ratioG, pvalG,df, JSD, pct.correct) 
   colnames(tab.res) <- c("VW Ratio", "VW p-value", "FT Ratio", "FT p-value",
                          "G Ratio", "G p-value", "df", "JSD", "% pred >50")
   out <- NULL
   if (!("VW" %in% utility.stats)) out <- 1:2             
   if (!("FT" %in% utility.stats)) out <- c(out,3:4)
   if (!("G" %in% utility.stats))  out <- c(out,5:6)
   if (!("JSD" %in% utility.stats)) out <- c(out,8)
   if (!("df" %in% utility.stats))  out <- c(out,7)
   if (!("pct.correct" %in% utility.stats))  out <- c(out,9)
   if (!is.null(out)) tab.res  <- tab.res[, -out]    

   # and save value for plotting
   if (plot.stat == "VW") val <- ratioVW
   else if (plot.stat == "G") val <- ratioG
   else if (plot.stat == "FT") val <- ratioFT
   else if (plot.stat == "JSD") val <- JSD
   else if (plot.stat == "pct.correct") val <- pct.correct
 }

 if (method == "gen"){
   pMSE <- utilR <- pMSEExp <- utilR <- pval <- 
     df <- pct.correct <- rep(NA, npairs)
   if (object$m > 1) cat("Fitting general utility models for all combinations of variables")
   for (i in 1:npairs) {
     if (object$m > 1) cat("\nCombination ", i, " ", sep = "")
     utility.list[[i]] <- utility.gen(object, data, vars = pairs[,i], maxit = maxit, 
                                      method = "logit", maxorder = 1, k.syn = k.syn)
     if (object$m == 1) {
       if(!utility.list[[i]]$fit$converged) any.not.converged <- TRUE
      } else {
        for (im in 1:object$m) {
          if(!utility.list[[i]]$fit[[im]]$converged) any.not.converged <- TRUE  
        }
      }
    
     # now make tab.res of results
     pMSE[i] <- mean(utility.list[[i]]$pMSE)
     pMSEExp[i] <- median(utility.list[[i]]$pMSEExp)
     utilR[i] <- mean(utility.list[[i]]$utilR)
     pval[i] <- median(utility.list[[i]]$pval)
     pct.correct[i] <- mean(utility.list[[i]]$pct.correct)
     df[i] <- median(utility.list[[i]]$df)
     #utilStd[i] <- (utilR[i]*df[i] - df[i]) / sqrt(df[i]*2)
        
     if (any.not.converged ) {
       if (convergence.action == "setNA")  utilR[i] <- NA
       # if (convergence.action == "setNA" ) utilStd[i] <- NA
       if (convergence.action == "setNA") pct.correct[i] <- NA
     }
   }
      tab.res <- cbind.data.frame( pMSE, pMSEExp, utilR, pval, df, pct.correct)
      colnames(tab.res) <- c( "pMSE", "pMSEExp", "utilR", "p-value", "df", "% pred >50")
      out <- NULL
      if (!("pMSE" %in% utility.stats)) out <- 1             
      if (!("pMSEExp" %in% utility.stats)) out <- c(out, 2)
      if (!("utilR" %in% utility.stats)) out <- c(out, 3)
      #if (!("utilStd" %in% utility.stats)) out <- c(out,4)
      if (!("pval" %in% utility.stats)) out <- c(out, 4)
      if (!("df" %in% utility.stats)) out <- c(out, 5)
      if (!("pct.correct" %in% utility.stats)) out <- c(out, 6)
      if (!is.null(out)) tab.res  <- tab.res[,-out]

      # and save value for plotting
      if (plot.stat == "utilR") val <- utilR
      else if (plot.stat == "pct.correct") val <- pct.correct
 }
 
 if (tables == "oneway") dimnames(tab.res)[[1]] <- paste(X1)
 else if (tables == "twoway") dimnames(tab.res)[[1]] <- paste(X1, X2)
 else if (tables == "threeway") dimnames(tab.res)[[1]] <- paste(X1, X2, X3)

 # find  worst set of variables
 worstn <- val[order(-val)][1:nworst]
 names(worstn) <- dimnames(tab.res)[[1]][order(-val)][1:nworst]
 worstnind <- (1:npairs)[order(-val)][1:nworst]
 if (method == "tab"){
   worsttabs <- as.list(1:nworst)
   for (i in 1:nworst) {
     worsttabs[[i]] <- list(tab.obs = utility.list[[tab.obs = worstnind[i]]]$tab.obs, tab.syn = utility.list[[tab.obs = worstnind[i]]]$tab.syn)
   }
   names(worsttabs) <- names(worstn)
 } else worsttabs <- NULL

 # calculate vars w highest scores 
 var.scores <- NULL
 nn <- dimnames(tab.res)[[1]]
 num <- score <- rep(NA, nv)
 names(num) <- names(score) <- names.vars
 for ( i in 1:nv) {  
   num[i] <- sum(grepl(names.vars[i], nn))
   score[i] <- sum(val[grepl(names.vars[i], nn)])
 }
 var.scores <- sort(score/num, decreasing = TRUE)
 if (tables == "threeway"){
   if (is.null(third.var)) third.var <- names(var.scores)[1] #select worst as third.var if not specified
   else third.var <- names.vars[vars == third.var]
 } 
 if (tables == "oneway") toplot <- data.frame(X1 = X1, X2 = X2, val = val)
 else if (tables %in% c("twoway")) toplot <- rbind(data.frame(X1 = X1, X2 = X2, val = val), 
                                                   data.frame(X1 = X2, X2 = X1, val = val))
 else if (tables %in% c("threeway")) {
   toplot <- data.frame(X1 = X1, X2 = X2, X3 = X3, val = val)
   toplot <- toplot[(third.var == toplot$X1 | third.var == toplot$X2 | third.var == toplot$X3), ]
   toplot[third.var == toplot$X1, 1:2] <- toplot[third.var == toplot$X1, 2:3]
   toplot[third.var == toplot$X2, 2] <- toplot[third.var == toplot$X2, 3]
   v2 <- toplot[,-3]
   v2$X1 <- toplot$X2; v2$X2 <- toplot$X1
   toplot <- rbind(toplot[, c(1, 2, 4)], v2)
 }
    
 if (is.null(plottitle)){  
   if (tables == "twoway") plottitle <- paste("Two-way utility\n", plot.stat, 
                                              " for each pair of variables.", sep = "")
   else if (tables == "oneway") plottitle <- paste("One-way utility\n", plot.stat, 
                                                   " for each variable.", sep = "")
   else if (tables == "threeway") plottitle <- paste("Three-way utility", plot.stat, 
                                                     " for each pair of variables\n along with ", 
                                                     third.var, ".", sep = "")
   plottitle <- gsub("pct.correct", "% correctly predicted above 50%", plottitle)
 }

 res <- list(tab.res = tab.res, method = method, 
             plot.stat = plot.stat, tables = tables, 
             third.var = third.var, toplot = toplot,  
             any.not.converged =  any.not.converged, 
             convergence.action = convergence.action, 
             var.scores = var.scores, plot = plot,
             print.tab.res = print.tab.res, 
             digits.tab.res = digits.tab.res, 
             plottitle = plottitle, max.scale = max.scale,
             ntabstoprint = ntabstoprint, nworst = nworst, 
             worstn = worstn, worsttabs = worsttabs, ...)

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
  res
}


###-----print.utility.tables-----------------------------------------------
print.utility.tables <- function(x, print.tab.res = NULL, digits.tab.res = NULL,  
                                 plot = NULL, plottitle = NULL, max.scale = NULL, 
                                 nworst = NULL, ntabstoprint = NULL, ...) {
  
  if (x$method == "gen" && x$any.not.converged) {
    cat("\nWarning: Some general utility models did not converge, check results.\n")
    if (x$convergence.action == "setNA") cat("\nUtility ratio set to NA for non-converged models.\n")
  }                    

  if (is.null(max.scale)) max.scale <- x$max.scale
  if (!is.null(max.scale)) {
    if (max(x$toplot$val, na.rm = TRUE) <  max.scale) {
      cat("Selected scale expanded from ", max(x$toplot$val, na.rm = TRUE),
          " to ", max.scale, ".\n", sep = "")
    } else {
      cat("\nMaximum of scale set to  ", max.scale, " lower than maximum in results ",
          max(x$toplot$val, na.rm = TRUE), ".\n", sep = "")
      x$toplot$val[x$toplot$val > max.scale] <- max.scale
    }
  } else max.scale = max(x$toplot$val, na.rm = TRUE) 
  
  if (is.null(print.tab.res)) print.tab.res <- x$print.tab.res
  if (is.null(digits.tab.res)) digits.tab.res <- x$digits.tab.res
  if (is.null(plot)) plot <- x$plot
  if (is.null(plottitle)) plottitle <- x$plottitle 
  if (is.null(ntabstoprint)) ntabstoprint <- x$ntabstoprint
  if (is.null(nworst)) nworst <- x$nworst
  
  if (x$plot.stat == "pct.correct") x$plot.stat <- "% correctly predicted above 50%"
  if (x$tables == "twoway") cat("Two-way utility by method ", x$method,
                                " utility measure (", x$plot.stat,") plotted for ", 
                                dim(x$tab.res)[[1]], " pairs of variables.\n", sep = "")
  if (x$tables == "oneway") cat("Univariate utility by method ", x$method,
                                " and value plotted ", x$plot.stat," for ", 
                                dim(x$tab.res)[[1]], "variables.\n", sep = "")
  if (x$tables == "threeway") {
    cat("Three-way utility (total ", dim(x$tab.res)[[1]],") by method ",
        x$method,".\nAverage of 3-way scores ", x$plot.stat,
        " (ordered) for 3 way tables including each variable.\n", sep = "")
    print(x$var.scores)
    cat("\nVariable with highest average score, ", x$third.var,
        " selected to make plots. To see others set parameter 'third.var'.\n", sep = "")
  }

  if (x$tables %in% c("threeway")) cat("\nThree-way utility for pairs of variables along with ", 
                                       x$third.var, " by method ", x$method, ".\n", sep = "")

  cat("Variable combinations with worst ", nworst ,
      " utility scores (", x$plot.stat,"):\n", sep = "")
  print(x$worstn)
  if (ntabstoprint > nworst) cat("Only ", nworst, 
    " tables can be printed. For more rerun with 'nworst' set to a larger number.\n", sep = "")

  if (ntabstoprint > 0 & x$method == "tab")  {
    cat("\nPrinting table of observed and synthetic data for the ", 
        ntabstoprint, " table with the worst utility\n", sep = "")
    for (i in 1:ntabstoprint) {
      cat("Tables of ", names(x$worsttabs[i]), "\nOriginal data\n", sep = "")
      print(x$worsttabs[[i]]$tab.obs)
      cat("Synthetic data\n")
      print(x$worsttabs[[i]]$tab.syn)
    }
  }
  
  if (plot){ 
    if (x$convergence.action == "setNA" & x$any.not.converged == TRUE) plottitle <- paste(plottitle, 
      "Grey indicates that logit model did not converge", sep = "\n")
    p <- ggplot(x$toplot, aes(x = x$toplot$X2, y = x$toplot$X1)) + 
      geom_raster(aes(fill = x$toplot$val)) + 
      scale_fill_gradient(low = "grey90", high = "red", limits = c(0,max.scale)) +
      labs(x = "", y = "", title = plottitle) +
      theme_bw() + theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.7),
                         axis.text.y = element_text(size = 7),
                         title = element_text(size = 11),
                         legend.title = element_blank())
    print(p)
  }

  if (!print.tab.res) {
    cat("\nMedian and maximum of selected utility measures for all tables compared\nMedian\n")
    tab <- x$tab.res
    medians <- apply(tab, 2, function(x) median(x, na.rm = TRUE))
    medians <- medians[!names(medians) == "df"]
    medians <- medians[!grepl("p-value", names(medians))]
    print(round(medians, digits.tab.res))
    
    maxima <- apply(tab, 2, function(x) max(x, na.rm = TRUE))
    maxima <- maxima[!names(maxima) == "df"]
    maxima <- maxima[!grepl("p-value", names(maxima))]
    cat("Maxima\n")
    print(round(maxima, digits.tab.res))
    cat("\nFor more details of all scores use print.tab.res = TRUE.\n")
  } else {
    cat("\nTable of selected utility measures\n")
    tab.res <- x$tab.res
    tab.res[, grepl("p-value", names(tab.res)) | grepl("pMSE", names(tab.res))]  <- 
      round(tab.res[, grepl("p-value", names(tab.res)) | grepl("pMSE", names(tab.res))], 4)
    tab.res[, !grepl("p-value", names(tab.res)) & !grepl("pMSE", names(tab.res))] <- 
      round(tab.res[, !grepl("p-value", names(tab.res)) & !grepl("pMSE", names(tab.res))], digits.tab.res)
    print(tab.res)
  }

  invisible(x)
}
