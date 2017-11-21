###-----compare------------------------------------------------------------
compare <- function(object, data, ...) UseMethod("compare")

###-----compare.default----------------------------------------------------
compare.default <- function(object, ...)
  stop("No compare method associated with class ", class(object), call. = FALSE)


###-----compare.synds------------------------------------------------------
compare.synds <- function(object, data, vars = NULL, 
  msel = NULL, breaks = 20, nrow = 2, ncol = 2, rel.size.x = 1,
  cols = c("#1A3C5A","#4187BF"), ...){ 
                                                                         ## GR6 & GR7 extra parameter nrow and drop plot.na
 if (is.null(data)) stop("Requires parameter data to give name of the real data\n", call.=FALSE)
 if (!is.data.frame(data)) stop("Argument data must be a data frame\n", call.=FALSE)          ## GR2 extra check
 #if (class(object)!="synds") stop("Object must have class synds\n", call.=FALSE )                                                                    
 if (!is.null(msel) & !all(msel %in% (1:object$m))) stop("Invalid synthesis number(s)", call.=FALSE)                                                                        

 # browser()
 # single / pooled synthetic data sets                                  
 if (object$m==1) {
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

 if (!is.null(vars)){
   if (!(all(vars %in% synnames))) stop("Variable(s) ", 
     paste0(vars[is.na(match(vars,synnames))], collapse=", "),
     " not in synthetic data \n", call.=FALSE)
   if (!(all(vars %in% realnames))) stop("Variable(s) ", 
     paste0(vars[is.na(match(vars,realnames))], collapse=", "),
     " not in observed data \n", call.=FALSE)
   commonnames <- commonnames[match(vars,commonnames)]
 }
 
 if (!(all(synnames %in% realnames))) cat("Warning: Variable(s)", 
   paste0(synnames[is.na(match(synnames,realnames))], collapse=", "),
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
 num <- sapply(df.synall, is.numeric) | sapply(df.synall, is.integer)  
 fac <- sapply(df.synall, function(x) is.factor(x) | is.logical(x))  

 # frequency tables for factors
 if (sum(fac) > 0) {
   any.fac.na <- unlist(apply(df.obs[,fac,drop=FALSE],2,function(x) any(is.na(x))))    
   per.obs.fac <- ggfac(df.obs[, fac, drop = FALSE], anyna = any.fac.na)

   if (length(msel) <= 1) per.syn.facall <- ggfac(df.synall[, fac, drop = FALSE], 
     name = "synthetic", anyna = any.fac.na)
   if (length(msel) > 1) {
     per.syn.fac <- vector("list",length(msel))
     for (i in 1:length(msel)) per.syn.fac[[i]] <- ggfac(df.syn[[i]][, fac, drop = FALSE], 
       name = paste0("syn=", msel[i]), anyna = any.fac.na)
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
   any.na <- unlist(apply(df.obs[,num,drop=FALSE],2,function(x) any(is.na(x))))  
   
   lbreaks    <- as.list(rep(breaks, length(na)))                   
  ## get limits(=breaks) from both observed and synthetic
  #--
   df.both  <- rbind.data.frame(df.obs, df.synall)                
   per.both <- ggnum(df.both[, num, drop = FALSE], na = na, 
     breaks = lbreaks, anyna = any.na) 
  #--  
   per.obs.num <- ggnum(df.obs[, num, drop = FALSE], na = na, 
     breaks = per.both$hbreaks, anyna = any.na)
   
   if (length(msel) <= 1) {
     per.syn.numall <- ggnum(df.synall[, num, drop = FALSE], name = "synthetic", 
       na = na, breaks = per.both$hbreaks, anyna = any.na) 
   } 
   if (length(msel) > 1) {
     per.syn.num <- vector("list",length(msel))
     for (i in 1:length(msel)) per.syn.num[[i]] <- ggnum(df.syn[[i]][, num, drop = FALSE], 
       name =  paste0("syn=",msel[i]), 
       na = na, breaks = per.both$hbreaks, anyna = any.na)
   } 
 } else {
   per.obs.num <- NULL
   per.syn.numall <- NULL
 }

 # data frame for ploting 
 if (length(msel) <= 1) per.fac <- rbind.data.frame(per.obs.fac$perdf, 
   per.obs.num$perdf, per.syn.facall$perdf, per.syn.numall$perdf )

 if (length(msel) > 1) {
   per.fac <- rbind.data.frame(per.obs.fac$perdf, per.obs.num$perdf)
   for (i in 1:length(msel)) {
     if (sum(fac)>0) temp.fac <- per.syn.fac[[i]]$perdf else temp.fac <- NULL
     if (sum(num)>0) temp.num <- per.syn.num[[i]]$perdf else temp.num <- NULL  
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
     if (sum(fac)>0) temp.fac <- per.syn.fac[[i]]$perlist else temp.fac <- NULL
     if (sum(num)>0) temp.num <- per.syn.num[[i]]$perlist else temp.num <- NULL
     os.table.fac <- mapply(rbind, os.table.fac, temp.fac, SIMPLIFY = FALSE)    
     os.table.num <- mapply(rbind, os.table.num, temp.num, SIMPLIFY = FALSE)     
   }
  }
 os.table <- c(os.table.fac, os.table.num)
 if (is.null(msel)){
   for (i in 1:length(os.table)) dimnames(os.table[[i]])[[1]] <- c("observed","synthetic")
 } else {
   for (i in 1:length(os.table)) dimnames(os.table[[i]])[[1]] <- c("observed",paste0("syn=", msel))
 }
 Value <- Percent <- Data <- NULL                               ## otherwise 'no visible binding for global variables'
 # sorts the factor labels in the right order for numeric vars
 per.fac$Value <- as.character(per.fac$Value)
 vals <- unique(per.fac$Value)
 valsnum <- unique(per.fac$Value[per.fac$Variable %in% names(num[num==TRUE])])
 valsnum.nonmiss <- sort(as.numeric(vals[vals %in% valsnum & substr(vals,1,4)!="miss"]))
 valsnum.miss <- sort(vals[vals %in% valsnum & substr(vals,1,4)=="miss"])
 vals[vals %in% valsnum] <- c(valsnum.nonmiss,valsnum.miss)
 per.fac$Value <- factor(as.character(per.fac$Value),levels=vals)

 # get different plots in order of data
 nperplot <- nrow*ncol
 nplots   <- ceiling(length(commonnames)/nperplot)
 plots    <- vector("list",nplots)
 tables   <- vector("list",nplots)
 
 for (i in 1:nplots) {
   min <-(i-1)*nperplot+1
   max <- min(length(commonnames),(i-1)*nperplot+nperplot)
 # tables
   ttables <- vector("list",(max-min+1))
   names(ttables) <- commonnames[min:max]
   for (j in commonnames[min:max]) {
     ttables[[j]] <- os.table[[j]]
   } 
   tables[[i]] <- ttables

  # plots
   per.fact <- per.fac[per.fac$Variable %in% commonnames[min:max],]
   # per.fact <- per.fact[order(match(per.fact$Variable,commonnames[min:max])),]
   p <- ggplot(data=per.fact, aes(x=Value,y=Percent,fill=Data))
   p <- p + geom_bar(position="dodge",colour=cols[1], stat="identity") + 
      facet_wrap(~ Variable, scales="free", ncol=ncol)  
   p <- p + guides(fill = guide_legend(override.aes = list(colour = NULL))) + 
        theme(axis.text.x=element_text(angle=-30, hjust=0, vjust=1,size=rel(rel.size.x)), 
              legend.position="top", 
              legend.key = element_rect(colour = cols[1])) 
   p <- p + theme(legend.title=element_blank())
   if (length(msel) > 1) p <- p + scale_fill_manual(values = c(cols[1],rep(cols[2],length(msel))))
   if (length(msel) <= 1) p <- p + scale_fill_manual(values = cols)
   plots[[i]] <- p
  }
 
 if (length(tables)==1) {
   tables <- tables[[1]]
   plots  <- plots[[1]]
 }  
 #browser()
 res <- list(tables = tables, plots = plots)
 class(res) <- "compare.synds"
 return(res)
}


###-----pertable-----------------------------------------------------------
# calculate percentages
pertable <- function(x) {
 res <- table(x, useNA = "always")*100/sum(table(x, useNA = "always")) 
 return(res)
}


###-----ggfac--------------------------------------------------------------
# calculate percentages for factors and store in a data frame (a long format)
ggfac <- function(data, name = "observed", anyna){ 
  data <- as.data.frame(data)
  perlist  <- lapply(data, pertable)
  for (i in 1:length(perlist)) {
    if (anyna[i]==FALSE) perlist[[i]] <- perlist[[i]][-length(perlist[[i]])]  
  }
  Percent  <- unlist(perlist, use.names = F)
  Value    <- unlist(lapply(perlist,names), use.names = F)
  Variable <- rep(names(perlist),sapply(perlist,length))
  perdf    <- cbind.data.frame(Percent, Value, Variable)
  perdf$Data <- name
  return(list(perdf = perdf, perlist = perlist))
} 


###-----ggnum--------------------------------------------------------------
# calculate percentages for numeric variables (non-missing values and 
# missing data categories seperately) and store in a data frame (a long format)
ggnum <- function(data, name = "observed", na = as.list(rep(NA,ncol(data))), 
                  breaks = as.list(rep(30,ncol(data))), anyna){ 
  data <- as.data.frame(data)

# non-missing values  
  nvar <- ncol(data)
  perlist <- vector("list", nvar)
  hbreaks <- vector("list", nvar)
  
  for (i in 1:nvar){
  ## counts for non-missing values  
    vardata <- data[!(data[,i] %in% na[[i]]),i]                  
    hh      <- hist(vardata, breaks = breaks[[i]], plot=F)
    counts  <- hh$counts
    names(counts) <- hh$breaks[-length(hh$breaks)]
    hbreaks[[i]]  <- hh$breaks
  ## counts for missing values  
    dataNA <- data[(data[,i] %in% na[[i]]),i] 
    
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
    if (!is.null(NAcounts)) names(NAcounts) <- paste("miss", names(NAcounts), sep=".") 
    counts <- c(counts, NAcounts)
    perlist[[i]] <- counts*100/length(data[,i])
  }
  names(perlist) <- colnames(data)

# create data frame in a long format  
  Percent  <- unlist(perlist, use.names = F)
  Value    <- unlist(lapply(perlist,names), use.names = F)
  Variable <- rep(names(perlist),sapply(perlist,length))
  perdf    <- cbind.data.frame(Percent, Value, Variable)
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
  dodge.height = .5, point.size = 2.5, incomplete = FALSE, 
  population.inference = FALSE, ci.level = 0.95, ...) {   # c("#132B43", "#56B1F7")

 # compares and plots fits to synthetic and original data
 # first parameter must be a fit to synthetic data from glm.synds()
 # or lm.synds()
 value <- "Value"
 coefficient <- c("Coefficient", "Model")
 m <- object$m
 n <- sum(object$n)
 if (is.list(object$k)) k <- sum(object$k[[1]]) else k <- sum(object$k)
 fitting.function <- object$fitting.function
 syn.coef    <- object$mcoefavg
 
 call <- match.call()
 # if (!class(object)=="fit.synds") stop("Object must have class fit.synds\n")
 if (!is.data.frame(data)) stop("Data must be a data frame\n")  # theoretically can be a matrix (?)
 if (ci.level <= 0 | ci.level > 1) stop("ci.level must be beteen 0 and 1.\n")
 if (incomplete == TRUE & m < length(syn.coef) + 1) stop("\nYou have selected incompletely synthesised data with m less than number of coefficients.\nLack-of-fit test cannot be calculated unless m > ", length(syn.coef),".\n", call. = FALSE)

 # get fit to real data
 if (fitting.function=="lm") {
 real.fit <- summary(do.call(object$fitting.function,
                     args=list(formula=formula(object),data=call$data)))
 } else {
 real.fit <- summary(do.call(object$fitting.function,
                     args=list(formula=formula(object),
                     family=object$call$family,data=call$data)))
 }
 
 real.varcov <- real.fit$cov.unscaled                  
 real.coef   <- real.fit$coefficients[,1] 
 real.se     <- real.fit$coefficients[,2]  


 if (population.inference==FALSE) { syn.fit  <- summary.fit.synds(object,real.varcov=real.varcov,incomplete=incomplete)}
 else syn.fit  <- summary.fit.synds(object,incomplete=incomplete,population.inference=TRUE,real.varcov=real.varcov)
 
 if (return.plot == TRUE){
   yvar <- as.character(formula(object)[[2]])

   # # prepare data for plotting confidence intervals (one data frame)
   # if (plot=="Z"){
   #   BetaCI <- dfCI(real.fit, Z = TRUE, ci.level = ci.level)
   #   BsynCI <- dfCI(syn.fit, Z = TRUE, name.Z = "Z.syn",
   #               model.name = "synthetic", ci.level = ci.level)
   #   xlab = "Z value"
   #   title = paste0("Z values for fit to ",yvar)
   # } else {
   #   BetaCI <- dfCI(real.fit)
   #   if (population.inference == FALSE) {
   #     BsynCI <- dfCI(syn.fit, names.est.se = c("B.syn","se(B).syn"),
   #                    model.name = "synthetic", ci.level = ci.level)  
   #   } else {
   #     BsynCI <- dfCI(syn.fit, names.est.se = c("B.syn","se(B.syn)"),
   #                    model.name = "synthetic", ci.level = ci.level)
   #   }
   #   xlab = "Value"
   #   title = paste0("Coefficients for fit to ",yvar)
   # }
   
   # prepare data for plotting confidence intervals (one data frame)
   if (plot=="Z"){
     BetaCI <- dfCI(real.fit, Z = TRUE, ci.level = ci.level)
     
     if (population.inference == FALSE) {  ## get interval from real var
       BsynCI <- BetaCI
       for (i in c(1,3,4)) BsynCI[,i] <- BsynCI[,i] + (syn.coef - real.coef)/real.se    #!GR ci by shifting
       BsynCI[,5] <- "synthetic"
     } else {
       BsynCI <- dfCI(syn.fit, Z = TRUE, name.Z = "Z.syn", model.name = "synthetic", ci.level = ci.level)
     }
     xlab = "Z value"
     title = paste0("Z values for fit to ",yvar)
   
   } else {
     BetaCI <- dfCI(real.fit, Z = FALSE, ci.level = ci.level)  #!GR 9/17  next line name.Z changed
     
     if (population.inference==FALSE) {  ## get interval from real var
       BsynCI <- BetaCI
       for (i in c(1,3,4)) BsynCI[,i] <- BsynCI[,i] + (syn.coef - real.coef)    #!GR ci by shifting
       BsynCI[,5] <- "synthetic"
     } else {
       BsynCI <- dfCI(syn.fit, Z = FALSE, name.Z = "syn.coef", model.name = "synthetic", ci.level = ci.level)
     }

     xlab = "Value"
     title = paste0("Coefficients for fit to ",yvar)
   }
   
   modelCI <- rbind.data.frame(BetaCI, BsynCI)
   rownames(modelCI) <- 1:nrow(modelCI)

   if(!plot.intercept) modelCI <- modelCI[modelCI$Coefficient!="(Intercept)",]

   CI.geom <- geom_errorbar(aes_string(ymin = "LowCI", ymax = "HighCI",
     color = "Model", linetype = "Model"), data = modelCI, width = 0, 
     lwd = lwd, lty = lty, position = position_dodge(width = dodge.height))

   point.geom <- geom_point(aes_string(#ymin = value, ymax = value,  #BN-03/02/2017 commented
     color = "Model", shape = "Model"), data = modelCI, 
     size = point.size, position = position_dodge(width = dodge.height))

   p <- ggplot(data = modelCI, aes_string(x = coefficient, y = value))
   p <- p + geom_hline(yintercept = 0, colour = "grey", linetype = 2, lwd = 1)
   p <- p + CI.geom + point.geom + labs(title = title, y = xlab)
   p <- p + scale_shape_manual(values = c(16:17), breaks = c("synthetic","observed")) +
            scale_colour_manual(values = lcol, breaks = c("synthetic","observed"))
   p <- p + coord_flip()
   #p <- p + theme_bw()
   # scale_colour_manual(values = rev(brewer.pal(3,"Blues")))
   # scale_colour_grey(start = 0, end = .6)
   p
 } else p <- NULL

 # detailed results
 res.obs <- real.fit$coefficients[,-4]
 colnames(res.obs) <- c("Beta","se(Beta)","Z")
 # res.syn  <- syn.fit$coefficients
 res.syn  <- syn.fit$coefficients[,-4] 
 res.syn  <- res.syn[order(match(rownames(res.syn), rownames(res.obs))), ]
 res.overlap <- compare.CI(res.syn, res.obs, ci.level = ci.level, intercept = TRUE)
 ncoef <- nrow(res.obs) 
 
 diff <- syn.coef - real.coef
 if(incomplete == TRUE){   ####  var cov for incomplete
   QB <- matrix(NA, m, length(object$mcoefavg))
     for(i in 1:m){
       QB[i,] <- object$mcoef[i,] - object$mcoefavg
     }
     lof.varcov <- t(QB) %*% QB/(m-1)/m
 }
 else if (object$proper == TRUE) lof.varcov <- real.varcov * (1+n/k)/m
 else lof.varcov <- real.varcov*n/k/m

 lack.of.fit <- t(diff) %*% solve(lof.varcov) %*% diff  #!GR multiply by m for combined test
 if (incomplete == FALSE) lof.pvalue  <- 1 - pchisq(lack.of.fit,ncoef)
 else {
   lack.of.fit <- lack.of.fit * (object$m-ncoef)/ncoef/(object$m-1) #!GR18/10/17 Hotellings T square
   lof.pvalue  <- 1 - pf(lack.of.fit, ncoef, object$m - ncoef)      #!GR18/10/17 changed degrees of freedom for denom
 }
 
 # res.diff <- compute.diff(res.syn, res.obs, intercept = TRUE, m = m, n = n, k = k)
 # res.diff <- compute.diff(syn.coef, res.obs, intercept = TRUE,m=object$m,k=k,n=dim(data)[1]) #gr
 res.diff <- compute.diff(syn.coef, res.obs, intercept = TRUE, lof.varcov = lof.varcov)
 
 # calculate summary measures
 mean.ci.overlap   <-  mean(res.overlap[,1])
 mean.abs.std.diff <-  mean(abs(res.diff[,1]))
 
  #mean.lof          <-  mean(lack.of.fit)
  #mean.lof.exp      <-  mean(lack.of.fit)/ncoef
  #std.lof           <- (mean(lack.of.fit) - ncoef)/sqrt(2*ncoef)

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
  
  CI <- qnorm(1- (1 - ci.level)/2)
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


###-----compare.CI---------------------------------------------------------
compare.CI <- function(synthetic, observed, ci.level, intercept, ...){
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
 obs.upper <- observed[i, 1] + (CI * observed[i, 2])
 obs.lower <- observed[i, 1] - (CI * observed[i, 2])
	  
## CI overlap
 overlap.lower <- max(obs.lower, syn.lower)       
 overlap.upper <- min(obs.upper, syn.upper)       
#   if(overlap.lower >= overlap.upper){            #! BN-06/01/17 to comment out ??
#	   CIoverlap[i, 1] <- 0                         #! BN-06/01/17 to comment out ??
#	 } else {                                       #! BN-06/01/17 to comment out ??
	   CIoverlap[i, 1] <- 0.5 * 
       (((overlap.upper - overlap.lower) / (obs.upper - obs.lower)) + 
		    ((overlap.upper - overlap.lower) / (syn.upper - syn.lower)))
#	 }                                              #! BN-06/01/17 to comment out ??
##Coef coverage
#coverage[i, 1] = (observed[i, 1] <= syn.upper & observed[i, 1] >= syn.lower)
##Z coverage
#coverage[i, 2] = (observed[i, 3] <= synthetic[i, 4] + (1.96 * synthetic[i, 5]) & 
#	 observed[i, 3] >= synthetic[i, 4] - (1.96 * synthetic[i, 5]))
 }

  if(intercept == FALSE) {
    CIoverlap = CIoverlap[-1, , drop = FALSE]
  }
  
	return(as.data.frame(CIoverlap))
}


###-----compute.diff-------------------------------------------------------
compute.diff <- function(synthetic, observed, intercept, lof.varcov, ...){  
  
  # tempStdDiff <- cbind((synthetic[,"B.syn"] - observed[,"B"]) / synthetic[,2]) 
  tempStdDiff <- (synthetic - observed[,1]) / observed[,2] 
  # pval <- round(2 * (1 - pnorm(sqrt(m*k/n)*abs(tempStdDiff))), 3)
  
  pval <- round(2 * (1 - pnorm(abs(synthetic - observed[,1])/sqrt(diag(lof.varcov)))), 3) 
  coefdiff <- data.frame("Std. coef diff" = tempStdDiff, 
                         "p value" = pval, check.names = FALSE)

  if(intercept == FALSE){
    coefdiff = coefdiff[-1, ]
  }
  return(coefdiff)
} 