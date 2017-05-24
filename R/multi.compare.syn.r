###-----multi.compare------------------------------------------------------
multi.compare <- function (object, data, var = NULL, by = NULL, msel = NULL, 
  barplot.position = "fill", cont.type = "hist", y.hist = "count", 
  boxplot.point = TRUE, binwidth = NULL, ...) {
 
# CHECKS
#---  
 if (is.null(data)) stop("Requires parameter data to give name of the real data.\n", call. = FALSE)
 if (!is.data.frame(data)) stop("Argument data must be a data frame.\n", call. = FALSE)
 if (!is.null(msel) & !all(msel %in% (1:object$m))) stop("Invalid synthesis number(s).", call. = FALSE)

 if (length(var)!=1){
   cat("\nParameter var set to ",var,", should be a single variable.\nOnly first variable used.\n") 
   var <- var[1]
 }
 if (!(var %in% names(data))) stop("\nParameter var set to ",var,
   ", should be the name of a variable in data\n", call. = FALSE)
 if (!all(by %in% names(data))) stop("\nParameter by set to ", by,
   ", should all be names of variables in data.\n", call.=FALSE)
 if (any(var %in% by)){  
   cat("\nParameter var: ",var,"included in parameter by: ", by,
   " now removed from by.\n")
   by <- by[-match(var,by)]
 }
#--- 
 
 for (i in match(by,names(data))) data[,i] <- addNA(data[,i], ifany = TRUE)     # numeric will be changed to factor (create seperate data for this?)
 bytable <- xtabs(as.formula(paste("~", paste(by, collapse = " + "))), 
   data = data, exclude = NULL, na.action = na.pass) 
 cat("\nPlots of ",var," by ",by,"\nNumbers in each plot (observed data):\n\n")
 print(bytable)
 nplots <- length(bytable) 
 if (nplots > 100) cat("\nCAUTION: You have ",nplots," sections in your plot.\n")
    
 add <- NULL 
 
 # data prep
 #----
 if (is.null(msel) & object$m > 1) msel <- 1:object$m
 if (object$m == 1) {
   synall <- cbind(object$syn, source = "syn")
 } else if (length(msel) == 1) {
   synall <- cbind(object$syn[[msel]], source = "syn")
 } else if (object$m > 1 & length(msel) > 1) {
   synall <- Map(cbind, object$syn[msel], source = paste0("syn=", msel))
   synall <- do.call(rbind, synall)
 }
 obssyn <- rbind(cbind(data, source = "obs"), synall)
 #----
 
 ..count.. <- ..density.. <- NULL
 
 if (is.numeric(data[, var])) {
   if (cont.type == "hist") {
     p <- ggplot(data = obssyn, aes(x = eval(parse(text = var))))
     plabs <- labs(x = var, fill = "")
     if (y.hist == "count") {
       ptype <- geom_histogram(aes(y = ..count.., fill = source), 
         position = "dodge", binwidth = binwidth)
     } else if (y.hist == "density"){
       ptype <- geom_histogram(aes(y = ..density.., fill = source), 
         position = "dodge", binwidth = binwidth)
     }
   } else if (cont.type == "boxplot") {
     p <- ggplot(data = obssyn, aes(x = source, y = eval(parse(text = var))))
     ptype <- geom_boxplot(aes(colour = source), alpha = 0.7)
     plabs <- labs(y = var, x = "", colour = "")
     if (boxplot.point) add <- geom_jitter(size = 0.2, alpha = 0.2)  #!BN-size was 0.5
   }
 } else {
   p <- ggplot(data = obssyn, aes(x = source))
   ptype <- geom_bar(aes(fill = eval(parse(text = var))), position = barplot.position)
   plabs <- labs(x = "", fill = var)
 }
 
 if (length(by) == 1){
   playaout <- facet_wrap(by)   # nrow = 1
 } else if (length(by) > 1) {
   form <- paste(by[1], "~", paste0(by[-1], collapse = "+"))
   playaout <- facet_grid(eval(parse(text = form)))
 }
 
 p <- p + add + ptype + playaout + plabs + scale_fill_brewer(palette="Set1") + 
      scale_colour_brewer(palette="Set1")
 return(p)
}
