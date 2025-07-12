###----------------------- multi.disclosure------------------------------
multi.disclosure <- function(object, data, ...) UseMethod("multi.disclosure")

###-----multi.disclosure.default-----------------------------
multi.disclosure.default <- function(object, ...)
  stop("No multi.disclosure method associated with class ", class(object), call. = FALSE)

###-----multi.disclosure.data.frame---multi.disclosure.list--------
multi.disclosure.data.frame <- multi.disclosure.list <- 
  function(object, data, cont.na = NULL, keys , targets = NULL, print.flag = TRUE, 
           denom_lim = 5, exclude_ov_denom_lim = FALSE,
           not.targetslev = NULL,  
           usetargetsNA = TRUE,  usekeysNA = TRUE, 
           exclude.keys = NULL, exclude.keylevs = NULL, exclude.targetlevs = NULL,
           ngroups_targets = NULL, ngroups_keys = NULL, 
           ident.meas = "repU", attrib.meas = "DiSCO",
           thresh_1way = c(50, 90),thresh_2way = c(4, 80), 
           digits = 2, plot = TRUE,  compare.synorig = TRUE, ...)
    
  {
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
if(compare.synorig){
      if (m ==1) adjust.data <- synorig.compare(object,data, print.flag = FALSE) else
      if (m > 1) adjust.data <- synorig.compare(object[[1]],data, print.flag = FALSE)
   
     if (!adjust.data$unchanged) {
        data <- adjust.data$orig
        object <- adjust.data$syn

        cat("Synthetic data or original or both adjusted with synorig.compare to try to make them comparable.\n")
        if (m > 1) cat("only first element of the list has been adjusted and will be used here\n")
        m <- 1 }
}

    object <- list(syn = object, m = 1, cont.na = cont.na) 
    class(object) <- "synds"

    res <- multi.disclosure.synds(object, data, keys = keys , targets = targets, 
           denom_lim = denom_lim,  exclude_ov_denom_lim = exclude_ov_denom_lim, 
           not.targetslev = not.targetslev, print.flag = print.flag ,  
           usetargetsNA = usetargetsNA, usekeysNA = usekeysNA, ngroups_targets = ngroups_targets,
           ngroups_keys = ngroups_keys, ident.meas = ident.meas, attrib.meas = attrib.meas, 
           thresh_1way = thresh_1way, thresh_2way = thresh_2way,
           digits = digits, plot = plot,...) 
    res$call <- match.call()
    return(res)
  }


###-----multi.disclosure.synds-------------------------
multi.disclosure.synds <-     function(object, data, keys , targets = NULL, print.flag = TRUE,  
                                         denom_lim = 5, exclude_ov_denom_lim = FALSE,
                                         not.targetslev = NULL, 
                                         usetargetsNA = TRUE,  usekeysNA = TRUE, 
                                         exclude.keys =NULL, exclude.keylevs = NULL, exclude.targetlevs = NULL,
                                         ngroups_targets = NULL, ngroups_keys = NULL, 
                                         ident.meas = "repU", attrib.meas = "DiSCO",
                                         thresh_1way = c(50, 90),thresh_2way = c(4, 80), 
                                         digits = 2, plot = TRUE, ...)
  

{
  ###----------------check input parameters ----
  if (!(is.data.frame(data)) )   stop("data  must be a data frame \n\n", call. = FALSE)
  data <- data.frame(data) ## in case it is table or a tibble
  if (!( inherits(object,"synds")))   stop(" object must be an object of class synds\n\n", call. = FALSE)

  if (object$m ==1) {
    names.syn <- names(object$syn)
  }
  else names.syn <- names(object$syn[[1]])

  if (is.numeric(keys) & !all(keys %in% 1:dim(data)[2])) stop("If keys are numeric they must be in range 1 to ",
                                                              dim(data)[2] , call. = FALSE  )
  if (is.numeric(keys)) keys <- names(data)[keys]
  
  if (is.null(targets) ) targets <- names.syn[!names.syn %in% keys]
  if (is.numeric(targets) & !all(targets %in% 1:dim(data)[2])) stop("If targets is numeric they must be in range 1 to ",
                                                                  dim(data)[2] , call. = FALSE  )
  if (is.numeric(targets)) targets <- names(data)[targets]

  # targets must be variables in  data and object$syn
  # keys must be a vector of variable names in data and in object$syn
  # target must not be in keys

  if (!(all(keys %in% names(data)) && all(keys %in% names.syn) 
        && all(targets %in% names(data)) & all(targets %in% names.syn )))
    stop("keys and targets must be variables in data and synthetic data. \n", call. = FALSE)
  if ( any(targets %in% keys) )
    stop("no targets can be in keys \n", call. = FALSE)
  if (any(names(data)  == "target"))
    stop("your data has a variables called 'target' please rename in original and synthetic data.\n\n", call. = FALSE) 
  
  if (is.null(not.targetslev) )  not.targetslev <- rep("", length(targets))
  else if (length(not.targetslev) != length(targets) || !is.character(not.targetslev)) 
    stop("not.targetslev must be a character vector of the same length ",length(targets)," as targets", call. = FALSE  )
  
  
  if (is.null(usekeysNA) ) usekeysNA <- TRUE
  
  # get keys in same order as in data
  keys <- names(data)[names(data) %in% keys]
  if (!(length(ident.meas) == 1  && length(attrib.meas) == 1 &&
        ident.meas %in% c("repU", "UiSiO","baseCAPd") && attrib.meas %in% c("DiSCO", "DiSDiO","DCAP")) )
        stop('ident.meas and attrib.meas must be single values
             from c("repU", "UiSiO""repU", "UiSiO","baseCAPd") and  c("DiSCO", "DiSDiO","DCAP") respectively\n\n', call. = FALSE)

  if (length(usekeysNA) == 1) usekeysNA <- rep(usekeysNA, length(keys))
  if (length(usetargetsNA) == 1) usetargetsNA <- rep(usetargetsNA, length(targets))
  
  if (length(usekeysNA) != length(keys)) stop("usekeysNA must be same length as keys", call. = FALSE)
  if (length(usetargetsNA) != length(targets)) stop("usetargetsNA must be same length as targets", call. = FALSE)

  if (!is.null(ngroups_targets)) {
  if (length(ngroups_targets) == 1) ngroups_targets <- rep(ngroups_targets, length(targets))
  if (!length(ngroups_targets) == length(targets))  stop("ngroups_targets must be same length as targets", call. = FALSE)
  }
  
Norig <- dim(data)[1]

 if (!is.null(exclude.keys)) {
   if (!is.list(exclude.keys) && ! length(exclude.keys) == length(targets))
   stop("exclude.keys must be a list of same length as targets", call. = FALSE)
   if (!is.list(exclude.keylevs) && ! length(exclude.keylevs) == length(targets))
     stop("exclude.keylevs must be a list of same length as targets", call. = FALSE)
   if (!is.list(exclude.targetlevs) && ! length(exclude.targetlevs) == length(targets))
     stop("exclude.targetlevs must be a list of same length as targets", call. = FALSE)

 }
###-------------------------- create output data frames-------------------
 ident.orig <- ident.syn <- attrib.orig <- attrib.syn <- n2way <- rep(NA, length(targets))
 check1 <- check2 <- rep("", length(targets))
 output.list <-  as.list(1:length(targets))
 names(output.list) <- targets



  for (i in 1:length(targets)) {
    if (print.flag) cat("------------------",i,targets[i],"-------------------","\n")
    if (not.targetslev[i] =="") not.targetlev = NULL
    else not.targetlev <-not.targetslev[i]
    ttt <-disclosure(object, data, keys = keys, target = targets[i],  denom_lim = denom_lim,
            exclude_ov_denom_lim = exclude_ov_denom_lim, print.flag = print.flag, digits =digits,
            usetargetNA = usetargetsNA[i], usekeysNA = usekeysNA, not.targetlev = not.targetslev[i],
            exclude.keys =exclude.keys[[i]], exclude.keylevs = exclude.targetlevs[[i]],
            exclude.targetlevs = exclude.targetlevs[[i]], ngroups_target = ngroups_targets[i],
            ngroups_keys = ngroups_keys, thresh_1way =thresh_1way ,thresh_2way = thresh_2way)
    
    class(ttt) <- "disclosure"
    output.list[[i]] <- ttt

    if (ident.meas == "repU")  ident.syn[i] = mean(ttt$ident$repU[1])
    if (ident.meas == "UiSiO") ident.syn[i] = mean(ttt$ident$UiSiO)
    ident.orig[i] <- mean(ttt$ident$UiO)

    if (attrib.meas == "DiSCO" ){
       attrib.orig[i] <- mean(ttt$attrib$Dorig)
       attrib.syn[i] <- mean(ttt$attrib$DiSCO )}
    if (attrib.meas == "DiSDiO" ){
      attrib.orig[i] <- mean(ttt$attrib$DiO)
        attrib.syn[i] <- mean(ttt$attrib$DiSDiO )}
    if (attrib.meas == "DCAP" ){
      attrib.orig[i] <- mean(ttt$allCAPs$CAPd)
      attrib.syn[i] <-  mean(ttt$allCAPs$DCAP)}

    check1[i] <- ttt$check1
    check2[i] <- ttt$check2[1]
    if (ttt$check2[1] == "") n2way[i] <- 0
    else n2way[i] <- length(ttt$check2)
  }

##------------------------------ end of i loop-------------------------
    result <- data.frame(attrib.orig , attrib.syn, check1 = check1,
                    Npairs = n2way, check2 = check2)
 
    dimnames(result)[[1]] <- targets
    names(output.list) <- targets
    
###-----ntoc---------------------------------------------------------------
    # to make labels for variables of constant length from an integer
    ntoc <- function(x)
    {
      nch <- nchar(max(x))
      res <- as.character(x)
      res <- str_pad(res, nch, "left", "0")
      return(res)
    }
    ##-------------------- end of ntoc-------------------------------
    result <- data.frame(result)
    result <- result[order(result[,2]),]

dimnames(result)[[1]] <- paste(ntoc(1:dim(result)[1]),dimnames(result)[[1]] )
    #changed <- rep(FALSE, dim(result)[1])
 
            dimnames(result)[[1]][result[,5] != "" & result[,3] == ""] <- 
      paste(dimnames(result)[[1]][result[,5] != "" & result[,3] == ""],"2way checks")
            dimnames(result)[[1]][result[,3] != "" & result[,5] == ""] <- 
      paste(dimnames(result)[[1]][result[,3] != "" & result[,5] == ""],"1way checks")
            dimnames(result)[[1]][result[,3] != "" & result[,5] != ""] <- 
      paste(dimnames(result)[[1]][result[,3] != "" & result[,5] != ""],"1&2way checks")
          
   ###---------------------- now create data for plot--------------------

      toplot <- result[rep( 1:dim(result)[[1]], rep(2,dim(result)[[1]]) ),]
      for ( i in 1:(dim(result)[1]) )  toplot[i*2,1] <- toplot[i*2,2]
   
      toplot$measure <- dimnames(result)[[2]][1]
      for ( i in 1:(dim(result)[1]) )  toplot$measure[i*2] <- dimnames(result)[[2]][2]
      toplot <- data.frame(name = rep(dimnames(result)[[1]],rep(2, dim(result)[1])), toplot)
      dimnames(toplot)[[1]] <- 1:dim(toplot)[1]
  
      names(toplot)[2] <- "VALUE"
      
      if (attrib.meas %in% c("DCAP","TCAP")) oratt = "CAPd"
      else oratt <- "Dorig"

      attrib.plot <- ggplot(toplot) + 
       geom_point(data = toplot, size=5, aes(colour = .data$measure, shape = .data$measure, x=.data$VALUE, y=.data$name)) +
       xlim(0,100)  +
        labs(x = "Disclosure measure", y = "", 
         title = "Comparison of attribute disclosure measures",
         subtitle= paste( attrib.meas,"for synthetic data  to",oratt ,"for original data.")) +
        geom_line(data = toplot, mapping = aes(x=.data$VALUE, y=.data$name), arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"))
#cat("line 207------------------------------------------------\n")

       res <- list( attrib.table = result, attrib.plot = attrib.plot, keys = keys,
                         ident.orig = ident.orig, ident.syn = ident.syn, Norig = Norig,
                         denom_lim = denom_lim, exclude_ov_denom_lim = exclude_ov_denom_lim,
                         digits = digits, usetargetsNA = usetargetsNA, usekeysNA = usekeysNA, 
                         ident.meas = ident.meas, attrib.meas = attrib.meas, m = object$m,
                         plot = plot, output.list = unclass(output.list))
       
       class(res) <- "multi.disclosure"
       return(res)
    }        
  
  
###-----print.multi.disclosure-----------------------------------------------
  print.multi.disclosure <- function(x,  digits = NULL,  
                                   plot = NULL, to.print = c("ident","attrib"), ...) {
    
     if (is.null(digits)) digits <- x$digits
    if (is.null(plot)) plot <- x$plot
  
    if ("ident" %in% to.print) {
      cat("Disclosure risk for",x$Norig,"records in the original data\n")
      if (x$m >1) cat("Results are averaged over" ,x$m,"syntheses")
      cat("\nIdentity disclosure measures\n")
      cat("from keys:", x$keys,"\n")
      cat("For original  ( UiO ) ", round(x$ident.orig[1],x$digits),"%\n")
      cat("For synthetic (", x$ident.meas,")",round(x$ident.syn[1],x$digits),"%.\n")
    }
    if ("attrib" %in% to.print) {
      cat("\nTable of attribute disclosure measures for",x$keys,"\n")
      if (x$attrib.meas %in% c("DCAP","TCAP")) oratt = "CAPd"
      else oratt <- "Dorig"

      cat("Original measure is ",oratt,"and synthetic measure is",x$attrib.meas,"\n")
      cat("Variables Ordered by synthetic disclosure measure\n\n")
      toprint <- x$attrib.table
      for (i in 1:2) toprint[,i] <- round(toprint[,i], x$digits)
      print(toprint)
  if (plot & !is.null(x$attrib.plot)) {
    print(x$attrib.plot)
  }
    }
 
  if (any(x$attrib$check1 != "")) cat("\n\nCheck check_1way in output from disclosure()\n",
   "for targets with check1 \n")
    if (any(x$attrib$check2 != "")){ cat("\n\nCheck check_2way in output from disclosure()\n",
                                        "for targets with check2 \n")
      cat("Only the first two-way pair is given here when npairs >1\n")
    }
  invisible(x)
}


  