###-----numtocat.syn-------------------------------------------------------
# group numeric variables in a data frame into categories

numtocat.syn <- function(data, numtocat = NULL,
                         print.flag = TRUE, cont.na = NULL, 
                         catgroups = 5, style.groups = "quantile")
{
 if (!is.data.frame(data)) stop("Data must be a data.frame.\n", call. = FALSE)
 varnames <- names(data)

 # checks on numtocat
 if (!is.null(numtocat)) {
   if (is.numeric(numtocat)) {
     if (!all(numtocat %in% 1:length(data))) stop("Column indices must be between 1 and ", 
                                                  length(data), ".", sep = "", call. = FALSE)  
     varnos <- numtocat
     numtocat <- names(data)[varnos]
   } else {
     if (!all(numtocat %in% varnames)) stop("Variable(s) ", 
                                            paste(numtocat[!numtocat %in% varnames], collapse = ", "),
                                            " not in data.\n", sep = "", call. = FALSE)
     varnos <- match(numtocat,varnames)
   }
   vclass <- sapply(data[, varnos, drop = FALSE], is.numeric)
   if (!all(vclass)) stop("Variable(s) in numtocat (", 
                          paste(numtocat[!vclass], collapse = ", "), 
                          ") not numeric.\n", sep = "", call. = FALSE)
 } else { 
   ## if NULL use all numeric variables
   varnos   <- (1:length(data))[sapply(data, is.numeric)]
   numtocat <- names(data)[varnos]
 }

 # checks on catgroups
 if (length(catgroups) == 1) catgroups <- rep(catgroups,length(numtocat))
 else if (length(catgroups) != length(numtocat)) stop("catgroups must be a single number or a vector of the same length as numtocat.\n", call. = FALSE)

 if (any(catgroups < 2)) stop("catgroups must all be > 1.", call. = FALSE)
 # checks on cont.na
 if (!is.null(cont.na)) {
   if (!is.list(cont.na)) stop("cont.na must be a  list.\n", call. = FALSE)
   if (!all(names(cont.na) %in% numtocat)) stop("Variable(s): ",
        paste(names(cont.na)[!names(cont.na) %in% numtocat],collapse = ", "),
        " in cont.na not in the variables being grouped.\n", call. = FALSE)
   cna <- as.list(rep(NA,length(numtocat)))
   for (i in 1:length(cont.na)) cna[[match(names(cont.na)[i],numtocat)]] <- cont.na[[i]]
 } else {
   cna <- as.list(rep(NA,length(numtocat)))
 }
 names(cna) <- numtocat

 if (print.flag == TRUE) cat("Variable(s) ", paste(numtocat, collapse = ", "),
                             " grouped into categories.\n", sep = "")
 breaks <- vector("list", length(varnos)); names(breaks) <- numtocat
 levels <- vector("list", length(varnos)); names(levels) <- numtocat
 orig <- data[, varnos, drop = FALSE]
 names(orig) <- paste("orig", names(orig), sep = ".")
 
 ###-----group_var----------------------------------------------------------
 
 group_var <- function(x,  n = 5, style = style.groups, cont.na = NA) {
   # categorise one continous variable into groups
   if (!is.numeric(x)) stop("x must be numeric.\n", call. = FALSE)
   # select non-missing(nm) values 
   xnm <- x[!(x %in% cont.na) & !is.na(x)]
   my_breaks   <- unique(classIntervals(xnm, n = n, style = style, breaks = breaks)$brks)
   xnm_grouped <- cut(xnm, breaks = my_breaks, dig.lab = 8, 
                      right = FALSE, include.lowest = TRUE)
   my_levels   <- c(levels(xnm_grouped), cont.na[!is.na(cont.na)])
   # assign groupings to non-missing data
   x[!(x %in% cont.na) & !is.na(x)] <- as.character(xnm_grouped)
   x <- as.factor(x)
   list(x = x, breaks = my_breaks, levels = my_levels)  
 } 
 
 ###------------------ end of group_var------------------------------------
 for (i in 1:length(varnos)) {
   grpd <- group_var(data[, varnos[i]], cont.na = cna[[i]], n = catgroups[i], style = style.groups)
   breaks[[i]] <- grpd$breaks
   nb <- length(breaks[[i]])

   if (length(breaks[[i]]) <= 3 ) {
     grpd <- group_var(data[, varnos[i]], cont.na = cna[[i]], n = catgroups[i], style = "equal")
     breaks[[i]] <- grpd$breaks

     if (length(breaks[[i]]) <= 3 ) {
       cat("Only",length(breaks[[i]]) - 1,"groups produced for", numtocat[i],"even after changing method.\n")
       cat("Check data\n\n")
     }
     else {
      cat("Grouping changed from 'quantile' to  'equal' in function numtocat.syn for",
               numtocat[i],"because only",nb - 1," groups produced\n")
     }
   }
   data[, varnos[i]] <- grpd$x
   breaks[[i]] <- grpd$breaks
   levels[[i]] <- grpd$levels
 } 
 return(list(data = data, breaks = breaks, 
             levels = levels, orig = orig, cont.na = cna, 
             numtocat = numtocat, ind = varnos))
}

 


###-----merge levels of factors.syn-------------------------------------------------------
# merge factor levels with fewer than minsize values

mergelevels.syn <- function(data, vars = NULL, newlabel = FALSE, addNA = FALSE,
                         print.flag = FALSE, minsize = 10, merge.byhand = FALSE,
                         merge.details = NULL)
{
  if (!is.data.frame(data)) stop("Data must be a data.frame.\n", call. = FALSE)
  varnames <- names(data)
 
  print(ls())
if (merge.byhand == FALSE) {
  # checks on levels
  if (!is.null(vars)) {
    if (is.numeric(vars)) {
      if (!all(vars %in% 1:length(data))) stop("Column indices must be between 1 and ", 
                                               length(data), ".", sep = "", call. = FALSE)  
      varnos <- vars
      vars <- names(data)[varnos]
    } else {
      if (!all(vars %in% varnames)) stop("Variable(s) ", 
                                         paste(vars[!vars %in% varnames], collapse = ", "),
                                         " not in data.\n", sep = "", call. = FALSE)
      varnos <- match(vars,varnames)
    }
    
    vclass <- sapply(data[, varnos, drop = FALSE], is.factor)
    if (!all(vclass)) stop("Variable(s) (", 
                           paste(vars[!vclass], collapse = ", "), 
                           ") in mergelevels.syn not factor(s).\n", sep = "", call. = FALSE)
  } else { 
    ## if NULL use all factor variables
    varnos   <- (1:length(data))[sapply(data, is.factor)]
    vars <- names(data)[varnos]
  }
  if (addNA) for (i in 1:length(varnos)) data[,varnos[i]] <- addNA(data[,varnos[i]])
  # find factors that need fixing
  tofix <- rep(FALSE, length(varnos))
  for (i in 1:length(varnos)) {
    tb <- table(data[,varnos[i]])
 
    if (min(tb) < minsize ){ 
      tofix[i] <- TRUE
      if (print.flag) cat("Categories merged for ",vars[i],"\n")
      if (print.flag) {cat("Original data\n") ;print(tb)}
      data[,varnos[i]] <- fct_lump_min(data[,varnos[i]],min = minsize)
      lenlev <- length(levels(data[,varnos[i]]))
      if ( table(data[,varnos[i]])[lenlev] < minsize ) {
        ########### agregate with next lowest
        levnxtlow <- levels(data[,varnos[i]])[lenlev -1 ]
        data[,varnos[i]] <- fct_lump_n(data[,varnos[i]],n = lenlev - 2)
        if (newlabel == TRUE) levels(data[,varnos[i]])[lenlev -1 ] <- paste(levnxtlow,names(tb)[tb < minsize], collapse = " ")
      } else {
        ######################## relabel combined category
        if (newlabel  == TRUE) levels(data[,varnos[i]])[lenlev] <- paste(names(tb)[tb < minsize], collapse = " ")
      }
      if (print.flag) {cat("Merged categories data\n") ;print(table(data[,varnos[i]]))}
      }
    }
    if (sum(tofix) >0) cat("The following" ,sum(tofix),"variables have categories merged:\n",vars[tofix],"\n")
    else cat("No variables had categories with fewer than ", minsize,"counts\n")
} else {
  if ( is.null(merge.details) ) stop("Need to supply merge.details if merge,byhand is TRUE\n", call. = FALSE)
  if (!all(names(merge.details) %in% names(data))) stop("Some names of the list merge.details is not a variable in data:\n",
        names(merge.details)[!(names(merge.details) %in% names(data))],"\n", call. = FALSE)
  vars <-  names(merge.details)
  varnos   <- (1:length(data))[match(names(merge.details),names(data))] 
  if (addNA) for (i in 1:length(varnos)) data[,varnos[i]] <- addNA(data[,varnos[i]])

  for ( i in 1:length(varnos)) { 
    if (print.flag) {
      cat("Original data\n")
      print(table(data[,varnos[i]], useNA = "ifany"))
      }
    templevs <- levels(data[,varnos[i]] )
                      
    if (addNA & any(is.na(templevs)))   templevs[is.na(templevs)] <- "NA"
    if (  !all( merge.details[[i]][-1] %in% templevs ) )   stop("Levels in merge.details for ", 
            vars[i], " are not in its levels problem entries shown here:\n",
            merge.details[[i]][-1][!(merge.details[[i]][-1] %in% replace)],"\n", call. = FALSE)
    levels(data[,varnos[i]])[match(merge.details[[i]][-1],templevs)] <-
                                                               merge.details[[i]][1]
    if (print.flag) {
      cat("Merged cateories data\n")
      print(table(data[,varnos[i]], useNA = "ifany"))
    }
  }
}
  return(data)
}


#
