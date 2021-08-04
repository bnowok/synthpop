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
 for (i in 1:length(varnos)) {
   grpd <- group_var(data[, varnos[i]], cont.na = cna[[i]], n = catgroups[i], style = style.groups)
   data[, varnos[i]] <- grpd$x
   breaks[[i]] <- grpd$breaks
   levels[[i]] <- grpd$levels
 } 

 return(list(data = data, breaks = breaks, 
             levels = levels, orig = orig, cont.na = cna, 
             numtocat = numtocat, ind = varnos))
}


###-----group_var----------------------------------------------------------

group_var <- function(x,  n = 5, style = "quantile", cont.na = NA) {
  # categorise one continous variable into groups
  if (!is.numeric(x)) stop("x must be numeric.\n", call. = FALSE)
  # select non-missing(nm) values 
  xnm <- x[!(x %in% cont.na) & !is.na(x)]
  my_breaks   <- unique(classIntervals(xnm, n = n, style = style)$brks)
  xnm_grouped <- cut(xnm, breaks = my_breaks, dig.lab = 8, 
                     right = FALSE, include.lowest = TRUE)
  my_levels   <- c(levels(xnm_grouped), cont.na[!is.na(cont.na)])
  # assign groupings to non-missing data
  x[!(x %in% cont.na) & !is.na(x)] <- as.character(xnm_grouped)
  x <- as.factor(x)
  list(x = x, breaks = my_breaks, levels = my_levels)  
}  

