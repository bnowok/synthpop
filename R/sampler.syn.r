#-----------------------------sampler.syn-------------------------------
# The sampler controls the generation of conditional distributions
# This function is called by syn()

sampler.syn <- function(p, data, m, syn, visit.sequence,
                        rules, rvalues, event, proper,
                        print.flag, k, pred.not.syn, 
                        models, numtocat,  ...)
{
 #--- Assign optional parameters (...) to appropriate synthesising function   
 dots  <- as.list(substitute(list(...)))[-1L]         
 meth.with.opt <- paste(c("cart", "cartbboot", "ctree", "survctree", "polyreg",
                          "polr", "rf", "ranger", "bag", "ipf", "catall"), collapse = "\\.|")
 meth.check <- grep(meth.with.opt, names(dots), value = TRUE)
 args.err <- !(names(dots) %in% meth.check)
 if (any(args.err)) stop("Unknown optional parameter(s): ", 
   paste(names(dots)[args.err], collapse = ", "),
   "\nNote that they have to be method specific, e.g. 'ctree.minbucket' and NOT 'minbucket'\n", 
   call. = FALSE)
 if (length(dots) == 0) {
   mth.args <- NULL
 } else {  
   #mth.args.dots <- strsplit(names(dots), "\\.")
   mth.args.dots <- regmatches(names(dots), regexpr("\\.", names(dots)), invert = TRUE)
   mth.dots  <- unique(lapply(mth.args.dots, "[[", 1))
   args.dots <- lapply(mth.args.dots, "[[", -1)
   mth.args  <- setNames(vector("list", length(mth.dots)), unlist(mth.dots))
  
   for (i in 1:length(mth.dots)) { 
     ind <- grep(mth.dots[[i]], names(dots))
     mth.args[[i]] <- setNames(dots[ind], args.dots[ind])
   } 
 } 
 #---
  
 fits <- NULL 
    
 if (m > 0) {
   if (models) fits <- rep(list(setNames(vector("list", length(p$method)),
                                names(p$method))), m) 
   for (i in 1:m) {  # Synthesising loop
     if (print.flag & m > 1) cat("Synthesis number ", i, 
                         "\n--------------------\n", sep = "")  
     if (print.flag & m == 1) cat("Synthesis\n-----------\n", sep = "")  
     
     # Code for methods that take more than one variable together: ipf & catall      
     #--------------------------------------------------------------------------
     rest.visit.sequence <- p$visit.sequence  # when no grouped methods used
          
     if (any(p$method %in% c("catall", "ipf"))) {
       ordmethod <- p$method[p$visit.sequence]
       grind <- (1:length(p$visit.sequence))[ordmethod %in% ordmethod[1]]
            
       ## to reorder any dummies for grouped variables
       if (any(names(p$visit.sequence) %in% 
           paste(names(p$visit.sequence[grind]), "1", sep = "."))) {  
         dumind <- (1:length(p$visit.sequence))[names(p$visit.sequence) %in% 
                    paste(names(p$visit.sequence[grind]), "1", sep = ".")]
         othind <- (1:length(p$visit.sequence))[-c(grind, dumind)]
         p$visit.sequence <- p$visit.sequence[c(grind, dumind, othind)]
         ordmethod <- p$method[p$visit.sequence]
       }
            
       grouped <- p$visit.sequence[ordmethod %in% ordmethod[1]]
       
       if (print.flag == TRUE) {
         if (length(rest.visit.sequence) > 0  && 
             ncol(data) - length(numtocat) > length(grouped)) {
           cat("First ", length(grouped), " variables (", 
               paste(names(grouped), collapse = ", "),
               ") synthesised together by method '", ordmethod[1], "'\n", sep = "")
         } else {
           cat("All ", length(grouped), 
               " variables in the data synthesised together by method '", 
               ordmethod[1], "'\n", sep = "")
         }   
       }   
       x <- p$data[, grouped]
       if (!(ordmethod[1] %in% names(mth.args))) fun.args <- NULL else
       fun.args  <- mth.args[[ordmethod[1]]]
       f <- paste("syn", ordmethod[1], sep = ".")
       synfun <- do.call(f, args = c(list(x = x, k = k, 
                                          proper = proper), fun.args))
       p$syn[, grouped] <- synfun$res
       if (models) {
         fits[[i]][[grouped[1]]] <- synfun$fit
         for (j in 2:length(grouped)) fits[[i]][[grouped[j]]] <- 
           paste("See first in group:", names(grouped)[1])
       }
       
       rest.visit.sequence <- p$visit.sequence[-(1:length(grouped))]
       if (length(rest.visit.sequence) > 0 & print.flag & 
           ncol(data) - length(numtocat) > length(grouped)) cat("\nRemaining variables:\n")
     }
     
     # Other variables 
     #--------------------------------------------------------------------------
     if (length(rest.visit.sequence) > 0) {           
       prcount <- 0 # to get new lines come out right
       for (j in rest.visit.sequence) {

         theMethod <- p$method[j]
         # get optional parameters for theMethod if they are provided
         if (!(theMethod %in% names(mth.args))) fun.args <- NULL else            
           fun.args  <- mth.args[[theMethod]]                                     
         
         vname <- dimnames(p$data)[[2]][j]
         if (print.flag & theMethod != "dummy"  
             & j <= (ncol(data) - length(numtocat))) { 
           cat(" ", vname, sep = "")
           prcount <- prcount + 1
         }  
         if (print.flag & prcount %% 10 == 0 & 
             j <= (ncol(data) - length(numtocat))) cat("\n")                                                  

         ya  <-  1:nrow(p$data) 
         ypa <- 1:k    
                     
         # ya = yavailable, ym = ymissing                                            
         if (any(p$rules[[j]] != "")) {
           com.rules  <- paste(p$rules[[j]], collapse = " | ")
           evalrul.y  <- with(p$data,eval(parse(text = com.rules)))
           ym         <- which(evalrul.y == TRUE & !is.na(evalrul.y))
           ya         <- setdiff(1:nrow(p$data), ym)                                  
           evalrul.yp <- with(p$syn,eval(parse(text = com.rules)))         
           ypm        <- which(evalrul.yp == TRUE & !is.na(evalrul.yp))        
           ypa        <- setdiff(1:nrow(p$syn), ypm)       
         }                                                                       
           
         # != "", != "dummy", != "passive"
         if (theMethod != "" & (!is.passive(theMethod)) & theMethod != "dummy" ) {
           
           if (theMethod %in% c("sample", "sample.proper", "constant")) {
             
             y   <- p$data[ya, j]
             if (is.factor(y)) y <- y[, drop = TRUE]
             xp  <- length(ypa)
             x   <- length(ya)
             nam <- vname
             f   <- paste("syn", theMethod, sep = ".")
             if (theMethod == "constant") {
               synfun <- do.call(f, args = list(y = y, xp = xp, ...))    
             } else if (is.numeric(y)) {
               synfun <- do.call(f, args = list(y = y, xp = xp,       
                 smoothing = p$smoothing[j], cont.na = p$cont.na[[j]], 
                 proper = proper, ...)) 
             } else {
               synfun <- do.call(f, args = list(y = y, xp = xp, 
                proper = proper, ...)) 
             }
             p$syn[ypa, j]  <- synfun$res
             if (models) fits[[i]][[j]] <- synfun$fit
            
           } else {
          
             x    <- p$data[ya, p$predictor.matrix[j, ] == 1, drop = FALSE]
             xp   <- p$syn[ypa, p$predictor.matrix[j, ] == 1, drop = FALSE]
             y    <- p$data[ya, j]
             if (is.factor(y)) y <- y[, drop = TRUE]
             nam  <- vname
             f    <- paste("syn", theMethod, sep = ".") 
             if (!theMethod %in% c("collinear", "nested")) {   # nested needs added to allow missing values 
             #if(theMethod!="collinear"){                  
               keep <- remove.lindep.syn(x, y, ...)
               x    <- x[, keep, drop = FALSE]
               xp   <- xp[, keep, drop = FALSE]
             }                                            
             if (theMethod == "survctree") {
               if (p$event[j] == -1) yevent <- rep(1,length(y))                   
               else yevent  <- p$data[ya,p$event[j]]
               survres      <- do.call(f, args = c(list(y = y, yevent = yevent,
                                       x = x, xp = xp, proper = proper), 
                                       fun.args))
               p$syn[ypa, j] <- survres[[1]]                                # synthetic data survival goes to p$syn
               if (p$event[j] != -1) p$syn[ypa,p$event[j]] <- survres[[2]]  # synthetic data event goes to p$syn
               if (models) fits[[i]][[j]] <- survres$fit 
             } else if (theMethod == "logreg" & p$denom[j] != 0) {                 
               synfun <- do.call(f, args = list(y = y, x = x, xp = xp,
                                 denom = p$data[ya,p$denom[j]], 
                                 denomp = p$syn[ypa, p$denom[j]],       
                                 proper = proper, ...))
               p$syn[ypa, j] <- synfun$res
               if (models) fits[[i]][[j]] <- synfun$fit           
             
             } else if (theMethod == "nested") {
               if (is.numeric(y)) {
                synfun <- do.call(f, args = c(list(y = y, x = x, xp = xp,
                  smoothing = p$smoothing[j], cont.na = p$cont.na[[j]], 
                  proper = proper), fun.args))
               } else {
                synfun <- do.call(f, args = c(list(y = y, x = x, xp = xp,
                  proper = proper), fun.args))
               }
               p$syn[ypa, j] <- synfun$res
               if (models) fits[[i]][[j]] <- synfun$fit
               
             } else {
               if (is.numeric(y)) {
                synfun <- do.call(f, args = c(list(y = y, x = x, xp = xp,
                  smoothing = p$smoothing[j],
                  proper = proper), fun.args))
               } else {
                synfun <- do.call(f, args = c(list(y = y, x = x, xp = xp,
                  proper = proper), fun.args))
               }
               p$syn[ypa, j] <- synfun$res
               if (models) fits[[i]][[j]] <- synfun$fit
             }
           }

           if (any(p$rules[[j]] != "")) {
             if (length(p$rules[[j]]) == 1 & length(ypm) > 0) {
               p$syn[ypm,j] <- p$rvalues[[j]] 
             } else {
               for (r in 1:length(p$rules[[j]])) {
                 revalrul.yp  <- with(p$syn,eval(parse(text = p$rules[[j]][r])))  
                 rypm <- which(revalrul.yp == TRUE & !is.na(revalrul.yp))
                 if (length(rypm) > 0) p$syn[rypm,j] <- p$rvalues[[j]][r]
               }
             }                 
           }  
         } # end of !="", !="dummy", !="passive"

         else if (is.passive(theMethod)) {
           class0 <- class(p$syn[,j])
           synfun <- syn.passive(data = p$syn, func = theMethod) 

           if (is.factor(synfun$res[[1]]) & any(is.na(synfun$res[[1]]))) {
             synfun$res[[1]] <- addNA(synfun$res[[1]], ifany = TRUE)
             levels(synfun$res[[1]])[is.na(levels(synfun$res[[1]]))] <- "NAtemp"
           }

           p$syn[, j] <- synfun$res
           class(p$syn[,j]) <- class0
           if (models) fits[[i]][[j]] <- synfun$fit 
         }

         else if (theMethod == "dummy") {    # replace dummy variables in p$syn
           # getting dummy values from a synthesised categorical variable
           cat.columns <- p$syn[, p$categories[j, 4]]  # this is the single column with the data for which this is the dummy
           model.frame(~cat.columns - 1, data = p$syn) 
           p$syn[, (j:(j + p$categories[p$categories[j, 4], 2] - 1))] <- # replaces all the dummies for this variable with
           matrix((model.matrix(~cat.columns - 1)[, -1]),               # dummies calculated from the synthesised data
                   ncol = p$categories[p$categories[j, 4], 2],
                   nrow = nrow(p$syn))
           p$syn[,j] <- as.numeric(p$syn[, j])
           remove("cat.columns")
           if (models) fits[[i]][[j]] <- "dummy"    
         }
       } # end j loop 
     } # end other variables
     if (print.flag) cat("\n")  
     
     #if (k==dim(data)[1]) syn[[i]] <- p$syn[,1:dim(data)[2]]
     #else syn[[i]] <- p$syn[sample(1:dim(data)[1],k),1:dim(data)[2]]

     syn[[i]] <- p$syn[, 1:dim(data)[2], drop = FALSE]
     nms <- names(data)
     # exclude unsynthesised if drop.pred.only set to true
     if (sum(pred.not.syn ) > 0) {
       syn[[i]] <- syn[[i]][, !pred.not.syn]
       nms <- nms[!pred.not.syn]  # GR save names to use below if data just one column
     }
     # Prevent a single character column being changed to a factor
     chgetochar <- (sum(!pred.not.syn) == 1 & any(class(syn[[i]][, 1]) == "character"))       
     syn[[i]] <- as.data.frame(syn[[i]])
     if (chgetochar) {
       syn[[i]][, 1] <- as.character(syn[[i]][, 1])
       names(syn[[i]]) <- nms
     }
   
     #turn NA level in factors / logical to missing NA's
     # and remove contrasts 
     for (j in (1:ncol(syn[[i]]))) {
       if (is.factor(syn[[i]][,j])) {                                    #!BN-20/04/16
         if ("NAlogical" %in% levels(syn[[i]][,j])) {
           levels(syn[[i]][,j])[levels(syn[[i]][,j]) == "NAlogical"] <- NA
           syn[[i]][,j] <- as.logical(syn[[i]][,j])
         } else {                                                            
        # syn[[i]][,j] <- factor(syn[[i]][,j],exclude=NA,levels=levels(syn[[i]][,j]))
           levels(syn[[i]][,j])[levels(syn[[i]][,j]) == "NAtemp"] <- NA  #!BN 10/08/15 
         }
        #! attributes(syn[[i]][,j])$contrasts <- NULL                    #!BN-28/04/16   UNCOMMENT???? 
       }                                                                       
     }
   } # end i loop (m)
 } # end synthesising (m > 0)

 return(list(syn = syn, fits = fits))
}


###-----remove.lindep.syn--------------------------------------------------

remove.lindep.syn <- function(x, y, eps = 0.00001, maxcor = 0.99999, 
                              allow.na = FALSE, ...) 
{
  if (ncol(x) == 0) return(NULL) 
  if (eps <= 0) stop("\n Argument 'eps' must be positive.", call. = FALSE)
  xobs <- sapply(x, as.numeric)                                       
  yobs <- as.numeric(y)
  keep <- unlist(apply(xobs, 2, var) > eps)
  keep[is.na(keep)] <- FALSE
  keep <- keep & suppressWarnings((unlist(apply(xobs, 2, cor, yobs)) < maxcor)) # if y includes NA -> NAs error
  if (all(!keep)) warning("\nAll predictors are constant or have too high correlation.\n")
  ksum <- sum(keep)
  cx   <- cor(xobs[, keep, drop = FALSE], use = "all.obs")
  eig  <- eigen(cx, symmetric = TRUE)
  ncx  <- cx
  while (eig$values[ksum]/eig$values[1] < eps) {
    j   <- (1:ksum)[order(abs(eig$vectors[, ksum]), decreasing = TRUE)[1]]
    keep[keep][j] <- FALSE
    ncx  <- cx[keep[keep], keep[keep], drop = FALSE]
    ksum <- ksum - 1
    eig  <- eigen(ncx)
  }
  # if (!all(keep)) cat("\tVariable(s): ", paste(dimnames(x)[[2]][!keep], collapse = ", "),
  #   " removed due to linear dependency",sep="")
  return(keep)
}
