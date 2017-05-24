#-----------------------------sampler.syn-------------------------------

sampler.syn <- function(p, data, m, syn, visit.sequence,
                        rules, rvalues, event, proper,
                        print.flag, k, pred.not.syn, models, ...){
# The sampler controls the generation of conditional distributions
# This function is called by syn

  #print(p[-c(1:3,6,7)])                                                         #--TEMP
  
  #--- Assign optional parameters (...) to appropriate synthesising function   
  #rpart.args   <- c(names(formals(rpart)), names(formals(rpart.control)))
  #ctree.args   <- c(names(formals(ctree)), names(formals(ctree_control)))
  #polyreg.args <- c(names(formals(multinom)),names(formals(nnet.default)))
  #polr.args    <- c(names(formals(polr)),names(formals(nnet.default)))
  #
  #method.args0 <- list(cart = rpart.args, cartboot = rpart.args, 
  #                     ctree = ctree.args, survctree = ctree.args,
  #                     polyreg = polyreg.args, polr = polr.args)
  #dots         <- as.list(substitute(list(...)))[-1L]  
  #dots.names   <- names(dots)
  #
  #method.args  <- sapply(method.args0,intersect,dots.names)
  #method.args  <- method.args[sapply(method.args,length)>0]
  #---


  #--- Assign optional parameters (...) to appropriate synthesising function   
  dots  <- as.list(substitute(list(...)))[-1L]         
  meth.with.opt <- paste(c("cart","cartbboot","ctree","survctree","polyreg","polr","rf","bag"), collapse="\\.|")
  meth.check <- grep(meth.with.opt,names(dots),value=TRUE)
  args.err <- !(names(dots) %in% meth.check)
  if (any(args.err)) stop("Unknown optional parameter(s): ", 
    paste(names(dots)[args.err],collapse=", "),
    "\nNote that they have to be method specific, e.g. 'ctree.minbucket' and NOT 'minbucket'\n", 
    call. = FALSE)
  if (length(dots)==0){
    mth.args <- NULL
  } else {  
    mth.args.dots <- strsplit(names(dots), "\\.")
    mth.dots  <- unique(lapply(mth.args.dots, "[[", 1))
    args.dots <- lapply(mth.args.dots, "[[", -1)
    mth.args  <- setNames(vector("list", length(mth.dots)),unlist(mth.dots))
  
    for (i in 1:length(mth.dots)) { 
      ind <- grep(mth.dots[[i]], names(dots))
      mth.args[[i]] <- setNames(dots[ind], args.dots[ind])
    } 
  } 
  #---
  
  fits <- NULL 
    
  if (m > 0){
	  if (print.flag) cat("syn  variables")
    
    if (models) fits <- rep(list(setNames(vector("list", length(p$method)),
      names(p$method))), m) 

    for (i in 1:m){  # begin i loop : repeated synthesising loop
      if (print.flag) cat("\n",i,"   ",sep="")
## augment the data with the actual dummy variables  maybe not needed now??  GR
 # This next code replaces the dummy variables for a synthesised variable with 
 # the new values 
 #     for (j in setdiff(p$visit.sequence,visit.sequence)){
 #       cat.columns <- p$syn[, p$categories[j, 4]]
 #       p$syn[,(j:(j+p$categories[p$categories[j,4],2]-1))] <- 
 #                        matrix((model.matrix(~cat.columns-1)[,-1]),
 #                        ncol=p$categories[p$categories[j,4],2],nrow=nrow(p$data)) 
 #       p$syn[,j]<-as.numeric(p$syn[,j])
 #     }

      for(j in p$visit.sequence) {
        #print(p$denom[j])
        theMethod <- p$method[j]
        # get optional parameters for theMethod if they are provided
        if (!(theMethod %in% names(mth.args))) fun.args <- NULL else            #BN--13/05/2015
        fun.args  <- mth.args[[theMethod]]                                      #BN--13/05/2015
        vname     <- dimnames(p$data)[[2]][j]
        
        if(print.flag & theMethod!="dummy"  & j<=ncol(data)) cat(" ",vname,sep="")
        if(j%%10==0 & j <= ncol(data)) cat("\n    ")
        
########*********
        ya  <-  1:nrow(p$data) 
        ypa <- 1:k    
                     
        # ya=yavailable, ym=ymissing                                            
        if(any(p$rules[[j]]!="")) {
          com.rules  <- paste(p$rules[[j]],collapse=" | ")
          evalrul.y  <- with(p$data,eval(parse(text=com.rules)))
          ym         <- which(evalrul.y==TRUE & !is.na(evalrul.y))
          ya         <- setdiff(1:nrow(p$data),ym)                                  
          evalrul.yp <- with(p$syn,eval(parse(text=com.rules)))         
          ypm        <- which(evalrul.yp==TRUE & !is.na(evalrul.yp))        
          ypa        <- setdiff(1:nrow(p$syn),ypm)       
        }                                                                       
           
        # !="", !="dummy", !="passive"
        if(theMethod!="" & (!is.passive(theMethod)) & theMethod!="dummy" ){
          if (theMethod %in% c("sample","sample.proper","constant")) {
            y   <- p$data[ya, j]
            if (is.factor(y)) y <- y[,drop=TRUE]
            xp  <- length(ypa)
            x   <- length(ya)
            nam <- vname
            f   <- paste("syn", theMethod, sep = ".")
            synfun <- do.call(f, args = list(y=y,xp=xp,
              smoothing=p$smoothing[j],cont.na=p$cont.na[[j]],proper=proper, ...)) 
            p$syn[ypa, j]  <- synfun$res
            if (models) fits[[i]][[j]] <- synfun$fit
            
          } else {
          
            x    <- p$data[ya, p$predictor.matrix[j, ] == 1, drop = FALSE]
            xp   <- p$syn [ypa, p$predictor.matrix[j, ] == 1, drop = FALSE]
            y    <- p$data[ya, j]
            if (is.factor(y)) y <- y[,drop=TRUE]
            nam  <- vname
            f    <- paste("syn", theMethod, sep = ".")
            if(theMethod!="collinear"){                  #! BN to check
              keep <- remove.lindep.syn(x, y, ...)
              x    <- x[, keep, drop = FALSE]
              xp   <- xp[, keep, drop = FALSE]
            }                                            #! BN to check
            if (theMethod=="survctree") {
              if (p$event[j]==-1) yevent <- rep(1,length(y))                    #BN--27/05/2015 'p$' added
              else yevent  <- p$data[ya,p$event[j]]
              survres      <- do.call(f, args = c(list(y=y,yevent=yevent,
                x=x,xp=xp,proper=proper), fun.args))
              p$syn[ypa, j] <- survres[[1]]# synthetic data survival goes to p$syn
              if (p$event[j]!=-1) p$syn[ypa,p$event[j]] <- survres[[2]] # synthetic data event goes to p$syn
              if (models) fits[[i]][[j]] <- survres$fit 
            } else if (theMethod=="logreg" & p$denom[j]!=0) {     #GR denom new   #BN--27/05/2015 '& j<=ncol(data)' removed           
              synfun <- do.call(f, args = list(y=y,x=x,xp=xp,
                denom=p$data[ya,p$denom[j]],denomp=p$syn[ypa,p$denom[j]],       #BN--27/05/2015 'ya','ypa' added
                proper=proper, ...))
              p$syn[ypa, j] <- synfun$res
              if (models) fits[[i]][[j]] <- synfun$fit           
            } else {
              synfun <- do.call(f, args = c(list(y=y, x=x, xp=xp,
                smoothing=p$smoothing[j], proper=proper), fun.args))
              p$syn[ypa, j] <- synfun$res
              if (models) fits[[i]][[j]] <- synfun$fit
            }
          
          }

          if(any(p$rules[[j]]!="")){
            if(length(p$rules[[j]])==1 & length(ypm)>0){
              p$syn[ypm,j] <- p$rvalues[[j]] 
            } else {
              for (r in 1:length(p$rules[[j]])){
                revalrul.yp  <- with(p$syn,eval(parse(text=p$rules[[j]][r])))  
                rypm <- which(revalrul.yp==TRUE & !is.na(revalrul.yp))
                if (length(rypm)>0) p$syn[rypm,j] <- p$rvalues[[j]][r]
              }
            }                 
          }  
        } # end of !="", !="dummy", !="passive"

        else if (is.passive(theMethod)) {
          class0 <- class(p$syn[,j])
          synfun <- syn.passive(data = p$syn, func = theMethod) 
          #p$syn[,j] <- suppressWarnings(model.frame(as.formula(theMethod), p$syn, na.action=na.pass))	#BN 25/08 added suppressWarnings to avoid NAs by coersion for NAtemp
          p$syn[, j] <- synfun$res
          class(p$syn[,j]) <- class0
          if (models) fits[[i]][[j]] <- synfun$fit 
        }

        else if (theMethod=="dummy") {    # replace dummy variables in p$syn
          # getting dummy values from a synthesised categorical variable
          cat.columns <- p$syn[,p$categories[j,4]]  # this is the single column with the data for which this is the dummy
          model.frame(~cat.columns-1,data=p$syn) 
          p$syn[,(j:(j+p$categories[p$categories[j,4],2]-1))] <-   # replaces all the dummies for this variable with
          matrix((model.matrix(~cat.columns-1)[,-1]),              # dummies calculated from the synthesised data
                  ncol=p$categories[p$categories[j,4],2],
                  nrow=nrow(p$syn))
          p$syn[,j] <- as.numeric(p$syn[,j])
          remove("cat.columns")
          if (models) fits[[i]][[j]] <- "dummy"    
        }

      } # end j loop 
    
      #if (k==dim(data)[1]) syn[[i]] <- p$syn[,1:dim(data)[2]]
      #else syn[[i]] <- p$syn[sample(1:dim(data)[1],k),1:dim(data)[2]]
      
########*********
      syn[[i]] <- p$syn[,1:dim(data)[2],drop=FALSE]
      nms <- names(data)
      # exclude unsynthesised if drop.pred.only set to true
      if (sum(pred.not.syn )>0) {
        syn[[i]] <- syn[[i]][,!pred.not.syn]
        nms <- nms[!pred.not.syn]         # GR save names to use below if data just one column
      }
                                          # GR changes extra lines needed # to prevent a single character column being changed to a factor
      chgetochar <- (sum(!pred.not.syn)==1 & class(syn[[i]][,1])=="character")       
  
      syn[[i]] <- as.data.frame(syn[[i]])
      if (chgetochar) {
         syn[[i]][,1] <- as.character(syn[[i]][,1])
         names(syn[[i]]) <- nms
      }
   

      #turn NA level in factors to missing NA's    # to delete - replaced by the code below
  #    for (j in (1:ncol(syn[[i]]))){
  #      if(is.factor(syn[[i]][,j])) {
  #        syn[[i]][is.na(as.character(syn[[i]][,j])),j] <- NA
  #        syn[[i]][,j] <- factor(syn[[i]][,j])
  #      }
  #    }

      #turn NA level in factors / logical to missing NA's
      # and remove contrasts 
      for (j in (1:ncol(syn[[i]]))){
        if(is.factor(syn[[i]][,j])) {                                    #!BN-20/04/16
          if("NAlogical" %in% levels(syn[[i]][,j])){
            levels(syn[[i]][,j])[levels(syn[[i]][,j])=="NAlogical"] <- NA
            syn[[i]][,j] <- as.logical(syn[[i]][,j])
          } else {                                                            
        # syn[[i]][,j] <- factor(syn[[i]][,j],exclude=NA,levels=levels(syn[[i]][,j]))
            levels(syn[[i]][,j])[levels(syn[[i]][,j])=="NAtemp"] <- NA     #!BN 10/08/15 
          }
        #! attributes(syn[[i]][,j])$contrasts <- NULL                           UNCOMMENT???? BN-28/04/16
        }                                                                       
      }

    } # end i loop (m)
 } # end synthesising
  
 if (print.flag) cat("\n")
 return(list(syn = syn, fits = fits))
}



###-----remove.lindep.syn--------------------------------------------------

remove.lindep.syn <- function(x, y, eps=0.00001, maxcor=0.99999, 
                              allow.na=FALSE, ...) {
  if (ncol(x)==0) return(NULL) 
  if (eps <= 0) stop("\n Argument 'eps' must be positive.")
  xobs <- sapply(x,as.numeric)                                       #!BN-1605 was xobs <- x
  yobs <- as.numeric(y)
  keep <- unlist(apply(xobs, 2, var) > eps)
  keep[is.na(keep)] <- FALSE
  keep <- keep & suppressWarnings((unlist(apply(xobs, 2, cor, yobs)) < maxcor)) # if y includes NA -> NAs error
  if (all(!keep)) warning("All predictors are constant or have too high correlation.")
  ksum <- sum(keep)
  cx   <- cor(xobs[, keep, drop=FALSE], use = "all.obs")
  eig  <- eigen(cx, symmetric = TRUE)
  ncx  <- cx
  while(eig$values[ksum]/eig$values[1] < eps) {
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
