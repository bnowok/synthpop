padMis.syn <- function(data, method, predictor.matrix, visit.sequence,
                       nvar, rules, rvalues, default.method, cont.na, 
                       smoothing, event, denom) {

 # Function called by syn to make dummy/factor variable for missing values
 # in continuous variables. Data is augmented by columns for dummy/factor 
 # variables when they are used in sythesis. 

 # Check presence of missing values not covered by missing rules
 # (missing values for non-numeric variables are not counted)   
 No.NA <- vector("list",nvar)
 yes.rules <- sapply(rules, function(x) any(x!=""))
 com.rules <- lapply(rules, paste, collapse=" | ")
 for (j in 1:nvar){
   if (yes.rules[j]){
     No.NA[j] <- with(data,sum(data[!eval(parse(text=com.rules[[j]])),j] %in% cont.na[[j]]))      
   } else {
     No.NA[j] <- sum(data[,j] %in% cont.na[[j]]) 
   }
 }
 No.NA    <- sapply(No.NA,function(x) x>0)
 inpred   <- apply(predictor.matrix!=0,1,any)|apply(predictor.matrix!=0,2,any)
 factorNA <- rep(FALSE,nvar)

 for(j in 1:nvar){
    if (No.NA[j] & is.numeric(data[,j]) & inpred[j]==TRUE){
 
    # augment the data with a column for the original continuous variable with 
    # missing values replaced by zeros and a column for a new factor for 
    # missing values 
      nonmiscode <- 10^(nchar(round(max(data[,j],na.rm=TRUE)))+1)-1               #BN13/11
      y.0  <- ifelse(data[,j] %in% c(cont.na[[j]],rvalues[[j]]),0,data[,j])
      y.NA <- ifelse(data[,j] %in% c(cont.na[[j]],rvalues[[j]]),data[,j],nonmiscode) #BN13/11 0 changed with nonmiscode
      y.NA <- addNA(y.NA,ifany=TRUE) 
      levels(y.NA)[is.na(levels(y.NA))] <- "NAtemp"                         #BN25/08 to allow random forest
      data <- cbind(data,y.0,y.NA)           
      name.0  <- paste(attr(data,"names")[j],0,sep=".")
      name.NA <- paste(attr(data,"names")[j],NA,sep=".")
      names(data)[(ncol(data)-1):ncol(data)] <- c(name.0,name.NA)
      factorNA[(ncol(data)-1):ncol(data)] <- c(FALSE,TRUE) 

    # predictor.matrix is given two extra rows and columns for the new variables
    # rows and columns are copied from an original variable j in predictor.matrix
      predictor.matrix <- rbind(predictor.matrix,  matrix(rep(predictor.matrix[j,], 
                               times=2),byrow=TRUE,nrow=2))
      predictor.matrix <- cbind(predictor.matrix, matrix(rep(predictor.matrix[,j], 
                               times=2),ncol=2))
    # the original variable is removed from predictors (=insert zeros) 
      predictor.matrix[,j] <- 0
    # the original variable is imputed passively so its predictors can be removed as well
      predictor.matrix[j,] <- 0

    # add methods for new variables
      method[ncol(data)-1] <- method[j]
      if (method[j] %in% c("ctree","ctree.proper","cart","cart.proper","rf","bag")) {
        method[ncol(data)] <- method[j] 
      #} else if (method[j] %in% c("rf","bag")) {   
      #  method[ncol(data)] <- "cart" 
      } else if (method[j]=="sample") {   
        method[ncol(data)] <- "sample"                                          
      } else {
        method[ncol(data)] <- ifelse(nlevels(data[,ncol(data)])==2,
                                     default.method[2],default.method[3])
      }   
     
    # pass smoothing to new variable and remove from original one
      smoothing[ncol(data)-1] <- smoothing[j]
      smoothing[ncol(data)]   <- ""
      smoothing[j]            <- ""                                             #BN20V
    
    # pass denom and event for new variable and remove from original one
      denom[ncol(data)-1] <- denom[j]    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! check if correct
      denom[ncol(data)]   <- 0 
      denom[j]            <- 0
      event[ncol(data)-1] <- event[j]    
      event[ncol(data)]   <- 0 
      event[j]            <- 0
    
    # insert the column numbers for the new variables into the visit sequence 
    # before the jth column
      if (any(visit.sequence==j)){                    
        newcols <- c(ncol(data),ncol(data)-1)
        idx <- (1:length(visit.sequence))[visit.sequence==j]-1
        visit.sequence <- append(visit.sequence,newcols,idx)
      # modify method for the original variable
        method[j] <- paste0("~(ifelse(",name.NA,"!=", nonmiscode," | is.na(",name.0,
            "),as.numeric(levels(",name.NA,"))[",name.NA,"],",name.0,"))")    
      }

    # update missing rules and values for the new variables
      if (any(rules[[j]]!="")) rules[[ncol(data)-1]] <-
        c(rules[[j]][rules[[j]]!=""],paste(name.NA,"!=",nonmiscode,sep=""))  #BN13/11
      else rules[[ncol(data)-1]] <- paste(name.NA,"!=",nonmiscode,sep="")    #BN13/11 
      rules[[ncol(data)]]        <- rules[[j]]                      
      rules[[j]]                 <- ""
      #!BN1513 rule "year_death.NA!=0" should have only one correspnding 
      #! rvalue equal to 0; before a vector c(NA,0) was assigned instead of 0 
      if (length(rules[[j]])==1) rvalues[[ncol(data)-1]] <- 0           
      else rvalues[[ncol(data)-1]] <- c(rvalues[[j]],0)                 
      rvalues[[ncol(data)]]        <- rvalues[[j]]
    }
  }
   
  varnames <- dimnames(data)[[2]]  # now includes new names
  dimnames(predictor.matrix) <- list(varnames,varnames)
  names(method) <- varnames
  names(visit.sequence) <- varnames[visit.sequence]
  return(list(data = as.data.frame(data), 
              nvar = ncol(data),
              predictor.matrix = predictor.matrix, 
              method = method, 
              visit.sequence = visit.sequence, 
              rules = rules,
              rvalues = rvalues,
              factorNA = factorNA,
              smoothing = smoothing,
              event = event,
              denom = denom))
}
