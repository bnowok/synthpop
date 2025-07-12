###-----new version with big denom checks---------------------
disclosure <- function(object, data, ...) UseMethod("disclosure")

###-----disclosure.default-----------------------------
disclosure.default <- function(object, ...)
  stop("No disclosure method associated with class ", class(object), call. = FALSE)

###-----disclosure.data.frame---disclosure.list--------
disclosure.data.frame <- disclosure.list <- 
  function(object, data, cont.na = NULL, keys , target , print.flag = TRUE,
           denom_lim = 5, exclude_ov_denom_lim = FALSE, 
           not.targetlev = NULL,
           usetargetNA = TRUE, usekeysNA = TRUE, 
           exclude.keys =NULL, exclude.keylevs = NULL, exclude.targetlevs = NULL,
           ngroups_target = NULL, ngroups_keys = NULL, 
           thresh_1way = c(50, 90),thresh_2way = c(4, 80),
           digits = 2, to.print =c("short"), compare.synorig = TRUE,...) 
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
 ##------------------------------- check with synorig.compare ---------------------------------   
 
    if (compare.synorig) {
    
        if (m ==1) adjust.data <- synorig.compare(object,data, print.flag = FALSE) else
        if (m > 1) adjust.data <- synorig.compare(object[[1]],data, print.flag = FALSE)
        
        if (!adjust.data$unchanged) {
          object <- adjust.data$syn
          data <- adjust.data$orig
          cat("Synthetic data or original or both adjusted with synorig.compare to try to make them comparable\n\n")
          if (m > 1) {cat("only first element of the list has been adjusted and will be used here\n\n")
          m <- 1 }
        }
    } 

    object <- list(syn = object, m = m, cont.na = cont.na) 
    class(object) <- "synds"
    
    res <- disclosure(object, data, keys , target = target , denom_lim = denom_lim,
                      exclude_ov_denom_lim = exclude_ov_denom_lim, 
                      print.flag = print.flag, digits = digits, 
                      usetargetNA = usetargetNA, usekeysNA = usekeysNA,  
                      not.targetlev = not.targetlev,
                      exclude.keys = exclude.keys, exclude.keylevs =  exclude.keylevs,
                      exclude.targetlevs =  exclude.targetlevs, ngroups_target = ngroups_target,
                      ngroups_keys = ngroups_keys, thresh_1way = thresh_1way,thresh_2way = thresh_2way,
                      to.print = to.print, ...) 
    
    res$call <- match.call()
    return(res)
  }
    
###-----disclosure.synds-------------------------
disclosure.synds <-  function(object, data, keys , target , print.flag = TRUE,
                              denom_lim = 5, exclude_ov_denom_lim = FALSE, 
                              not.targetlev = NULL,
                              usetargetNA = TRUE, usekeysNA = TRUE, 
                              exclude.keys =NULL, exclude.keylevs = NULL, exclude.targetlevs = NULL,
                              ngroups_target = NULL, ngroups_keys = NULL, 
                              thresh_1way = c(50, 90),thresh_2way = c(4, 80),
                              digits = 2, to.print =c("short"),...) 
{ 

###-----------------------check input parameters ----
    if (!(is.data.frame(data)) )   stop("data  must be a data frame \n\n", call. = FALSE)
   data <- data.frame(data) ## in case it is table or a tibble
    if (!( inherits(object,"synds")))   stop(" object must be an object of class synds\n\n", call. = FALSE)
   
   if (is.numeric(keys) & !all(keys %in% 1:dim(data)[2])) stop("If keys are numeric they must be in range 1 to ",
                                      dim(data)[2] , call. = FALSE  )
   if (is.numeric(keys)) keys <- names(data)[keys]
   if (is.numeric(target) & !all(target %in% 1:dim(data)[2])) stop("If target is numeric it must be in range 1 to ",
                                                                   dim(data)[2] , call. = FALSE  )
   if (is.numeric(target)) target <- names(data)[target]
  
  Norig <- dim(data)[1]
  if (object$m ==1 ){ names.syn <- names(object$syn) 
  } else {names.syn <- names(object$syn[[1]])}
  # target must be a single variable in  data and object$syn
  # keys must be a vector of variable names in data and in s
  # target must not be in keys
if (!(all(keys %in% names(data)) && all(keys %in% names.syn) 
      && all(target %in% names(data)) & all(target %in% names.syn )))
  stop("keys and target must be variables in data and synthetic data. \n", call. = FALSE)
   if (any(duplicated(keys)))
     stop("keys cannot include duplicated values. \n", call. = FALSE)
  if (!(length(target)==1))
    stop("target must be a single variable \n", call. = FALSE)
   if (target %in% keys)
     stop("target cannot be in keys \n", call. = FALSE)
   if (any(names(data)  == "target"))
     stop("your data have a variables called 'target' please rename in original and synthetic data.\n\n", call. = FALSE) 

   if (!(length(usetargetNA)==1))
     stop("usetargetNA must be a single logical value \n", call. = FALSE)
   if (length(usekeysNA) ==1 ) usekeysNA <- rep(usekeysNA, length(keys))
   if (!(length(usekeysNA)==length(keys)))
     stop("usekeysNA must be a logical value of same length as keys\n", call. = FALSE)
  # get keys and usekeysNA in same order as in data
   
   #cat(keys,"keys\n", usekeysNA,"usekeysNA\n")
   #cat((1:length(names(data)))[names(data) %in% keys],"(1:length(names(data)))[names(data) %in% keys]\n")
   oldkeys <- keys
   keys <- names(data)[(1:length(names(data)))[names(data) %in% keys]]
   usekeysNA <- usekeysNA[match(oldkeys,keys)]

     # check excluded combinations
  if (!is.null(exclude.keys)) {
    if (! length(exclude.keylevs) == length(exclude.keys) & 
        length(exclude.keylevs) == length(exclude.targetlevs) ) stop("All excludes must be the same length\n" , call. = FALSE)
    if (!all(exclude.keys %in% keys)) stop("exclude.keys must be the name of one of your keys", call.= FALSE)
  }
   if(!is.null(denom_lim)){
     if (!(round(denom_lim) == denom_lim && denom_lim > 0 )){
       cat(denom_lim ,"denom_lim\n")
       stop("\ndenom_lim must be an integer >0\n", call.= FALSE)
     }
   }
   
    if (!is.null(ngroups_target)) {
      if(!(length(ngroups_target) ==1 )) stop("\nngroups_target must be a single value", call.= FALSE)
      if( ngroups_target  == 1 )
        stop("\nTarget ngroups cannot be set to 1", call.= FALSE)
    }
   
   if (!is.null(ngroups_keys)) {
     if (length(ngroups_keys) == 1) ngroups_keys <- rep(ngroups_keys, length(keys))
       if(!(length(ngroups_keys) == length(keys) ))
       stop("\nngroups_keys must be an vector of length 1 or same as keys", call.= FALSE)
     if( any( ngroups_keys  ==1) )
       stop("\nElements of ngroups cannot be set to 1", call.= FALSE)
   }
###---------------------- define output items-------------
   
  allCAPs         <- matrix(NA,object$m,5)
  attrib   <- matrix(NA,object$m,8)
  ident  <- matrix(NA,object$m ,4)

  dimnames(allCAPs) <- list(1:object$m, c("baseCAPd","CAPd", "CAPs", "DCAP","TCAP"))
  dimnames(attrib) <- list(1:object$m,c("Dorig","Dsyn","iS", "DiS","DiSCO",  "DiSDiO", "max_denom","mean_denom"))
  dimnames(ident) <- list(1:object$m,c("UiO", "UiS","UiOiS", "repU"))
  Nexclusions <- list(1:object$m)
  Nexclusions <- check_1way <- check_2way <-list(1:object$m)

###-----------------------restrict data sets to targets and keys------------------------------
m <- object$m
if (m == 1) syndata <- list(object$syn) else syndata <- object$syn
  dd <- data ## rename target variable to target
  targetx =dd[, names(dd) == target]
  dd$target <- targetx
 # cat(names(dd),"names(dd)\n")
  dd <- as.data.frame(dd[,c( "target",keys)])
  for (jj in 1:m) {
      targetx = syndata[[jj]][, names(syndata[[jj]]) == target]
      syndata[[jj]]$target <- targetx
      syndata[[jj]] <- as.data.frame(syndata[[jj]][, c( "target",keys)])
  }
###----- get cont.na parameters for stratified synthesis----------------------------------------
  # --------
  if (!is.null(object$strata.syn)) cna <- object$cont.na[1, ] 
  else   cna <- object$cont.na
###------ get cna levels for target and keys  ------------------------

  cna <- cna[c(target,keys)]

  for ( i in 1:length(cna)) {
    nm <- names(cna)[i]
    vals <- unique(cna[[i]][!is.na(cna[[i]])])  # get variables with cont.na other than missing
    if (length(vals) > 0){
      for (j in 1:length(vals))
        n_cna <- sum(vals[j] == data[,nm] & !is.na(data[,nm]))
      if (n_cna == 0) stop("\nValue ", vals[j], " identified as denoting a special or missing in cont.na for ",nm, " is not in data.\n",sep = "", call. = FALSE)
      else if (n_cna < 10 & print.flag) cat ("\nWarning: Only ",n_cna ," record(s) in data with value ",vals[j]," identified as denoting a missing value in cont.na for ",nm, "\n\n", sep = "")
    }
  }
  ###----------------------- group any continuous variables if ngroups not NULL----------------------
  ### needs a seperate loop to fix all syn data sets with the same breaks-----------------
  # Numeric variables
  if ( is.null(ngroups_target) ) ngroups_target <- 0
  if ( is.null(ngroups_keys) )    ngroups_keys  <- 0
  if (length(ngroups_keys) == 1)  ngroups_keys  <- rep(ngroups_keys,length(keys))
  ngroups <- c(ngroups_target,  ngroups_keys)

  if (print.flag) if (any (ngroups >0 & !sapply(dd,is.numeric))) cat("\n\nWith target",target,"variable(s) you have asked to group are not numeric:\n",
                   names(dd)[ngroups >0 & !sapply(dd,is.numeric)]," no grouping done for them.\n")
  if (any(ngroups > 0)){
   togroup <- (1:dim(dd)[2])[ ngroups >0 & sapply(dd,is.numeric) ]
   for (i in togroup) {
      syn0 <- c(sapply(syndata, '[[', i)) ## all synthetic values for ith var
      for (j in 1:m) {
            grpd <- group_num(dd[,i], syndata[[j]][,i], syn0,
                            ngroups[i], cont.na = cna[[i]], ...)
           if (length(table(grpd[[1]])) < 3 ) {
             grpd <- group_num(dd[,i],  syndata[[j]][,i], syn0,
                            cont.na = cna[[i]], n = ngroups[i], style = "equal")
             if (length(table(grpd[[1]])) < 3 ) cat("Only",length(table(grpd[[1]])),"groups produced for", names(data)[j],"even after changing method.\n")
             else if (print.flag) cat("Grouping changed from 'quantile' to  'equal' in function numtocat.syn for",names(data)[j],"because only",length(table(grpd[[1]]))," groups produced\n")
           }
       syndata[[j]][,i] <- grpd[[2]]
      }
      dd[, i] <- grpd[[1]]
      if (print.flag) {
        if (i  == 1) cat("Numeric values of",names(dd)[i],target,"grouped into ", length(table(dd[,i])),"groups\n with levels", names(table(dd[,i])),"\n")
        else cat("Numeric values of key",names(dd)[i],"grouped into ", length(table(dd[,i])),"groups\nwith levels", names(table(dd[,i])),"\n")
      }
   }
}
###--- make any remaining numeric values into factors---------------------------------------

  if (any(sapply(dd,is.numeric))){
    numvars <- (1:dim(dd)[2])[sapply(dd,is.numeric)]
    for ( i in numvars) {
      dd[,i] <- factor(dd[,i])
      for (jj in 1:object$m){ syndata[[jj]][,i] <- factor(syndata[[jj]][,i])
     }
    }
  }
  check1 <- check2 <- ""

###------------------- loop over object$m syntheses---------------

  for ( jj in 1:object$m) {
###-------------------------------- PRINT WHERE AT ----------------------------------   
    if (print.flag) cat("-------------------Synthesis",jj,"--------------------\n") 
    
    ss <- syndata[[jj]]
###----------------------- sort missings-----------------     
###  replace Missing values with factor value of "Missing" to make tables easier
    tomissing <- function(x){
      if (!is.factor(x)) stop("x must be a factor\n",.call =FALSE)
      x <- as.character(x)
      x[is.na(x)] <- "Missing"
      x <- factor(x)
    }
    if ( any(is.na(dd$target)) ) dd$target <- tomissing(dd$target)
    if ( any(is.na(ss$target)) ) ss$target <- tomissing(ss$target)

    for ( i in 1:length(keys)) {
      if ( any(is.na(dd[,names(dd) == keys[i]])) ) dd[,names(dd) == keys[i]] <-
                                      tomissing(dd[,names(dd) == keys[i]])
      if ( any(is.na(ss[,names(ss) == keys[i]])) ) ss[,names(ss) == keys[i]] <-
          tomissing(ss[,names(ss) == keys[i]])

 }

     Nd <- dim(dd)[1]
     Ns <- dim(ss)[1]
###------------------------ make composite variable for keys --------------------------
if (length(keys) >1) {
   ss$keys <- apply(ss[, names(ss) %in% keys],1,function(x) paste(x, collapse = " | "))
   dd$keys <-   apply(dd[, names(dd)  %in% keys],1,function(x) paste(x, collapse = " | "))
} else {
  ss$keys <- ss[, names(ss) == keys]
  dd$keys <- dd[, names(dd) == keys]
}


  
###--------------------- make tables ---------------------------

  NKd <- length(table(dd$keys))
  NKs <- length(table(ss$keys))

  tab_kts <- table(ss$target,ss$keys)
  tab_ktd <- table(dd$target,dd$keys)   ## two way target and keys table orig
  if (print.flag) cat("Table for target",target, "from GT alone with keys has",
                     dim(tab_ktd)[1], "rows", dim(tab_ktd)[2], "colums.\n")
         
  
  ###  get keys in s and d
  #
  Kd <- names(table(dd$keys))
  Ks <- names(table(ss$keys))
  Kboth <- Kd[Kd %in% Ks]
  Kall <- c(Kd,Ks)
  Kall <- Kall[!duplicated(Kall)]

  ### same thing for target
  
  Td <- names(table(dd$target))
  Ts <- names(table(ss$target))
  Tboth <- Td[Td %in% Ts]
  Tall <- c(Td,Ts)
  Tall <- Tall[!duplicated(Tall)]

  ### augment keys tables to match

   if (!(all(Kd %in% Ks))) { ## some original keys not found in synthetic
    extraKd <- Kd[!(Kd %in% Ks) ]
    extra_tab <- matrix(0,dim(tab_kts)[1],length(extraKd))
    dimnames(extra_tab) <- list(dimnames(tab_kts)[[1]],extraKd)
    tab_kts <- cbind(tab_kts,extra_tab)
    tab_kts  <- tab_kts[, order(dimnames(tab_kts)[[2]]), drop = FALSE] 
   }


  if (!(all(Ks %in% Kd))) {  ## extra synthetic keys not in original
    extraKs <- Ks[!(Ks %in% Kd) ]
    extra_tab <- matrix(0,dim(tab_ktd)[1],length(extraKs))
    dimnames(extra_tab) <- list(dimnames(tab_ktd)[[1]],extraKs)
    tab_ktd <- cbind(tab_ktd,extra_tab)
    tab_ktd <- tab_ktd[,order(dimnames(tab_ktd)[[2]]), drop = FALSE]
  }
  if (!(all(Td %in% Ts))) { ## some original target levels not found in synthetic
    extraTd <- Td[!(Td %in% Ts) ]
    if (is.null(dim(tab_kts)))  extra_tab <- matrix(0,length(extraTd),1)
    else extra_tab <- matrix(0,length(extraTd),dim(tab_kts)[2])
    dimnames(extra_tab) <- list(extraTd , dimnames(tab_kts)[[2]])
    tab_kts <- rbind(tab_kts,extra_tab)
    tab_kts  <- tab_kts[order(dimnames(tab_kts)[[1]]), , drop = FALSE] 
  }   else extraTd <- NULL

  if (!(all(Ts %in% Td))) {  ## extra synthetic target levels not in original  ############### 
    extraTs <- Ts[!(Ts %in% Td) ]
    extra_tab <- matrix(0,length(extraTs),dim(tab_ktd)[2],)
    dimnames(extra_tab) <- list(extraTs,dimnames(tab_ktd)[[2]])
    tab_ktd <- rbind(tab_ktd,extra_tab)
    tab_ktd <- tab_ktd[order(dimnames(tab_ktd)[[1]]),] 
  }   else extraTs <- NULL
 
  if(print.flag) cat("Table for target ",target, "from GT & SD with all key combinations has",
                     dim(tab_ktd)[1], "rows", dim(tab_ktd)[2], "colums.\n")  
###------------------------- calculate proportions and margins ---------
  tab_ktd_p <- sweep(tab_ktd,2,apply(tab_ktd,2,sum),"/")
  tab_ktd_p[is.na(tab_ktd_p)] <- 0
  tab_kts_p <- sweep(tab_kts,2,apply(tab_kts,2,sum),"/")
  tab_kts_p[is.na(tab_kts_p)] <- 0

  tab_kd <- apply(tab_ktd,2,sum)
  tab_td<- apply(tab_ktd,1,sum)
  tab_ks <- apply(tab_kts,2,sum)
  tab_ts <- apply(tab_kts,1,sum)

  NKall <- length(tab_kd)
  NKboth <- length(Kboth)
  NTd <- length(Td)
  NTs <- length(Ts)
  Nboth <-  sum(tab_kd[names(tab_kd) %in% Kboth])
  Nd_ins <- sum(tab_kd[names(tab_kd) %in% Ks])
###------------------------get tables for  calculating attribute disclosure measures-----------------------------------------------------

did <- tab_ktd ; did[tab_ktd_p != 1] <- 0
dis <- tab_kts ; dis[tab_kts_p != 1] <- 0
keys_syn <- apply(tab_kts,2,sum)
tab_iS <- tab_ktd
tab_iS[,keys_syn ==0 ] <- 0
tab_DiS <- tab_ktd
anydis <- apply(tab_kts_p,2,function(x) any(x ==1))
tab_DiS[,!anydis] <- 0
tab_DiSCO <- tab_iS
tab_DiSCO[tab_kts_p != 1] <- 0
tab_DiSDiO <- tab_DiSCO
tab_DiSDiO[tab_ktd_p != 1] <- 0


Nout <- rep(0, 6)
names(Nout) = c("excluded target","missing target","missing in  keys", "set to exclude","over denom_lim","remaining")
Nexcludes <- matrix(0,8,6)
dimnames(Nexcludes) <- list(c("original","synthetic","Dorig","Dsyn","iS", "DiS", "DiSCO", "DiSDiO"),
                            c("excluded target","missing target","missing in  keys", "set to exclude","over denom_lim","remaining"))
###----------------------------- now exclusions----------------------------------

tab_exclude <- function(xx,col,Nexcludes) {
 total <- sum(xx)
 ###------------------- drop all target records with not.targetlev---------------------------------
 if (!is.null(not.targetlev) && !(not.targetlev == ""))
  {
   Nout[1] <- sum(xx[dimnames(xx)[[1]] %in% not.targetlev,])
   xx[dimnames(xx)[[1]] %in% not.targetlev,] <- 0
 }
  ###------------------------ items excluded fOr missing values -------------------------
  #
##### missings excluded values from tables if set

if (!usetargetNA && any(dd$target == "Missing")) {
  Nout[2] <- sum(xx[dimnames(xx)[[1]] == "Missing",])
  xx[dimnames(xx)[[1]] == "Missing",] <- 0
}

 for (i in 1:length(keys)){
   if (!usekeysNA[i]) {   ## do not use NA values for ith key
     key_levs <- dimnames(xx)[[2]]
     drop_d <-  word(key_levs,i, sep = fixed(" | ")) == "Missing"
     Nout[3] <- Nout[3] + sum(xx[,drop_d])
     xx[,drop_d] <- 0
  }
}
###-------------------------- remove any excluded two way  combinations----------------------------------
 if (!is.null(exclude.keys)) {

  if (!all(exclude.targetlevs %in% levels(dd$target))) stop("exclude.targetlevs must be one of levels of ",target,"\n", call. =FALSE)
  
   for (i in 1:length(exclude.keys)){
    vout <- (1:(dim(xx)[1]))[dimnames(xx)[[1]] == exclude.targetlevs[i]]
    klev <- levels(dd[, names(dd) == exclude.keys[i]])
    if (!all(exclude.keylevs[i] %in% klev)) stop("exclude.keylevs position ",i, " must be one of levels of ",keys[i],"\n", call. =FALSE)#
    kind <- (1:length(keys))[keys == exclude.keys[i]]
    wordk <- word(dimnames(xx)[[2]], start = rep(kind, dim(xx)[2]), sep = fixed(" | "))
    kout <- (1:dim(tab_ktd)[2])[wordk == exclude.keylevs[i]]
    Nout[4] <- Nout[4] + sum(xx[vout, kout])
    xx[rep(vout,length(kout)), kout] <- 0
   }
 }  

###------------------ exclude if over denom_lim ----------------------------------
 if (exclude_ov_denom_lim)  {
   Nout[5] <- sum(xx[xx > denom_lim ])
   xx[xx > denom_lim ] <- 0
 }
 Nout[6] <-  total -sum(Nout[1:5])
 
###---------------------------- copy Nout into rows of Nexcludes ------------------------------
 Nexcludes[col,] <- Nout
   return(list(tab = xx, Nexcludes = Nexcludes))
}

yy <- tab_exclude(tab_ktd,1,Nexcludes)
tab_ktd <- yy$tab; Nexcludes <- yy$Nexcludes
yy <- tab_exclude(tab_kts,2,Nexcludes)
tab_kts <- yy$tab; Nexcludes <- yy$Nexcludes
yy <- tab_exclude(did,3,Nexcludes)
did <- yy$tab; Nexcludes <- yy$Nexcludes
yy <- tab_exclude(dis,4,Nexcludes)
dis <- yy$tab; Nexcludes <- yy$Nexcludes
yy <- tab_exclude(tab_iS,5,Nexcludes)
tab_iS<- yy$tab; Nexcludes <- yy$Nexcludes
yy <- tab_exclude(tab_DiS,6,Nexcludes)
tab_DiS<- yy$tab; Nexcludes <- yy$Nexcludes
yy <- tab_exclude(tab_DiSCO,7,Nexcludes)
tab_DiSCO <- yy$tab; Nexcludes <- yy$Nexcludes
yy <- tab_exclude(tab_DiSDiO,8,Nexcludes)
tab_DiSDiO <- yy$tab; Nexcludes <- yy$Nexcludes
###--------------------- exclusions done-----------------------------------------
###-----------------------get  identity disclosure measures -------------

tab_ks <- apply(tab_kts,2,sum)
tab_kd <- apply(tab_ktd,2,sum)
tab_ks1 <- tab_ks
tab_ks1[tab_ks1>1] <- 0
tab_kd1 <- tab_kd 
tab_kd1[tab_kd1>1] <- 0
tab_kd1_s <- tab_kd1[names(tab_kd1) %in% Ks]
tab_ksd1 <-tab_kd[tab_ks == 1 & tab_kd == 1] ## repU

tab_ktd_p <- sweep(tab_ktd,2,apply(tab_ktd,2,sum),"/")
tab_ktd_p[is.na(tab_ktd_p)] <- 0
tab_kts_p <- sweep(tab_kts,2,apply(tab_kts,2,sum),"/")
tab_kts_p[is.na(tab_kts_p)] <- 0

UiS<- sum(tab_ks1)/Ns*100
UiO<- sum(tab_kd1)/Nd*100
UiOiS<- sum(tab_kd1_s)/Nd*100
repU <- sum(tab_ksd1)/Nd*100

ident[jj,] <- c( UiO,UiS, UiOiS,repU )

###----------------------------- attrib dis measures-------------------------

Dorig <- sum(did)/Nd*100
Dsyn <- sum(dis)/Ns*100
iS <- sum(tab_iS)/Nd*100
DiS <- sum(tab_DiS)/Nd*100
DiSCO <- sum(tab_DiSCO)/Nd*100
DiSDiO <- sum(tab_DiSDiO)/Nd*100

attrib[jj,] <- c( Dorig,Dsyn,iS,DiS,DiSCO,DiSDiO,max(tab_DiSCO),mean(tab_DiSCO[tab_DiSCO>0]))
Nexclusions[[jj]] <- Nexcludes

###----------------- get  CAP and DCAP measures---------------------

baseCAPd <- sum((tab_td/Nd)**2)*100 ## just from univariates
CAPd <-   sum(apply(tab_ktd_p^2,2,sum)*tab_kd)/Nd*100 ## from tables

CAPs <-   sum(apply(tab_kts_p^2,2,sum)*tab_ks)/Ns*100
DCAP <-   sum(apply(tab_kts_p*tab_ktd,2,sum))/Nd*100
##restrict to uniques in d
tab_kts_pU <- tab_kts_p
tab_kts_pU[,tab_kd != 1] <- 0
TCAP_denom <- sum(tab_kd[tab_ks >0])
#TCAP0 <-   sum(apply(tab_kts_pU*tab_ktd,2,sum))/Nd*100
TCAP <-   sum(tab_DiSCO)/TCAP_denom*100

#cat(baseCAPd,"baseCAPd",TCAP0,"TCAP0",TCAP,"TCAP\n")
allCAPs[jj,] <- c( baseCAPd,  CAPd,  CAPs, DCAP,TCAP)

#if (print.flag) cat("-------------------Disclosure measures completed\n")

#if (print.flag) cat("------------------Now 1 and 2 way checks\n")
###----------------------- checks for most_dis_lev 1 way ---------------------------
  
tab_target <- apply(tab_ktd,1,sum)
tab_dis_target <- apply(tab_DiSCO,1,sum)
pctdisLev <- max(tab_dis_target)/sum(tab_dis_target)*100
totalDisclosive <- sum(tab_dis_target)
PctDisAll <- totalDisclosive/sum(tab_target)*100
tab_dis_target <- sort(tab_dis_target, decreasing = TRUE)
most_dis_lev <- names(tab_dis_target)[1]
PctLevelAll <- tab_target[names(tab_target) == most_dis_lev]/sum(tab_target)*100
PctLevelDis <- tab_dis_target[names(tab_dis_target) == most_dis_lev]/sum(tab_dis_target)*100
nLevelDis <-tab_dis_target[names(tab_dis_target) == most_dis_lev]

if (nLevelDis > thresh_1way[1] &  PctLevelDis > thresh_1way[2]) {
    check_1way[[jj]] <- data.frame(Level = most_dis_lev, All = sum(tab_target) ,PctLevelAll =  PctLevelAll,
                                totalDisclosive =  totalDisclosive, nLevelDis = nLevelDis,PctLevelDis =  PctLevelDis)
    check1 <- paste("Check " ,target," level ",most_dis_lev)
}
else check_1way[[jj]] <- ""
#if (print.flag) cat("------------------One-way checks completed\n")
###------------------------get details for check_2way-----------------
### only implemented for DiSCO  though could be changed

xx <- tab_DiSCO
if (any(xx > thresh_2way[1]))  {
  
  denoms <- as.vector(xx[xx > thresh_2way[1]])
  rows <- rep(1:dim(xx)[1],dim(xx)[2])
  rows <- rows[xx > thresh_2way[1]]
  cols <- rep(1:dim(xx)[2],rep(dim(xx)[1],dim(xx)[2]))
  cols <- cols[xx > thresh_2way[1]]
  
  target_levs = dimnames(xx)[[1]][rows]
  key_levs = dimnames(xx)[[2]][cols]
  temp <- data.frame(denoms,target_levs,key_levs)
  for (i in 1:length(keys)) {
    temp[,3+i] <- paste(target_levs, word(key_levs,i, sep = fixed(" | ")), sep="|")
    if (i ==1) allpairs <- data.frame(npairs = table(rep(temp[,3+i],denoms)),keys = keys[i])
    else allpairs <- rbind(allpairs, data.frame(npairs = table(rep(temp[,3+i],denoms)),keys = keys[i]))
  }
  names(allpairs)<- c("target_key_levs","npairs","key")
  allpairs <- allpairs[order(-allpairs$npairs),]
  allpairs$target_lev <- word(allpairs$target_key_levs, 1,sep=fixed("|"))
  allpairs$key_lev <- word(allpairs$target_key_levs, 2,sep=fixed("|"))
  allpairs$key_total <- allpairs$key_target_total <- 0

  for (ii in 1:dim(allpairs)[1]) {

    tttab <- table(dd$target[dd[,names(dd)==allpairs$key[ii]] == allpairs$key_lev[ii] ])
    allpairs$key_total[ii] <- sum(tttab)
    allpairs$key_target_total[ii] <- tttab[names(tttab) == allpairs$target_lev[ii]]
     }

  allpairs$PctTargetKeyLevel <- allpairs$key_target_total*100/allpairs$key_total
 #


  details <- allpairs[allpairs$PctTargetKeyLevel > thresh_2way[2],c(1:3,6:8)]
  if (dim(details)[1] >0 ) {
    check_2way[[jj]] <- details 
    check2 <- details[1,]
    check2 <- paste(dim(details)[1],"pairs need checks")
  }
  else  check_2way[[jj]] <- ""
} 
else  check_2way[[jj]] <- ""

#if (print.flag) cat("----------------------------Two-way checks completed.\n")
###-------------------------- end of jj loop -----------------------------------

} ######## end of jj loop
################################################ save ###############################
if (object$m ==1){
  if (!check1 == "") check_1way = check_1way[[1]]
  if (!check2 == "") check_2way = check_2way[[1]]
}
    
res <- list(ident = data.frame(ident), attrib = data.frame(attrib), 
            allCAPs = data.frame(allCAPs), check_1way = check_1way, 
            check1 = check1, check_2way = check_2way,
            check2 =check2, Nexclusions = Nexclusions,
            keys = keys, target = target, digits = digits, Norig = Norig,
            to.print = to.print, call = match.call())
class(res) <- "disclosure"
return(res)
}
###---------------------------print.disclosure---------------------------- 
print.disclosure <- function(x, to.print = NULL, digits = NULL,   ...)
{
  if (is.null(to.print)) to.print <- x$to.print
  if (is.null(digits)) digits <- x$digits

  if (!all(to.print %in% c("short","allCAPs","ident","attrib",
                           "check_2way","check_1way","all","exclusions")))   
  stop("to.print must be  choices from 'short','allCAPs',
        'ident','attrib', 'check_2way','check_1way','all','exclusions'\n", call. = FALSE)
  
  if (length(to.print) >1 & ("short" %in% to.print)) stop("A 'Short' entry in to.print should not be combined with other options.\n", call. = FALSE)
  
  
  cat("Disclosure measures from synthesis for",x$Norig, "records in original data.\n")

  
  nexcluded <- sapply(x$Nexclusions,function(x) sum(x[,-6]))
  if  (any(nexcluded>0)) {
    cat("\nSome records excluded from the evaluation of disclosure risks")
    cat("\nSee details by adding 'exclusions' to the parameter to.print of disclosure,
    or by printing the details with to,print including 'exclusions'.")
  }
  
  if (length(to.print) == 1 && to.print == "short") {
    cat("\nIdentity  measures for keys", x$keys,
        "\nand attribute measures for",x$target,"from the same keys\n" )
    short <- rbind(c(x$ident[1,1],x$attrib[1,1]),
                   cbind(x$ident[,4],x$attrib[,5]))
    dimnames(short) <- list(c("Original",paste("Synthesis", 1:dim(x$attrib)[1])),c("Identity (UiO/repU)","Attrib (Dorig/DiSCO)"))
    print(round(short, digits))

  }
  if ("all" %in% to.print) {
    cat("\nOutput from the following call to function disclosure()")
    print(x$call)
  }
    
    if (any(c("ident","all") %in% to.print))  { 
      cat("\nIdentity disclosure measures for", dim(x$ident)[[1]] ,"synthetic data set(s) from keys:\n", x$keys,"\n\n")
      print(round(x$ident, digits))
      }
    
    if (any(c("attrib","all") %in% to.print)) {
      cat("\nAttribute disclosure measures for",x$target,"from keys:",x$keys,"\n")
       print(round(x$attrib, digits))
    }
if (any(x$Nidentexclude >0)) {
  cat("\nNumber of records excluded from identity disclosure for missing keys\n")
  print(x$Nidentexclude)  
}
###---------------- print CAP measures-----------------------------------------
  if (any(c("allCAPs","all") %in% to.print)) {
    cat("\nCAP measures for the target", x$target, "with keys\n")
    cat(x$keys,"\n\n")
    print(round(x$allCAPs, digits))
  }   
  
  ###------------------------ print check messages  
  if (x$check1 != "" & !("check_1way" %in% to.print)){
    cat("\nThe 1 way distributions of",x$target,"has a large contribution to disclosure from level",x$check1,"\n")
    cat("Please add 'check_1way' to to.print (e.g. print(disclosure_result, to.print = 'check_1way'),\n")
    cat("and look at original data to decide if backround knowledge would make this disclosure likely\n")
    cat("Consider excluding this level with the not.targetlev parameter to the disclosure function.\n\n")
  }
       
  if ("check_1way" %in% to.print){
      if (!is.null(x$check_1way) ) {
      cat("\nDetails of target level contributing disproportionately to disclosure\n")
      print(x$check_1way)
      }
      else {cat("No 1 way checks for your target are flagged at your settings for thresh_1way\n\n ")}
  }
  
  if (!(length(x$check2) == 1 && x$check2 == "") & !("check_2way" %in% to.print)){
    cat("\nDetails of target-key pairs contributing disproportionately to disclosure of",x$target,"\n")
    if (length(x$check2) > 3) cat("Note only the first 3 printed here\n")
    for (i in 1:(min(length(x$check2),3))) cat(x$check2[i],"\n")
    cat("\nPlease examine component $check_2way of the disclosure object and look at original data.\n")
    cat("Consider excluding these key-target pairs with some the following parameters to disclosure:\n")
    cat("exclude_ov_denom_lim = TRUE or defining key-target combinations from exclude.targetlevs,\n")
    cat("exclude.keys and exclude.keylevs\n\n" )
  }
  
  if ("check_2way" %in% to.print){
      cat("\nDetails of target-key combinations contributing disproportionately to disclosure\n")
      if (!is.null(x$check_2way)) {
        cat("This is a list with one for each of",dim(x$check_2way),"syntheses\n")
        print(x$check_2way)
       }
      else {cat("No 2 way checks of key-target combinations are flagged at your settings for thresh_2way\n\n ")}
    }
  ###----------------------  print exclusions -----------------
  if ("exclusions" %in% to.print){
    for (i in 1:length(x$Nexclusions)){
      if (sum(x$Nexclusions[[i]]) > 0 ) {
        cat("\nRecords excluded from attribute exclusion measures for synthesis",i,"\n") 
        print(x$Nexclusions[[i]])}
      else cat("\nNo records excluded from attribute exclusion measures for synthesis",i,"\n")
    }
  }

  
  invisible(x)
}
