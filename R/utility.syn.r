###-----utility.gen--------------------------------------------------------

utility.gen<-function (object, data, method="logit",tree.method="rpart",maxorder = 1,
 vars =NULL,aggregate=FALSE,maxit=200,groups=FALSE,ngroups=5,mincriterion=0,nperms=50,
   cp=1e-8,minbucket=5,digits=2,print.zscores = FALSE,
 usethresh =TRUE,zthresh=1.6, print.variable.importance=FALSE, ...)
{
 cna<-object$cont.na
 ################################### check input parameters ###################################
  if (any(is.na(match(method, c("cart", "logit"))))) {
    stop("Invalid method type must be logit or cart.", call. = FALSE) ############# just 2 methods allowed now
  }
  if (aggregate==TRUE & !method=="logit") {
    stop("Aggregation only works for logit models.", call. = FALSE) 
  }

  if (is.null(data)) stop("Requires parameter 'data' to give name of the Sreal data.\n",  call. = FALSE)
  if (!class(object) == "synds") stop("Object must have class 'synds'.\n", call. = FALSE)
  m <- object$m

  if (!(is.null(vars))) {               ############# checks selected select variables and remove from object$syn
 	if (!(all(vars %in% names(data))))	stop("some vars specified not in names of data\n")
    names(object$cont.na)<-names(object$syn)
    data<-data[,vars]
 	  if (m==1) {
 	   names(object$cont.na)<-names(object$syn) ####  check this always right
     object$syn<-object$syn[,vars]
    }
    else for (i in 1:m)  {
     names(object$cont.na)<-names(object$syn[[1]])
     object$syn[[i]]<-object$syn[[i]][,vars]
    }
   }
   else if (m==1) vars=names(object$syn)
   else vars=names(object$syn[[1]])
   if  (!all(vars %in% names(data))) stop("Some variables in synthetic data not in data.\n", call. = FALSE)
   else data<-data[,vars] ##  make data match synthetic
   unsyn.vars <- names(object$method)[object$method == ""]################## removes unsynthesised
   if (any(vars %in% unsyn.vars) & !is.null(unsyn.vars)) {
	  remove<- (vars %in% unsyn.vars)
	  cat("The following unsynthesised variables removed from utility calculations\n",vars[remove],"\n")
	  vars<-vars[!remove]
	  if (length(vars)==0) stop ("Nothing left \n")
	  data<-data[,vars]
	  if (m==1) object$syn<-object$syn[,vars]
	  else for (i in 1:m) object$syn[[i]]<-object$syn[[i]][,vars]
   } 
   ################ remove single value columns
   #####  may not be needed now as their method may have been set to "" in the synthesis
  leneq1<-function(x) length(table(as.numeric(x))) %in% (0:1) 
  dout<-sapply(data,leneq1)
  if (m==1) sout<-sapply(object$syn,leneq1)
  else  sout<-sapply(object$syn[[1]],leneq1)
  if (any(dout)==T | any(sout==T)) {
   cat("Some columns with single values or all missing values excluded from utility comparisons\n")
   cat("In original  data: ",names(data)[dout],"\n")
   if (m==1) cat("In synthetic data: ",names(object$syn[[1]])[sout],"\n")
   else cat("In synthetic data: ",names(object$syn[[1]])[sout],"\n")
   if (length(names(data)[dout])!=length(names(object$syn)[sout]) || any(names(dout[dout==TRUE])!=names(sout[sout==TRUE]))) 
             	stop("Variables don't match - check.\n", call. = FALSE)
   data<-data[,!dout];
   if (m==1) object$syn<-object$syn[,!sout]
   else for (i in 1:m) object$syn[[i]]<-object$syn[[i]][,!sout]
  }
  #####################categoricasl var##############################################
  #################################### make missing data part of factor
  #
  catvars<-(1:dim(data)[2])[sapply(data,is.factor)]
  for (i in catvars){
   data[,i]<-factor(data[,i])
   if (m==1) object$syn[,i]<-factor(object$syn[,i])
   else for (j in 1:m) object$syn[[j]][,i]<-factor(object$syn[[j]][,i])  
   if (any(is.na(data[,i]))) {
    #cat(names(data)[i],"missings \n")
    data[,i]<-addNA(data[,i])
    if (m==1) object$syn[,i]<-addNA(object$syn[,i])
	  else for (j in 1:m)   object$syn[[j]][,i]<-addNA(object$syn[[j]][,i])
   }
  }
  ######################################################################
  ###########  now numeric variables####################################################
  #
  numvars<-(1:dim(data)[2])[sapply(data,is.numeric)]
  names(numvars)<-names(data)[numvars]
  #
  # if groups==TRUE divide numeric variables into groups
  #
  data0<-data ## to save if m>1
  if  (groups==TRUE){
     for (i in numvars){
     if (m==1) {
   	 groups<-group_num(data[,i],object$syn[,i],ngroups,cont.na=cna[[i]])
   	 data[,i]<-groups[[1]]
   	 object$syn[,i]<-groups[[2]]
	  }
 	  else for (j in 1:m) {
   	 groups<-group_num(data0[,i],object$syn[[j]][,i],ngroups)
   	 data[,i]<-groups[[1]]
   	 object$syn[[j]][,i]<-groups[[2]]
	  }
   }
  }
   for (i in numvars){
    if (anyNA(data[,i])){
     #cat(names(data)[i],"numeric w missings \n")
     newname<-paste(names(data)[i],"NA",sep="_")

     data<-data.frame(data,1*(is.na(data[,i])))
     names(data)[length(data)]<-newname
     data[is.na(data[,i]),i]<-0
	   if (m==1){
	    object$syn<-data.frame(object$syn,1*(is.na(object$syn[,i])))
      names(object$syn)[length(object$syn)]<-newname
      object$syn[is.na(object$syn[,i]),i]<-0
	   }
	   else{
      for (j in 1:m) {
	     object$syn[[j]]<-data.frame(object$syn[[j]],1*(is.na(object$syn[[j]][,i])))
       names(object$syn[[j]])[length(object$syn[[j]])]<-newname
	     object$syn[[j]][is.na(object$syn[[j]][,i]),i]<-0

		  }
     }
    }
   }
  ###########################################################################################################
   propcalcs<-function(syndata){
    df.prop<-rbind(syndata,data) ############# make data frame for calculating propensity score ################
    n1<-dim(data)[1]
    n2<-dim(syndata)[1]
    df.prop<-data.frame(df.prop,t=c(rep(1,n2),rep(0,n1)))
    if (aggregate==TRUE){
	   aggdat=aggregate(df.prop[,1],by=df.prop,FUN=length)
	   wt=aggdat$x
	   aggdat=aggdat[,-dim(aggdat)[2]]
	  }
    if (method=="logit"){
      if (!(maxorder %in% 0:4)) 
      stop("Invalid maximum order of interactions must be in 0:4.\n", call. = FALSE)
     if (maxorder >= 1) logit.int <- as.formula(paste("t ~ .^", maxorder + 1))
	   else logit.int <- as.formula(paste("t ~ ."))
     if (aggregate==TRUE)  propLogit<-glm(logit.int,data = aggdat,  family = "binomial",control=list(maxit=maxit),weights=wt)
	   else               propLogit<-glm(logit.int,data = df.prop, family = "binomial",control=list(maxit=maxit))
     if(propLogit$converged==FALSE) cat("Warning: logistic model did not converge in ",maxit," iterations, \n   you should increase parameter maxit? \n")
     km1=length(propLogit$coef)-1
     logitPhat <- predict(propLogit, type = "response")
     if (aggregate==TRUE){ 
      utilVal <- (sum(wt*(logitPhat - n2/(n1+n2))^2, na.rm = T))*(n1+n2)^3/n1^2/n2
     }
     else{
      utilVal = (sum((logitPhat - n2/(n1+n2))^2, na.rm = T))*(n1+n2)^3/n1^2/n2
     } 
     utilR<-utilVal/km1
     utilStd=(utilVal-km1)/sqrt(km1*2)
     fit=propLogit
     pscore=logitPhat
          res.ind <- list(utilVal=utilVal,utilExp=km1,utilR=utilR,utilStd=utilStd,
     fit=fit)
    }
    else if (method=="cart"){
      if (tree.method=="rpart") {
        #cat("cp for real",cp,"\n")
        propCART = rpart(t ~ ., data = df.prop, method = 'class', control = rpart.control(cp = cp,minbucket=minbucket))
        cartPhat = predict(propCART)[,2]
      } 
      else if (tree.method=="ctree") {
        propCART = ctree(t ~ ., data = df.prop, 
          controls = ctree_control(mincriterion=mincriterion,minbucket=minbucket))
            cartPhat = predict(propCART)
      }
      utilVal = (sum((cartPhat - n2/(n1+n2))^2, na.rm = T))*(n1+n2)^3/n1^2/n2
      ################################
      # now permutation test
      # 
      simutil<-rep(0,nperms)
      if (m==1) j<-1
      if (j==1 ) cat("now running ",nperms," permutation tests printing every 10th\n")
      cat("synthesis ",j,": ")
      for (i in 1:nperms){
        if (floor(i/100)==i/100) cat("\n")
        if (floor(i/10)==i/10)  cat(i," ",sep="")

        pdata<-df.prop
        pdata$t<-sample(pdata$t)
        if (tree.method=="rpart") {
         propCART = rpart(t ~ ., data = pdata, method = 'class', control = rpart.control(cp = cp,minbucket=minbucket))
         cartPhat = predict(propCART)[,2]
        }
        else if (tree.method=="ctree") {
          propCART = ctree(t ~ ., data = pdata, 
                         controls = ctree_control(mincriterion=mincriterion,minbucket=minbucket))
        }
        simutil[i] = (sum((cartPhat - n2/(n1+n2))^2, na.rm = T))*(n1+n2)^3/n1^2/n2/2 ### reduced by factor of 2 for lack of conditioning
      }
      cat("\n")
      utilExp<-mean(simutil)
      utilR<-utilVal/utilExp
      utilStd<-(utilVal-utilExp)/sqrt(var(simutil))
      res.ind <- list(utilVal=utilVal,utilExp=utilExp,utilR=utilR,utilStd=utilStd,
                    fit=propCART)
    }
    return(res.ind)
    ######################################## end of propclac function######################
   }
   if (m==1){
    res.ind<-propcalcs(object$syn)
    res<-list(call=match.call(),m=m,method=method,tree.method=tree.method,nperms=nperms,
          utilVal=res.ind$utilVal, utilExp=res.ind$utilExp,utilR=res.ind$utilR,
               utilStd=res.ind$utilStd,fit=res.ind$fit,digits=digits,print.zscores = print.zscores,
          usethresh =usethresh,zthresh=zthresh, print.variable.importance=print.variable.importance)
   }
   else {
     utilVal<-utilExp<-utilR<-utilStd<-rep(NA,m)
     fit<-as.list(1:m)
     for (j in 1:m){
       res.ind<-propcalcs(object$syn[[j]])
       utilVal[j]<-res.ind$utilVal; utilExp[j]<-res.ind$utilExp; utilR[j]<-res.ind$utilR 
       utilStd[j]<-res.ind$utilStd; fit[[j]]<-res.ind$fit
     }
     res<-list(call=match.call(),m=m,method=method,tree.method=tree.method,nperms=nperms,
            utilVal=utilVal,utilExp=utilExp,utilR=utilR,
               utilStd=utilStd,fit=fit,digits=digits,print.zscores = print.zscores,
            usethresh =usethresh,zthresh=zthresh, print.variable.importance=print.variable.importance)
   }
 class(res) <- "utility.gen"
 return(res)
}


###-----utility.tab--------------------------------------------------------

utility.tab <- function (object, data, vars = NULL, ngroups = 5, 
                         with.missings = TRUE, digits = 2, 
                         print.tables = FALSE, print.zdiff = FALSE, ...) 
{
  # CHECKS #
  #---
  if (is.null(data)) stop("Requires parameter 'data' to give name of the real data.\n", call. = FALSE)
  if (!class(data) == "data.frame") stop("Data must have class 'data.frame'.\n", call. = FALSE)
  if (!class(object) == "synds") stop("Object must have class 'synds'.\n", call. = FALSE)
  
  m <- object$m
  if (m==1) syndata <- object$syn
  else syndata <- object$syn[[1]]
  
  if (dim(data)[1] != dim(syndata)[1]) stop("\nThis function is not for case when sizes of original and synthetic data differ.\n", call.=FALSE)
  if (is.null(vars)) stop("Need to set variables with vars parameter.\n", call. = FALSE)
  else if (!(all(vars %in% names(data)))) stop("Unrecognized variable(s) in vars parameter: ",
    paste(vars[!vars %in% names(data)], collapse = ", "),"\n", call. = FALSE)
  #---
  
  data<-data[,vars]
  nvars<-ncol(data)
   data<-data[,vars]
  if (m==1){
    # make all into factors
    syndata<-syndata[,vars]
    for (i in 1:nvars) {
      if (is.numeric(data[,i])){
        grpd<-group_num(data[,i],syndata[,i],ngroups)
        data[,i]<-grpd[[1]];syndata[,i]<-grpd[[2]]
      }
      else if(is.character(data[,i])){
        data[,i]<-factor(data[,i])
        syndata[,i]<-factor(syndata[,i],levels<-levels(data[,i]))
      }
      if (with.missings==TRUE){
        nobs.missings<-sum(apply(is.na(data),1,sum)>0)
        nsyn.missings<-sum(apply(is.na(syndata),1,sum)>0)
        if (any(is.na(data[,i])))  syndata[,i]<-addNA(syndata[,i]) ### makes missings into part of factors
        if (any(is.na(data[,i]))) data[,i]<-addNA(data[,i]) ### makes missings into part of factors 

      }
    }
    nobs.missings<-sum(apply(is.na(data),1,sum)>0)
    nsyn.missings<-sum(apply(is.na(syndata),1,sum)>0)
    tab.syn<-table(syndata)
    tabd<-table(data)
    totcells<-length(tabd)
    df<-totcells - sum(tabd+tab.syn==0)-1
    diff<-(tab.syn-tabd)
    expect<-(tab.syn+tabd)/2
    UtabFT<-4*sum((tab.syn^(0.5)-tabd^(0.5))^2)
    tabsq<-diff^2/expect
    tabsq[expect==0]<-0
    UtabVW<-sum(tabsq)
    ratioFT=UtabFT/df
    stdFT=(UtabFT-df)/sqrt(2*df)
    ratioVW=UtabVW/df
    stdVW=(UtabVW-df)/sqrt(2*df)
    tab.zdiff<-diff/sqrt(expect)
    nempty= sum(tabd+tab.syn==0)
  }
  else{############################################# now for m>1
    df<-ratioFT<-ratioVW<-stdVW<-stdFT<-UtabFT<-UtabVW<-nempty<-nsyn.missings<-rep(NA,m)
    tab.syn<-tab.zdiff<-as.list(1:m)
    data0<-data
    nobs.missings<-sum(apply(is.na(data),1,sum)>0)
    for (ii in 1:m){
      #
      # make all into factors
      #
      syndata<-object$syn[[ii]][,vars]
      nsyn.missings[ii]<-sum(apply(is.na(syndata),1,sum)>0)
      for (i in 1:nvars) {
        if (is.numeric(syndata[,i])){
          grpd<-group_num(data0[,i],syndata[,i],ngroups)
          data[,i]<-grpd[[1]];syndata[,i]<-grpd[[2]]
        }
        else if(is.character(data[,i])){
          data[,i]<-factor(data[,i])
          syndata[,i]<-factor(syndata[,i],levels<-levels(data[,i]))
        }
        if (with.missings==TRUE){
          if (any(is.na(data0[,i]))){syndata[,i]<-addNA(syndata[,i])} ### makes missings into part of factors
          if (any(is.na(data0[,i]))) {data[,i]<-addNA(data[,i])}
        }
      }
      tab.syn[[ii]]<-table(syndata)
      tabd<-table(data)

      totcells<-length(tabd)
      df[[ii]]<-totcells - sum(tabd+tab.syn[[ii]]==0)-1
      diff<-(tab.syn[[ii]]-tabd)
      expect<-(tab.syn[[ii]]+tabd)/2
      UtabFT[[ii]]<-4*sum((tab.syn[[ii]]^(0.5)-tabd^(0.5))^2)
      tabsq<-diff^2/expect
      tabsq[expect==0]<-0
      UtabVW[[ii]]<-sum(tabsq)
      ratioFT[[ii]]=UtabFT[[ii]]/df[[ii]]
      stdFT[[ii]]=(UtabFT[[ii]]-df[[ii]])/sqrt(2*df[[ii]])
      ratioVW[[ii]]=UtabVW[[ii]]/df[[ii]]
      stdVW[[ii]]=(UtabVW[[ii]]-df[[ii]])/sqrt(2*df[[ii]])
      tab.zdiff[[ii]]<-diff/sqrt(expect)
      nempty[[ii]]= sum(tabd+tab.syn[[ii]]==0)
    }
  }    
  res <- list(m = m, UtabFT = unlist(UtabFT), UtabVW = unlist(UtabVW), 
              df = unlist(df), ratioFT = unlist(ratioFT), stdFT = unlist(stdFT),
              ratioVW = unlist(ratioVW), stdVW = unlist(stdVW),
              nempty = unlist(nempty), 
              with.missings = with.missings,
              nobs.missings = nobs.missings,
              nsyn.missings = nsyn.missings,
              tab.obs = tabd, tab.syn = tab.syn, tab.zdiff = tab.zdiff,
              digits = digits, print.zdiff = print.zdiff, print.tables = print.tables)
  
  class(res) <- "utility.tab"
  return(res)
}


#-------------------------------group_num---------------------------
#
# function to make groups of any continuous variables
#
group_num<-function(x1,x2,n=5,cont.na=NULL){
  #
  # makes 2 cont variable into a factor of n groups with same groupings
  # determined by first one
  # groups will often be of uneven size because of missing values
  # there may even be fewer than n groups if there are big sets of repeats
  #
  #
  if (!is.numeric(x1)) stop ("x must be numeric \n")
  xn<-x1[!(x1 %in% cont.na | is.na(x1))  ]
  xr<-rank(xn)
  ptiles=round(c(length(xr)*(1:(n-1))/n,max(xr)))
  #ptiles
  ptiles<-xn[order(xn)][ptiles]
  #ptiles
  ptiles<-ptiles[!duplicated(ptiles)]##  to allow for repeated values
  xnew=rep(1,length(xn))
  for (i in 2:length(ptiles)) xnew[xn>ptiles[i-1] & xn<=ptiles[i]]<-i
  #xnew<-factor(xnew)
  res1<-x1
  res1[!(x1 %in% cont.na | is.na(x1))]<-xnew
  #
  # now next one
  #
  xn<-x2[!(x2 %in% cont.na | is.na(x2))  ]
  xnew=rep(1,length(xn))
  for (i in 2:length(ptiles)) xnew[xn>ptiles[i-1] & xn<=ptiles[i]]<-i
  #xnew<-factor(xnew)
  res2<-x2
  res2[!(x2 %in% cont.na | is.na(x2))]<-xnew
  #
  # now make into factors including missing and cont.na data
  #
  nlev<-length(table(c(res1[!(x1 %in% cont.na | is.na(x1))],res2[!(x2 %in% cont.na | is.na(x2))])))
  labels=paste("Gp",1:nlev)
  if (any(is.na(c(res1,res2)))) {
    res1[is.na(res1)]<-nlev+1
    res2[is.na(res2)]<-nlev+1
    labels<-c(labels,"missing")
  }
  nlev<-length(labels)
  if (!is.null(cont.na)) {
    exlevs<-as.numeric(cont.na[!is.na(cont.na)])
    for ( i in 1:length(exlevs)){
      res1[x1==exlevs[i]]<-nlev+i
      res2[x2==exlevs[i]]<-nlev+i
    }
    labels<-c(labels,paste("cont.na",as.character(exlevs)))
    #print(table(res1))
  }
  res1<-factor(res1,labels=labels)
  res2<-factor(res2,labels=labels)
  list(res1,res2)
}          

