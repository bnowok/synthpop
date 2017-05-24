###-----replicated.uniques-------------------------------------------------
# unique units in the synthesised data that replicates unique real units
# (+number +percent of all observations)

replicated.uniques <- function(object, data, exclude = NULL){
  # check names of vars to be excluded
  if (!is.null(exclude)) {
    exclude.cols <- match(exclude, colnames(data))
    if (any(is.na(exclude.cols))) stop("Unrecognized variable(s) in exclude parameter for uniques: ",
      paste(exclude[is.na(exclude.cols)],collapse=", "), call. = FALSE)
  }
  
  # exclude vars from the original data
  if (!is.null(exclude)) data <- data[,-exclude.cols, drop = FALSE] else data <- data
  # extract uniques from the original data
  uReal <- data[!(duplicated(data) | duplicated(data,fromLast=TRUE)),]
  no.uniques <- nrow(uReal)

  if (is.null(no.uniques) || no.uniques == 0) {
    no.uniques <- 0
    no.duplicates <- per.duplicates <- rep(0,object$m)
    if (object$m == 1) rm.Syn <- rep(FALSE,nrow(object$syn))
    if (object$m > 1) rm.Syn <- matrix(FALSE,nrow=nrow(object$syn[[1]]),ncol=object$m)

  } else {

    if (object$m == 1){
      if (!is.null(exclude)) Syn <- object$syn[,-exclude.cols, drop = FALSE] else Syn <- object$syn
      rm.Syn <- rep(FALSE,nrow(Syn))
      i.unique.Syn <- which(!(duplicated(Syn) | duplicated(Syn,fromLast=TRUE)))
      uSyn <- Syn[i.unique.Syn,]
      uAll <- rbind.data.frame(uReal,uSyn)
      dup.of.unique <- duplicated(uAll)[(nrow(uReal)+1):nrow(uAll)]
      rm.Syn[i.unique.Syn] <- dup.of.unique
      no.duplicates <- sum(dup.of.unique)
    }

    if (object$m > 1){
    rm.Syn <- matrix(FALSE,nrow=nrow(object$syn[[1]]),ncol=object$m)
      for (i in 1:object$m){
        if (!is.null(exclude)) Syn <- object$syn[[i]][,-exclude.cols,drop = FALSE] else Syn <- object$syn[[i]]
        i.unique.Syn <- which(!(duplicated(Syn) | duplicated(Syn,fromLast=TRUE)))
        uSyn <- Syn[i.unique.Syn,]
        uAll <- rbind.data.frame(uReal,uSyn)
        dup.of.unique <- duplicated(uAll)[(nrow(uReal)+1):nrow(uAll)]
        rm.Syn[i.unique.Syn,i] <- dup.of.unique
      }
      no.duplicates <- colSums(rm.Syn)
    }
    per.duplicates <- no.duplicates/nrow(data)*100
  }

  return(list(replications = rm.Syn, no.uniques = no.uniques,
    no.replications = no.duplicates, per.replications = per.duplicates))
}


###-----sdc----------------------------------------------------------------
# sdc - statistical disclosure control:
# labeling, removing unique replicates of unique real individuals

sdc <- function(object, data, label = NULL, rm.replicated.uniques = FALSE, 
 uniques.exclude = NULL, recode.vars = NULL, bottom.top.coding = NULL,
 recode.exclude = NULL){

 if (!is.null(recode.vars)){
   if (!is.null(bottom.top.coding) && !is.list(bottom.top.coding)) 
       bottom.top.coding <- list(bottom.top.coding)
   if (!is.null(recode.exclude) && !is.list(recode.exclude)) 
       recode.exclude <- list(recode.exclude)
   if (length(bottom.top.coding)!=length(recode.vars) | 
       any(sapply(bottom.top.coding,length)!=2)) 
       stop("Bottom and top codes have to be provided for each variable in recode.vars.\nUse NA if there is no need for bottom or top recoding.\nFor more than one variable to be recoded provide a list of two-element vectors, e.g. list(c(0,60),c(NA,5000))",
       call. = FALSE)
   if (!is.null(recode.exclude) && length(bottom.top.coding)!=length(recode.exclude))
       stop("recode.exclude have to include codes for each variable in recode.vars.\nUse NA if all values should be considered for recoding.\nFor more than one variable to be recoded provide a list, e.g. list(NA,c(NA,-8)).",
       call. = FALSE)
 }

 if (object$m==1){
   if (!is.null(recode.vars)) {
     cols <- match(recode.vars,colnames(object$syn)) 
     for(i in cols){
       j <- match(i,cols) 
       recoded <- bottom.top.recoding(object$syn[,i],bottom.top.coding[[j]][1],
         bottom.top.coding[[j]][2],recode.exclude[[j]])
       object$syn[,i] <- recoded$x
       cat("\n",recode.vars[j],": no. of bottom-coded values - ",
         recoded$no.recoded.bottom,", no. of top-coded values - ",
         recoded$no.recoded.top, sep="")
     }
   cat("\n")
   }
   if (rm.replicated.uniques) {
     du <- replicated.uniques(object, data, exclude = uniques.exclude) 
     object$syn <- object$syn[!du$replications,]
     cat("no. of replicated uniques: ", du$no.replications, "\n", sep="")
   }
   if (!is.null(label)) object$syn <- cbind.data.frame(flag=label, object$syn)
 }
  
 if (object$m>1){
   if (!is.null(recode.vars)) {
     cols <- match(recode.vars,colnames(object$syn[[1]])) 
     for(k in 1:object$m){
       cat("\nm =",k)
       for(i in cols){
         j <- match(i,cols) 
         recoded <- bottom.top.recoding(object$syn[[k]][,i],bottom.top.coding[[j]][1],
           bottom.top.coding[[j]][2],recode.exclude[[j]])
         object$syn[[k]][,i] <- recoded$x
         cat("\n",recode.vars[j], ": no. of bottom-coded values - ",
         recoded$no.recoded.bottom, ", no. of top-coded values - ",
         recoded$no.recoded.top, sep="")
       }
     }
   cat("\n")
   }
   if (rm.replicated.uniques) {
     du <- replicated.uniques(object, data, exclude = uniques.exclude) 
     for (i in 1:object$m){ 
       object$syn[[i]] <- object$syn[[i]][!du$replications[,i],]
     }
   cat("no. of replicated uniques: ", 
     paste0(du$no.replications, collapse=", "),"\n", sep="")
   }
   if (!is.null(label)) object$syn <- mapply(cbind.data.frame, flag=label,
     object$syn, SIMPLIFY=FALSE, USE.NAMES=FALSE)
 }
 return(object) 
}


###---- bottom.top.recoding -----------------------------------------------

bottom.top.recoding <- function(x,bottom,top,exclude=NULL){
  below <- which(x < bottom & !x%in%exclude); no.below <- length(below)
  above <- which(x > top & !x%in%exclude); no.above <- length(above)
  x[below] <- bottom
  x[above] <- top
  return(list(x=x,no.recoded.bottom=no.below,no.recoded.top=no.above))
}


###---- read.obs ----------------------------------------------------------

read.obs <- function(file, convert.factors = TRUE, lab.factors = FALSE, 
                     export.lab = FALSE, ...){

 pos <- regexpr("\\.([[:alnum:]]+)$", file)
 ext <- ifelse(pos > -1L, substring(file, pos + 1L), "")

 if (ext=="sav") {
   real.data <- read.spss(file, to.data.frame = FALSE, 
                  use.value.labels = convert.factors, 
                  trim.factor.names = TRUE, ...)
  # trim.factor.names=T - trim trailing spaces from factor levels
  # use.value.labels=F -> to prevent combining factor levels with missing labels
  # for read.spss -> cbind(value=attributes(data)$label.table$...)
  # {Hmisc} real.data <- spss.get(file,use.value.labels=TRUE,max.value.labels=20)
   varlab <- attr(real.data, "variable.labels")
   vallab <- attr(real.data, "label.table")
   vallab <- lapply(vallab,rev)
   codepage <- attr(real.data,"codepage")
   real.data <- as.data.frame(real.data)
   attr(real.data, "variable.labels") <- varlab
   attr(real.data, "label.table") <- vallab
   if (codepage > 500) attr(real.data, "codepage") <- codepage
   
   labs <- list(var.lab=varlab,val.lab=vallab)
   if (export.lab) dput(labs,"SPSSlabels")
   
   if (convert.factors == FALSE & lab.factors == TRUE){  
   # convert completly labeled variables into factors
     ff <- !sapply(vallab,is.null)
     suppressWarnings(llall <-
       sapply(vallab, function(x) !(sum(!is.na(as.numeric(names(x))))!=0)))
     factors <- which(ff & llall)
     for (i in factors){
       real.data[,i] <- factor(real.data[,i],levels=vallab[[i]])
     }
   }
 
 } else if (ext=="dta") {
   real.data <- read.dta(file, convert.factors = convert.factors, ...)
   
   varlab <- attr(real.data, "val.labels")
   vallab <- attr(real.data, "label.table")
   labs <- list(var=varlab,val=vallab)
   if (export.lab) dput(labs,"Statalabels")
   
   if (convert.factors == FALSE & lab.factors == TRUE){  
   # convert completly labeled variables into factors
      ff <- which(varlab != "")
      suppressWarnings(varlaball <-
       sapply(vallab, function(x) !(sum(!is.na(as.numeric(names(x))))!=0)))
      factors <- ff[varlaball[varlab[ff]]]
      for (i in factors){
        real.data[,i] <- factor(real.data[,i],levels=vallab[[varlab[i]]])
      }
   }

  # for read.dta -> cbind(value=attributes(data)$label.table$...)
  # convert.factors=T, convert.underscore=F, warn.missing.labels=T, missing.type=T
  # {Hmisc} real.data <- stata.get(file,...)

 } else if (ext=="xpt") {
   real.data <- read.xport(file)
  # {Hmisc} real.data <- sasxport.get(file, ...)

 } else if (ext=="csv") {
   real.data <- read.csv(file, header=TRUE, ...)

 } else if (ext=="txt") {
   real.data <- read.table(file, header=TRUE, ...)

 } else if (ext=="tab") {
   real.data <- read.table(file, header=TRUE, sep="\t")

 } else {
   stop(".",ext," is an unrecognized data format",call.=FALSE)
 }

 attr(real.data,"filetype") <- ext
 return(real.data)
}
# R files (*.RData, *.rda)
# load(".rda") - don't assign to an object!


###---- write.syn ---------------------------------------------------------

write.syn <- function(object, filename,
  filetype = c("SPSS","Stata","SAS","csv","tab","rda","RData","txt"),
  convert.factors = "numeric", data.labels = NULL,
  save.complete = TRUE, extended.info = TRUE, ...){

# pos <- regexpr("\\.([[:alnum:]]+)$", file)
# ext <- ifelse(pos > -1L, substring(file, pos + 1L), "")
# without.ext <- strsplit(file,"\\.")[[1]][1]      #! will include path

 m <- object$m
 if (m == 1) object$syn <- list(object$syn)
 filetype <- match.arg(filetype)
 if (is.null(data.labels)) data.labels <- list(var.lab = object$var.lab,
                                               val.lab = object$val.lab)

 if (filetype=="SPSS"){
   if (m==1) {
     f1 <- paste0(filename,".sps"); f2 <- paste0(filename,".txt")
   } else {
     f1 <- paste0(filename,"_",1:m,".sps"); f2 <- paste0(filename,"_",1:m,".txt")
   }
   for (i in 1:m) {
     #write.foreign(object$syn[[i]], codefile=f1, datafile=f2, package=filetype, ...)
     write.syn.SPSS(object$syn[[i]], codefile = f1[i], datafile = f2[i],
       varnames = names(object$syn[[i]]), data.labels = data.labels, ...)
   }
 } else if (filetype == "Stata"){
   if (m==1) f1 <- paste0(filename,".dta") else f1 <- paste0(filename,"_",1:m,".dta")
   for (i in 1:m){
     write.dta(object$syn[[i]], file = f1[i], convert.factors = convert.factors, ...)
       #!### check why default convert.factors="labels" cuts the names
   }
 } else if (filetype == "SAS"){
   if (m==1) {
     f1 <- paste0(filename,".sas"); f2 <- paste0(filename,".txt")
   } else {
     f1 <- paste0(filename,"_",1:m,".sas"); f2 <- paste0(filename,"_",1:m,".txt")
   }
   for (i in 1:m){
     write.foreign(object$syn[[i]], codefile = f1[i], datafile = f2[i], package = filetype, ...)
   }
 } else if (filetype == "rda" | filetype == "RData"){
   if (m==1) f1 <- paste(filename,filetype,sep=".") else f1 <- paste0(filename,"_",1:m,".",filetype)
   for (i in 1:m){
     syn <- object$syn[[i]]
     save(syn, file = f1[i],...)
   }
 } else if (filetype == "csv"){
   if (m==1) f1 <- paste0(filename,".csv") else f1 <- paste0(filename,"_",1:m,".csv")
   for (i in 1:m){
     write.csv(object$syn[[i]], file = f1[i], row.names = FALSE, ...)
   }
 } else if (filetype == "txt"){
   if (m==1) f1 <- paste0(filename,".txt") else f1 <- paste0(filename,"_",1:m,".txt")
   for (i in 1:m){
     write.table(object$syn[[i]], file = f1[i], row.names = FALSE, ...)
   }
 } else if (filetype == "tab"){
   if (m==1) f1 <- paste0(filename,".tab") else f1 <- paste0(filename,"_",1:m,".tab")
   for (i in 1:m){
     write.table(object$syn[[i]],file=f1[i], row.names=FALSE, sep="\t", ...)
   }
 }

 call <- match.call()
 infofile <- paste("synthesis_info_",filename,".txt",sep="")
  #---
 sink(infofile)
 cat("Date saved:",format(Sys.time(), "%d %b %Y, %H:%M"), "\n")
 cat("Data frame with original data:", deparse(object$call$data), "\n")
 cat("Number of synthetic data sets:", m, "\n")
 cat("Output file(s):")
 if (filetype=="SPSS" | filetype=="SAS"){
   cat(paste0("(",filetype,") "))
   cat(paste(f1, f2, collapse=" "))
 } else {
   cat(paste0("(",filetype,") ",f1))
 }
 if (save.complete) {                            
   addname <- paste0("synobject_",filename,".RData")
   save(object,file=addname)
   cat("\nAdditional file: ", addname, sep="")
 }
 if (extended.info) {
   cat("\nMethods used:\n")
   print(object$method)
   cat("Seed used:",object$seed,"\n")
 }
 sink()
 #---

 cat("Synthetic data exported as ",filetype," file(s).",sep="")
 cat("\nInformation on synthetic data written to\n ",
     paste(getwd(),"/",infofile,sep=""),"\n")

}


###---- write.syn.SPSS ----------------------------------------------------

write.syn.SPSS <- function (df, datafile, codefile, varnames = names(df),
  data.labels = NULL, ...)
{
  varlabels <- data.labels$var.lab
  vallabels <- data.labels$val.lab
  if (is.null(varnames)) varnames <- names(df)
  
  dfn <- df
  for (i in 1:ncol(dfn)){
    if (!is.null(vallabels[[varnames[i]]])){
      dfn[,i] <- mapvalues(dfn[,i], from = names(vallabels[[varnames[i]]]), 
        to = vallabels[[varnames[i]]])
    }
  }
  write.table(dfn, file = datafile, row.names = FALSE, col.names = FALSE,
    sep = ",", quote = FALSE, na = "", eol = ",\n")
  varnames  <- gsub("[^[:alnum:]_\\$@#]", "\\.", varnames)
  if (is.null(varlabels)) varlabels <- varnames 
  dl.varnames <- varnames
  if (any(chv <- sapply(df, is.character))) {
    lengths <- sapply(df[chv], function(v) max(nchar(v)))
    if (any(lengths > 255L))
      stop("Cannot handle character variables longer than 255")
    lengths <- paste0("(A", lengths, ")")
    star <- ifelse(c(TRUE, diff(which(chv) > 1L)), " *"," ")
    dl.varnames[chv] <- paste(star, dl.varnames[chv], lengths)
  }
  cat("DATA LIST FILE=", adQuote(paste(getwd(), datafile, sep = "/")),
    " free (\",\")\n", file = codefile)
  cat("/", dl.varnames, " .\n\n", file = codefile, append = TRUE)
  cat("VARIABLE LABELS\n", file = codefile, append = TRUE)
  cat(paste(varnames, adQuote(varlabels[varnames]), "\n"), ".\n", file = codefile,
    append = TRUE)
  factors <- sapply(df, is.factor) & !sapply(vallabels[varnames], is.null)
  if (any(factors)) {
    cat("\nVALUE LABELS\n", file = codefile, append = TRUE)
    for (v in which(factors)) {
      cat("/\n", file = codefile, append = TRUE)
      cat(varnames[v], " \n", file = codefile, append = TRUE, sep = "")
      levs     <- vallabels[[varnames[v]]]
      levslabs <- names(vallabels[[varnames[v]]])
      cat(paste((levs), adQuote(levslabs), "\n", sep = " "),
        file = codefile, append = TRUE)
    }
    cat(".\n", file = codefile, append = TRUE)
  }
  cat("\nEXECUTE.\n", file = codefile, append = TRUE)
}


###---- adQuote -----------------------------------------------------------

adQuote <- function (x) paste("\"", x, "\"", sep = "")


###---- maplabs -----------------------------------------------------------
maplabs <- function (x, from, to) 
{
  if (is.factor(x)) {
    levels(x) <- maplabs(levels(x), from, to)
    return(x)
  }
  mapidx   <- match(x, from)
  mapidxNA <- is.na(mapidx)
  x[!mapidxNA] <- to[mapidx[!mapidxNA]]
  return(x)
}


###---- mapvalues -----------------------------------------------------------
mapvalues <- function (x, from, to, warn_missing = TRUE) 
{
    if (length(from) != length(to)) {
        stop("`from` and `to` vectors are not the same length.")
    }
    if (!is.atomic(x)) {
        stop("`x` must be an atomic vector.")
    }
    if (is.factor(x)) {
        levels(x) <- mapvalues(levels(x), from, to, warn_missing)
        return(x)
    }
    mapidx <- match(x, from)
    mapidxNA <- is.na(mapidx)
    from_found <- sort(unique(mapidx))
    if (warn_missing && length(from_found) != length(from)) {
        message("The following `from` values were not present in `x`: ", 
            paste(from[!(1:length(from) %in% from_found)], collapse = ", "))
    }
    x[!mapidxNA] <- to[mapidx[!mapidxNA]]
    x
}
