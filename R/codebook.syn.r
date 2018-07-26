###-----codebook.syn-------------------------------------------------------

codebook.syn <- function(data, maxlevs = 3)
{
# function to decribe variables in a data frame 
 
 if (!(is.data.frame(data))) stop("codebook.syn() requires a data frame as a parameter.\n", call. = FALSE)
 n <- dim(data)[[1]]
 p <- dim(data)[[2]]
 
 # calculate number and % of missing and non-missing values
 nmiss <- sapply(data, function(x) length(x[is.na(x)]))
 perctmiss <- round(nmiss/n * 100, 2)
 nok <- sapply(data, function(x) length(x[!is.na(x)]))
 ndistinct <- sapply(data, function(x) length(table(x)))
 dfclass <- sapply(data, class)

 fortab2 <- details <- rep("", length(nmiss))
 
 for (i in 1:p) {
   if (class(data[,i]) == "character") details[i] <- paste("Max length: ", 
                                            max(nchar(data[,i])), sep = "")
   if (class(data[,i]) == "numeric") details[i] <- paste("Range: ", 
     min(data[!is.na(data[,i]), i]), " - ", max(data[!is.na(data[,i]),i]), 
     sep = "")
   if (class(data[,i]) == "factor" & ndistinct[i] > maxlevs ) { 
     details[i] <- "See table in labs"
     fortab2[i] <- paste("'", paste(names(table(data[,i])), collapse = "' '"), 
                         "'", sep = "")
   }
   if (class(data[,i]) == "factor" & ndistinct[i] <= maxlevs ) details[i] <-
     paste("'", paste(names(table(data[,i])), collapse = "' '"),"'", sep = "")
 }

 if (any(sapply(data, class) == "factor" & ndistinct >= maxlevs )) {
   vnum <- (1:p)[sapply(data, class) == "factor" & ndistinct >= maxlevs]
   tabs2 <- vector("list", length(vnum))
   names(tabs2) <- names(data)[vnum]
   for (i in 1:length(vnum)) {
     tabs2[[i]] <- data.frame(label = names(table(data[, vnum[i]])))
   }
 } else tabs2 <- NULL

 result <- data.frame(variable = names(data), class = dfclass, 
                      nmiss = nmiss, perctmiss = perctmiss, 
                      ndistinct = ndistinct, details = details)
 rownames(result) <- 1:p
 list(tab = result, labs = tabs2)
}


