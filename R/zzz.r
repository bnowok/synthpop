.onAttach <- function(...) {
  info <- "Find out more at https://www.synthpop.org.uk/"
    packageStartupMessage(paste(strwrap(info), collapse = "\n"))
}
