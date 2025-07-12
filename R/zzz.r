.onAttach <- function(...) {
  info <- c("New version of synthpop (1.9-0) with disclosure functions",
            "see disclosure.pdf for details and NEWS file for other changes",
           "","Find out more at https://www.synthpop.org.uk/")
    packageStartupMessage(paste(strwrap(info), collapse = "\n"))
}
