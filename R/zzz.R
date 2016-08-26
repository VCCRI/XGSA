
#.onLoad <- function(libname, pkgname) {


#  biomart_ID <<- "ENSEMBL_MART_ENSEMBL"
#  host_ID <<-  "www.ensembl.org"
#  supported_species <<- find_supported_datasets()


#  packageStartupMessage("Initialised Ensembl connection")

#  invisible()
#}

.onLoad <- function(libname, pkgname) {

  op <- options()
  op.xgsa <- list(
    xgsa.name = "xgsa name",
    xgsa.date = "220816"
  )
  toset <- !(names(op.xgsa) %in% names(op))
  if(any(toset)) options(op.xgsa[toset])


  biomartID <<- "ENSEMBL_MART_ENSEMBL"
  hostID <<-  "www.ensembl.org"
  supportedSpecies <<- find_supported_datasets()

  invisible()
}

