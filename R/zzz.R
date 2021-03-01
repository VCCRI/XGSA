
.onLoad <- function(libname, pkgname) {

  op <- options()
  op.xgsa <- list(
    xgsa.name = "xgsa name",
    xgsa.date = "220816"
  )
  toset <- !(names(op.xgsa) %in% names(op))
  if(any(toset)) options(op.xgsa[toset])


  biomartID <- "ENSEMBL_MART_ENSEMBL"

  hostID <-  "www.ensembl.org"
  #hostID <-  "asia.ensembl.org"
  #hostID <-  "aug2017.archive.ensembl.org"

  supportedSpecies <<- find_supported_datasets()

  invisible()
}

