dllunload <- function(){
  dyn.unload(
    system.file("libs", "x64", "RCGAL.dll", package = "RCGAL")
  )
}

makedoc <- function(){
  roxygen2::roxygenise(load_code = roxygen2::load_installed)
}
