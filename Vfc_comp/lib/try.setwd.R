try.setwd <- function(x) {
  out <- tryCatch( {
    message("*** Setting working directory in [script_path]/results.")
    setwd(x)
  },
  error = function(cond){
    message("Cannot set working directory correctly. Please check directory rights.")
    message(cond)
    return(NA)
  },
  warning = function(cond){
    message("Warning issued:")
    message(cond)
    return(NULL)
  }, 
  finally = {message("*** Working directory is:")}
  )
  
  return(getwd())
}