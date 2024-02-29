hicaggrStartupMessage <- function()
{
# Startup message obtained as 
# > figlet HicAggR
  msg <- c(paste0(
"    _/    _/  _/              _/_/                        _/_/_/    
   _/    _/        _/_/_/  _/    _/    _/_/_/    _/_/_/  _/    _/   
  _/_/_/_/  _/  _/        _/_/_/_/  _/    _/  _/    _/  _/_/_/      
 _/    _/  _/  _/        _/    _/  _/    _/  _/    _/  _/    _/     
_/    _/  _/    _/_/_/  _/    _/    _/_/_/    _/_/_/  _/    _/      
                                       _/        _/                 
                                  _/_/      _/_/                    ",
".......version ",
packageVersion("HicAggR")),
"\nType 'citation(\"HicAggR\")' for citing this R package in publications.",
"\n\nImported---------------------",
paste0("\nInteractionSet--------------",
ifelse("package:InteractionSet"%in%search(),"\u{2713}","\u{0058}"))
)
  return(msg)
}

.onAttach <- function(lib, pkg)
{
    suppressMessages(library(InteractionSet))
    msg <- hicaggrStartupMessage()
    if(!interactive())
        msg[1] <- paste("Package 'HicAggR' version", packageVersion("HicAggR"))
    packageStartupMessage(msg)
    invisible()
}
