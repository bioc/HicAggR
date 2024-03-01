.onAttach <- function(lib, pkg)
{
        # suppressPackageStartupMessages(requireNamespace(InteractionSet))
msg <- c(paste0(
"_/    _/  _/              _/_/                        _/_/_/   
_/    _/        _/_/_/  _/    _/    _/_/_/    _/_/_/  _/    _/   
_/_/_/_/  _/  _/        _/_/_/_/  _/    _/  _/    _/  _/_/_/      
_/    _/  _/  _/        _/    _/  _/    _/  _/    _/  _/    _/     
_/    _/  _/    _/_/_/  _/    _/    _/_/_/    _/_/_/  _/    _/      
                                    _/        _/                 
                                _/_/      _/_/     ",
"version: ",
utils::packageVersion("HicAggR")),
"\nType 'citation(\"HicAggR\")' 
for citing this R package in publications."
# ,
# "\n\nImported---------------------",
# paste0("\nInteractionSet              ",
# ifelse("package:InteractionSet" %in% search(),"\u{2713}","\u{0058}"))
)
    if(!interactive())
        msg[1] <- paste("Package 'HicAggR' version", 
        utils::packageVersion("HicAggR"))
    packageStartupMessage(msg)
    invisible()
}
