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
utils::packageVersion("HicAggR"))
)
    if(!interactive())
        msg[1] <- paste("Package 'HicAggR' version", 
        utils::packageVersion("HicAggR"))
    packageStartupMessage(msg)
    invisible()
}
