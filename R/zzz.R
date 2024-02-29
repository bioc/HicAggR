hicaggrStartupMessage <- function()
{
# Startup message obtained as 
# > figlet HicAggR
  msg <- c(paste0(
" __  __              ______                  ____       
 /\ \/\ \  __        /\  _  \                /\  _`\     
 \ \ \_\ \/\_\    ___\ \ \L\ \     __      __\ \ \L\ \   
  \ \  _  \/\ \  /'___\ \  __ \  /'_ `\  /'_ `\ \ ,  /   
   \ \ \ \ \ \ \/\ \__/\ \ \/\ \/\ \L\ \/\ \L\ \ \ \\ \  
    \ \_\ \_\ \_\ \____\\ \_\ \_\ \____ \ \____ \ \_\ \_\
     \/_/\/_/\/_/\/____/ \/_/\/_/\/___L\ \/___L\ \/_/\/ /
                                   /\____/ /\____/       
                                   \_/__/  \_/__/        .......version ", 
packageVersion("HicAggR")),
"\nType 'citation(\"HicAggR\")' for citing this R package in publications.")
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  msg <- hicaggrStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'HicAggR' version", packageVersion("HicAggR"))
  packageStartupMessage(msg)
  invisible()
}
