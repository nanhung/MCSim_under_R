# Download all files form this repository (the current version of MCSim is 6.0.1)
# Check the chek whether the compiler is in the PATH by using
# Sys.getenv("PATH") 

if(!require(pkgbuild)) install.packages("pkgbuild")
library(pkgbuild)

set_PATH <- function(PATH = "c:/Rtools/mingw_32/bin"){

    if (Sys.info()[['sysname']] == "Windows") {
    if(!pkgbuild::find_rtools()){
      warning("Please make sure you had installed Rtools. Or, assign the PATH to Rtools.")
    }
    if(Sys.which("gcc") == ""){
      Sys.setenv(PATH = paste(PATH, Sys.getenv("PATH"), sep=";"))
    }
  }
  
  # The macos used clang as default, the following command is used to switch to GCC
  # Sys.setenv(PATH = paste("/usr/local/bin", Sys.getenv("PATH"), sep=";"))
  # port select --list gcc
  # sudo port select --set gcc mp-gcc8
  
  # Check the GNU compiler 
  Sys.which("gcc")
  system('gcc -v')
}

makemod <- function(){
  if(Sys.which("gcc") == ""){
    stop("Please setting the PATH of compiler")
  }
  system(paste("gcc -o ./MCSim/mod.exe ./MCSim/mod/*.c ", sep = "")) 
  
  if(file.exists("MCSim/mod.exe")){
    message("The mod.exe had been created.")
  }
}
set_PATH()
makemod()

makemcsim <- function(model, deSolve = F){
  exe_file <- paste0("mcsim.", model, ".exe")
  #if(file.exists(exe_file)) stop(paste0("* '", exe_file, "' had been created."))
  if(file.exists(model)) {
    invisible(file.copy(from = paste0(getwd(),"/", model), to = paste0(getwd(),"/modeling/", model)))
    invisible(file.remove(model))
  }
  if (deSolve == T){
    system(paste("./MCSim/mod.exe -R modeling/", model, " ", model, ".c", sep = "")) 
    system (paste0("R CMD SHLIB ", model, ".c")) # create *.dll files
    dyn.load(paste(model, .Platform$dynlib.ext, sep="")) # load *.dll
    source(paste0(model,"_inits.R"))
  } else {
    system(paste("./MCSim/mod.exe modeling/", model, " ", model, ".c", sep = "")) 
    system(paste("gcc -O3 -I.. -I./MCSim/sim -o mcsim.", model, ".exe ", model, ".c ./MCSim/sim/*.c -lm ", sep = ""))  
    invisible(file.remove(paste0(model, ".c")))
    if(file.exists(exe_file)) message(paste0("* Created executable program '", exe_file, "'.")) 
  }
}

mcsim <- function(model, input){
  exc = paste0("mcsim.", model, ".exe")
  if (file.exists(exc) == F) {
    makemcsim(model)
    if (file.exists(exc) == F) {
    stop("* '", exc, "' is not exist .")
    }
  }
  
  if(file.exists(input)) {
    invisible(file.copy(from = paste0(getwd(),"/", input), to = paste0(getwd(),"/modeling/", input)))
    invisible(file.remove(input))
  }
  tx  <- readLines(paste0("modeling/", input))
  MCMC_line <- grep("MCMC", x=tx)
  if (length(MCMC_line) != 0){
    #file_defore <- list.files()
    RandomSeed <- runif(1, 0, 2147483646)
    tx2 <- gsub(pattern = "10101010", replace = paste(RandomSeed), x = tx)
    writeLines(tx2, con=paste0("modeling/", input))
    system(paste("./mcsim.", model, ".exe ", "modeling/", input, sep = ""))
    #file_after <- list.files() # exist a bug if file had been created
    #outfile <- setdiff(file_after,file_defore)[1]
    outfile <- "sim.out"
    tx2 <- gsub(pattern = ",0,", replace = ",1,", x = tx)
    checkfile <- "chk.out"
    tx3 <- gsub(pattern = paste0("\"", outfile, "\",\"\""), 
                replace = paste0("\"", checkfile, "\",\"", outfile, "\""), 
                x = tx2)
    writeLines(tx3, con=paste0("modeling/", input))
    system(paste("./mcsim.", model, ".exe ", "modeling/", input, sep = ""))
    writeLines(tx, con=paste0("modeling/", input))
    message(paste0("* Create '", checkfile, "' from the last iteration."))
    #invisible(file.remove(paste0(outfile, ".kernel")))
    df <- read.delim("sim.out")
  } else {
    system(paste("./mcsim.", model, ".exe ", "modeling/", input, sep = ""))
    df <- read.delim("sim.out", skip = 1)
  }
  return(df)
}

clear <- function(){
  files <- c(dir(pattern = c("*.out")), 
             dir(pattern = c("*.R.exe")),
             dir(pattern = c("*.R.so")),
             dir(pattern = c("*.R.o")),
             dir(pattern = c("*.R.dll")),
             dir(pattern = c("*.R.c")),
             dir(pattern = c("*.R_inits.R")),
             dir(pattern = c("*.perks")))
  invisible(file.remove(files))
}

report <- function(){
  cat("\n\n-----Report started line-----\n\n")
  cat(Sys.getenv("PATH"), "\n")
  print(Sys.which("gcc"))
  system('gcc -v')
  cat("\n-----Report ended line-----\n\n")
}

## The following is test function ####
plotmcsim <- function(filename, sim = 1, ...){
  ncols <- max(count.fields(filename, sep = "\t"))
  tmp <- read.delim(filename, col.names=1:ncols, sep="\t", fill=T, header=F)
  index <- which(tmp[,1] == "Time")
  tmp$simulation <- ""
  for(i in seq(index)){
    str <- index[i]+1
    end <- ifelse(i == length(index), nrow(tmp),index[i+1]-2)
    tmp$simulation[c(str:end)] <- i
  }
  df <- (subset(tmp, simulation == sim))
  x <- as.numeric(as.character(df$X1))
  y <- as.numeric(as.character(df$X2))
  plot(x, y, ...)
}


