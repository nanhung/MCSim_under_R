# Download all files form this repository (the current version of MCSim is 6.0.1)
# Check the chek whether the compiler is in the PATH by using
# Sys.getenv("PATH") 

if(!require(pkgbuild)) install.packages("pkgbuild")

set_PATH <- function(){
  PATH <- "c:/Rtools/mingw_32/bin"
  
  if (Sys.info()[['sysname']] == "Windows") {
    if(!pkgbuild::find_rtools()){
      warning("Please make sure you had installed Rtools.")
    }
    if(Sys.which("gcc") == ""){
      Sys.setenv(PATH = paste(PATH, Sys.getenv("PATH"), sep=";"))
    }
  }
  
  # You have two options to use GNU compiler:
  # If you have installed MinGW in your PC you can use
  # Sys.setenv(PATH = paste("c:\\MinGW\\bin", Sys.getenv("PATH"), sep=";"))  
  # Otherwise, if you have Rtools installed, you can assign the bin location manually,
  # Sys.setenv(PATH = paste("c:\\Rtools\\mingw_32/bin", Sys.getenv("PATH"), sep=";"))
  
  # The macos used clang as default, the following command is used to switch to GCC
  # Sys.setenv(PATH = paste("/usr/local/bin", Sys.getenv("PATH"), sep=";"))
  # port select --list gcc
  # sudo port select --set gcc mp-gcc8
  
  # Check the GNU compiler 
  Sys.which("gcc")
  system('g++ -v')
}

compile_mcsim <- function(){
  if(Sys.which("gcc") == ""){
    stop("Please setting the PATH of compiler")
  }
  system(paste("gcc -o ./MCSim/mod.exe ./MCSim/mod/*.c ", sep = "")) 
  
  if(file.exists("mod/mod.exe")){
    message("The mod.exe had been created.")
  }
}
set_PATH()
compile_mcsim()

clear <- function(){
  files <- c(dir(pattern = c("*.out")), dir(pattern = c("*.exe")), dir(pattern = c("*.exe")))
  invisible(file.remove(files))
}

compile_mod <- function(mName){
  exe_file <- paste0("mcsim.", mName, ".exe")
  if(file.exists(exe_file)) stop(paste0("* '", exe_file, "' had been created."))
  if(file.exists(mName)) {
    invisible(file.copy(from = paste0(getwd(),"/", mName), to = paste0(getwd(),"/input/", mName)))
    invisible(file.remove(mName))
  }
    
  # Compile the "simple.model" to "simple.c"   
  system(paste("./MCSim/mod.exe input/", mName, " ", mName, ".c", sep = "")) 
  # Compile the "simple.model.c" to the executable program named "mcsim.simple.model.exe"
  system(paste("gcc -O3 -I.. -I./MCSim/sim -o mcsim.", mName, ".exe ", mName, ".c ./MCSim/sim/*.c -lm ", sep = ""))
  
  if(file.exists(exe_file)) message(paste0("* Created executable file '", exe_file, "'."))
  invisible(file.remove(paste0(mName, ".c")))
}

run_mcsim <- function(mName, inName){
  if(file.exists(inName)) {
    invisible(file.copy(from = paste0(getwd(),"/", inName), to = paste0(getwd(),"/input/", inName)))
    invisible(file.remove(inName))
  }
  tx  <- readLines(paste0("input/", inName))
  MCMC_line <- grep("MCMC", x=tx)
  if (length(MCMC_line) != 0){
    #file_defore <- list.files()
    RandomSeed <- runif(1, 0, 2147483646)
    tx2 <- gsub(pattern = "10101010", replace = paste(RandomSeed), x = tx)
    writeLines(tx2, con=paste0("input/", inName))
    system(paste("./mcsim.", mName, ".exe ", "input/", inName, sep = ""))
    #file_after <- list.files() # exist a bug if file had been created
    #outfile <- setdiff(file_after,file_defore)[1]
    outfile <- "sim.out"
    tx2 <- gsub(pattern = ",0,", replace = ",1,", x = tx)
    checkfile <- "chk.out"
    tx3 <- gsub(pattern = paste0("\"", outfile, "\",\"\""), 
                replace = paste0("\"", checkfile, "\",\"", outfile, "\""), 
                x = tx2)
    writeLines(tx3, con=paste0("input/", inName))
    system(paste("./mcsim.", mName, ".exe ", "input/", inName, sep = ""))
    writeLines(tx, con=paste0("input/", inName))
    message(paste0("* Create '", checkfile, "' from the last iteration."))
    invisible(file.remove(paste0(outfile, ".kernel")))
    df <- read.delim("sim.out")
  } else {
    system(paste("./mcsim.", mName, ".exe ", "input/", inName, sep = ""))
    df <- read.delim("sim.out", skip = 1)
  }
  return(df)
}

plot_sim <- function(filename, sim = 1, ...){
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

pkplot <- function(out, var = i, ...){
  plot(out[,1], out[,var], main = var,...)
}



