# Download all files form this repository (the current version of MCSim is 6.0.1)
# Check the chek whether the compiler is in the PATH by using
# Sys.getenv("PATH") 

set_PATH <- function(PATH = "c:\\MinGW\\bin"){
  if(Sys.which("gcc") == "")
    Sys.which("gcc")
  # You have two options to use GNU compiler:
  # If you have installed MinGW in your PC you can use
  Sys.setenv(PATH = paste(PATH, Sys.getenv("PATH"), sep=";")) # Recommend
  
  # Otherwise, if you have Rtools installed, you can assign the bin location manually,
  # Sys.setenv(PATH = paste("c:\\Rtools\\bin", Sys.getenv("PATH"), sep=";"))
  # Sys.setenv(PATH = paste("c:\\Rtools\\mingw_32/bin", Sys.getenv("PATH"), sep=";"))
  # Sys.setenv(PATH = paste("c:\\Rtools\\mingw_64/bin", Sys.getenv("PATH"), sep=";"))
  # Sys.setenv(BINPREF = "c:\\Rtools\\mingw_64/bin/") # Danger zone
  
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

compile_mod <- function(mName){
  exe_file <- paste0("mcsim_", mName, ".exe")
  if(file.exists(exe_file)) stop(paste0("* '", exe_file, "' had been created."))
    
  # Compile the "simple.model" to "simple.c"   
  system(paste("./MCSim/mod.exe input/", mName, " ", mName, ".c", sep = "")) 
  # Compile the "simple.model.c" to the executable program named "mcsim.simple.model.exe"
  system(paste("gcc -O3 -I.. -I./MCSim/sim -o mcsim_", mName, ".exe ", mName, ".c ./MCSim/sim/*.c -lm ", sep = ""))
  
  if(file.exists(exe_file)) message(paste0("* Created executable file '", exe_file, "'."))
  invisible(file.remove(paste0(mName, ".c")))
}

run_mcsim <- function(mName, inName){
  tx  <- readLines(paste0("input/", inName))
  MCMC_line <- grep("MCMC", x=tx)
  if (length(MCMC_line) != 0){
    #file_defore <- list.files()
    RandomSeed <- runif(1, 0, 2147483646)
    tx2 <- gsub(pattern = "10101010", replace = paste(RandomSeed), x = tx)
    writeLines(tx2, con=paste0("input/", inName))
    system(paste("./mcsim_", mName, ".exe ", "input/", inName, sep = ""))
    #file_after <- list.files() # exist a bug if file had been created
    #outfile <- setdiff(file_after,file_defore)[1]
    outfile <- "sim.out"
    tx2 <- gsub(pattern = ",0,", replace = ",1,", x = tx)
    checkfile <- "chk.out"
    tx3 <- gsub(pattern = paste0("\"", outfile, "\",\"\""), 
                replace = paste0("\"", checkfile, "\",\"", outfile, "\""), 
                x = tx2)
    writeLines(tx3, con=paste0("input/", inName))
    system(paste("./mcsim_", mName, ".exe ", "input/", inName, sep = ""))
    writeLines(tx, con=paste0("input/", inName))
    message(paste0("* Create '", checkfile, "' from the last iteration."))
    invisible(file.remove(paste0(outfile, ".kernel")))
  } else {
    system(paste("./mcsim_", mName, ".exe ", "input/", inName, sep = ""))
    df <- read.delim("sim.out", skip = 1)
  }
  return(df)
}

pkplot <- function(out, var = i, ...){
  plot(out[,1], out[,var], main = var,...)
}




