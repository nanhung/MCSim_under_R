compile_mcsim <- function(){
  if(Sys.which("gcc") == ""){
    stop("Please setting the PATH of compiler")
  }
    system(paste("gcc -o ./MCSim/mod.exe ./MCSim/mod/*.c ", sep = "")) 
    
  if(file.exists("mod/mod.exe")){
    message("The mod.exe had been created.")
  }
}

compile_mod <- function(mName){
  # Compile the "simple.model" to "simple.c" 
  system(paste("./MCSim/mod.exe input/", mName, " ", mName, ".c", sep = "")) 
  # Compile the "simple.model.c" to the executable program named "mcsim.simple.model.exe"
  system(paste("gcc -O3 -I.. -I./MCSim/sim -o mcsim_", mName, ".exe ", mName, ".c ./MCSim/sim/*.c -lm ", sep = ""))
  invisible(file.remove(paste0(mName, ".c")))
}

run_mcsim <- function(mName, inName){
  tx  <- readLines(paste0("input/", inName))
  MCMC_line <- grep("MCMC", x=tx)
  if (length(MCMC_line) != 0){
    #file_defore <- list.files()
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
    message(paste0("*Create '", checkfile, "' from the last iteration."))
    invisible(file.remove(paste0(outfile, ".kernel")))
  } else system(paste("./mcsim_", mName, ".exe ", "input/", inName, sep = ""))
}






