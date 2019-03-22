compile_mcsim <- function(){
  if(Sys.which("gcc") == ""){
    stop("Please setting the PATH of compiler")
  }
    system(paste("gcc -o ./mod/mod.exe ./mod/*.c ", sep = "")) 
  if(file.exists("mod/mod.exe")){
    message("The mod.exe had been created.")
  }
}

compile_mod <- function(mName){
  # Compile the "simple.model" to "simple.c" 
  system(paste("./mod/mod.exe model/", mName, " ", mName, ".c", sep = "")) 
  # Compile the "simple.model.c" to the executable program named "mcsim.simple.model.exe"
  system(paste("gcc -O3 -I.. -I./sim -o mcsim_", mName, ".exe ", mName, ".c ./sim/*.c -lm ", sep = ""))
  invisible(file.remove(paste0(mName, ".c")))
}

run_mcsim <- function(mName, inName){
  tx  <- readLines(paste0("infile/", inName))
  MCMC_line <- grep("MCMC", x=tx)
  if (class(MCMC_line) == "integer"){
    file_defore <- list.files()
    system(paste("./mcsim_", mName, ".exe ", "infile/", inName, sep = ""))
    file_after <- list.files() # exist a bug if file had been created
    outfile <- setdiff(file_after,file_defore)[1]
    tx2 <- gsub(pattern = ",0,", replace = ",1,", x = tx)
    checkfile <- "check.out"
    tx3 <- gsub(pattern = paste0("\"", outfile, "\",\"\""), 
                replace = paste0("\"", checkfile, "\",\"", outfile, "\""), 
                x = tx2)
    writeLines(tx3, con=paste0("infile/", inName))
    system(paste("./mcsim_", mName, ".exe ", "infile/", inName, sep = ""))
    writeLines(tx, con="infile/digoxin.mcmc.in.R")
    message(paste0("*Create '", checkfile, "' from the last iteration."))
  } else system(paste("./mcsim_", mName, ".exe ", "infile/", inName, sep = ""))
}






