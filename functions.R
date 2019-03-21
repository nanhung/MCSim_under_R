compile_mcsim <- function(){
  if(Sys.which("gcc") == ""){
    stop("Please setting the PATH of compiler")
  }
  if(!file.exists("mod/mod.exe")){
    system(paste("gcc -o ./mod/mod.exe ./mod/*.c ", sep = "")) 
  }
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
  system(paste("./mcsim_", mName, ".exe ", "infile/", inName, sep = ""))
}
