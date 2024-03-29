# Download all files form this repository (the current version of MCSim is 6.0.1)
# Check the chek whether the compiler is in the PATH by using
# Sys.getenv("PATH") 

set_PATH <- function(PATH = "c:/rtools40/mingw32/bin;c:/Rtools/mingw_32/bin"){
  
  if (Sys.info()[['sysname']] == "Windows") {
    if(Sys.which("gcc") == ""){ # echo $PATH
      Sys.setenv(PATH = paste(PATH, Sys.getenv("PATH"), sep=";"))
    } # PATH=$PATH:/c/Rtools/mingw_32/bin; export PATH
  } # PATH=$PATH:/c/MinGW/msys/1.0/local/bin
  
  # The macos used clang as default, the following command is used to switch to GCC
  # Sys.setenv(PATH = paste("/usr/local/bin", Sys.getenv("PATH"), sep=";"))
  # port select --list gcc
  # sudo port select --set gcc mp-gcc8
  
  # Check the GCC compiler 
  Sys.which("gcc")
  system('gcc -v')
}

makemod <- function(){
  if(Sys.which("gcc") == ""){
    stop("Please setting the PATH of compiler")
  }
  mod_c_files <- paste(paste0("./MCSim/mod/", list.files("MCSim/mod", "*\\.c$")), collapse = ' ')
  system(paste0("gcc -o ./MCSim/mod.exe ", mod_c_files))
  
  if(file.exists("MCSim/mod.exe")){
    message("The mod.exe had been created.")
  }
}

makemcsim <- function(model, deSolve = F, dir = "modeling"){
  exe_file <- paste0("mcsim.", model, ".exe")
  #if(file.exists(exe_file)) stop(paste0("* '", exe_file, "' had been created."))
  
  if(file.exists("modeling") == F){
    dir.create("modeling")
  }
  
  if(file.exists(model)) { # Move model file from working directory to modeling folder
    invisible(file.copy(from = paste0(getwd(),"/", model), to = paste0(getwd(),"/modeling/", model)))
    invisible(file.remove(model))
    message(paste0("* The '", model, "' has been moved to modeling folder."))
  }
  
  if (deSolve == T){
    system(paste("./MCSim/mod.exe -R ", dir, "/", model, " ", model, ".c", sep = "")) 
    system (paste0("R CMD SHLIB ", model, ".c")) # create *.dll files
    dyn.load(paste(model, .Platform$dynlib.ext, sep="")) # load *.dll
    source(paste0(model,"_inits.R"))
  } else {
    system(paste("./MCSim/mod.exe ", dir, "/", model, " ", model, ".c", sep = "")) 
    
    sim_c_files <- paste(paste0("./MCSim/sim/", list.files("MCSim/sim", pattern = "*\\.c$")), collapse = ' ')
    system(paste("gcc -O3 -I.. -I./MCSim/sim -o mcsim.", model, ".exe ", model, ".c ", sim_c_files,
                 " -lm ", sep = ""))  
    invisible(file.remove(paste0(model, ".c")))
    if(file.exists(exe_file)) message(paste0("* Created executable program '", exe_file, "'.")) 
  }
}

mcsim <- function(model, input, dir = "modeling", parallel = F){
  
  if(file.exists("modeling") == F){
    dir.create("modeling")
  }
  
  exc = paste0("mcsim.", model, ".exe")
  if (file.exists(exc) == F) {
    makemcsim(model, dir = dir)
    if (file.exists(exc) == F) {
      stop("* '", exc, "' is not exist .")
    }
  }
  
  if(file.exists(input)) {
    invisible(file.copy(from = paste0(getwd(),"/", input), to = paste0(getwd(),"/modeling/", input)))
    invisible(file.remove(input))
  }
  
  tx  <- readLines(paste0(dir, "/", input))
  MCMC_line <- grep("MCMC \\(", x=tx)
  MonteCarlo_line <- grep("MonteCarlo \\(", x=tx)
  SetPoints_line <- grep("SetPoints \\(", x=tx)
  
  if (length(MCMC_line) != 0){
    #file_defore <- list.files()
    RandomSeed <- exp(runif(1, min = 0, max = log(2147483646.0)))
    tx2 <- gsub(pattern = "10101010", replace = paste(RandomSeed), x = tx)
    checkfile <- "MCMC.check.out"
    
    if(file.exists(checkfile)){
      file.remove(checkfile)
    }
    
    if (parallel == T){ 
      i <- sample(1111:9999, 1)
      name <- gsub("\\..*", "", input)
      mcmc_input <- paste0(name, "_", i, ".in")
      mcmc_output <- paste0(name, "_", i, ".out")
      tx3 <- gsub(pattern = "MCMC.default.out", replace = mcmc_output, x = tx2)
      writeLines(tx3, con = mcmc_input)
      system(paste("./mcsim.", model, ".exe ", mcmc_input, sep = ""))
      
    } else{ 
      tmp <- "tmp_mcmc.in.R"
      writeLines(tx, con=paste0(dir, "/", input))
      writeLines(tx2, con=paste0(dir, "/", tmp))
      system(paste("./mcsim.", model, ".exe ", dir, "/", tmp, sep = ""))
      outfile <- "MCMC.default.out"
      tx2 <- gsub(pattern = ",0,", replace = ",1,", x = tx)
      tx3 <- gsub(pattern = paste0("\"", outfile, "\",\"\""), 
                  replace = paste0("\"", checkfile, "\",\"", outfile, "\""), 
                  x = tx2)
      writeLines(tx3, con=paste0(dir, "/", tmp))
      
      system(paste("./mcsim.", model, ".exe ", dir, "/", tmp, sep = ""))
      file.remove(paste0(dir, "/", tmp))
    }
    
    if(file.exists(checkfile)){
      message(paste0("* Create '", checkfile, "' from the last iteration."))
    }
    
    if (parallel == T){ 
      df <- read.delim(mcmc_output)
    } else {
      df <- read.delim("MCMC.default.out")
    }
    
  } else if (length(MonteCarlo_line) != 0){
    RandomSeed <- runif(1, 0, 2147483646)
    tx2 <- gsub(pattern = "10101010", replace = paste(RandomSeed), x = tx)
    writeLines(tx2, con=paste0(dir, "/", input))
    message(paste("Execute:", " ./mcsim.", model, ".exe ", dir, "/", input, sep = ""))
    system(paste("./mcsim.", model, ".exe ", dir, "/", input, sep = ""))
    writeLines(tx, con=paste0(dir, "/", input))
    df <- read.delim("simmc.out")
  } else if (length(SetPoints_line) != 0){
    message(paste("Execute:", " ./mcsim.", model, ".exe ", dir, "/", input, sep = ""))
    system(paste("./mcsim.", model, ".exe ", dir, "/", input, sep = ""))
    df <- read.delim("simmc.out")
  } else {
    message(paste("Execute:", " ./mcsim.", model, ".exe ", dir, "/", input, sep = ""))
    system(paste("./mcsim.", model, ".exe ", dir, "/", input, sep = ""))
    df <- read.delim("sim.out", skip = 1)
  }
  return(df)
}

clear <- function(){
  files <- c(dir(pattern = c("*.out")),
             dir(pattern = c("sim.in")),
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
}

readsims <- function(x, exp = 1){
  ncols <- ncol(x)
  index <- which(x[,1] == "Time")
  str <- ifelse(exp == 1, 1, index[exp-1]+1)
  end <- ifelse(exp == length(index)+1, nrow(x), index[exp]-2)
  X <- x[c(str:end),]
  ncolX <- ncol(X) 
  X <- as.data.frame(matrix(as.numeric(as.matrix(X)), ncol = ncolX))
  if (exp > 1) names(X) <- as.matrix(x[index[exp-1],])[1:ncols] else names(X) <- names(x)
  X <- X[, colSums(is.na(X)) != nrow(X)]
  return(X)  
}

mcmc_array <- function(data, start_sampling = 0){
  n_chains <- length(data)
  sample_number <- dim(data[[1]])[1] - start_sampling
  dim <- c(sample_number, n_chains, dim(data[[1]])[2])
  n_iter <- dim(data[[1]])[1]
  n_param <- dim(data[[1]])[2]
  x <- array(sample_number:(n_iter * n_chains * n_param), dim = dim)
  for (i in 1:n_chains) {
    x[, i, ] <- as.matrix(data[[i]][(start_sampling + 1):n_iter, ])
  }
  dimnames(x)[[3]] <- names(data[[1]])
  x
}