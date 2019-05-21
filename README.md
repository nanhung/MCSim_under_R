# MCSim under R

This MCSim sandbox aim to help the beginner (especially Windows user) run GNU MCSim (current version 6.1.0) in R. 

## Prerequisites
- R (<https://cran.r-project.org/>)  
- RStudio (<https://www.rstudio.com/>)  
- Rtools (<https://cran.r-project.org/bin/windows/Rtools/>)  

## Instruction 

(1) Download all files from this repository.

(2) Open `"MCSim_under_R.Rproj"`.

(3) Open the R script in examples folder and follow the guidance to do the simple test run.

- **Note:** Use `getwd()` in R to make sure your working directory is in `MCSim_under_R`, such as `C:/Users/nanhung/MCSim_under_R`. The default install location of Rtools is `c:/Rtools`. 

### Workflow

The workflow of MCSim under R can separate into three levels as following diagram,

![](https://raw.githubusercontent.com/nanhung/MCSim_under_R/master/doc/fig/flowchart.png)

## Functions

Here are the R functions that can help you run MCSim in R environment more easily. All R functions are put in `functions.R` in MCSim folder.

- `makemcsim(model, deSolve)`:  Preprocessing and compiling the model-file to the executable file as  [makemcsim](https://www.gnu.org/software/mcsim/mcsim.html#Using-makemcsim) in GNU MCSim. The `model` assignment is a string giving the name of the model-file (e.g., `"pbpk.model.R"`). The `deSolve` assignment is a logical factor to use **deSolve** package as an ODE solver. 
- `mcsim(model, input)`: Using the compiled program with the input-file to run simulation. See [Running Simulations](https://www.gnu.org/software/mcsim/mcsim.html#Running-Simulations). The `input` assignment is a string giving the name of the input-file (e.g., `"linear.in.R"`).

### Help 

- Welcome to submit your problem in [issues](https://github.com/nanhung/MCSim_under_R/issues)

- For more detail, please see [this tutorial slide](https://nanhung.rbind.io/slide/190418_tutorial.html#1)

### Reference

- [GNU MCSim website](https://www.gnu.org/software/mcsim/)
- R. Woodrow Setzer. [Dynamic Modeling using MCSim and R](https://www.toxicology.org/groups/ss/BMSS/DynamicModelingwith%20MCsimandR.pdf)
- [MCSim under R (Windows)](https://nanhung.rbind.io/post/mcsim-under-r-windows/)
- [Using MCSim on Windows with RStudio](https://rpubs.com/Nanhung/MCSim_with_RStudio)
