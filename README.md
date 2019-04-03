# MCSim under R

This MCSim sandbox aim to help the beginner (especially Windows user) run GNU MCSim (current version 6.1.0) in R. 

Note: Please install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) in your system first.

## Instruction 

(1) Download all files from this repository.

(2) If you use RStudio, open "MCSim_under_R.Rproj".

(3) Open the R script in example folder and follow the guidance to do the simple test run.

- Remember to use `source("function.R")` to call the additional R functions.
- The models and inputs file have to put into the `input` folder.
- Use `getwd()` in R to make sure your working directory is in `MCSim_under_R`, such as `C:/Users/nanhung/MCSim_under_R`.

## Functions

Here are the R functions that can help you run MCSim in R environment more easily. All R functions are put in `functions.R`.

- `makemcsim(model-file)`:  Preprocessing and compiling the model file to the executable file as  [makemcsim](https://www.gnu.org/software/mcsim/mcsim.html#Using-makemcsim) in GNU MCSim.
- `mcsim(input-file, model-file)`: Using the compiled program to run simulation. See [Running Simulations](https://www.gnu.org/software/mcsim/mcsim.html#Running-Simulations).

## Examples

Some example R scripts are put into the `example` folder. 

### Help 

- Welcome to submit your problem in [issues](https://github.com/nanhung/MCSim_under_R/issues)

### Reference

- [GNU MCSim website](https://www.gnu.org/software/mcsim/)
- R. Woodrow Setzer. [Dynamic Modeling using MCSim and R](https://www.toxicology.org/groups/ss/BMSS/DynamicModelingwith%20MCsimandR.pdf)
- [MCSim under R (Windows)](https://nanhung.rbind.io/post/mcsim-under-r-windows/)