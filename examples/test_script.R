## Digoxin MCMC ####
# Define the input variable
mName <- "digoxin.model.R" # the model file put in the model folder
inName <- "digoxin.mcmc.in.R" # the input file put in the infile folder

# Create the executable file
makemcsim(mName)

# Run!!
set.seed(1111)
out <- mcsim(mName, inName)
out

par(mfrow= c(2,2))
names(out)
plot(out$k_12.1., type = "l")
plot(out$k_21.1., type = "l")
plot(out$k_10.1., type = "l")
plot(out$V_central.1., type = "l")

check_df <- read.delim("chk.out")
par(mfrow= c(1,1))
plot(check_df$Time, check_df$Data)
lines(check_df$Time, check_df$Prediction)
plot(check_df$Data, check_df$Prediction)
abline(0,1)


## simple model #####
# Define the name of model and input files
mName <- "simple.model.R" # the model file put in the model folder
inName <- "simple.in.R" # the input file put in the infile folder

# Create the executable file
makemcsim(mName, dir = "modeling/simple") # the files are located in modeling/simple 
# file.exists("mcsim_simple.model.R.exe") # check if you create the "mcsim_simple.model.exe" file successfully

# Run!!
out <- mcsim(mName, inName, dir = "modeling/simple")

# Print result
out

# Plot 
out <- out[2:13,] # omit 0
par(mfrow=c(2,2))
plot(out$Time, out$y0, log = "x", type = "b")
plot(out$Time, out$y1, log = "x", type = "b")
plot(out$Time, out$y2, log = "x", type = "b")


#
library(httk)
library(tidyverse)
library(sensitivity)
library(pksensi)

params <- parameterize_1comp(chem.name = "Bisphenol A")
MW <- params$MW
Vdist <- params$Vdist
kelim <- params$kelim
kgutabs <- params$kgutabs
Fgutabs <- params$Fgutabs * params$hepatic.bioavailability

#
t <- seq(0, 20, 0.1) # d
out <- solve_1comp(chem.name = "Bisphenol A", days = 50, doses.per.day = 1, daily.dose = 1, times = t)
data <- as.data.frame(out)
head(data)
dose <- 1 / 1000 * MW # mg -> uM
css <- dose * Fgutabs / kelim / Vdist 
httk_css <- calc_analytic_css(chem.name = "Bisphenol A", model = "1compartment")
plot(data$time, data$Ccompartment, type = "l")
abline(h = css)
abline(h = httk_css)

#
out_mcsim <- mcsim("pbtk1cpt.model.R", "pbtk1cpt.in.R", dir = "modeling/pbtk1cpt")
head(out_mcsim)
lines(out_mcsim$Time / 24, out_mcsim$Ccompartment, type = "l", col="red")

# 
Vdist_dist <- rnorm(n = 1000, Vdist, Vdist * 0.2)
kelim_dist <- rnorm(n = 1000, kelim, kelim * 0.2)
Fgutabs_dist <- runif(n = 1000, 0.8, 1)

css_dist <- dose * Fgutabs_dist / kelim_dist / Vdist_dist 
hist(css_dist)
summary(css_dist)

out <- mcsim("pbtk1cpt_css.model.R", "pbtk1cpt_css.mtc.in.R", dir = "modeling/pbtk1cpt")
head(out)
hist(out$css_1.1)
summary(out$css_1.1)

plot(density(css_dist))
lines(density(out$css_1.1), col = "red")

# Reverse
# serum BPA levels were 2.84 μg/L (arithmetic mean) and 0.18 μg/L (geometric mean) (He et al. 2009)
#  https://doi.org/10.1016/j.envres.2009.04.003
# The detectable rate was 17% for serum samples
# the detection limit of 0.31 μg/L for urine and 0.39 μg/L for serum
Css <- rlnorm(1000, log(0.18/1000), log(7)) # μg/L to mg/l
exp(sd(log(Css))) # Geometric Standard deviation

oral_equiv_dist <- Css * kelim_dist * Vdist_dist / Fgutabs_dist # mg/kg-d
summary(oral_equiv_dist)
hist(oral_equiv_dist)
plot(density(oral_equiv_dist))
boxplot(oral_equiv_dist, log = "y")


plot(Css, oral_equiv_dist, log = "xy", xlab = "Css (uM)", ylab = "oral equivalent dose (mg/kg-d)")
abline(v = 0.31/1000) # detection limit
x <- subset(Css, Css > 0.31/1000)
length(x)/1000 # Detection rate
mean(x) * 1000 # arithmetic mean

set.seed(1234)
out <- mcsim("pbtk1cpt_rtk.model.R", "pbtk1cpt_rtk.in.R", dir = "modeling/pbtk1cpt")
plot(out$Css, out$Dose_1.1, log = "xy", xlab = "Css (uM)", ylab = "oral equivalent dose (mg/kg-d)")
summary(out$Dose_1.1)
abline(v = 0.31/1000) # detection limit
DL <- quantile(x, prob = 0.13)
abline(v = DL, col = "red") # 13 %

x <- subset(out$Css, out$Css > 0.31/1000)
length(x)/1000 # Detection rate
mean(x) * 1000 # arithmetic mean
hist(out$Dose_1.1)
boxplot(out$Dose_1.1, oral_equiv_dist, log = "y", names = c("R","MCSim"))


#
for (i in 1:10)
{
  out <- mcsim("pbtk1cpt_rtk.model.R", "pbtk1cpt_rtk.in.R", dir = "modeling/pbtk1cpt")
  x <- subset(out$Css, out$Css > 0.31/1000)
  Detect.rate <- length(x)/1000 # Detection rate
  Observ.mean <- mean(x) * 1000 # arithmetic mean
  plot(out$Css, out$Dose_1.1, log = "xy", 
       xlab = "Css (uM)", ylab = "oral equivalent dose (mg/kg-d)",
       main = paste("Detection rate: ", round(Detect.rate, digit = 2), 
                    "Observation mean: ", round(Observ.mean, digit = 2)))
  abline(v = 0.31/1000) # detection limit
  DL <- quantile(x, prob = 0.13)
  abline(v = DL, col = "red") # 13 %
  date_time<-Sys.time()
  while((as.numeric(Sys.time()) - as.numeric(date_time))<1){} 
}

#
easy_1 <- function (X) 
{
  X[, 1] + X[, 2] + X[, 3]
}
easy_2 <- function (X) 
{
  X[, 1] * X[, 2] / X[, 3]
}

set.seed(123)
x <- morris(model = easy_1, factors = 3, r = 18, 
            design = list(type = "oat", levels = 6, grid.jump = 3),  # grid.jump = levels/2
            binf = c(1, 1, 1), bsup = c(2, 4, 7))

par(mfrow = c(1,3))
for(i in 1:3){
  hist(x$X[,i], main = colnames(x$X)[i])
}

par(mfrow = c(1,3))
for(i in 1:3){
  cor <- cor(x$X[,i], x$y)
  plot(x$X[,i], x$y, 
       xlab = colnames(x$X)[i],
       main = paste("r = ", round(cor, digits = 2)))
}

par(mfrow = c(1,1))
plot(x)

x

# Analytical solution (first order)
2 - 1 # X1
4 - 1 # X2
7 - 1 # X3

q <- "qunif"
q.arg <- list(list(min = 1,  max = 2),
              list(min = 1,  max = 4),
              list(min = 1,  max = 7))
x <- fast99(model = easy_1, factors = 3, n = 192,
            q = q, q.arg = q.arg)

par(mfrow = c(1,3))
for(i in 1:3){
  hist(x$X[,i], main = colnames(x$X)[i])
}

par(mfrow = c(1,3))
for(i in 1:3){
  cor <- cor(x$X[,i], x$y)
  plot(x$X[,i], x$y, 
       xlab = colnames(x$X)[i],
       main = paste("r = ", round(cor, digits = 2)))
}

par(mfrow = c(1,1))
plot(x)
x

# Analytical solution
Total <- 1^2 + 3^2 + 5^2
S1 <- 1^2 / Total
S1
S2 <- 3^2 / Total
S2
S3 <- 5^2 / Total
S3

# easy 2
set.seed(123)
x <- morris(model = easy_2, factors = 3, r = 18, 
            design = list(type = "oat", levels = 6, grid.jump = 3),
            binf = c(1,1,1), bsup = c(2,6,7))
plot(x, xlim = c(0,5), ylim = c(0,5))
x

x <- fast99(model = easy_2, factors = 3, n = 512,
            q = q, q.arg = q.arg)
plot(x)
x

#
Css_fun <- function (X) 
{
  Dose <- 0.228291 # Ingestion dose (uM)
  Dose * X[, 1] / X[, 2] / X[, 3] # Fgutabs_dist / kelim_dist / Vdist_dist 
}

binf <- c(min(Fgutabs_dist), min(kelim_dist), min(Vdist_dist))
bsup <- c(max(Fgutabs_dist), max(kelim_dist), max(Vdist_dist))
binf
bsup

set.seed(1234)
x <- morris(model = Css_fun, factors = c("Fgutabs", "kelim", "Vdist"), r = 32, 
            design = list(type = "oat", levels = 6, grid.jump = 3), # grid.jump = levels/2
            binf = binf, bsup = bsup)
head(x$X)

par(mfrow = c(1,3))
for(i in 1:3){
  cor <- cor(x$X[,i], x$y)
  plot(x$X[,i], x$y, 
       xlab = colnames(x$X)[i],
       main = paste("r = ", round(cor, digits = 2)))
}

par(mfrow = c(1,1))
plot(x, xlim = c(0, 5), ylim = c(0, 5))
abline(0,1) # non-linear and/or non-monotonic
abline(0,0.5, lty = 2) # monotonic
abline(0,0.1, lty = 3) # almost linear
legend("topleft", legend = c("non-linear and/or non-monotonic",
                             "monotonic", "linear"), lty = c(1:3))

for (i in 1:10)
{
  x <- morris(model = Css_fun, factors = c("Fgutabs", "kelim", "Vdist"), r = 32, # test r = 32
              design = list(type = "oat", levels = 6, grid.jump = 3), 
              binf = binf, bsup = bsup)
  plot(x, xlim = c(0, 6), ylim = c(0, 6))
  abline(0,1) # non-linear and/or non-monotonic
  abline(0,0.5, lty = 2) # monotonic
  abline(0,0.1, lty = 3) # almost linear
  legend("topleft", legend = c("non-linear and/or non-monotonic",
                               "monotonic", "linear"), lty = c(1:3))
  date_time<-Sys.time()
  while((as.numeric(Sys.time()) - as.numeric(date_time))<1){} 
}


q <- "qunif"
q.arg <- list(list(min = min(Fgutabs_dist),  max = max(Fgutabs_dist)),
              list(min = min(kelim_dist),  max = max(kelim_dist)),
              list(min = min(Vdist_dist),  max = max(Vdist_dist)))
x <- fast99(model = Css_fun, factors = c("Fgutabs", "kelim", "Vdist"), n = 512,
            q = q, q.arg = q.arg)
x

par(mfrow = c(1,3))
for(i in 1:3){
  cor <- cor(x$X[,i], x$y)
  plot(x$X[,i], x$y, 
       xlab = colnames(x$X)[i],
       main = paste("r = ", round(cor, digits = 2)))
}

par(mfrow = c(1,1))
plot(x)

# Monte Carlo
model <- "pbtk1cpt.model.R"
out <- mcsim(model, "pbtk1cpt.mtc.in.R", dir = "modeling/pbtk1cpt")
head(out)
index <- which(names(out) == "Ccompartment_1.1" | names(out) == "Ccompartment_1.24")

X <- apply(out[,index[1]:index[2]], 2, quantile,  c(0.5, 0.025, 0.975))
dat <- t(X)
colnames(dat) <- c("median", "LCL", "UCL")
df <- as.data.frame(dat)
df$time <- seq(0, 23, 1)
ggplot(df, aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), fill = "grey70", alpha = 0.5) + 
  geom_line()

# 
n <- 1000
parameters <- c("Vdist","kelim", "kgutabs", "Fgutabs")  
outputs <- c("Ccompartment")
times <- seq(0.1, 480.1, 1)
dist<-c("Normal_cv", "Normal_cv", "Normal_cv", "Uniform") # MCSim definition
q.arg<-list(list(6.137241, 0.2),
            list(0.02283233, 0.2),
            list(2.18, 0.2),
            list(0.8, 1.0))
condition <- c("MW = 228.291", "Period = 24", "IngDose = 1.0")

set.seed(2222)
y<-solve_mcsim(mName = model, params = parameters, vars = outputs, monte_carlo = n,
               dist = dist, q.arg = q.arg, time = times, condition = condition)
pksim(y, log = T)
summary(y)

# Sensitivity
q <- "qunif"
q.arg <- list(list(min =  Vdist * 0.8, max = Vdist * 1.2),
              list(min = kelim * 0.8, max = kelim * 1.2),
              list(min = kgutabs * 0.8, max = kgutabs * 1.2),
              list(min = 0.8, max = 1)) 
set.seed(1234)
params <- parameters
x <- rfast99(params, n = 400, q = q, q.arg = q.arg, replicate = 10)
times <- seq(400, 480, 1)
x$a

y <- solve_mcsim(x, mName = model,  params = parameters, time = times,  vars = outputs, condition = condition)
#tell2(x,y$y)

plot(y)

check(y)

pksim(y)


## Monte Carlo ####
model <- "perc.model.R"
input <- "perc.mtc.in.R" 
out <- mcsim(model, input, dir = "modeling/perc")
head(out)

last.param <- which(names(out) == "Flow_pul")

par(mfrow = c(4, 4), mar = c(2,2,4,1))
for (i in 2:last.param){
  hist(out[,i], main = names(out[i]), xlab = "")
}

vars <- names(out)
sim1.1 <- which(vars == "C_exh_ug_1.1" | vars == "C_exh_ug_1.10")
sim1.2 <- which(vars == "C_ven_1.1" | vars == "C_ven_1.8") 
sim2.1 <- which(vars == "C_exh_ug_2.1" | vars == "C_exh_ug_2.10")
sim2.2 <- which(vars == "C_ven_2.1" | vars == "C_ven_2.8")

ggPK <- function(data, index, time){
  X <- apply(data[,index[1]:index[2]], 2, quantile,  c(0.5, 0.025, 0.975))
  dat <- t(X)
  colnames(dat) <- c("median", "LCL", "UCL")
  df <- as.data.frame(dat)
  df$time <- time
  ggplot(df, aes(x = time, y = median)) +
    geom_ribbon(aes(ymin = LCL, ymax = UCL), fill = "grey70", alpha = 0.5) + 
    geom_line() +
    scale_y_log10()
}

t_exh_ug <- c(239.9, 245, 270, 360, 1320, 2760, 4260, 5700, 8580, 10020)
t_Ven <- c(239.9, 360, 1320, 2760, 4260, 5700, 8580, 10020)
ggPK(data = out, index = sim1.1, time = t_exh_ug)
ggPK(data = out, index = sim1.2, time = t_Ven)
ggPK(data = out, index = sim2.1, time = t_exh_ug)
ggPK(data = out, index = sim2.2, time = t_Ven)

# pksensi (PERC) - Monte Carlo
n <- 1000
parameters <- c("LeanBodyWt","Pct_M_fat", "Pct_LM_liv", "Pct_LM_wp",
                "Pct_Flow_fat", "Pct_Flow_liv", "Pct_Flow_pp",
                "PC_fat", "PC_liv", "PC_wp", "PC_pp", "PC_art",
                "Vent_Perf", "sc_Vmax", "Km", "Flow_pul")  
outputs <- c("C_exh_ug", "C_ven")
times <- c(239.9, 245, 270, 360, 1320, 2760, 4260, 5700, 8580, 10020)
dist <- rep("Uniform", 16) # MCSim definition
q.arg<-list(list(50, 70),
            list(0.2, 0.3),
            list(0.03, 0.04),
            list(0.25, 0.3),
            list(0.06, 0.08),
            list(0.2, 0.3),
            list(0.2, 0.25),
            list(110, 150),
            list(5.0, 8.0),
            list(5.0, 8.5),
            list(1.6, 1.8),
            list(12, 15),
            list(0.8, 1.3),
            list(0.04, 0.06),
            list(7, 13),
            list(7.4, 7.6))
condition <- c("InhMag = 72", "Period = 1e10", "Exposure = 240")

set.seed(2222)
y<-solve_mcsim(mName = model, params = parameters, vars = outputs, monte_carlo = n,
               dist = dist, q.arg = q.arg, time = times, condition = condition)

par(mfrow = c(1, 2), mar = c(2,3,1,1))
pksim(y, vars = "C_exh_ug", log = T)
pksim(y, vars = "C_ven", log = T)

# pksensi (PERC) - Sensitivity
q <- rep("qunif", 16)
x <- rfast99(params = parameters, n = 1024, q = q, q.arg = q.arg, replicate = 10)
dim(x$a)

par(mfrow = c(4, 4), mar = c(2,2,4,1))
for(i in 1:16){
  hist(x$a[,,i])  
}

y <- solve_mcsim(x, mName = model, params = parameters, vars = outputs,
                 time = times, condition = condition)

dim(y)

check(y)

plot(y)

heat_check(y) 
heat_check(y, index = "CI") 
