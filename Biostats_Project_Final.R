#---Loading Yi et al. Table II data — same approach as Jan 25 lecture where we learnt how to load data in---#
#data <- read.csv(/Users/rj/Desktop/Prof_Hernando/yi_2023_GI50.csv) #or use file.choose()
data <- read.csv(file.choose()) #This will open a file dialog where you can navigate to the location of your data file and select it. This is a convenient way to load data without hardcoding the file path.
head(data) #Simply displays first 6 rows of my data frame [Sanity check] & should always be run after loading data in. 
str(data) #10 observations means 10 rows, 9 variables means 9 columns: Patient | age | location ...

# Data is now Loaded and ready to go!

# ============================================================
# Statistical In Silico Critique of Yi et al. (2023)
# A Monte Carlo Power Analysis of an Organoid-Based Drug Study
# Author: Reda Jalal
# ============================================================

# === Quantitative Concern 1: Sample Size and Statistical Power ===
# Linked to: January 28 (sampling distributions) and March 1 (sample size)

# Split into the two groups older than 69, or younger <= 69
# From the data frame called data, give me the column called cisplatin_uM & from that I'll take either [data$age <= 69] or
# [data$age > 69] and give me the cisplatin_uM values for those two groups separately
# data[data$age <= 69, c("patient", "age")] & same for > 69 is what I need, remember c means combine, otherwise I can only 
#obtain one column at a time

group_under69 <- data$cisplatin_uM[data$age <= 69]
group_over69 <- data$cisplatin_uM[data$age > 69]

# ---- PART A: Sampling Variability at n=5 ----
# Parameters taken directly from my constructed table of Yi et al. data (cvs file)
# Cisplatin sensitivity split by age (median split at 69 years)
# Means and SDs obtained from the values reported in the paper
mu_under69 <- 24.81
mean(group_under69) #For confirmation
sd_under69 <- 11.65
sd(group_under69) #For confirmation
mu_over69 <- 49.10
mean(group_over69) #For confirmation
sd_over69 <- 5.07
sd(group_over69) #For confirmation

n <- 5 #This is the sample size per group in Yi et al. (2023)
B <- 1000 #Number of Monte Carlo simulations
means_under69 <- numeric(B) #Empty vector to store means for under 69 group
means_over69 <- numeric(B) 

set.seed(123) #keep before any loop from here on out

for (b in 1:B){
  samp_under69 <- rnorm(n, mu_under69, sd_under69)
  samp_over69  <- rnorm(n, mu_over69, sd_over69)
  means_under69[b] <- mean(samp_under69)
  means_over69[b]  <- mean(samp_over69)
}

# Plotting sampling distributions to show how unstable sample means are at n=5
par(mfrow = c(1,2)) # plot windoww layout: sets up 2 plots side by side
hist(means_under69,
     main = "Age ≤69: Cisplatin GI50",
     xlab = "Sample Mean (µM)",
     col  = "lightblue",
     border = "white")

hist(means_over69,
     main = "Age >69: Cisplatin GI50",
     xlab = "Sample Mean (µM)",
     col  = "lightcoral",
     border = "white")

# My Visual analysis:
#Firstly the plugged in values are: 
# Under 69: SD = 11.65 — large spread in the population
# Over 69: SD = 5.07 — small spread in the population
# I can see a bell shaped curve: Normal distribution demonstration of dentral limit theorem. 
# Under 69 group (blue) — wider, more spread out, ranging roughly from 10 to 50 µM
# Over 69 group (red) — narrower, more concentrated, ranging roughly from 44 to 57 µM
# Conclusion: Both histograms are wide. The under 69 group spans roughly 40 µM of range. That means if Yi et al. had recruited a 
# slightly different set of 5 patients, their reported mean of 24.81 could easily have come out anywhere between 10 and 50. That 
# single number they reported is just one possible outcome from a highly unstable experiment.
# This is  visual proof that n=5 is too small to reliably estimate the true group mean.

# ---- PART B: Power curve ----
# Step 1: Calculate effect size from Yi et al. reported values
# Pooled SD: weighted average of both group SDs
# I've used the same formula from Feb-04 lecture t-test code (See lesson)

w1 <- (n - 1) / ((n - 1) + (n - 1))
w2 <- (n - 1) / ((n - 1) + (n - 1))
pooled_sd <- sqrt(w1 * sd_under69^2 + w2 * sd_over69^2)
paste("the pooled SD is:",pooled_sd)

# Under 69 group had SD = 11.65 — patients in this group vary quite a lot from each other
# Over 69 group had SD = 5.07 — patients in this group are more similar to each other
# Pooled SD = 8.98 (this is the value obtained from 'pooled_sd', this = the average spread across both groups combined
# Effect size = difference between means / pooled SD

effect_size <- (mu_over69 - mu_under69) / pooled_sd
paste("the effect side is:",effect_size)

# I got 2.7, this is a huge value; a value around 0.5 would be considered 'medium'
# Due to this large effect size, because it's so huge, even a small sample size of n=5 per group can yield high power. This is because the difference between the groups is so large relative to the variability within the groups that it becomes easier to detect a statistically significant difference, even with a small number of samples.'

# I'll now test a range of sample sizes — say n = 3, 5, 8, 10, 15, 20, 25, 30 — and for each one run 1000 simulations, calculate 
# the t-statistic each time, and count how often it exceeds the threshold. That proportion is the power at that sample size.
# Store the power for each n, then plot them all together = Power curve . In order to find what sample size gives 80% power 

alpha <- 0.05
sample_sizes <- c(3, 5, 8, 10, 15, 20, 25, 30)
power_values <- numeric(length(sample_sizes)) #has 8 slots for the 8 sample sizes I'm trying
B <- 1000

# Loop through each sample size
# For each n, simulate 1000 experiments and calculate power
# There will be 2 for loops nested in one another: =
# The outer loop iterates through sample sizes, and the inner loop simulates B experiments for each sample size, 
# calculating the t-statistic and determining whether it exceeds the critical threshold for significance. The proportion of 
#simulations where the null hypothesis is rejected gives us the estimated power for that sample size.

for(i in 1:length(sample_sizes)){
  
  n_i <- sample_sizes[i] # current sample size
  df <- n_i + n_i - 2 # degrees of freedom for t-test
  thresh <- qt(1 - alpha/2, df) # critical value threshold
  
  decisions <- numeric(B)
  
  for(b in 1:B){
    s1 <- rnorm(n_i, mu_under69, sd_under69)
    s2 <- rnorm(n_i, mu_over69, sd_over69)
    
    m1 <- mean(s1); m2 <- mean(s2)
    v1 <- var(s1);  v2 <- var(s2)
    
    w1 <- (n_i-1)/((n_i-1)+(n_i-1))
    w2 <- (n_i-1)/((n_i-1)+(n_i-1))
    S2 <- w1*v1 + w2*v2
    
    t.stat <- (m1 - m2) / sqrt(S2*(1/n_i + 1/n_i)) #Calculate the t statstic, if exceeds threshold store 1, otherwise store 0
    decisions[b] <- 1*(abs(t.stat) > thresh)
  }

  power_values[i] <- sum(decisions)/B #Count how many /1000 sims detected the effect size and divide by 1000 to get me the empirical power for that sample size, stored in power_values vector
}

# now I'll print power at each sample size, I'll also use 'paste' for clearer for formatting answers

for(i in 1:length(sample_sizes)){
  cat(paste("n =", sample_sizes[i], "| Power =", round(power_values[i], 3), "\n"))
}

# Now before I even plot, the values that I got, even if they used only n=5, they had 0.95 power, which is very high. 
# This is because of the huge effect size of 2.7. Even at n=3, they had 0.87 power, which is still very high. This means
# that even with a very small sample size, they had a very high chance of detecting the true effect if it exists, due to the
# large difference between the groups relative to the variability within the groups.
# Anything 10 & above is 1 meaning they would have detected the effect in all 1000 simulations, which is perfect power.
# If this comparison was well powered, why run 16 ANOVAs? The other comparisons — oxaliplatin by location with groups of n=7 
# and n=3 — will tell a completely different story." teh fact that they ran 19 ANOVAS makes me question did they continue 
# with teh experiment knowing they'd be fine with their 5 sampels per group...OR was it pure chance...

# Plotting the power curve
# Red dashed line marks 80% power threshold
# Blue dashed line marks n=5 used by Yi et al.

par(mfrow = c(1,1)) # reset to single plot

plot(sample_sizes, power_values,
     type = "b",
     pch  = 19,
     xlab = "Sample size per group (n)",
     ylab = "Estimated Power",
     main = "Power Curve: Cisplatin by Age\nYi et al. (2023)",
     ylim = c(0.5, 1.0),
     yaxt = "n",
     xlim = c(0, 32))

axis(2, at = seq(0.5, 1.0, by = 0.1),
     labels = seq(0.5, 1.0, by = 0.1))

abline(h = 0.80, col = "red",  lty = 2, lwd = 2)
abline(v = 5,    col = "blue", lty = 2, lwd = 2)

text(7,  0.76, "n=5 (Yi et al.)", col = "blue", cex = 0.8)
text(20, 0.83, "80% threshold",   col = "red",  cex = 0.8)

# ---Post-meeting---: Plotting multiple power curves on the same graph for different effect sizes: 
# 2.7, 1.0, and 0.8 — to show that the result only looks well-powered because the
# effect happened to be enormous. My point 1 Analysis will be completed with this. 

# We already have power_values for ES=2.7
# Now simulate for ES=1.0 and ES=0.8 to show contrast
# Linked to: March 1 — power curves and sample size

# New means to achieve target effect sizes
# Formula: mu_over69 = mu_under69 + (effect_size x pooled_sd)

mu_over69_ES10 <- mu_under69 + 1.0 * pooled_sd
mu_over69_ES08 <- mu_under69 + 0.8 * pooled_sd

# Empty containers for new power values

power_ES10 <- numeric(length(sample_sizes))
power_ES08 <- numeric(length(sample_sizes))

set.seed(123)

# Loop for effect size 1.0

for(i in 1:length(sample_sizes)){
  n_i    <- sample_sizes[i]
  df     <- n_i + n_i - 2
  thresh <- qt(1 - alpha/2, df)
  decisions <- numeric(B)
  
  for(b in 1:B){
    s1 <- rnorm(n_i, mu_under69, sd_under69)
    s2 <- rnorm(n_i, mu_over69_ES10, sd_over69)
    m1 <- mean(s1); m2 <- mean(s2)
    v1 <- var(s1);  v2 <- var(s2)
    w1 <- (n_i-1)/((n_i-1)+(n_i-1))
    w2 <- (n_i-1)/((n_i-1)+(n_i-1))
    S2 <- w1*v1 + w2*v2
    t.stat <- (m1-m2)/sqrt(S2*(1/n_i+1/n_i))
    decisions[b] <- 1*(abs(t.stat) > thresh)
  }
  power_ES10[i] <- sum(decisions)/B
}

# Loop for effect size 0.8

for(i in 1:length(sample_sizes)){
  n_i    <- sample_sizes[i]
  df     <- n_i + n_i - 2
  thresh <- qt(1 - alpha/2, df)
  decisions <- numeric(B)
  
  for(b in 1:B){
    s1 <- rnorm(n_i, mu_under69, sd_under69)
    s2 <- rnorm(n_i, mu_over69_ES08, sd_over69)
    m1 <- mean(s1); m2 <- mean(s2)
    v1 <- var(s1);  v2 <- var(s2)
    w1 <- (n_i-1)/((n_i-1)+(n_i-1))
    w2 <- (n_i-1)/((n_i-1)+(n_i-1))
    S2 <- w1*v1 + w2*v2
    t.stat <- (m1-m2)/sqrt(S2*(1/n_i+1/n_i))
    decisions[b] <- 1*(abs(t.stat) > thresh)
  }
  power_ES08[i] <- sum(decisions)/B
}

# Plot all three power curves on same graph

par(mfrow = c(1,1))

plot(sample_sizes, power_values,
     type = "b", pch = 19, col = "black",
     xlab = "Sample size per group (n)",
     ylab = "Estimated Power",
     main = "Power Curves at Different Effect Sizes\nYi et al. (2023) Cisplatin by Age",
     ylim = c(0, 1),
     xlim = c(0, 32))

lines(sample_sizes, power_ES10,
      type = "b", pch = 17, col = "blue")

lines(sample_sizes, power_ES08,
      type = "b", pch = 15, col = "red")

abline(h = 0.80, lty = 2, col = "grey50")
abline(v = 5,    lty = 2, col = "grey50")

legend("bottomright",
       legend = c("ES = 2.70 (reported)",
                  "ES = 1.00",
                  "ES = 0.80"),
       col    = c("black", "blue", "red"),
       pch    = c(19, 17, 15),
       lty    = 1)

#Analysis: "Yi et al.'s study only appears adequately powered because their cisplatin by age comparison happened to have an enormous
#effect size of 2.70. Had the true effect been a more modest 1.0 — which is still considered large — n=5 would have given only 30% power.
#For their other comparisons where effect sizes are smaller, the study was almost certainly severely underpowered."



# ============================================================
# POINT 2: Multiple Comparisons Without Correction
# Linked to: February 4 — Type I error and false positive rates
# ============================================================

# Yi et al. ran 16 ANOVAs on the same dataset at alpha = 0.05
# No correction for multiple comparisons was applied
# Bonferroni correction would require p < 0.05/16 = 0.003125 per test as mentioned by professor

# I'll simulate 16 hypothesis tests under the null hypothesis (no true effect — both groups identical)
# and count how often at least one returns p < 0.05 by chance

n_tests <- 16
alpha <- 0.05 
B <- 10000 # Number of simulations, each of these 10k times will represent one hypothetical version of this study, where a researcher runs 16 tests ona. dataset where no true effect exists
# Make 2 x for loops, one that runs the simulation 10k times  and the inner loop to generate two identical groups drawn from same distribution , run a t test for each of the 16, and store the p value
# I then need to score 1 if any of these p values are below 0.05, and 0 otherwise
# I can then count out of these 10k studies, how many contained at least one false positive

set.seed(123)
at_least_one_FP <- numeric(B)

for(b in 1:B){
  p_values <- numeric(n_tests)
  
  for(t in 1:n_tests){
    group1 <- rnorm(5, 0, 1)
    group2 <- rnorm(5, 0, 1)
    p_values[t] <- t.test(group1, group2)$p.value
  }
  
  at_least_one_FP[b] <- 1*(min(p_values) < alpha)
}

false_positive_rate <- sum(at_least_one_FP)/B

cat(paste("Empirical false positive rate across 16 tests:", 
          round(false_positive_rate, 3), "\n"))
cat(paste("Analytical expectation: 1-(0.95^16) =", 
          round(1-0.95^16, 3), "\n"))
analytical_expectation <- round(1-0.95^16, 3)
cat(paste("The difference between my calculated empirical false positive rate and the analytical expectation = ",analytical_expectation - false_positive_rate, ",is due to random sampling variability in the simulations. With a large number of simulations (B=10,000), the empirical false positive rate to converge towards the analytical expectation of approximately 0.56."))

#---Let's try now using with the Benferroni correction, caclulate it:---

bonferroni_thresh <- alpha / n_tests 
set.seed(123)
at_least_one_FP_bonf <- numeric(B)

for(b in 1:B){
  p_values <- numeric(n_tests)
  
  for(t in 1:n_tests){
    group1 <- rnorm(5, 0, 1)
    group2 <- rnorm(5, 0, 1)
    p_values[t] <- t.test(group1, group2)$p.value
  }
  
  at_least_one_FP_bonf[b] <- 1*(min(p_values) < bonferroni_thresh)
}

false_positive_rate_bonf <- sum(at_least_one_FP_bonf)/B

cat(paste("Empirical false positive rate across 16 tests:", 
          round(false_positive_rate, 3), "\n"))
cat(paste("With Bonferroni correction (alpha = 0.003):", 
          round(false_positive_rate_bonf, 3), "\n"))

s# Yi et al. ran 16 uncorrected tests at alpha = 0.05, giving a 53% chance of at least one false positive. 
# Had they applied Bonferroni correction — requiring p < 0.003 per test — that risk would have dropped to 3.7%, consistent
# with the intended 5% error rate, what this would mean is that their Irinotecan results that they considered "signnificant" is likely not significant
# using this new corrected threshold, for Cisplatin, which was the other 'significant' finding, we know p < 0.01, what we don't know the actual value...hence we can't comment!
# Hence the statistical basis for the clincal recommednation that patients characteristics can guide chemotherapy decisions fades...

# ---- POINT 2: Bar chart — False positive rates ----
# Visual comparison of uncorrected vs Bonferroni corrected
# false positive rates across 16 simultaneous tests

par(mfrow = c(1,1))

bp <- barplot(fp_rates,
              names.arg = labels,
              col       = colours,
              border    = "white",
              ylim      = c(0, 0.70),
              ylab      = "False Positive Rate",
              main      = "Effect of Bonferroni Correction\non Family-wise False Positive Rate",
              cex.names = 0.85,
              width     = 0.5,
              space     = 1)

abline(h = 0.05, lty = 2, col = "red", lwd = 2)

text(bp[1], false_positive_rate + 0.10,
     paste0(round(false_positive_rate * 100, 1), "%"),
     col = "darkred", cex = 1.2, font = 2)

text(bp[2], false_positive_rate_bonf + 0.10,
     paste0(round(false_positive_rate_bonf * 100, 1), "%"),
     col = "darkblue", cex = 1.2, font = 2)

text(bp[2] + 0.4, 0.07, "5% target",
     col = "red", cex = 0.8)



#==== Point_3: A rough demo done in meeting to show when n = 1 v1 returns NA, and the only way
# to obtain is the statiscically INVALID workaround
# of S2 = V2 in order to even get a t statisctic...but it's meaningless:

n1 <- 1
n2 <- 10

group_1 <- rnorm(n1, 0, 1)
group_2 <- rnorm(n2, 5, 1)
m1 <- mean(group_1); m2 <- mean(group_2)
v1 <- var(group_1);  v2 <- var(group_2)


w1 <- (n1-1)/((n1-1)+(n2-1))
w2 <- (n2-1)/((n2-1)+(n1-1))
S2 <- v2

t.stat <- (m1 - m2) / sqrt(S2*(1/n1 + 1/n2))
t.stat
v1
v2
S2

# I'll re write this cleanly, and also run an instability simulation, running same setup 1000x, for eahch time 'll draw 1 observation for 
# group 1 and record the t statistic. I'll plot the distribution of those 1000x t stats to show how much they vary for visual prood
# of our discussed unreliability of this: 

# ---- PART A: Demonstrating the NA problem ----
# Reconstructing Yi et al.'s gross type 1 subgroup (n=1)
# Using simplified groups to demonstrate the statistical problem

set.seed(123)

n1 <- 1   # gross type 1 — single patient S123
n2 <- 6   # gross type 2 — six patients

# Generate groups using paper values for 5-FU gross type
# Gross type 1: S123 GI50 = 3.096 mM
# Gross type 2: mean = 0.41, SD = 0.25
group1 <- rnorm(n1, mean = 3.096, sd = 1)
group2 <- rnorm(n2, mean = 0.41,  sd = 0.25)

# Compute means and variances
m1 <- mean(group1); m2 <- mean(group2)
v1 <- var(group1);  v2 <- var(group2)

# Print to show the problem
cat(paste("Mean group 1:", round(m1, 3), "\n"))
cat(paste("Mean group 2:", round(m2, 3), "\n"))
cat(paste("Variance group 1 (V1):", v1, "\n"))
cat(paste("Variance group 2 (V2):", round(v2, 4), "\n"))

# How we suspect they g0t around it:  
# V1 is NA so we cannot compute a proper pooled variance
# The only option is to use V2 alone as the common variance

S2 <- v2  # using only group 2 variance — group 1 contributes nothing

# Compute degrees of freedom and threshold
df     <- n1 + n2 - 2
thresh <- qt(1 - alpha/2, df)

# Compute t-statistic using the workaround
t.stat <- (m1 - m2) / sqrt(S2 * (1/n1 + 1/n2))

cat(paste("Pooled variance (S2 = V2 workaround):", round(S2, 4), "\n"))
cat(paste("t-statistic:", round(t.stat, 3), "\n"))
cat(paste("Threshold at alpha=0.05:", round(thresh, 3), "\n"))
cat(paste("Decision:", ifelse(abs(t.stat) > thresh, 
                              "Reject H0 — significant", 
                              "Fail to reject H0 — not significant"), "\n"))
#See now if I remove set seed, I'd get a completely dfferent value for the t stat...

# ---- PART B: Instability simulation ----
# Repeat the same setup 1000 times
# Show how wildly the t-statistic varies when n1=1
# This proves the test is unreliable

B <- 1000
t_stats_n1 <- numeric(B)
decisions_n1 <- numeric(B)

set.seed(123)

for(b in 1:B){
  g1 <- rnorm(1,  mean = 3.096, sd = 1)
  g2 <- rnorm(6,  mean = 0.41,  sd = 0.25)
  
  m1 <- mean(g1); m2 <- mean(g2)
  v2 <- var(g2)
  
  # Forced workaround — S2 = V2 only
  S2 <- v2
  
  t_stats_n1[b]   <- (m1 - m2) / sqrt(S2 * (1/1 + 1/6))
  decisions_n1[b] <- 1*(abs(t_stats_n1[b]) > thresh)
}

cat(paste("Proportion of simulations rejecting H0:", 
          round(sum(decisions_n1)/B, 3), "\n"))
cat(paste("Range of t-statistics: from", 
          round(min(t_stats_n1), 2), "to", 
          round(max(t_stats_n1), 2), "\n"))

# Plot distribution of t-statistics
# Shows wild instability when n1=1
hist(t_stats_n1,
     main   = "Distribution of t-statistics\nwhen n=1 in one group",
     xlab   = "t-statistic",
     col    = "lightcoral",
     border = "white",
     breaks = 50,
     xlim   = c(-5, 30))

abline(v =  thresh, col = "red", lty = 2, lwd = 2)
abline(v = -thresh, col = "red", lty = 2, lwd = 2)

text(20, 50, 
     paste("Threshold = ", round(thresh, 2)), 
     col = "red", cex = 0.8)

# Issue with this graph is g1 and g2's means are too widely apart 3.096 vs 0.41 making all t stat values large and +ve
# I'll try making means identucal, so there's no true effect and show t stat: 
# Revised instability simulation — equal means, no true effect
# This better demonstrates the instability caused by n=1
# With equal means the t-statistic should hover around 0
# but with n=1 it bounces wildly

set.seed(123)
t_stats_null <- numeric(B)

for(b in 1:B){
  g1 <- rnorm(1, mean = 0, sd = 1)
  g2 <- rnorm(6, mean = 0, sd = 1)
  
  m1 <- mean(g1); m2 <- mean(g2)
  v2 <- var(g2)
  S2 <- v2
  
  t_stats_null[b] <- (m1 - m2) / sqrt(S2 * (1/1 + 1/6))
}

hist(t_stats_null,
     main   = "t-statistic instability when n=1\n(equal group means: No true effect)",
     xlab   = "t-statistic",
     col    = "lightcoral",
     border = "white",
     breaks = 40)

abline(v =  thresh, col = "red", lty = 2, lwd = 2)
abline(v = -thresh, col = "red", lty = 2, lwd = 2)

# SUMMARY — Point 3 Finding:
# V1 = NA confirmed — variance cannot be computed from n=1
# The only workaround is statistically invalid (S2 = V2 only)
# Under the null hypothesis (equal means, no true effect)
# the t-statistic ranges wildly — showing the test is completely
# unreliable when one group contains a single observation
# Yi et al.'s F(2,6) = 4.701 from this analysis cannot be trusted



# ============================================================
# POINT 5: Unreliable GI50 Estimates From Poor Curve Fitting
# Linked to: February 7 and March 1 — power analysis
# ============================================================

# BACKGROUND:
# For every patient-drug combination, Yi et al. fitted a sigmoidal
# dose-response curve to 7 viability measurements taken at 7 drug
# concentrations. The quality of each curve fit is reported as R²
# in Table II of the paper. This gives 40 R² values in total
# (10 patients x 4 drugs). The two worst fits in the entire dataset
# both occur in the irinotecan column:
# S123 irinotecan: R² = 0.366 — curve explains only 36.6% of data
# S165 irinotecan: R² = 0.397 — curve explains only 39.7% of data

# WHY WE FOCUS ON IRINOTECAN BY STAGE:
# The irinotecan by stage comparison was one of only two results
# Yi et al. reported as statistically significant (p < 0.05).
# The authors split 10 patients into two clinical stage groups:
# Stage I-IIA (early disease): S115, S123, S165, S176, S199 — n=5
# Stage IIB-IV (advanced disease): S126, S137, S150, S222, S236 — n=5
# This grouping was based on Korean clinical guidelines for adjuvant
# chemotherapy. The individual stage codes (IIA, IIIB, IV etc.) from
# the CSV are collapsed into these two broader clinical categories.
# Crucially, S123 — one of the two patients with a poor R² — falls
# in the Stage I-IIA group, meaning an uncertain GI50 estimate
# directly influences the group mean used in the ANOVA.

# THE PROBLEM:
# The GI50 values fed into this ANOVA are not equally reliable.
# Some are derived from excellent curve fits (R² > 0.90) while
# others are derived from poor fits (R² < 0.40). The ANOVA treats
# every GI50 identically regardless of how well it was estimated.
# This inflates within-group variance and reduces statistical power.

# NOTE ON THE NOISE FORMULA/# WHY THIS IS TRICKY:
# R² alone does not directly give us a standard deviation we can plug
# into rnorm(). We need to translate curve fit quality into measurement
# noise. Our approach: noise_SD = reported_SD x sqrt(1 - R²)
# This means:
# noise_SD = reported_SD x sqrt(1 - R²)
# Perfect fit (R²=1.0) → sqrt(0) = 0 → no noise added
# Poor fit (R²=0.366)  → sqrt(0.634) = 0.796 → substantial noise added
# This formula was proposed based on statistical reasoning and
# will be confirmed with Prof. Hernando before finalising.


# Where I'm getting values from (quoted from paper):
# Patients with stages I-IIA showed a lower
# GI50 to irinotecan (Mstage I-IIA=12.52, SD=11.53 vs. Mstage
# IIB-IV=30.95, SD=11.64).

# ---- PART A: Clean power simulation ----
mu_stage1 <- 12.52
mu_stage2 <- 30.95
sd_stage1 <- 11.53
sd_stage2 <- 11.64
n_stages <- 5

# These SD Values already tell me that there's similar variabilities in each group..
# But the SD values themselves relevant to their means respectively are large = meaning the 
# Patients within these groups vary a lot from each other in their response to irinotecan, some very sensitive & some the exact opposite

# Add the Simulation parameters: 
Alpha <- 0.05
B <- 1000 # Number of simulations
sample_sizes <- c(3, 5, 8, 10, 15, 20, 25, 30)
power_clean  <- numeric(length(sample_sizes))

#Now for the double for loop like what I've done before: 
set.seed(123)
power_clean <- numeric(length(sample_sizes))

for(i in 1:length(sample_sizes)){
  
  n_i      <- sample_sizes[i]
  df_i     <- 2 * n_i - 2
  thresh_i <- qt(1 - 0.05/2, df_i)
  
  decisions_i <- numeric(1000)
  
  for(b in 1:1000){
    s1 <- rnorm(n_i, 12.52, 11.53)
    s2 <- rnorm(n_i, 30.95, 11.64)
    
    m1 <- mean(s1); m2 <- mean(s2)
    v1 <- var(s1);  v2 <- var(s2)
    S2 <- 0.5*v1 + 0.5*v2
    
    t.stat <- (m1-m2) / sqrt(S2*(2/n_i))
    decisions_i[b] <- 1*(abs(t.stat) > thresh_i)
  }
  
  power_clean[i] <- sum(decisions_i)/1000
}

# Print results
for(i in 1:length(sample_sizes)){
  cat(paste("n =", sample_sizes[i],
            "| Clean Power =", round(power_clean[i], 3), "\n"))
}

# SUMMARY — Point 5 Part A (Baseline — Clean Simulation):
# This is the OPTIMISTIC estimate — treating every GI50 as perfectly
# reliable with no measurement uncertainty whatsoever.
# We asked a completely different question:
# If the true population means really are 12.52 and 30.95 — how often would an experiment 
# with n=5 per group actually detect that difference?
# The answer was 57.6% of the time.
# That means even if the effect is completely real — even if younger stage patients genuinely
# do respond better to irinotecan across the whole population — an experiment with only 5 
# patients per group would miss that effect 42% of the time....
# 80% power requires approximately n=8 per group.
# Part B will add measurement noise from poor R² curve fits —
# inflating within-group variance further and reducing power
# even more. If the baseline is already insufficient, the true
# situation accounting for GI50 measurement uncertainty is worse.

#---- PART B: Adding noise simulation ----
# Adding measurement noise based on R² values from Table II
# Each patient's GI50 uncertainty reflected by their curve fit quality

# R² values for irinotecan by stage (from Table II): 
r2_irinotecan <- c(0.8261,  # S115
                  0.3663,  # S123 — poor fit
                  0.5561,  # S126
                  0.6901,  # S137
                  0.5437,  # S150
                  0.3972,  # S165 — poor fit
                  0.8376,  # S176
                  0.8114,  # S199
                  0.8853,  # S222
                  0.6369)  # S236
# Now to split these values by stage group: Same as they did: 
# Stage I-IIA: S115, S123, S165, S176, S199
r2_stage1 <- c(0.8261, 0.3663, 0.3972, 0.8376, 0.8114)
# Stage IIB-IV: S126, S137, S150, S222, S236
r2_stage2 <- c(0.5561, 0.6901, 0.5437, 0.8853, 0.6369)
# obtain means for each: 
mean_r2_stage1 <- mean(r2_stage1)
mean_r2_stage2 <- mean(r2_stage2)
cat(paste("Mean R² for Stage I-IIA group:", round(mean_r2_stage1, 3), "\n"))
cat(paste("Mean R² for Stage IIB-IV group:", round(mean_r2_stage2, 3), "\n"))
# Something to note here: stage I-IIA has tow very bad fits whereas the other group is more consistently mediocre with its results...
# This means the measurement reliability will differ between these two groups: Stage I-IIA will have an inflated within group variance. 
# Calculate noise SD for each group
# Formula: noise_SD = reported_SD x sqrt(1 - mean_R²)

noise_sd_stage1 <- sd_stage1 * sqrt(1 - mean_r2_stage1)
noise_sd_stage2 <- sd_stage2 * sqrt(1 - mean_r2_stage2)

# Effective SD — reported SD plus measurement noise
eff_sd_stage1 <- sd_stage1 + noise_sd_stage1
eff_sd_stage2 <- sd_stage2 + noise_sd_stage2

cat(paste("Original SD Stage I-IIA:", sd_stage1, "\n"))
cat(paste("Effective SD Stage I-IIA:", round(eff_sd_stage1, 3), "\n"))
cat(paste("Original SD Stage IIB-IV:", sd_stage2, "\n"))
cat(paste("Effective SD Stage IIB-IV:", round(eff_sd_stage2, 3), "\n"))
cat(paste("So theres a difference of", round(eff_sd_stage1 - sd_stage1, 3),"for SD Stage I-IIA, and a difference of", round(eff_sd_stage2 - sd_stage2, 3), "for SD Stage IIB-IV, showing that the measurement noise inflates the variability in both groups."))
# Now to re-run the power simulation with these effective SDs, Use same loop as before but replacing the SD used with these eff SDs now.
set.seed(123)
power_noisy <- numeric(length(sample_sizes))

for(i in 1:length(sample_sizes)){
  
  n_i <- sample_sizes[i]
  df_i <- 2 * n_i - 2
  thresh_i <- qt(1 - 0.05/2, df_i)
  
  decisions_i <- numeric(1000)
  
  for(b in 1:1000){
    s1 <- rnorm(n_i, 12.52, eff_sd_stage1)
    s2 <- rnorm(n_i, 30.95, eff_sd_stage2)
    
    m1 <- mean(s1); m2 <- mean(s2)
    v1 <- var(s1);  v2 <- var(s2)
    S2 <- 0.5*v1 + 0.5*v2
    
    t.stat <- (m1-m2) / sqrt(S2*(2/n_i))
    decisions_i[b] <- 1*(abs(t.stat) > thresh_i)
  }
  
  power_noisy[i] <- sum(decisions_i)/1000
}

# Print results
for(i in 1:length(sample_sizes)){
  cat(paste("n =", sample_sizes[i],
            "| Clean Power =", round(power_clean[i], 3),
            "| Noisy Power =", round(power_noisy[i], 3), "\n"))
}

# Visualize results: Both curves
# Plot clean vs noisy power curves
par(mfrow = c(1,1))

plot(sample_sizes, power_clean,
     type = "b", pch = 19, col = "blue",
     xlab = "Sample size per group (n)",
     ylab = "Estimated Power",
     main = "Power Curve: Clean vs Noisy\nIrinotecan by Stage (Yi et al. 2023)",
     ylim = c(0, 1),
     xlim = c(0, 32))

lines(sample_sizes, power_noisy,
      type = "b", pch = 17, col = "red")

abline(h = 0.80, lty = 2, col = "grey50", lwd = 2)
abline(v = 5,    lty = 2, col = "grey50", lwd = 2)

legend("bottomright",
       legend = c("Clean — no noise",
                  "Noisy — R² adjusted"),
       col    = c("blue", "red"),
       pch    = c(19, 17),
       lty    = 1)

text(7, 0.85, col = "grey30", cex = 0.8)
text(30, 0.83, "80% threshold",  col = "grey30", cex = 0.8)

# Summary: 
# SUMMARY — Point 5 Complete Finding:
# Irinotecan by stage — Yi et al.'s second significant result
# Clean power at n=5 = 0.576 (57.6%) — already below 80%
# Noisy power at n=5 = 0.289 (28.9%) — drops to less than 1 in 3
# To reach 80% power:
# Without noise: need n=8 per group
# With R² noise: need n=20 per group
# Yi et al. used n=5 — under realistic conditions they had
# less than 30% power for their second significant result
# This finding should be confirmed with Prof. Hernando
# regarding the noise formula: noise_SD = SD x sqrt(1-R²)
# combined as: effective_SD = sqrt(SD² + noise_SD²)


#===Extension of PART B, displaying multiple noise levels===
# Hernando's suggestion: show how power degrades 
# progressively as R² worsens
# Using Stage I-IIA SD as reference (sd_stage1 = 11.53)
# and Stage IIB-IV SD (sd_stage2 = 11.64)

r2_levels <- c(0.90, 0.70, 0.50, 0.366)
r2_labels <- c("R² = 0.90 (good fit)",
               "R² = 0.70 (moderate fit)", 
               "R² = 0.50 (poor fit)",
               "R² = 0.366 (worst case — S123)")

# Empty matrix to store power values for each R² level
# Rows = sample sizes, Columns = R² levels
power_matrix <- matrix(NA, 
                       nrow = length(sample_sizes), 
                       ncol = length(r2_levels))

set.seed(123)

for(j in 1:length(r2_levels)){
  
  r2 <- r2_levels[j]
  
  # Calculate effective SD for this R² level
  noise_sd1 <- sd_stage1 * sqrt(1 - r2)
  noise_sd2 <- sd_stage2 * sqrt(1 - r2)
  eff_sd1   <- sqrt(sd_stage1^2 + noise_sd1^2)
  eff_sd2   <- sqrt(sd_stage2^2 + noise_sd2^2)
  
  for(i in 1:length(sample_sizes)){
    
    n_i      <- sample_sizes[i]
    df_i     <- 2 * n_i - 2
    thresh_i <- qt(1 - 0.05/2, df_i)
    
    decisions_i <- numeric(1000)
    
    for(b in 1:1000){
      s1 <- rnorm(n_i, 12.52, eff_sd1)
      s2 <- rnorm(n_i, 30.95, eff_sd2)
      m1 <- mean(s1); m2 <- mean(s2)
      v1 <- var(s1);  v2 <- var(s2)
      S2 <- 0.5*v1 + 0.5*v2
      t.stat <- (m1-m2) / sqrt(S2*(2/n_i))
      decisions_i[b] <- 1*(abs(t.stat) > thresh_i)
    }
    
    power_matrix[i, j] <- sum(decisions_i)/1000
  }
}
# Plotting all curves on the same graph
# Print results for each R² level
for(j in 1:length(r2_levels)){
  cat(paste("\n", r2_labels[j], "\n"))
  for(i in 1:length(sample_sizes)){
    cat(paste("n =", sample_sizes[i], 
              "| Power =", round(power_matrix[i,j], 3), "\n"))
  }
}

# Plot all power curves together
par(mfrow = c(1,1))

colours <- c("steelblue", "orange", "red3", "darkred")

plot(sample_sizes, power_clean,
     type = "b", pch = 19, col = "blue", lwd = 2,
     xlab = "Sample size per group (n)",
     ylab = "Estimated Power",
     main = "Power Degradation Across R² Levels\nIrinotecan by Stage (Yi et al. 2023)",
     ylim = c(0, 1),
     xlim = c(0, 32))

for(j in 1:length(r2_levels)){
  lines(sample_sizes, power_matrix[, j],
        type = "b", pch = 17,
        col  = colours[j],
        lwd  = 2)
}

abline(h = 0.80, lty = 2, col = "grey50", lwd = 2)
abline(v = 5,    lty = 2, col = "grey50", lwd = 2)

legend("bottomright",
       legend = c("Clean — no noise", r2_labels),
       col    = c("blue", colours),
       pch    = c(19, 17, 17, 17, 17),
       lty    = 1,
       lwd    = 2,
       cex    = 0.8)

