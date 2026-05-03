data <- read.csv(file.choose())
dim(data)
colnames(data)
unique(data$CANCER_TYPE)

# Isolate Colorectaal Carcinoma & all the 4 drugs used in yo et al if they're also here!

data_coread <- data[data$CANCER_TYPE == "Colorectal Carcinoma", ]
cat(paste("Number of colorectal carcinoma entries:", nrow(data_coread), "\n"))
unique(data_coread$DRUG_NAME)
length(unique(data_coread$DRUG_NAME))

grep("luorou", unique(data_coread$DRUG_NAME), value = TRUE)
grep("isplatin", unique(data_coread$DRUG_NAME), value = TRUE)
grep("xaliplatin", unique(data_coread$DRUG_NAME), value = TRUE)
grep("rinotecan", unique(data_coread$DRUG_NAME), value = TRUE)

# I want to check how many entres for each drug in colorectal carcinoma
yi_drugs <- c("5-Fluorouracil", "Cisplatin", "Oxaliplatin", "Irinotecan")

for(drug in yi_drugs){
  n <- nrow(data_coread[data_coread$DRUG_NAME == drug, ])
  cat(paste(drug, ":", n, "cell lines\n"))
}

# Must take into account,  here we have IC50 values & not GI50 values as in yi et al. We're also not working with organoids but with established cancer cell lines. 
# We have IC50s but I think I can also get the raw data needed to obtain those reported IC50s and hence check the r^2  curveis  fits, This was also an issue on Yi et al.


# Now I can continue as I am, OR do I look for an organoid based dataset instead?