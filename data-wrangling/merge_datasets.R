# Load in sorted datasets

setwd(paste(getwd(), "/B419", sep = ""));

d_tbhiv <- read.csv('data/d_tbhiv.csv')
d_tbresis <- read.csv ('data/d_tbresis.csv')
d_tbtreat <- read.csv('data/d_tbtreat.csv')
d_tobacco <- read.csv('data/d_tobacco.csv')
d_alcohol <-read.csv('data/d_alcohol.csv')
d_diabetes <- read.csv('data/d_diabetes.csv')
d_healthcare <- read.csv('data/d_healthcare.csv')
d_water <- read.csv('data/d_water.csv')
d_lifestats <- read.csv('data/d_lifestats.csv')
d_popden <- read.csv('data/d_popden.csv')
d_chnutri <- read.csv('data/d_chnutri.csv')

# Merge datasets

merged_data1 <- merge(d_tbhiv, d_tbresis, by = c("Country", "Year"), 
                      all = TRUE, sort = TRUE)
merged_data2 <- merge(merged_data1, d_tbtreat, by = c("Country", "Year"), 
                     all = TRUE, sort = TRUE)
merged_data3 <- merge(merged_data2, d_tobacco, by = c("Country", "Year"), 
                      all = TRUE, sort = TRUE)
merged_data4 <- merge(merged_data3, d_alcohol, by = c("Country", "Year"), 
                      all = TRUE, sort = TRUE)
merged_data5 <- merge(merged_data4, d_diabetes, by = c("Country", "Year"), 
                      all = TRUE, sort = TRUE)
merged_data6 <- merge(merged_data5, d_healthcare, by = c("Country", "Year"), 
                      all = TRUE, sort = TRUE)
merged_data7 <- merge(merged_data6, d_water, by = c("Country", "Year"), 
                      all = TRUE, sort = TRUE)
merged_data8 <- merge(merged_data7, d_lifestats, by = c("Country", "Year"), 
                      all = TRUE, sort = TRUE)
merged_data9 <- merge(merged_data8, d_popden, by = c("Country", "Year"), 
                      all = TRUE, sort = TRUE)
final_dataset <- merge(merged_data9, d_chnutri, by = c("Country", "Year"), 
                      all = TRUE, sort = TRUE)

write.csv(final_dataset, 'data/merged_all_datasets.csv')

