#MULTIPLE REGRESSION- GEM-MEDOME#
#UPDATED 2019-08-02
#DRAFTED BY KRISH BILIMORIA

#Calling libraries from other sourecs
#need (data.table)
library(magrittr)
library(glmnet)

getwd() #confirm that working directory is: "/~/essential-meds"

#Importing input variables, initial transformation
gem <- data.table::fread("input/gem.csv",data.table = F) #each medicine is one row, each country is one column
gem$sum <- rowSums(gem[,4:140]) #this calculates the number of times medicines show up on a country list, columns 1-4 are titles/headers
gem$maf <- gem$sum/137
rownames(gem) <- gem$medicine
drug.names <- gem$medicine

tgem <- gem[,4:140] %>% #creating a transposed variable of medicines for regression, each medicine is one column, each country is one row
  t %>%
  as.data.frame()

haq <- data.table::fread("input/haq.csv", data.table = F) #each country is one column, and the first row is "overall HAQ", which is our outcome of interest

thaq <- haq[1,2:138] %>% #subsetting for "HAQ_overall" from all of the different hHAQ measurements
  t %>%
  `colnames<-` (c("haq")) %>%
  as.data.frame()
thaq$country <- rownames(thaq)
           
covar <- data.table::fread("input/covar.csv", data.table = F) #covariates
cmeds <- data.table::fread("input/cleanmeds.csv", data.table = F) #medicines on cleanmeds list

#Variable subsetting- missing data summary

#MISSING COVARIATES (4 countries)
#Cook Islands, covar row 30- no health-exp, pop, gdp
#Niue, covar row 89- no halth-exp, pop, gdp
#Somalia, 112- no health-exp, gdp
#Syrian arab republic, row 118- no gdp

#MISSING HAQ INDICES: (7 countries)
#Botswana, thaq
#Cook Islands, thaq
#Democratic Republic of Congo, thaq
#Niue, thaq
#Palau, thaq
#Saint Kitts & Nevis, thaq
#Tuvalu, thaq

#Combining HAQ score with covariates
haq.covar <- merge(thaq,covar[,c('country',"health_expenditure", "population", "gdp")],by.x='country',by.y='country',all.x=T)

#Logging covariates
haq.covar[c("health_expenditure","population","gdp")] <- lapply(haq.covar[c("health_expenditure","population","gdp")],log)

#Combining predictors with covariates
tgem$country <- rownames(tgem) #137x2069 at this time
data.full <- merge(haq.covar,tgem, by.x='country', by.y='country', all.x=T) %>% #137x2073 
  as.data.frame()
data.full[c("haq","health_expenditure","population","gdp")] <- lapply(data.full[c("haq", "health_expenditure","population","gdp")],as.numeric)

#Parallel removal of NAs
data.full.na <- na.omit(data.full) #127x2073

#LINEAR MODEL WITH COVARIATES

#Test linear model with just covariates
my.lms.covar <- lm(haq ~ health_expenditure + population + gdp, data = data.full) 

#Running 2068 regressions, including rows with incomplete data, allowing regression model to account for differences
my.lms <- lapply(drug.names, function(x) lm(formula(haq ~ data.full[,x] + health_expenditure + population + gdp), data = data.full)) %>%
  setNames(drug.names)

#Running 2068 regressions, on a dataframe that has removed incomplete data
my.lms.na <- lapply(drug.names, function(x) lm(formula(haq ~ data.full.na[,x] + health_expenditure + population + gdp), data = data.full.na)) %>%
  setNames(drug.names)

#extracting summaries
my.lms.summaries <- lapply(my.lms, summary)
my.lms.na.summaries <- lapply(my.lms.na, summary)

#P-values, rsquared values
pval <- lapply(my.lms.summaries, function(x) c(x$coefficients[2,c(1:4)], r_sq=x$r.squared, adj_r_sq= x$adj.r.squared))
pval.na <- lapply(my.lms.na.summaries,function(x) c(x$coefficients[2,c(1:4)], r_sq=x$r.squared, adj_r_sq= x$adj.r.squared))

#Subsetting the list variable, turning into dataframe
pval.df=dplyr::bind_rows(pval) %>% t
colnames(pval.df) <- c("estimate", "standard-error", "tvalue", "pvalue", "r_sq", "adj_r_sq")
pval.df <- as.data.frame(pval.df)


pval.df.na=dplyr::bind_rows(pval.na) %>% t
colnames(pval.df.na) <- c("estimate", "standard-error", "tvalue", "pvalue", "r_sq", "adj_r_sq")
pval.df.na <- as.data.frame(pval.df.na)

#Adding ATC codes, sum, maf
pval.df$drug.names <- row.names(pval.df) 
pval.df <- merge(pval.df, gem[,c("medicine","ATC code primary", "ATC code secondary", "sum", "maf")], by.x = "drug.names", by.y = "medicine")
pval.df <- pval.df[order(pval.df$`ATC code primary`),]

pval.df.na$drug.names <- row.names(pval.df.na) 
pval.df.na <- merge(pval.df.na, gem[,c("medicine","ATC code primary", "ATC code secondary", "sum", "maf")], by.x = "drug.names", by.y = "medicine")
pval.df.na <- pval.df[order(pval.df$`ATC code primary`),]

#Removing rare medicines, and common medicines from the list of medicines
#dim(pval.df) = 2068x11
pval.df.maf <- pval.df[-which(pval.df$maf <0.05),] #remove rare medicines, 1137x11
pval.df.maf <- pval.df.maf[-which(pval.df.maf$maf > 0.95),] #remove very common medicines, 1114x11

#Annotating CLEANMEDS List with significant findings

#make med names lowercased
pval.df.lower <- pval.df
pval.df.lower$drug.names <- tolower(pval.df.lower$drug.names)

cleanmeds.comparison <- merge(cmeds, pval.df.lower, by.x="Medication", by.y="drug.names")

#The following names are not consistent between the two databases, will change later
#[1] "acetaminophen"                                          
#[2] "alendronate"                                            
#[3] "amoxicillin/clavulanic acid"                            
#[4] "beclomethasone"                                         
#[5] "benztropine"                                            
#[6] "cephalexin"                                             
#[7] "chlorthalidone"                                         
#[8] "epinephrine"                                            
#[9] "ethinylestradiol/levonorgestrel"                        
#[10] "ferrous fumarate"                                       
#[11] "levodopa/carbidopa"                                     
#[12] "levonorgestrel \xd0 releasing intrauterine system"      
#[13] "nitroglycerin"                                          
#[14] "nortriptyline"                                          
#[15] "polyethylene glycol 3350"                               
#[16] "polymyxin b"                                            
#[17] "senna"                                                  
#[18] "sulfamethoxazole/ trimethoprim"                         
#[19] "sulfasalazine"                                          
#[20] "thiamine"                                               
#[21] "urea"                                                   
#[22] "vaginal ring eluting etonogestrel and ethinyl estradiol"
#[23] "valacyclovir"                                           
#[24] "valproic acid"                                          
#[25] "vitamin B12"                                            
#[26] "vitamin D" 

#outputs
write.csv(pval.df, file="output/results-full.csv") #2068x11 (all countries, even with missing HAQ/covariate data)
write.csv(pval.df.na, file="output/results-full-na.csv") #2068x11 (removed countries with missing HAQ/covariate data)
write.csv(pval.df.maf, file="output/results-full-maf.csv") #1137x11
write.csv(cleanmeds.comparison, file="output/results-cleanmeds.csv") #102x20, 26 medicines removed due to inconsistent naming

#end