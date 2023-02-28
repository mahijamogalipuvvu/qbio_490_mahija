#knitr::opts_knit$set(root.dir = normalizePath("/Users/mahijamogalipuvvu/Desktop/QBIO490/qbio_490_mahija/analysis_data")) 
setwd("/Users/mahijamogalipuvvu/Desktop/QBIO490/qbio_490_mahija/analysis_data") #setting working directory
getwd()

#installing and loading all necessary packages
if(!require(survival)) {
  install.packages("survival")
}

if(!require(survminer)) {
  install.packages("survminer")
}

if(!require(ggplot2)) {
  install.packages("ggplot2")
}
library(survival)
library(survminer)
library(ggplot2)

library(BiocManager)
library(TCGAbiolinks)
library(maftools)

clin_query <- GDCquery(project = "TCGA-BRCA", 
                       data.category = "Clinical", 
                       file.type = "xml")
GDCdownload(clin_query)

#reading in clinical data
clinic <- read.csv("/Users/mahijamogalipuvvu/Desktop/QBIO490/qbio_490_mahija/analysis_data/brca_clinical_data.csv")

clinical_drug <- GDCprepare_clinic(query = clin_query,
                                   clinical.info = "drug")

clinical_rad <- GDCprepare_clinic(query = clin_query,
                                  clinical.info = "radiation")

colnames(clinic)                                  #printing column names
print(sum(!is.na(clinic$days_to_birth)))          #printing number of true values
#i chose days to birth as my variable which is a discrete variable.

colnames(clinical_drug)
print(sum(!is.na(clinical_drug$days_to_drug_therapy_end)))
#i chose the number of days till their treatment ends which is a discrete varaiable
# this variable provides the number of days till their treatment with the drug ends 

colnames(clinical_rad)
print(sum(!is.na(clinical_rad$radiation_dosage)))
#i chose the radiation dosage which is a discrete variable
#this variable tells the severity and amount of radiation given to the patient

#hypothesis 1: Patients with more days till their drug therapy end may have higher radiation dosage. 
#hypothesis 2: As the number of days till their drug therapy ends decreases, their chances of survival increase. 
#hypothesis 3: The higher the radiation dosage, the higher the chances of survival. 

#typeof(clinical_rad$radiation_dosage)

#cleaning data for plotting
cleanNaMaskDrug <- ifelse(is.na(clinical_drug$days_to_drug_therapy_end), FALSE, TRUE)
cleanClinicalDrug <- clinical_drug[cleanNaMaskDrug, ]
head(cleanClinicalDrug)

cleanNaMaskRad <- ifelse(!is.na(clinical_rad$radiation_dosage), TRUE, FALSE)
cleanClinicalRad <- clinical_rad[cleanNaMaskRad, ]
head(cleanClinicalRad)

#merging data
mergeClinic <- merge(cleanClinicalDrug, cleanClinicalRad, by = "bcr_patient_barcode")

#creating a scatter plot
plot(x = mergeClinic$days_to_drug_therapy_end, 
     y = mergeClinic$radiation_dosage,
     main = "The Relationship Between Number of Days Till the End of Drug Therapy and Radiation Dosage",
     xlab = "Number of Days Till the End of Drug Therapy",
     ylab = "Radiation Dosage")

#i created a scatter plot bc they're both discrete/continuous variables
#from the plot, i gleaned that there may be no relationship between the two variables

#VARIABLE 1: DAYS TO DRUG THERAPY END

#creating a column w/ same name to merge both datasets
clinical_drug$Tumor_Sample_Barcode <- clinical_drug$bcr_patient_barcode

#merging based off tumor sample barcode
mergeClinicDrug <- merge(clinic, clinical_drug, by = "Tumor_Sample_Barcode")

#cleaning for NAs
cleanNaMaskMergeDrug <- ifelse(is.na(mergeClinicDrug$Tumor_Sample_Barcode), FALSE, TRUE)
mergeClinicDrug <- mergeClinicDrug[cleanNaMaskMergeDrug, ]

#no categorical variable so creating + cleaning one for the KM plot
lowMask <- ifelse(mergeClinicDrug$days_to_drug_therapy_end <= 2000, T, F)
midMask <- ifelse(mergeClinicDrug$days_to_drug_therapy_end > 2000 & mergeClinicDrug$days_to_drug_therapy_end <= 4000, T, F)
mergeClinicDrug$daysStatus <- ifelse(lowMask, "Low", ifelse(midMask, "Middle", "High"))

#creating a survival time column to plot by
mergeClinicDrug$survivalTime <- ifelse (is.na(mergeClinicDrug$days_to_death),
                                        mergeClinicDrug$survivalTime <- mergeClinicDrug$days_to_last_followup,
                                        mergeClinicDrug$survivalTime <- mergeClinicDrug$days_to_death)



#remove -Inf values in the survivalTime column
inf_maskD <- ifelse(mergeClinicDrug$survivalTime == "-Inf", F, T)
mergeClinicDrug <- mergeClinicDrug[inf_maskD, ]

#create a death event column
mergeClinicDrug$death_event <- ifelse(mergeClinicDrug$vital_status == "Alive", mergeClinicDrug$death_event <- FALSE, mergeClinicDrug$death_event <- TRUE)

#creating a survival status column to plot by
mergeClinicDrug$survivalStatus <- ifelse(mergeClinicDrug$vital_status == "Dead", T, F)

#create survminer object
survivalTime <- ifelse(mergeClinicDrug$vital_status == "Dead", mergeClinicDrug$days_to_death, mergeClinicDrug$days_to_last_followup)
survivalStatus <- ifelse(mergeClinicDrug$vital_status == "Dead", T, F)

survivalObjectDays <- Surv(time = survivalTime, event = survivalStatus)

fitObject <- survfit(survivalObjectDays ~ daysStatus, data = mergeClinicDrug)

survplot <- ggsurvplot(fitObject, 
                       pval=TRUE,
                       ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                       legend = "right")

#saving plots to folder on local computer
jpeg("/Users/mahijamogalipuvvu/Desktop/QBIO490/KM_plotD.jpg")

#creating the KM plot
KM_plotD <- survplot$plot + theme_bw() + theme(axis.title =
                                                element_text(size=20), axis.text = element_text(size=16),
                                              legend.title = element_text(size=14), legend.text =
                                                element_text(size=12))

#the KM plots suggest that the more number of days till drug therapy ends may mean a lesser survival probability
#the p-value is 0.11 which means this data is not significant
#the differences seem slightly signfiicant

KM_plotD
dev.off()

#VARIABLE 2: RADIATION DOSAGE

#creating a column w/ same name to merge both datasets
clinical_rad$Tumor_Sample_Barcode <- clinical_rad$bcr_patient_barcode

#merging based off tumor sample barcode
mergeClinicRad <- merge(clinic, clinical_rad, by = "Tumor_Sample_Barcode")

cleanNaMaskMergeRad <- ifelse(is.na(mergeClinicRad$Tumor_Sample_Barcode), FALSE, TRUE)
mergeClinicRad <- mergeClinicRad[cleanNaMaskMergeRad, ]

#no categorical variable so creating + cleaning one for the KM plot
lowRMask <- ifelse(mergeClinicRad$radiation_dosage == 60, T, F)
midRMask <- ifelse(mergeClinicRad$radiation_dosage == 6040, T, F)
highRMask <- ifelse(mergeClinicRad$radiation_dosage == 10000, T, F)
mergeClinicRad$radiationStatus <- ifelse(lowRMask, "Low", ifelse(midRMask, "Middle", "High"))

#creating a survival time column to plot by
mergeClinicRad$survivalTime <- ifelse (is.na(mergeClinicRad$days_to_death),
                                       mergeClinicRad$survivalTime <- mergeClinicRad$days_to_last_followup,
                                       mergeClinicRad$survivalTime <- mergeClinicRad$days_to_death)



#remove -Inf values in the survivalTime column
inf_maskR <- ifelse(mergeClinicRad$survivalTime == "-Inf", F, T)
mergeClinicRad <- mergeClinicRad[inf_maskR, ]

#create a death event column
mergeClinicRad$death_event <- ifelse(mergeClinicRad$vital_status == "Alive", mergeClinicRad$death_event <- FALSE, mergeClinicRad$death_event <- TRUE)

#creating a survival status column to plot by
mergeClinicRad$survivalStatus <- ifelse(mergeClinicRad$vital_status == "Dead", T, F)

#create survminer object
survivalTimeR <- ifelse(mergeClinicRad$vital_status == "Dead", mergeClinicRad$days_to_death, mergeClinicRad$days_to_last_followup)
survivalStatusR <- ifelse(mergeClinicRad$vital_status == "Dead", T, F)

survivalObjectRad <- Surv(time = survivalTimeR, event = survivalStatusR)

fitObjectR <- survfit(survivalObjectRad ~ radiationStatus, data = mergeClinicRad)

survplotR <- ggsurvplot(fitObjectR, 
                       pval=TRUE,
                       ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                       legend = "right")

#saving plots to folder on local computer
jpeg("/Users/mahijamogalipuvvu/Desktop/QBIO490/KM_plotR.jpg")

#creating KM plot
KM_plotR <- survplotR$plot + theme_bw() + theme(axis.title =
                                                element_text(size=20), axis.text = element_text(size=16),
                                              legend.title = element_text(size=14), legend.text =
                                                element_text(size=12))

#the KM plots don't suggest a pattern of any sort 
#the p-value is 0.17 which means this data is not significant
#the differences do not signfiicant
KM_plotR
dev.off()

