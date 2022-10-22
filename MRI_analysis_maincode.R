## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)


## --------------------------------------------------------------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(kableExtra)
source("MRI_analysis_funs.R")
library(neuroCombat)
library(ggpubr)
library(DescTools)
library(flextable)
library(table1)
library(lme4)
library(caret)
library(pROC)
library(gridExtra)
library(mice)
library(purrr)
library(DescTools)
library(OptimalCutpoints)
library(grid)



## --------------------------------------------------------------------------------------------------------------------------------
##########Clean and format data##########################

#Use latest GCS but do not impute unknown values
gcs <- read.csv("../Data/CENTER_injury_data.csv", na.strings = c("", " ", "NA"))
gcs$Latest_GCS <- ifelse(is.na(gcs$InjuryHx.GCSEDDischScore)==FALSE, gcs$InjuryHx.GCSEDDischScore, 
                         ifelse(is.na(gcs$InjuryHx.GCSEDArrScore)==FALSE, gcs$InjuryHx.GCSEDArrScore,
                                ifelse(is.na(gcs$InjuryHx.GCSFirstHospScore)==FALSE, gcs$InjuryHx.GCSFirstHospScore,
                                       ifelse(is.na(gcs$InjuryHx.GCSPreHospBestScore)==FALSE, gcs$InjuryHx.GCSPreHospBestScore,
                                       gcs$InjuryHx.GCSPreHospBestReportedTotalScore))))
gcs$Latest_GCS_cat <- ifelse(gcs$Latest_GCS %in% c(15, 14, 13) | gcs$Subject.PatientType == "1", "mild",
                             ifelse(gcs$Latest_GCS %in% c(12, 11, 10, 9, 8, 7, 6, 5, 4, 3), "mod-sev", "untestable"))
gcs <- 
  gcs %>%
  select(Master_subject_ID = subjectId,
         Latest_GCS,
         Latest_GCS_cat)

gcs$Latest_GCS <- factor(gcs$Latest_GCS, 
                         ordered = TRUE,
                         levels = c("3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"))


data <- gcs

#Marshall score
marshall <- read.csv("../Data/CENTER_imaging_data.csv", na.strings = c("", " ", "NA"))
marshall <- 
  marshall %>%
  select(Master_subject_ID = subjectId,
         Marshall_score_central = Imaging.MarshallCTClassification) %>%
  filter(is.na(Marshall_score_central) == FALSE) %>%
  group_by(Master_subject_ID) %>%
  slice(1)
data <- merge(data, marshall, by = "Master_subject_ID", all.x = TRUE, all.y = TRUE)

#supplement missing Marshall scores from Virginias review
clinical <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Legacy_data/Most_uptodate_curated_data/Sophies_clinical_database_20220822.csv", na.strings = c("", " ", "NA"))
clinical <- 
  clinical %>%
  select(Master_subject_ID, 
         Marshall_score_vfjn = Marshall_score) %>%
  filter(is.na(Master_subject_ID)==FALSE)
data <- merge(clinical, data, by = "Master_subject_ID", all.x = FALSE, all.y = TRUE)
data$Marshall_score <- ifelse(is.na(data$Marshall_score_central)== FALSE, data$Marshall_score_central, data$Marshall_score_vfjn)
data <- 
  data %>%
  select(-c(Marshall_score_central,
            Marshall_score_vfjn))
data$Marshall_score <- ifelse(is.na(data$Marshall_score)==TRUE, "unknown", data$Marshall_score)

imdate <- read.csv("../Data/CENTER_imaging_data.csv", na.strings = c("", " ", "NA"))
imdate$Imaging.ExperimentDate <- as.Date(imdate$Imaging.ExperimentDate, format = "%d/%m/%Y")
imdate$Imaging.ExperimentId <- gsub("CTBI_", "", imdate$Imaging.ExperimentId)
imdate <-
  imdate %>% 
  filter(str_detect(Imaging.ExperimentLabel, "MR")) %>%
  group_by(subjectId,Imaging.ExperimentDate) %>%
  slice(1) %>% #pick one line per scan-date
  ungroup() %>%
  group_by(subjectId) %>%
  arrange(Imaging.ExperimentDate, .by_group = TRUE) %>%
  slice(c(1, 2)) %>% #keep first two scan-dates
  mutate(Scan_number = row_number()) %>%
  ungroup()
#calculate days to scan for first and second scan
#in case the first scan misses images or quantitative data
imdate$Days_to_MR <- (imdate$Imaging.ExperimentDate - as.Date("1970-01-01")) %>% as.numeric()
imdate1 <-
  imdate %>%
  select(subjectId, Scan_number, Days_to_MR) %>%
  spread(key = Scan_number, value = Days_to_MR)
colnames(imdate1) <- c("subjectId", "Days_to_MR1", "Days_to_MR2")
imdate2 <-
  imdate %>%
  select(subjectId, Scan_number, Imaging.ExperimentId) %>%
  spread(key = Scan_number, value = Imaging.ExperimentId) 
colnames(imdate2) <- c("subjectId", "Scan_ID_MR1", "Scan_ID_MR2")
imdate <-
  merge(imdate1, imdate2, by = "subjectId", all.x = TRUE, all.y = TRUE)
data <- 
  merge(data, imdate, by.x = "Master_subject_ID", by.y = "subjectId", all.x = TRUE, all.y = TRUE)

#Add of time of injury to calculate how long after the 
#serum biomarker was sampled
time <- read.csv("../Data/CENTER_timeofinj_data.csv", na.strings = c("", " ", "NA"))

time$Date_inj <- "1970-01-01"
time$Date_Time_inj <- as.POSIXct(paste(time$Date_inj, time$Subject.TimeInj), format="%Y-%m-%d %H:%M:%S")
time <- 
  time %>% 
  select(Master_subject_ID = subjectId,
         Date_Time_inj)

data <- merge(data, time, by = "Master_subject_ID", all.x = TRUE, all.y = TRUE)

#add serum biomarker data
bio <- read.csv("../Data/CENTER_biomarker_data.csv", na.strings = c("", " ", "NA"))
bio$Date_Time_bio <- as.POSIXct(paste(bio$Biomarkers.CollectionDate, bio$Biomarkers.CollectionTime), format="%d/%m/%Y %H:%M:%S")
colnames(bio) <- gsub("Biomarkers.", "", colnames(bio), fixed = TRUE)
bio <- 
  bio %>%
  select(Master_subject_ID = subjectId,
         Date_Time_bio,
         GFAP, 
         NFL,
         NSE,
         S100B,
         Tau,
         UCH.L1) %>%
  group_by(Master_subject_ID) %>% 
  slice(which.min(Date_Time_bio))

#Select only patients with biomarker samples within 24h of injury
data <- merge(data, bio, by = "Master_subject_ID", all.x = TRUE, all.y = TRUE)
data$Time_to_bio <- difftime(data$Date_Time_bio, data$Date_Time_inj, units = "hours")
data$Time_to_bio <- abs(as.numeric(data$Time_to_bio))

#Log transform biomarker concentrations
data <-
  data %>%
  mutate(across(c("GFAP", "NFL", "NSE", "S100B", "Tau", "UCH.L1"), 
         log, 
         .names = "{.col}"))

#add auxillary data
aux <- read.csv("../Data/CENTER_aux_data.csv", na.strings = c("", " ", "NA", "88"))
colnames(aux) <- gsub("InjuryHx.", "", colnames(aux))
colnames(aux) <- gsub("Subject.", "", colnames(aux))
aux$MajorEI <- ifelse(aux$PelvicGirdleAIS >= 3 |
                          aux$UpperExtremitiesAIS >= 3 |
                          aux$ThoracicSpineAIS >= 3 |
                          aux$AbdomenPelvicContentsAIS >= 3 |
                          aux$ExternaAIS >= 3 |
                          aux$ThoraxChestAIS >= 3 |
                          aux$LumbarSpineAIS >= 3 |
                          aux$LowerExtremitiesAIS >= 3, "present", "absent")

aux$MI_HeadNeck <- ifelse(is.na(aux$HeadNeckAIS) == TRUE, NA,
                          ifelse(aux$HeadNeckAIS >= 3, "present", "absent"))
aux$MI_Face <- ifelse(is.na(aux$FaceAIS) == TRUE, NA,
                          ifelse(aux$FaceAIS >= 3, "present", "absent"))
aux$MI_ThoracicSpine <- ifelse(is.na(aux$ThoracicSpineAIS) == TRUE, NA,
                          ifelse(aux$ThoracicSpineAIS >= 3, "present", "absent"))
aux$MI_Chest <- ifelse(is.na(aux$ThoraxChestAIS) == TRUE, NA,
                          ifelse(aux$ThoraxChestAIS >= 3, "present", "absent"))
aux$MI_AbdomenPelvis <- ifelse(is.na(aux$AbdomenPelvicContentsAIS) == TRUE, NA,
                          ifelse(aux$AbdomenPelvicContentsAIS >= 3, "present", "absent"))
aux$MI_LumbarSpine <- ifelse(is.na(aux$LumbarSpineAIS) == TRUE, NA,
                          ifelse(aux$LumbarSpineAIS >= 3, "present", "absent"))
aux$MI_UpperEx <- ifelse(is.na(aux$UpperExtremitiesAIS) == TRUE, NA,
                          ifelse(aux$UpperExtremitiesAIS >= 3, "present", "absent"))
aux$MI_LowerEx <- ifelse(is.na(aux$LowerExtremitiesAIS) == TRUE, NA,
                          ifelse(aux$LowerExtremitiesAIS >= 3, "present", "absent"))
aux$MI_Externa <- ifelse(is.na(aux$ExternaAIS) == TRUE, NA,
                          ifelse(aux$ExternaAIS >= 3, "present", "absent"))

#quantify the number of major EI
mycols <- c("PelvicGirdleAIS", "UpperExtremitiesAIS", "ThoracicSpineAIS", "AbdomenPelvicContentsAIS", "ExternaAIS", "ThoracicSpineAIS", "AbdomenPelvicContentsAIS", "ExternaAIS", "ThoraxChestAIS", "LumbarSpineAIS", "LowerExtremitiesAIS")
aux$MajorEI_count <- rowSums(aux[, mycols] > 2)
aux$MajorEI_count <- as.numeric(aux$MajorEI_count)
aux <- 
  aux %>%
  select(-c(contains("AIS")))
aux <-
  aux %>%
  select(Master_subject_ID = subjectId, everything())
aux <-
  aux %>%
  rename(GOSE_12months = GOSE12monthEndpointDerived,
         GOSE_6months = GOSE6monthEndpointDerived)
aux$GOSE_12months <- factor(aux$GOSE_12months, ordered = TRUE)
aux$GOSE_6months <- factor(aux$GOSE_6months, ordered = TRUE)
aux$Outcome_alive <- ifelse(is.na(aux$GOSE_12months) == TRUE, NA,
                      ifelse(aux$GOSE_12months %in% c("1"), "0", "1"))
aux$Outcome_5plus <- ifelse(is.na(aux$GOSE_12months) == TRUE, NA,
                      ifelse(aux$GOSE_12months %in% c("1", "2_or_3", "4"), "0", "1"))
aux$Outcome_7plus <- ifelse(is.na(aux$GOSE_12months) == TRUE, NA,
                      ifelse(aux$GOSE_12months %in% c("1", "2_or_3", "4", "5", "6"), "0", "1"))
aux$Outcome_8 <- ifelse(is.na(aux$GOSE_12months) == TRUE, NA,
                      ifelse(aux$GOSE_12months %in% c("1", "2_or_3", "4", "5", "6", "7"), "0", "1"))

#recode outcome for presentation in table1
aux$FavUnfav <- factor(aux$Outcome_5plus,
                       levels = c(0, 1),
                       labels = c("unfavourable", "favourable"))
aux$Recovery <- factor(aux$Outcome_8,
                       levels = c(0, 1),
                       labels = c("incomplete", "complete"))

aux$InjViolenceVictimAlcohol <- factor(aux$InjViolenceVictimAlcohol, ordered = TRUE,
                                       levels = c(0, 2, 1),
                                       labels = c("absent", "suspected", "present"))
aux$Hospital.ICUAdmisStatusIntubated <- factor(aux$Hospital.ICUAdmisStatusIntubated,
                                               levels = c(0, 1),
                                               labels = c("self-ventilating", "intubated"))
aux$PatientType <- factor(aux$PatientType,
                          levels = c(1, 2, 3),
                          labels = c("ED", "ward", "ICU"))

aux <- 
  aux %>%
  rename(Alcohol_intoxication = InjViolenceVictimAlcohol,
         ICU_admission_status = Hospital.ICUAdmisStatusIntubated)

aux$Admission_status <- ifelse(aux$PatientType == "ward", "ward",
                          ifelse(aux$PatientType == "ICU" & 
                                   is.na(aux$ICU_admission_status) == FALSE & 
                                   aux$ICU_admission_status == "self-ventilating", "ICU (self-ventilating)",
                                 ifelse(aux$PatientType == "ICU" & 
                                          is.na(aux$ICU_admission_status) == FALSE & 
                                          aux$ICU_admission_status == "intubated", "ICU (intubated)",
                                        ifelse(aux$PatientType == "ICU" & 
                                                 is.na(aux$ICU_admission_status) == TRUE, "ICU (airway unknown)",
                                               ifelse(aux$PatientType == "ED", "not admitted", NA)))))
aux$Admission_status <- factor(aux$Admission_status,
                               levels = c("not admitted",
                                          "ward",
                                          "ICU (self-ventilating)",
                                          "ICU (intubated)",
                                          "ICU (airway unknown)"))

aux <- 
  aux %>%
  mutate(across(c(contains("Outcome")),
                  as.numeric))
data <- 
  merge(data, aux, by = "Master_subject_ID", all.x = TRUE, all.y = TRUE)

#add pupil data
pup <- read.csv("../Data/CENTER_ImpactNoLab_data.csv", na.strings = c("88", " ", "", "NA"))
colnames(pup) <- c("Master_subject_ID", "Motor_score", "Unreactive_pupils", "Hypoxia", "Hypotension")
pup$Unreactive_pupils <- factor(pup$Unreactive_pupils,
                                ordered = TRUE,
                                levels = c(0, 1, 2))
pup$Hypoxia <- factor(pup$Hypoxia, ordered = TRUE,
                                       levels = c(0, 2, 1),
                                       labels = c("absent", "suspected", "present"))
pup$Hypotension <- factor(pup$Hypotension, ordered = TRUE,
                                       levels = c(0, 2, 1),
                                       labels = c("absent", "suspected", "present"))
data <- merge(data, pup, by = "Master_subject_ID", all.x = TRUE, all.y = TRUE)

## --------------------------------------------------------------------------------------------------------------------------------
#Add a column to say if patient has quantative MRI data
jhu_fa <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/Trajectory_paper/Raw_data/JHU_updated/fa_jhu_all_20220823.csv", na.strings = c("", " ", "NA")) %>% select(-X)
data$DTI_MR1 <- ifelse(data$Scan_ID_MR1 %in% jhu_fa$Scan_ID, "yes", "no")
data$DTI_MR2 <- ifelse(data$Scan_ID_MR2 %in% jhu_fa$Scan_ID, "yes", "no")


#Add column to say if MR1 has TAI as per
#central reporting
reports <- read.csv("C:/Users/sophie/Documents/Documents/Medicine/PhD/CTnegMild/Data/MRI_reports.csv", na.strings = c("", " ", "NA"))
reports <- 
  reports %>%
  select(Imaging.ExperimentId,
         Imaging.AnyIntracranTraumaticAbnormality, 
         Imaging.TAI) %>%
  group_by(Imaging.ExperimentId) %>%
  slice(1)
reports$Imaging.ExperimentId <- gsub("CTBI_", "", reports$Imaging.ExperimentId)
data <- 
  merge(data, reports, by.x = "Scan_ID_MR1", by.y = "Imaging.ExperimentId", all.x = TRUE, all.y = FALSE)

#rename variable
data <-
  data %>%
  rename(Cohort = Latest_GCS_cat)

#get CT reports
ct <- read.csv("../Data/CENTER_CT_reports.csv", na.strings = c("", " ", "uninterpretable"))

ct$CT_TAI <- ifelse(is.na(ct$Imaging.TAI)==TRUE, NA,
                 ifelse(ct$Imaging.TAI == "present", "present", "absent"))
ct$CT_Cistern <- ifelse(is.na(ct$Imaging.CisternalCompression)==TRUE, NA,
                 ifelse(ct$Imaging.CisternalCompression == "present", "present", "absent"))
ct$CT_SAH <- ifelse(is.na(ct$Imaging.TraumaticSubarachnoidHemorrhage)==TRUE, NA,
                 ifelse(ct$Imaging.TraumaticSubarachnoidHemorrhage == "present", "present", "absent"))
ct$CT_Mids <- ifelse(is.na(ct$Imaging.MidlineShift)==TRUE, NA,
                 ifelse(ct$Imaging.MidlineShift == "present", "present", "absent"))
ct$CT_EDH <- ifelse(is.na(ct$Imaging.EpiduralHematoma)==TRUE, NA,
                 ifelse(ct$Imaging.EpiduralHematoma == "present", "present", "absent"))
ct$CT_SDH <- ifelse(is.na(ct$Imaging.SubduralCollectionMixedDensity) == TRUE &
                   is.na(ct$Imaging.SubduralHematomaAcute)==TRUE &
                   is.na(ct$Imaging.SubduralHematomaSubacuteChronic) == TRUE, NA,
                 ifelse(ct$Imaging.SubduralCollectionMixedDensity == "present" |
                          ct$Imaging.SubduralHematomaAcute == "present" |
                          ct$Imaging.SubduralHematomaSubacuteChronic == "present", "present", "absent"))
ct$CT_Haematoma <- ifelse(is.na(ct$CT_EDH) == TRUE &
                         is.na(ct$CT_SDH) == TRUE, NA,
                       ifelse(ct$CT_EDH == "present" | ct$CT_SDH == "present", "present", "absent"))
cols <- ct  %>% select(contains("CT_")) %>% colnames()
ct[cols] <- lapply(ct[cols], factor) 
ct <- 
  ct %>%
  filter(Imaging.Timepoint == "CT Early") %>%
  select(Master_subject_ID = subjectId,
         contains("CT_")) %>%
  group_by(Master_subject_ID) %>%
  slice(1)
data <- 
  merge(data, ct, by = "Master_subject_ID", all.x = TRUE, all.y = FALSE)


#IMPACT lab variables
lab <- read.csv("../Data/CENTER_ImpactLab.csv", na.strings = c("NA", "", " "))
lab$Date_time <- as.POSIXct(paste(lab$Labs.DLDate, lab$Labs.DLTime), format="%Y-%m-%d %H:%M:%S")
lab <-
  lab %>%
  filter(is.na(Labs.DLGlucosemmolL)==FALSE | is.na(Labs.DLHemoglobingdL) == FALSE) %>%
  group_by(subjectId) %>%
  slice(which.min(Date_time)) %>%
  select(Master_subject_ID = subjectId,
         Glucose = Labs.DLGlucosemmolL,
         Hb = Labs.DLHemoglobingdL)
data <- 
  merge(data, lab, by = "Master_subject_ID", all.x = TRUE, all.y = FALSE)

#Record whether patient has an MRI used in the subgroup analysis i.e. an MRI compatible with Marshall 1 or 2
#Add TAI grade where available
#and record if TAi cannot be reported
adam1 <- read.csv("../Data/Grade_these_scans_filled_vfjn.csv", na.strings = c("", "NA", " "))
adam1 <-
  adam1 %>%
  select(SiteCode = Site, 
         Master_subject_ID, 
         Scan_ID_MR1 = Scan_ID, 
         Comment:BS_pons_right)
adam2 <- read.csv("../Data/More_scans_to_grade_filled.csv", na.strings = c("", "NA", " "))
adam2 <-
  adam2 %>%
  select(SiteCode,
         Master_subject_ID,
         Scan_ID_MR1,
         Comment:BS_pons_right)
adam <-
  bind_rows(adam1, adam2)

temp <- 
  adam %>%
  filter(Any_TAI %in% c("yes", "no"))
data$MRI_used <- ifelse(data$Master_subject_ID %in% temp$Master_subject_ID, "yes","no")


## --------------------------------------------------------------------------------------------------------------------------------
#Add information on intracranial surgery
surg <- read.csv("../Data/IntraCranSurgery.csv", na.strings = c("", " ", "88", NA))
surg <-
  surg %>%
  select(Master_subject_ID = subjectId,
         IC_surgery_any = Surgeries.CranialSurgDone,
         IC_surgery_emergency = InjuryHx.EmergSurgInterventionsIntraCran,
         IC_surgery_decomcran = Surgeries.DecompressiveCran)

surg$IC_surgery_any <- ifelse(is.na(surg$IC_surgery_any)==TRUE, NA, 
                              ifelse(surg$IC_surgery_any == "1", "performed", "not performed"))
surg$IC_surgery_decomcran <- ifelse(is.na(surg$IC_surgery_decomcran)==TRUE, NA, 
                              ifelse(surg$IC_surgery_decomcran == "1", "performed", "not performed"))
surg$IC_surgery_emergency <- ifelse(is.na(surg$IC_surgery_emergency)==TRUE, NA, 
                              ifelse(surg$IC_surgery_emergency == "1", "performed", "not performed"))

data <- 
  merge(data, surg, by = "Master_subject_ID", all.x = TRUE, all.y = FALSE)


#Add information on preinjury ASA and early seizure
asa <- read.csv("../Data/MedicalHistory.csv", na.strings = c("", " ", "88"))
asa <-
  asa %>%
  select(Master_subject_ID = subjectId,
         ASA = MedHx.MedHxPreInjASAPSClass,
         Seizure = InjuryHx.EDComplEventSeizures,
         Hx_neuro = MedHx.MedHxNeuro,
         Hx_psych = MedHx.MedHxPsychiatric,
         Hx_renal = MedHx.MedHxRenal,
         Hx_hepa = MedHx.MedHxHepatic)
asa$ASA <- factor(asa$ASA)
asa$Hx_hepa <- factor(asa$Hx_hepa,
                      levels = c(0, 1),
                      labels = c("absent", "present"))
asa$Hx_neuro <- factor(asa$Hx_neur,
                      levels = c(0, 1),
                      labels = c("absent", "present"))
asa$Hx_renal <- factor(as.character(asa$Hx_renal),
                      levels = c(0, 1),
                      labels = c("absent", "present"))
asa$Hx_psych <- factor(as.character(asa$Hx_psych),
                      levels = c(0, 1),
                      labels = c("absent", "present"))
asa$Seizure <- factor(asa$Seizure,
                      levels = c(0, 1, 2, 3),
                      labels = c("none", "partial", "generalised", "status epilepticus"))

data <- 
  merge(data, asa, by = "Master_subject_ID", all.x = TRUE, all.y = FALSE)

#Add information about when patient was intubated

intub <- read.csv("../Data/Intubation.csv", na.strings = c("", " ", "88"))
intub <-
  intub %>%
  select(Master_subject_ID = subjectId,
         Intub_prehosp = InjuryHx.EDArrivalAirway,
         Intub_onscene = InjuryHx.PresEmergencyCareIntubation)
intub$Intub_prehosp <- ifelse(is.na(intub$Intub_prehosp) == TRUE, NA,
                              ifelse(intub$Intub_prehosp == "4", "yes", "no"))
intub$Intub_prehosp <- factor(intub$Intub_prehosp)
intub$Intub_onscene <- ifelse(is.na(intub$Intub_onscene) == TRUE, NA,
                              ifelse(intub$Intub_onscene == "1", "yes", "no"))
intub$Intub_onscene <- factor(intub$Intub_onscene)

data <- 
  merge(data, intub, by = "Master_subject_ID", all.x = TRUE, all.y = FALSE)




## --------------------------------------------------------------------------------------------------------------------------------

#Generate labels that can later be used for Table 1
label(data$Admission_status) <-  "Admission status" 
label(data$MajorEI) <-  "Major extra-cranial injury" 
label(data$Marshall_score) <-  "Marshall CT score"
label(data$Days_to_MR1) <-  "Time to MRI"
units(data$Days_to_MR1) <- "days"
label(data$Time_to_bio) <-  "Time to blood protein sample"
units(data$Time_to_bio) <- "hours"
label(data$Unreactive_pupils) <-  "Unreactive pupils"
label(data$Alcohol_intoxication) <-  "Alcohol intoxication"
label(data$Outcome_5plus) <-  "Outcome"
label(data$Latest_GCS) <- "Glasgow Coma Scale"
label(data$Motor_score) <- "Motor score"
units(data$Glucose) <- "mmol/L"
label(data$Hb) <- "Haemoglobin"
units(data$Hb) <- "g/dL"
label(data$CT_TAI) <- "CT - TAI"
label(data$CT_Cistern) <- "CT - cisternal compression"
label(data$CT_Mids) <- "CT - midline shift"
label(data$CT_EDH) <- "CT - EDH"
label(data$CT_SAH) <- "CT - SAH"
label(data$CT_Haematoma) <- "CT - haematoma"
label(data$GOSE_6months) <- "GOSE at 6 months"
label(data$MajorEI_count) <- "Number of major extracranial injuries"
label(data$IC_surgery_any) <- "Any neurosurgery"
label(data$IC_surgery_decomcran) <- "Decompressive craniectomy"
label(data$IC_surgery_emergency) <- "Emergency neurosurgery"
label(data$Seizure) <- "Early seizure"
label(data$ASA) <- "Pre-injury ASA grade"
label(data$MI_AbdomenPelvis) <- "Major abdomen/pelvis injury"
label(data$MI_Chest) <- "Major chest injury"
label(data$MI_Externa) <- "Major external injury"
label(data$MI_Face) <- "Major face injury"
label(data$MI_HeadNeck) <- "Major head/neck injury"
label(data$MI_LowerEx) <- "Major lower extremity injury"
label(data$MI_LumbarSpine) <- "Major lumbar spine injury"
label(data$MI_ThoracicSpine) <- "Major thoracic spine injury"
label(data$MI_UpperEx) <- "Major upper extremity injury"
label(data$Hx_hepa) <- "Prior hepatic disease"
label(data$Hx_renal) <- "Prior renal disease"
label(data$Hx_neuro) <- "Prior neurological disease"
label(data$Hx_psych) <- "Prior psychiatric disease"
label(data$Intub_prehosp) <- "Intubated before study hospital"
label(data$Intub_onscene) <- "Intubated on scene"


#########start weeding out patients
# save overall CENTER_TBI population for comparison
center <- data
print(paste0("Number of patients in CENTER-TBI: ", n_distinct(data$Master_subject_ID)))

num <- 
  data %>%
  filter(Age <16) %>%
  nrow()
print(paste0("Of those, patients who were aged less than 16: ", num))
data <-
  data %>%
  filter(Age >= 16)
print(paste0("Number of CENTER-TBI patients aged >=16: ", n_distinct(data$Master_subject_ID)))

m <- 
  data %>%
  filter(Cohort == "mild") %>%
  nrow()
print(paste0("Of those, patients with mild TBI: ", m))
data <- 
  data %>% 
  filter(Cohort != "mild")
print(paste0("Number of adults in CENTER-TBI with moderate/severe/untestable GCS: ", n_distinct(data$Master_subject_ID)))

#only keep untestable if they were intubated (else GCS may simply not have been recorded)
not_icu <- 
  data %>%
  filter(Cohort ==  "untestable" & PatientType != "ICU") %>%
  nrow()
icu_self <- 
  data %>%
  filter(Cohort ==  "untestable" & PatientType == "ICU" & ICU_admission_status == "self-ventilating") %>%
  nrow()
icu_unknown <- 
  data %>%
  filter(Cohort ==  "untestable" & PatientType == "ICU" & is.na(ICU_admission_status)==TRUE) %>%
  nrow()
print(paste0("Number with unknown GCS who were on the ward/ on ICU but self-ventilating/on ICU but intubation status unclear: ", 
             not_icu, "/",
             icu_self, "/",
             icu_unknown))
data <- 
  data %>%
  filter(Cohort == "mod-sev" | (Cohort == "untestable" & ICU_admission_status == "intubated"))

#at this point I can label untestable GCS untestable (same for motor score)
data$Latest_GCS <- ifelse(is.na(data$Latest_GCS)==TRUE, "untestable", as.character(data$Latest_GCS))
data$Latest_GCS <- factor(data$Latest_GCS,
                          ordered = TRUE,
                          levels = c("untestable", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
data$Motor_score <- ifelse(data$Cohort== "untestable", "untestable", as.character(data$Motor_score))
data$Motor_score <- factor(data$Motor_score,
                          ordered = TRUE,
                          levels = c("untestable","1", "2", "3", "4", "5", "6"))
print(paste0("Patients with moderate-severe TBI or GCS untestable because of intubation: ", nrow(data)))

m3 <- summary(as.factor(data$Marshall_score))[3]
m4 <- summary(as.factor(data$Marshall_score))[4]
m5 <- summary(as.factor(data$Marshall_score))[5]
m6 <- summary(as.factor(data$Marshall_score))[6]
msum <- m3 + m4 + m5 +m6

print(paste0("Number of those patients with Marshall 3, 4, 5 or 6: ", msum))
data <- 
  data %>%
  filter(Marshall_score %in% c("1", "2", "unknown"))
print(paste0("Number of non-mild patients with Marshall 1, 2, or unknown: ", n_distinct(data$Master_subject_ID)))

temp <- 
  data %>%
  filter(Marshall_score == "unknown" & MRI_used == "no" & Master_subject_ID != "7bZC357")
print(paste0("Of those patients with unknown Marshall without an MRI to prove Marshall 1/2 equivalence: ", nrow(temp)))

data <- 
  data %>%
  filter(!Master_subject_ID %in% temp$Master_subject_ID)
print(paste0("Patients with low/untestable GCS, 24h biomarker and Marshall 1/2/equivalent: ", nrow(data)))


data$Marshall_score <- factor(data$Marshall_score,
                                       levels = c("1", "2", "unknown"),
                                       labels = c("1", "2", "MRI equivalent to 1 or 2"))

num <- 
  data %>%
  filter(Time_to_bio > 24) %>%
  nrow()
print(paste0("Of those, patients who did not have serum biomarker samples within 24h: ", num))
data <- 
  data %>%
  filter(Time_to_bio <= 24)
print(paste0("Number of non-mild patients with unimpressive CT, aged >=16 and serum biomarkers within 24h: ", n_distinct(data$Master_subject_ID)))

print(paste0("Of those, patients who never had an MRI: ", summary(data$Days_to_MR1)[7]))
data <-
  data %>%
  filter(is.na(Days_to_MR1)==FALSE)
print(paste0("Number of non-mild patients with unimpressive CT, aged >=16, 24h biomarker and at least one MR: ", n_distinct(data$Master_subject_ID)))
m <- 
  data %>%
  filter(Days_to_MR1 > 30) %>%
  nrow()
print(paste0("Number of patients with their first MRI more than 30 days since injury: ", m))
data <-
  data %>%
  filter(Days_to_MR1 <= 30)
print(paste0("Number of non-mild patients with unimpressive CT, aged >= 16, 24h biomarker and MRI within 30 days: ", n_distinct(data$Master_subject_ID)))

#Remove patients later identified as inelegible after
#inspection of scans

ex1 <- c("7bZC357", #amyloidosis
        "7gPQ888", "8yyG376", "7qrW488") #Marshall > 2 excluded earlier in this flowchart
ex2 <- c("2Exj679" , #no FLAIR/SWI
        "2sCG537" , #motion artefat
        "2UBy682" , "4xHA258" ,"7LSK389") #no FLAIR/SWI
ex3 <- c("7bZC357", #amyloidosis
         "2Exj679" , #no FLAIR/SWI
        "2sCG537" , #motion artefat
        "2UBy682" , "4xHA258" ,"7LSK389") #no FLAIR/SWI
ex4 <- c("8yyG376", "7qrW488") #CT unknown and MRI Marshall > 2 excluded earlier in this flowchart
ex5 <- "7gPQ888" # Marshall reported as 2 but MRI much worse than a Marshall 2 equivalent
print(paste0("Number of elegible patients excluded after visual inspection of MRI: ", length(ex3)+length(ex4)+length(ex5)))
data <- 
  data %>% 
  filter(!Master_subject_ID %in% c(ex3, ex5))
print(paste0("Total number of elegible patients (both cohorts combined): ", n_distinct(data$Master_subject_ID)))




## --------------------------------------------------------------------------------------------------------------------------------
#Make 2 versions of Table 1
# version 1 comparing mod-sev and untestable cohort
# version 2 comparing this study population with the overall CENTER-TBI cohort
#Summary Table
data$Admission_status <- droplevels(data$Admission_status)
label(data$Admission_status) <-  "Admission status"

table1 <- table1(~ Age + 
                   Sex + 
                   ASA +
                   Hx_neuro +
                   Hx_psych +
                   Hx_renal +
                   Hx_hepa +
                   Admission_status + 
                   MajorEI +
                   MajorEI_count +
                   MI_HeadNeck +
                   MI_Face +
                   MI_Chest +
                   MI_ThoracicSpine +
                   MI_LumbarSpine +
                   MI_AbdomenPelvis +
                   MI_UpperEx +
                   MI_LowerEx +
                   MI_Externa +
                   Latest_GCS +
                   Unreactive_pupils +
                   Seizure +
                   Hypoxia +
                   Hypotension +
                   Alcohol_intoxication +
                   Time_to_bio +
                   Marshall_score +
                   CT_TAI +
                   CT_Mids +
                   CT_Cistern +
                   CT_SAH +
                   CT_EDH +
                   CT_Haematoma +
                   IC_surgery_any +
                   IC_surgery_emergency +
                   IC_surgery_decomcran +
                   GOSE_6months | Cohort,
         data = data, 
         droplevels=F,
         render.continuous = my.render.cont, 
         render.categorical = my.render.cat)
table1
ft1 <- t1flex(table1)
save_as_docx(ft1, path = "Figures/Table1.docx")




## --------------------------------------------------------------------------------------------------------------------------------
center$Study <- ifelse(center$Master_subject_ID %in% data$Master_subject_ID, "included", "not included")
table1b <- table1(~ Age + 
                   Sex + 
                    ASA +
                   Admission_status + 
                   MajorEI +
                   Latest_GCS +
                   Unreactive_pupils +
                    Seizure +
                   Hypoxia +
                   Hypotension +
                   Alcohol_intoxication +
                   Time_to_bio +
                   Marshall_score +
                   CT_TAI +
                   CT_Mids +
                   CT_Cistern +
                   CT_SAH +
                   CT_EDH +
                   CT_Haematoma +
                   IC_surgery_any +
                   IC_surgery_emergency +
                   IC_surgery_decomcran +
                   GOSE_6months | Study,
         data = center, 
         droplevels=F,
         render.continuous = my.render.cont, 
         render.categorical = my.render.cat)
table1b
ft1b <- t1flex(table1b)
save_as_docx(ft1b, path = "Figures/Table1b.docx")


## --------------------------------------------------------------------------------------------------------------------------------
#Add TAI grade where available
#and record if TAI cannot be reported
adam1 <- read.csv("../Data/Grade_these_scans_filled_vfjn.csv", na.strings = c("", "NA", " "))
adam1 <-
  adam1 %>%
  select(SiteCode = Site, 
         Master_subject_ID, 
         Scan_ID_MR1 = Scan_ID, 
         Comment:BS_pons_right,
         BS_dorsal,
         BS_DuretCont)
adam2 <- read.csv("../Data/More_scans_to_grade_filled.csv", na.strings = c("", "NA", " "))
adam2 <-
  adam2 %>%
  select(SiteCode,
         Master_subject_ID,
         Scan_ID_MR1,
         Comment:BS_pons_right,
         BS_dorsal,
         BS_DuretCont)
adam <-
  bind_rows(adam1, adam2)

#Adam-Gentry classification
adam$AG_stage <- ifelse(is.na(adam$Any_TAI)==TRUE | adam$Any_TAI == "exclude", NA,
                        ifelse((adam$BS_midbrain_left == "yes" | 
                                               adam$BS_midbrain_right == "yes" |
                                               adam$BS_pons_left == "yes" |
                                               adam$BS_pons_right == "yes"), 3,
                               ifelse((adam$CC_body == "yes" | 
                                        adam$CC_genu == "yes" |
                                        adam$CC_splenium == "yes"), 2,
                                      ifelse(adam$Hemispheres == "yes", 1, 0))))
#Firsching classification
adam$F_grade <- ifelse(is.na(adam$Any_TAI)==TRUE | adam$Any_TAI == "exclude", NA,
                       ifelse(adam$BS_pons_left == "yes" & adam$BS_pons_right == "yes", 4,
                              ifelse(adam$BS_midbrain_left == "yes" & adam$BS_midbrain_right == "yes", 3,
                                     ifelse(adam$BS_pons_left == "yes" | 
                                              adam$BS_midbrain_left == "yes" |
                                              adam$BS_pons_right == "yes" | 
                                              adam$BS_midbrain_right == "yes", 2,
                                            ifelse(adam$Hemispheres == "yes" |
                                                     adam$CC_body == "yes" |
                                                     adam$CC_genu == "yes" |
                                                     adam$CC_splenium == "yes", 1, 0)))))
#Burden
adam$H_any <- ifelse(is.na(adam$Any_TAI)==TRUE | adam$Any_TAI == "exclude", NA,
                           ifelse(adam$Hemispheres == "yes", 1, 0))
adam$CC_any <- ifelse(is.na(adam$Any_TAI)==TRUE | adam$Any_TAI == "exclude", NA,
                      ifelse((adam$CC_body == "yes" | 
                                        adam$CC_genu == "yes" |
                                        adam$CC_splenium == "yes"), 1, 0))
adam$BS_any <- ifelse(is.na(adam$Any_TAI)==TRUE | adam$Any_TAI == "exclude", NA,
                      ifelse((adam$BS_midbrain_left == "yes" | 
                                adam$BS_midbrain_right == "yes" |
                                adam$BS_pons_left == "yes" |
                                adam$BS_pons_right == "yes"), 1, 0))
adam$Burden_visual <- ifelse(is.na(adam$Any_TAI)==TRUE | adam$Any_TAI == "exclude", NA,
                             rowSums(adam[, c("H_any", "CC_any", "BS_any")]))




## --------------------------------------------------------------------------------------------------------------------------------

#Pontine involvement
adam$BSI_pontine <- ifelse(is.na(adam$Any_TAI)==TRUE | adam$Any_TAI == "exclude", NA,
                           ifelse(adam$AG_stage %in% c("0", "1", "2"), "grade <3",
                       ifelse(adam$BS_pons_left == "yes" | adam$BS_pons_right == "yes", "present", "absent")))
adam$BSI_pontine <- factor(adam$BSI_pontine,ordered = TRUE,
                           levels = c("grade <3", "absent", "present"))


#Bilateral BSI
adam$BSI_bilateral <- ifelse(is.na(adam$Any_TAI)==TRUE | adam$Any_TAI == "exclude", NA,
                             ifelse(adam$AG_stage %in% c("0", "1", "2"), "grade <3",
                       ifelse((adam$BS_pons_left == "yes" | adam$BS_midbrain_left == "yes") &
                                (adam$BS_pons_right == "yes" | adam$BS_midbrain_right == "yes"), "present", "absent")))
adam$BSI_bilateral <- factor(adam$BSI_bilateral,ordered = TRUE,
                           levels = c("grade <3", "absent", "present"))

#Dorsal BSI
adam$BSI_dorsal <- ifelse(is.na(adam$Any_TAI)==TRUE | adam$Any_TAI == "exclude", NA,
                          ifelse(adam$AG_stage %in% c("0", "1", "2"), "grade <3",
                       ifelse(adam$BS_dorsal == "yes", "present", "absent")))
adam$BSI_dorsal <- factor(adam$BSI_dorsal, ordered = TRUE,
                           levels = c("grade <3", "absent", "present"))

#BSI with Duret or contusion
adam$BSI_DuretCont <- ifelse(is.na(adam$Any_TAI)==TRUE | adam$Any_TAI == "exclude", NA,
                             ifelse(adam$AG_stage %in% c("0", "1", "2"), "grade <3",
                       ifelse(adam$BS_DuretCont == "yes", "present", "absent")))
adam$BSI_DuretCont <- factor(adam$BSI_DuretCont,ordered = TRUE,
                           levels = c("grade <3", "absent", "present"))
                             
#merge with data
graded <- adam %>% select(Master_subject_ID, AG_stage, Burden_visual, BSI_pontine, BSI_bilateral, BSI_dorsal, BSI_DuretCont)
data <- merge(data, graded, by = "Master_subject_ID", all.x = TRUE, all.y = FALSE)




## --------------------------------------------------------------------------------------------------------------------------------
##Is the brainstem affected?
data$BS_visual <- ifelse(is.na(data$AG_stage)==TRUE, NA,
                         ifelse(data$AG_stage == 3, "present", "absent"))

data$Isolated_BSI <- ifelse(data$BS_visual == "absent", "no BSI",
                            ifelse(data$BS_visual == "present" & as.numeric(data$Burden_visual < 2), "isolated BSI", "non-isolated BSI"))
data$Isolated_BSI <- factor(data$Isolated_BSI,
                            levels = c("no BSI", "isolated BSI", "non-isolated BSI"))



## --------------------------------------------------------------------------------------------------------------------------------
catVars <- c("AG_stage",
             "BS_visual",
             "BSI_bilateral",
             "BSI_dorsal",
             "BSI_pontine", 
             "BSI_DuretCont")

data[catVars] <- lapply(data[catVars], factor) 

label(data$Days_to_MR1) <- "Time to MRI"
units(data$Days_to_MR1) <- "days"
label(data$AG_stage) <- "Adams-Gentry stage"
label(data$BS_visual) <- "Brainstem injury - any type"
label(data$BSI_bilateral) <- "Brainstem injury - bilateral"
label(data$BSI_dorsal) <- "Brainstem injury - dorsal"
label(data$BSI_pontine) <- "Brainstem injury - pontine"
label(data$BSI_DuretCont) <- "Brainstem injury - contusion or Duret haemorrhage"

#Summary Table

table2 <- table1(~ Days_to_MR1 +
                   AG_stage +
                   BS_visual +
                   BSI_bilateral +
                   BSI_dorsal +
                   BSI_pontine +
                   BSI_DuretCont| Cohort,
         data = data, 
         droplevels=F,
         render.continuous = my.render.cont, 
         render.categorical = my.render.cat)
table2
ft2 <- t1flex(table2)
save_as_docx(ft2, path = "Figures/Table_staging.docx")

## --------------------------------------------------------------------------------------------------------------------------------

#grade of TAI
compare_bio_stages(data = data,
                   mycohort = "all",
                   myvar = "AG_stage",
                   mytitle = "Adam-Gentry stage")

compare_bio_stages(data = data,
                   mycohort = "mod-sev",
                   myvar = "AG_stage",
                   mytitle = "Adam-Gentry stage derived visually")

compare_bio_stages(data = data,
                   mycohort = "untestable",
                   myvar = "AG_stage",
                   mytitle = "Adam-Gentry stage derived from visually")


#presence of BSI

compare_bio_levels(data = data,
                   mycohort = "all",
                   myvar = "BS_visual",
                   mytitle = "Brainstem injury identified")

compare_bio_levels(data = data,
                   mycohort = "mod-sev",
                   myvar = "BS_visual",
                   mytitle = "Brainstem injury identified visually") 

compare_bio_levels(data = data,
                   mycohort = "untestable",
                   myvar = "BS_visual",
                   mytitle = "Brainstem injury identified visually") 




## --------------------------------------------------------------------------------------------------------------------------------
#make a table showing the association of biomarkers with certain BSI features
features <- c("BSI_bilateral", "BSI_dorsal", "BSI_pontine", "BSI_DuretCont")
res_list <- vector(mode = "list", length = length(features))

for (f in 1:length(features)){
  myvar <- features[f]

mydata <- 
    data %>% 
    filter(Cohort == "mod-sev") %>%
    select(Master_subject_ID, myvar, GFAP:UCH.L1) %>% 
    gather(key = "Protein", value = "Log_conc", GFAP:UCH.L1) 
  
  proteins <- c("GFAP", "NFL", "NSE", "S100B", "Tau",  "UCH.L1")
  mycols   <- c("Feature", "Protein", "grade <3", "absent","present", "P_raw")
  res <- as.data.frame(matrix(,
                              nrow = length(proteins),
                              ncol = length(mycols)))
  colnames(res) <- mycols


  for (p in 1:length(proteins)){
    df <- 
      mydata %>% 
      spread(key = myvar, value = "Log_conc") %>%
      filter(Protein == proteins[p])
    
    res[p, "Feature"] <- features[f]
    res[p, "Protein"] <- proteins[p]
    
    med <- median(df$`grade <3`, na.rm = TRUE) %>% exp() 
    med <- format(round(as.numeric(med), 2), nsmall = 2)
    Q1 <- quantile(exp(df$`grade <3`), na.rm = TRUE)[2] 
    Q1 <- format(round(Q1,2), nsmall = 2)
    Q3 <- quantile(exp(df$`grade <3`), na.rm = TRUE)[4]
    Q3 <- format(round(Q3,2), nsmall = 2)
    res[p, "grade <3"] <- paste0(med, " (", Q1, "-", Q3, ")")
    
    med <- median(df$absent, na.rm = TRUE) %>% exp() 
    med <- format(round(as.numeric(med), 2), nsmall = 2)
    Q1 <- quantile(exp(df$absent), na.rm = TRUE)[2] 
    Q1 <- format(round(Q1,2), nsmall = 2)
    Q3 <- quantile(exp(df$absent), na.rm = TRUE)[4]
    Q3 <- format(round(Q3,2), nsmall = 2)
    res[p, "absent"] <- paste0(med, " (", Q1, "-", Q3, ")")
    
    med <- median(df$present, na.rm = TRUE) %>% exp() 
    med <- format(round(as.numeric(med), 2), nsmall = 2)
    Q1 <- quantile(exp(df$present), na.rm = TRUE)[2] 
    Q1 <- format(round(Q1,2), nsmall = 2)
    Q3 <- quantile(exp(df$present), na.rm = TRUE)[4]
    Q3 <- format(round(Q3,2), nsmall = 2)
    res[p, "present"] <- paste0(med, " (", Q1, "-", Q3, ")")
      
      df <- 
        mydata %>% 
        filter(Protein == proteins[p])
      df[, myvar] <- factor(df[, myvar], ordered = TRUE)
      res[p, "P_raw"] <- JonckheereTerpstraTest(df[,"Log_conc"], df[,myvar])$p.value
    
    
  }
  
res_list[[f]] <- res

}

temp <- bind_rows(res_list) 
temp <- temp %>% select(Protein, Feature, everything()) %>% arrange(Protein)
temp <- temp %>% group_by(Protein) %>% mutate(P_adj = p.adjust(P_raw, method = "fdr"))
temp <- temp %>% select(-P_raw)
temp$P_adj <- format(round(temp$P_adj, 3), nsmall = 3)
temp$P_adj <- as.numeric(temp$P_adj)
temp$Feature <- gsub("BSI_", "", temp$Feature)
temp$Feature <- gsub("DuretCont", "Duret/cont.", temp$Feature)

ft_BSI <- 
  temp %>%
  regulartable() %>%
  set_header_labels("grade <3" = "No BSI",
                    "absent" = "BSI without this feature",
                    "present" = "BSI with this feature",
                    P_adj = "Adj. p-value") %>%
  theme_vanilla() %>%
    vline_left() %>%
    vline_right() %>%
    vline(j = c(2,5)) %>%
    merge_v(j = 1) %>%
  rotate(j = 1, align = "center", rotation = "btlr", part = "body") %>%
  bold(part = "header") %>%
  bold(i = ~ P_adj <0.05, 
  j = ~ P_adj, bold = TRUE) %>%
  merge_v(j = 1) %>%
    align(j = c(3,4,5,6), align = "right") %>%
    padding(padding = 1, padding.top = 3) %>%
    fontsize(size = 10, part = "all") %>%
    valign(valign = "top", part = "body", j = c(2, 3, 4, 5, 6))

#add horizontal line after merged cell
  row_loc <- rle(cumsum(ft_BSI$body$spans$columns[,1] ))$values
  ft_BSI <- 
    ft_BSI %>% 
    border(border.bottom = fp_border_default(),
           i=row_loc, 
           j = 1:6, 
           part="body") 
  ft_BSI <- 
    ft_BSI %>% 
    border(border.bottom = fp_border_default(),
           i = ft_BSI$body$spans$columns[,1] > 1, 
           j = 1, 
           part="body") %>% 
    border(border.bottom = fp_border_default(), 
           border.top = fp_border_default(),
           part = "header") 

  #Bold formatting
  ft_BSI <-
    ft_BSI    %>%
  fix_border_issues()
  
ft_BSI
save_as_docx(ft_BSI, path = "Figures/BSI_modsev.docx")


## --------------------------------------------------------------------------------------------------------------------------------
# See if proteins can act as screening tools for BSI
plot_roc_curves(data, "mod-sev", "BS_visual")


## --------------------------------------------------------------------------------------------------------------------------------
#Check if there is a trend between biomarkers and visual TAI burden
mydata <- 
    data %>% 
    select(Master_subject_ID, Cohort, Burden_visual, GFAP:UCH.L1) %>% 
    gather(key = "Protein", value = "Log_conc", GFAP:UCH.L1)
mydata$Burden_visual <- factor(mydata$Burden_visual, ordered = TRUE)

proteins <- c("GFAP", "NFL", "NSE", "S100B", "Tau", "UCH.L1")
mycohorts <- c("mod-sev", "untestable")

res <- as.data.frame(matrix(,
                            nrow = length(proteins)*length(mycohorts),
                            ncol = 3))
colnames(res) <- c("Protein", "Cohort", "Raw_P")

for (c in 1:length(mycohorts)){
  for (p in 1:length(proteins)){
    r <- (c-1)*length(proteins) + p
    res[r, "Protein"] <- proteins[p]
    res[r, "Cohort"] <- mycohorts[c]
    j <- JonckheereTerpstraTest(Log_conc ~ Burden_visual, 
                              data = mydata %>% filter(Cohort == mycohorts[c] & Protein == proteins[p]))
    res[r, "Raw_P"] <- j$p.value
    }
}
res$Adj_P <- p.adjust(res$Raw_P, method = "fdr")
res$Adj_P <- format(round(res$Adj_P, 3), nsmall = 3)
res$Sig <- ifelse(res$Adj_P < 0.05, "*", " ")
res$label <- paste0("adj. p-value = ", res$Adj_P, res$Sig)
res <- merge(res, mydata, by = c("Cohort", "Protein"))

# New facet label names for protein variable
protein.labs <- c("GFAP", "NFL", "NSE", "S100B", "Tau", "UCH-L1")
names(protein.labs) <- c("GFAP", "NFL", "NSE", "S100B", "Tau", "UCH.L1")


p3 <- ggplot(data = mydata,
       aes(x = Burden_visual, y = Log_conc,  color = Burden_visual)) +
  geom_boxplot(color = "black", outlier.shape = NA) +
  geom_jitter(width=0.15, alpha=0.3, color = "firebrick3")+
  facet_grid(Protein ~ Cohort, scales = "free", labeller = labeller(Protein = protein.labs)) + 
  theme_bw()+
  theme(legend.position = "none",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  xlab("Areas with AI") +
  ylab("")+
   geom_text(
    size    = 2.5,
    color = "black",
    data    = res,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -0.61,
    vjust   = -0.3
  )+ 
 theme(strip.text.x = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.1)))
p3


## --------------------------------------------------------------------------------------------------------------------------------
#Check if there is a trend between protein and visual burden FOR BOTH COHORTS COMBINED
mydata <- 
    data %>% 
    select(Master_subject_ID, Burden_visual, GFAP:UCH.L1) %>% 
    gather(key = "Protein", value = "Log_conc", GFAP:UCH.L1)
mydata$Burden_visual <- factor(mydata$Burden_visual, ordered = TRUE)

proteins <- c("GFAP", "NFL", "NSE", "S100B", "Tau", "UCH.L1")

res <- as.data.frame(matrix(,
                            nrow = length(proteins),
                            ncol = 2))
colnames(res) <- c("Protein",  "Raw_P")


  for (p in 1:length(proteins)){
    res[p, "Protein"] <- proteins[p]
    j <- JonckheereTerpstraTest(Log_conc ~ Burden_visual, 
                              data = mydata %>% filter(Protein == proteins[p]))
    res[p, "Raw_P"] <- j$p.value
    }

res$Adj_P <- p.adjust(res$Raw_P, method = "fdr")
res$Adj_P <- format(round(res$Adj_P, 3), nsmall = 3)
res$Sig <- ifelse(res$Adj_P < 0.05, "*", " ")
res$label <- paste0("adj. p-value = ", res$Adj_P, res$Sig)
res <- merge(res, mydata, by = "Protein")


p3b <- ggplot(data = mydata,
       aes(x = Burden_visual, y = Log_conc,  color = Burden_visual)) +
  geom_boxplot(color = "black", outlier.shape = NA) +
  geom_jitter(width=0.15, alpha=0.3, color = "firebrick3")+
  facet_grid(vars(Protein), scales = "free") + 
  theme_bw()+
  theme(legend.position = "none") +
  ylab("Log serum protein concentration") +
  xlab("Areas with AI")+
   geom_text(
    size    = 2.5,
    color = "black",
    data    = res,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -0.82,
    vjust   = -0.3
  )+ 
 theme(strip.text.y = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.1)))


t1 <- textGrob("    Overall cohort")
t2 <- textGrob("            Moderate-severe")
t4 <- textGrob("Unrecorded")
lay <- rbind(c(1, 1, 1,  2, 2, 4, 4, 4),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6),
             c(5, 5, 5,  6, 6, 6, 6, 6))
mygrid <- grid.arrange(t1, t2, t4, p3b, p3, layout_matrix = lay)
ggsave(plot = mygrid, 
      filename = "Figures/plot_TAI_burden.tiff",
      width = 15, height = 15, dpi = 300,
      units = "cm")



## --------------------------------------------------------------------------------------------------------------------------------

cutdata <- 
  data %>%
  filter(Cohort == "mod-sev")%>%
  mutate(across(c("GFAP", "NFL", "NSE", "S100B", "Tau", "UCH.L1"), 
         exp, 
         .names = "{.col}"))
cutdata$BS_visual <- ifelse(cutdata$BS_visual == "present", 1, 0)



## --------------------------------------------------------------------------------------------------------------------------------
##########Calculate savings per patient#########

#Collect Raw costs
p_cost_uk   <- c(12.82,  18.45,  9.68,     17.56,     17.06,   12.03)
mri_cost_uk <- 385.8
phleb_uk    <- 4

p_cost_usa   <- c(50,  0,  0,     0,     0,   50)
mri_cost_usa <- 2758.56
phleb_usa     <- 15.16

#Set up results table
proteins <- c("GFAP", "NFL", "NSE", "S100B", "Tau", "UCH.L1")
resnames <- c("Protein", 
              "Threshold", 
              "Sensitivity", 
              "Specificity", 
              "Above_t_num", 
              "Above_t_perc", 
              "Protein_cost_all_uk",
              "MRI_cost_all_uk",
              "MRI_cost_above_t_uk",
              "Cost_protein_plus_selected_MRI_uk",
              "Savings_uk",
              "Protein_cost_all_usa",
              "MRI_cost_all_usa",
              "MRI_cost_above_t_usa",
              "Cost_protein_plus_selected_MRI_usa",
              "Savings_usa")
res <- as.data.frame(matrix(,
                            nrow = length(proteins),
                            ncol = length(resnames)))
colnames(res) <- resnames

#Calculate costs per patients
for (i in 1:length(proteins)){
  cutdata$Myprotein <- cutdata[, proteins[i]]
  
  #calculate optimal threshold for biomarker concentrations
  mycut <- optimal.cutpoints(Myprotein ~ BS_visual, data = cutdata, 
                                 methods = "MinValueSe",
                                 tag.healthy = 0, 
                                 direction = "<",
                           control = control.cutpoints(valueSe = 0.90 ))
  sum <- summary(mycut)[["MinValueSe"]][["Global"]][["optimal.cutoff"]]
  res[i, "Protein"]             <- proteins[i]
  res[i, "Threshold"]           <- sum$cutoff
  res[i, "Sensitivity"]         <- sum$Se[[1]]
  res[i, "Specificity"]         <- sum$Sp[[1]]
  num                           <- cutdata %>% filter(get(proteins[i]) >= sum$cutoff) %>% nrow()
  res[i, "Above_t_num"]         <- num
  res[i, "Above_t_perc"]        <- num/nrow(cutdata)*100
  
  #calculate UK costs
  res[i, "Protein_cost_all_uk"] <- p_cost_uk[i] + phleb_uk
  res[i, "MRI_cost_all_uk"]     <- mri_cost_uk
  res[i, "MRI_cost_above_t_uk"] <- (num*mri_cost_uk)/nrow(cutdata)
  res[i, "Cost_protein_plus_selected_MRI_uk"] <- res[i, "Protein_cost_all_uk"] + res[i, "MRI_cost_above_t_uk"]
  res[i, "Savings_uk"]          <- res[i, "MRI_cost_all_uk"] - res[i, "Cost_protein_plus_selected_MRI_uk"]
  
  #calculate US costs
  res[i, "Protein_cost_all_usa"] <- p_cost_usa[i] + phleb_usa
  res[i, "MRI_cost_all_usa"]     <- mri_cost_usa
  res[i, "MRI_cost_above_t_usa"] <- (num*mri_cost_usa)/nrow(cutdata)
  res[i, "Cost_protein_plus_selected_MRI_usa"] <- res[i, "Protein_cost_all_usa"] + res[i, "MRI_cost_above_t_usa"]
  res[i, "Savings_usa"]          <- res[i, "MRI_cost_all_usa"] - res[i, "Cost_protein_plus_selected_MRI_usa"]
  
}

#Format displayed names and numbers
res$Protein             <- gsub("UCH.L1", "UCH-L1", res$Protein)
##performance
res$Sensitivity         <- format(round(res$Sensitivity, 2), nsmall = 2)
res$Specificity         <- format(round(res$Specificity, 2), nsmall = 2)
res$Above_t_perc        <- format(round(res$Above_t_perc, 0), nsmall = 0)
res$Above_t_num         <- paste0(res$Above_t_num, " (", res$Above_t_perc, "%)")
##UK costs
res$Cost_protein_plus_selected_MRI_uk <- format(round(res$Cost_protein_plus_selected_MRI_uk, 2), nsmall = 2)
res$MRI_cost_all_uk     <- format(round(res$MRI_cost_all_uk, 2), nsmall = 2)
res$Savings_uk          <- format(round(res$Savings_uk, 2), nsmall = 2)
#USA costs
res$Cost_protein_plus_selected_MRI_usa <- format(round(res$Cost_protein_plus_selected_MRI_usa, 2), nsmall = 2)
res$MRI_cost_all_usa     <- format(round(res$MRI_cost_all_usa, 2), nsmall = 2)
res$Savings_usa          <- format(round(res$Savings_usa, 2), nsmall = 2)

#select columns for display
res <- 
  res %>%
  select(Protein, 
         Threshold, 
         Sensitivity, 
         Specificity, 
         Above_t_num, 
         MRI_cost_all_uk, 
         Cost_protein_plus_selected_MRI_uk, 
         Savings_uk,
         MRI_cost_all_usa, 
         Cost_protein_plus_selected_MRI_usa, 
         Savings_usa)

#USA results for biomarkers other than GFAP and UCHL1 need to be blank
res[c(2:5), c("Cost_protein_plus_selected_MRI_usa", "Savings_usa")] <- "-"

#Save results as pretty table
ft <- 
    res %>% 
    regulartable() %>%
    set_header_labels(Above_t_num = "Patients above threshold",   
                      MRI_cost_all_uk = "MRI for all",
                      Cost_protein_plus_selected_MRI_uk = "Protein plus selected MRI",
                      Savings_uk = "Savings",
                      MRI_cost_all_usa = "MRI for all",
                      Cost_protein_plus_selected_MRI_usa = "Protein plus selected MRI",
                      Savings_usa = "Savings") %>%
    autofit() %>%
    add_header_row(values = c(" "," "," ", " ", " ", 
                              "Costs in the United Kingdom\n(GBP per patient)", 
                              "Costs in the United States\n(USD per patient)"),
                   colwidths = c(1,1, 1, 1, 1, 
                                 3,
                                 3))%>%
    theme_vanilla() %>%
    padding(padding = 0, part = "all", padding.right = 3, padding.left = 3) %>%
    align(align = "right", part = "all")%>%
    align(align = "left", j = c(1), part = "all") %>%
    align(align = "center", part = "header", i = c(1))
ft

save_as_docx(ft, path = "Figures/Savings_90perc.docx")

