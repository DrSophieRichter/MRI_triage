##############################################
#       Functions used in the paper          #
# Serum biomarkers identify critically ill   #
# traumatic brain injury patients for MRI    #               
##############################################


###########################################
#rendering functions for table1
###########################################

#set up display of continuous variables
my.render.cont <- function(x) 
{
  #display 2 sig figs:
  with(stats.apply.rounding(stats.default(x), digits=2), 
       #in first column display variable name and "Median (IQR)":
       c("", "Median (Min-Max)" = 
           #in second column calculate and display:
           sprintf("%s (%s %s %s)", 
                   MEDIAN,  
                   round(min(x, na.rm = T),1) , 
                   "-", 
                   round(max(x, na.rm = T),1))))   #calculations requested
}


#set up display of categorical variables
my.render.cat <- function(x) 
{
  #in first column display variable name:
  c("", 
    #in second column calculate and display:
    sapply(stats.default(x), function(y) with(y, sprintf("%d (%0.0f %%)", 
                                                         FREQ, 
                                                         round(PCT, 0)))))
}


###################################################
#
# Compare biomarker levels between two groups
# with an imaging feature present vs absent
#
####################################################

compare_bio_levels <- function(data, mycohort, myvar, mytitle){
  
  if (mycohort == "all"){
    data <- 
      data 
  }else {
    data <- 
      data %>% 
      filter(Cohort == mycohort)
  }
  
  
  mydata <-
    data %>%
    select(Master_subject_ID, myvar, GFAP:UCH.L1) %>% 
    gather(key = "Protein", value = "Log_conc", GFAP:UCH.L1) %>%
    spread(key = myvar, value = "Log_conc")
  
  proteins <- c("GFAP", "NFL", "NSE", "S100B", "Tau",  "UCH.L1")
  res <- as.data.frame(matrix(,
                              nrow = length(proteins)+4,
                              ncol = 4))
  colnames(res) <- c("Protein", "No_BS", "With_BS", "P_raw")
  for (p in 1:length(proteins)){
    df <- mydata %>% filter(Protein == proteins[p])
    
    res[p+4, "Protein"] <- proteins[p]
    
    med <- median(df$absent, na.rm = TRUE) %>% exp() 
    med <- format(round(as.numeric(med), 2), nsmall = 2)
    Q1 <- quantile(exp(df$absent), na.rm = TRUE)[2] 
    Q1 <- format(round(Q1,2), nsmall = 2)
    Q3 <- quantile(exp(df$absent), na.rm = TRUE)[4]
    Q3 <- format(round(Q3,2), nsmall = 2)
    res[p+4, "No_BS"] <- paste0(med, " (", Q1, "-", Q3, ")")
    
    med <- median(df$present, na.rm = TRUE) %>% exp() 
    med <- format(round(as.numeric(med), 2), nsmall = 2)
    Q1 <- quantile(exp(df$present), na.rm = TRUE)[2] 
    Q1 <- format(round(Q1,2), nsmall = 2)
    Q3 <- quantile(exp(df$present), na.rm = TRUE)[4]
    Q3 <- format(round(Q3,2), nsmall = 2)
    res[p+4, "With_BS"] <- paste0(med, " (", Q1, "-", Q3, ")")
    
    res[p+4, "P_raw"] <- wilcox.test(df$present, df$absent)$p.value
  }
  
  
  res[1, "Protein"] <- "Number of patients"
  res[2, "Protein"] <- "Sample time in hours"
  res[3, "Protein"] <- "GOSE < 8"
  res[4, "Protein"] <- "GOSE < 5"
  
  temp <- 
    data %>%
    filter(is.na(get(myvar))==FALSE)
  if (mycohort == "all"){
    temp <- temp
  } else {
    temp <- 
      temp %>% 
      filter(Cohort == mycohort)
  }
  
  res[2, "P_raw"] <- 
    wilcox.test(temp[,"Time_to_bio"] ~ temp[,myvar])$p.value %>%
    as.numeric()
  
  tab <- table(temp[,"Outcome_8"],temp[,myvar])
  res[3, "P_raw"] <- 
    fisher.test(tab)$p.value %>%
    as.numeric()
  res[3, "No_BS"] <- paste0(tab["0", "absent"], " vs ", tab["1", "absent"])
  res[3, "With_BS"] <- paste0(tab["0", "present"], " vs ", tab["1", "present"])
  
  tab <- table(temp[,"Outcome_5plus"],temp[,myvar])
  res[4, "P_raw"] <- 
    fisher.test(tab)$p.value %>%
    as.numeric()
  res[4, "No_BS"] <- paste0(tab["0", "absent"], " vs ", tab["1", "absent"])
  res[4, "With_BS"] <- paste0(tab["0", "present"], " vs ", tab["1", "present"])
  
  temp <-
    temp %>%
    group_by(get(myvar))%>%
    summarise(num=n(),
              freq = round(n()/nrow(temp)*100, 0), 
              time = round(median(Time_to_bio),0),
              Q1 = round(quantile(Time_to_bio)[2],0),
              Q3 = round(quantile(Time_to_bio)[4], 0))
  temp$num <- paste0(temp$num, " (", temp$freq, "%)")
  temp$time <- paste0(temp$time, " (", temp$Q1, "-", temp$Q3, ")")
  temp <- 
    temp %>%
    select(myvar = "get(myvar)", 
           num, 
           time)
  temp <- 
    t(temp) %>%
    as.data.frame()
  #make first row title
  colnames(temp) <- temp[1,]
  temp <- temp[-1,]
  res[1, "No_BS"] <- temp[1, "absent"]
  res[1, "With_BS"] <- temp[1, "present"]
  res[2, "No_BS"] <- temp[2, "absent"]
  res[2, "With_BS"] <- temp[2, "present"]
  
  res$P_adj <- p.adjust(res$P_raw, method = "fdr") 
  res$P_adj <- format(round(res$P_adj,3), nsmall = 3)
  res$P_raw <- format(round(as.numeric(res$P_raw),3), nsmall = 3)
  res$P_adj <- ifelse(is.na(res$P_raw)==TRUE, NA, as.numeric(res$P_adj))
  res$Sig <- ifelse(is.na(res$P_adj)==TRUE, NA,
                    ifelse(res$P_adj < 0.05, "signif", "ns"))
  res[1, "P_raw"] <- NA
  
  colnames(res)[2:3] <- c("Absent", "Present")
  
  res %>% 
    regulartable() %>% 
    autofit %>%
    align(align = "right", part = "all") %>%
    theme_zebra() %>% 
    add_header_lines(paste0(mytitle, " in ", mycohort, " cohort")) %>%
    save_as_docx(path = paste0("Figures/", mycohort, "_", myvar,".docx"))
}


#################################################
#
# Similar function to compare_bio_levels
# but for a 4-level factor such as AG_stage
# rather than a binary variable
#
#################################################

compare_bio_stages <- function(data, mycohort, myvar, mytitle){
  
  if (mycohort == "all"){
    data <- 
      data 
  }else {
    data <- 
      data %>% 
      filter(Cohort == mycohort)
  }
  
  mydata <- 
    data %>% 
    filter((is.na(get(myvar)) == FALSE)) %>%
    select(Master_subject_ID, myvar, GFAP:UCH.L1) %>% 
    gather(key = "Protein", value = "Log_conc", GFAP:UCH.L1) 
  
  proteins <- c("GFAP", "NFL", "NSE", "S100B", "Tau",  "UCH.L1")
  res <- as.data.frame(matrix(,
                              nrow = length(proteins)+4,
                              ncol = 6))
  colnames(res) <- c("Protein", "Stage_0", "Stage_1","Stage_2", "Stage_3", "P_raw")
  
  for (p in 1:length(proteins)){
    df <- 
      mydata %>% 
      spread(key = myvar, value = "Log_conc") %>%
      filter(Protein == proteins[p])
    
    
    res[p+4, "Protein"] <- proteins[p]
    
    med <- median(df$`0`, na.rm = TRUE) %>% exp() 
    med <- format(round(as.numeric(med), 2), nsmall = 2)
    Q1 <- quantile(exp(df$`0`), na.rm = TRUE)[2] 
    Q1 <- format(round(Q1,2), nsmall = 2)
    Q3 <- quantile(exp(df$`0`), na.rm = TRUE)[4]
    Q3 <- format(round(Q3,2), nsmall = 2)
    res[p+4, "Stage_0"] <- paste0(med, " (", Q1, "-", Q3, ")")
    
    med <- median(df$`1`, na.rm = TRUE) %>% exp() 
    med <- format(round(as.numeric(med), 2), nsmall = 2)
    Q1 <- quantile(exp(df$`1`), na.rm = TRUE)[2] 
    Q1 <- format(round(Q1,2), nsmall = 2)
    Q3 <- quantile(exp(df$`1`), na.rm = TRUE)[4]
    Q3 <- format(round(Q3,2), nsmall = 2)
    res[p+4, "Stage_1"] <- paste0(med, " (", Q1, "-", Q3, ")")
    
    med <- median(df$`2`, na.rm = TRUE) %>% exp() 
    med <- format(round(as.numeric(med), 2), nsmall = 2)
    Q1 <- quantile(exp(df$`2`), na.rm = TRUE)[2] 
    Q1 <- format(round(Q1,2), nsmall = 2)
    Q3 <- quantile(exp(df$`2`), na.rm = TRUE)[4]
    Q3 <- format(round(Q3,2), nsmall = 2)
    res[p+4, "Stage_2"] <- paste0(med, " (", Q1, "-", Q3, ")")
    
    if (ncol(df) < 6) {
      res[p+4, "Stage_3"] <- NA
    } else {
      med <- median(df$`3`, na.rm = TRUE) %>% exp() 
      med <- format(round(as.numeric(med), 2), nsmall = 2)
      Q1 <- quantile(exp(df$`3`), na.rm = TRUE)[2] 
      Q1 <- format(round(Q1,2), nsmall = 2)
      Q3 <- quantile(exp(df$`3`), na.rm = TRUE)[4]
      Q3 <- format(round(Q3,2), nsmall = 2)
      res[p+4, "Stage_3"] <- paste0(med, " (", Q1, "-", Q3, ")")
      
      df <- 
        mydata %>% 
        filter(Protein == proteins[p])
      df[, myvar] <- factor(df[, myvar], ordered = TRUE)
      res[p+4, "P_raw"] <- JonckheereTerpstraTest(Log_conc ~ get(myvar), data = df)$p.value
    }
    
  }
  
  
  res[1, "Protein"] <- "Number of patients"
  res[2, "Protein"] <- "Sample time in hours"
  res[3, "Protein"] <- "GOSE < 8"
  res[4, "Protein"] <- "GOSE < 5"
  
  temp <- 
    data %>%
    filter((is.na(get(myvar))==FALSE))
  
  if (mycohort == "all"){
    temp <- temp
  } else {
    temp <- 
      temp %>% 
      filter(Cohort == mycohort)
  }
  
  temp[, myvar] <- factor(temp[, myvar], ordered = TRUE)
  
  res[2, "P_raw"] <- 
    JonckheereTerpstraTest(temp[,"Time_to_bio"] ~ temp[,myvar])$p.value %>%
    as.numeric()
  
  tab <- table(temp[,"Outcome_8"],temp[,myvar])
  res[3, "P_raw"] <- 
    fisher.test(tab)$p.value %>%
    as.numeric()
  res[3, "Stage_0"] <- paste0(tab["0", "0"], " vs ", tab["1", "0"])
  res[3, "Stage_1"] <- paste0(tab["0", "1"], " vs ", tab["1", "1"])
  res[3, "Stage_2"] <- paste0(tab["0", "2"], " vs ", tab["1", "2"])
  df <- 
    mydata %>% 
    spread(key = myvar, value = "Log_conc") %>%
    filter(Protein == proteins[p])
  if (ncol(df) < 6) {
    res[3, "Stage_3"] <- NA
  } else {
    res[3, "Stage_3"] <- paste0(tab["0", "3"], " vs ", tab["1", "3"])
  }
  
  
  tab <- table(temp[,"Outcome_5plus"],temp[,myvar])
  res[4, "P_raw"] <- 
    fisher.test(tab)$p.value %>%
    as.numeric()
  res[4, "Stage_0"] <- paste0(tab["0", "0"], " vs ", tab["1", "0"])
  res[4, "Stage_1"] <- paste0(tab["0", "1"], " vs ", tab["1", "1"])
  res[4, "Stage_2"] <- paste0(tab["0", "2"], " vs ", tab["1", "2"])
  if (ncol(df) < 6) {
    res[4, "Stage_3"] <- NA
  } else {
    res[4, "Stage_3"] <- paste0(tab["0", "3"], " vs ", tab["1", "3"])
  }
  
  temp <-
    temp %>%
    group_by(get(myvar))%>%
    summarise(num=n(),
              freq = round(n()/nrow(temp)*100, 0), 
              time = round(median(Time_to_bio),0),
              Q1 = round(quantile(Time_to_bio)[2],0),
              Q3 = round(quantile(Time_to_bio)[4], 0))
  temp$num <- paste0(temp$num, " (", temp$freq, "%)")
  temp$time <- paste0(temp$time, " (", temp$Q1, "-", temp$Q3, ")")
  temp <- 
    temp %>%
    select(myvar = "get(myvar)", 
           num, 
           time)
  temp <- 
    t(temp) %>%
    as.data.frame()
  #make first row title
  colnames(temp) <- temp[1,]
  temp <- temp[-1,]
  res[1, "Stage_0"] <- temp[1, "0"]
  res[1, "Stage_1"] <- temp[1, "1"]
  res[1, "Stage_2"] <- temp[1, "2"]
  
  res[2, "Stage_0"] <- temp[2, "0"]
  res[2, "Stage_1"] <- temp[2, "1"]
  res[2, "Stage_2"] <- temp[2, "2"]
  
  if (ncol(df) < 6) {
    res[1, "Stage_3"] <- NA
    res[2, "Stage_3"] <- NA
  } else {
    res[1, "Stage_3"] <- temp[1, "3"]
    res[2, "Stage_3"] <- temp[2, "3"]
  }
  
  
  res$P_adj <- p.adjust(res$P_raw, method = "fdr") 
  res$P_adj <- format(round(res$P_adj,3), nsmall = 3)
  res$P_raw <- format(round(as.numeric(res$P_raw),3), nsmall = 3)
  res$P_adj <- ifelse(is.na(res$P_raw)==TRUE, NA, as.numeric(res$P_adj))
  res$Sig <- ifelse(is.na(res$P_adj)==TRUE, NA,
                    ifelse(res$P_adj < 0.05, "signif", "ns"))
  res[1, "P_raw"] <- NA
  
  res %>% 
    regulartable() %>% 
    autofit %>%
    align(align = "right", part = "all") %>%
    theme_zebra() %>% 
    add_header_lines(paste0(mytitle, " in ", mycohort, " cohort")) %>%
    save_as_docx(path = paste0("Figures/", mycohort, "_", myvar,".docx"))
}

################################################################
#
# Fit models for proteins to predict BSI
# and plot roc curves
# displaying AUC
#
################################################################

plot_roc_curves <- function(data, mycohort, myvar){
  
  if (mycohort == "all"){
    temp <- data
  } else {
    temp <- data %>% filter(Cohort == mycohort)
  }
  
  temp[,myvar] <- ifelse(is.na(temp[, myvar])==TRUE, NA,
                        ifelse(temp[, myvar] == "absent", 0, 1))
  
  proteins <- c("GFAP", "NFL", "NSE", "S100B", "Tau", "UCH.L1")
  protein_list <- vector("list", length = length(proteins))
  names(protein_list) <- proteins
  elements <- c("fit", "probs", "g","auc", "plot")
  roc_list <- vector("list", length = length(proteins))
  names(roc_list) <- proteins
  
  for (p in 1:length(proteins)){
    protein_list[[p]]            <- vector("list", length = length(elements))
    names(protein_list[[p]])     <- elements
    #fit the model
    protein_list[[p]][["fit"]]     <- glm(get(myvar) ~ get(proteins[p]), data = temp, family = "binomial")
    #make predictions
    protein_list[[p]][["probs"]]   <- predict(object = protein_list[[p]][["fit"]], 
                                              newdata = temp, 
                                              type= "response")
    #compare with observations to create ROC curve
    protein_list[[p]][["g"]]           <- roc(get(myvar) ~ protein_list[[p]][["probs"]], data = temp)
    
    auc <- protein_list[[p]][["g"]]$auc[[1]] %>% as.numeric
    auc <- format(round(auc, 2), nsmall = 2)
    LL <- ci.auc(protein_list[[p]][["g"]])[[1]] %>% as.numeric
    LL <- format(round(LL, 2), nsmall = 2)
    UL <- ci.auc(protein_list[[p]][["g"]])[[3]] %>% as.numeric
    UL <- format(round(UL, 2), nsmall = 2)
    
    protein_list[[p]][["auc"]] <- paste0(auc, " (", LL, "-", UL, ")")
    
    #plot roc curve
    protein_list[[p]][["plot"]] <- 
      ggroc(protein_list[[p]][["g"]], legacy.axes= F, size=0.6) +
      theme_bw(base_size = 6) +
      annotate("rect", 
               xmin = 1, 
               xmax = 0, 
               ymin = 0.90, 
               ymax = 1,
               alpha = .2, 
               fill = "darkolivegreen4") +
      geom_abline(intercept = 1, 
                  slope = 1, 
                  color = "darkgrey", 
                  linetype = "dashed") +
      annotate("text", 
               x = 0.2, 
               y = 0.1, 
               size = 2.1,
               label = paste(names(protein_list)[p], protein_list[[p]][["auc"]])
      )
    
  }
  
  
  mygrid <- 
    grid.arrange(protein_list[[1]][["plot"]], 
                 protein_list[[2]][["plot"]],
                 protein_list[[3]][["plot"]], 
                 protein_list[[4]][["plot"]],
                 protein_list[[5]][["plot"]], 
                 protein_list[[6]][["plot"]])
  
  
  ggsave(filename = paste0("Figures/", mycohort, "_", myvar, "_roc_bio.tiff"), 
         plot = mygrid,
         dpi = 300,
         width = 13,
         height = 13,
         units = "cm")
}





