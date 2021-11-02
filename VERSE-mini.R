# Install necessary supporting programs
install.packages("ggrepel")
install.packages("prevR")
install.packages("sf")
install.packages("matchmaker")
install.packages("rmapshaper")
install.packages("PHEindicatormethods")
install.packages("radiant")
install.packages("StatMeasures")
#devtools::install_github("ropensci/rdhs", ref = "issue33_path")
devtools::install_github("ropensci/rdhs")

#Library needed programs
rm(list=ls()) 
library(usethis)
library(rmapshaper)
library(PHEindicatormethods)
library(radiant)
library(StatMeasures)
library(grid)
library(Matrix)
library(dplyr)
library(tidyverse)
library(devtools)
library(rdhs)
library(survey)
library(haven)
library(margins)
library(tibble)
library(rineq)
library(labelled)
library(ggplot2)
library(ggforce)
library(ggrepel)
library(prevR)
library(sf)
library(matchmaker)
library(RColorBrewer)
library(foreign)

# Set Working Directory
setwd("~/Downloads")

#Create Input Vectors
VACCINES<-c("HIB1")
SCHEDULE<-c("DEFAULT")
GEO = "Region"
COUNTRY<-"China"
FACTORS<-c("region","rural","education","wealth","sex","insurance")
DATA<-read_csv("VERSE-China-data format2.csv")


VERSE <- function(DATA,COUNTRY,VACCINES,SCHEDULE,FACTORS,GEO){

# Create Storage Vectors & Dataframes for national-level results
  output <- data.frame(matrix(NA, nrow = length(FACTORS)+3, ncol = 1))
  output_list <- list()
  efficiency_list <- list()
  CI_Results<- c()
  CI_E_Results<- c()
  AEG_Composite_Results<- c()
  Coverage_Results<- c()
  
# Rename DATA to correspond with VERSE programming
    dhs_data <- DATA

#Create Vaccine Schedule dataframe from SCHEDULE vector
    schedule<-as.data.frame(matrix(as.numeric(SCHEDULE),nrow=1))
    schedule_names<-VACCINES                          
    colnames(schedule) <- schedule_names
  
  #Defining fair inequity based upon Vaccine Schedule by creating vaccine-specific underage cutoffs
  for (i in VACCINES){
    dhs_data$q<- ifelse(dhs_data$hw1 < schedule[,i][1], 1, 
                        ifelse(dhs_data$hw1>=schedule[,i][1], 0, 0))
    dhs_data$q<- replace(dhs_data$q,is.na(dhs_data$q),0)
    names(dhs_data)[names(dhs_data) == "q"] <- paste("underage_",i, sep="")
  }
  
  
# Create Reference Values Vector
  FACT<-c()
  
  for (l in 1:length(FACTORS)){
    ifelse(FACTORS[l]=="region", FACT[l]<-"v101",
           ifelse(FACTORS[l]=="rural", FACT[l]<-"v025",
                  ifelse(FACTORS[l]=="education", FACT[l]<-"v106",
                         ifelse(FACTORS[l]=="wealth", FACT[l]<-"v190",
                                ifelse(FACTORS[l]=="sex", FACT[l]<-"b4",
                                       ifelse(FACTORS[l]=="insurance", FACT[l]<-"v481","NA"))))))
  }
  

  # Create Reference Levels Based on Fully Immunized for Age Outcome
  REF<- c()
  
  for(l in FACT){
    REF_DATA<-dhs_data
    REF_DATA[,"FULL"] <- as.numeric(unlist(REF_DATA[,"FULL"]))
    REF_DATA<-subset(dhs_data, dhs_data[,"FULL"]==1|dhs_data[,"FULL"]==0)
    REF_DATA[,"FULL"] <- as.numeric(unlist(REF_DATA[,"FULL"]))
    invisible(capture.output(REF_DATA<- REF_DATA %>% 
                               group_by(REF_DATA[,l]) %>% 
                               summarise(coverage=weighted.mean(FULL,v005))))
    REF_DATA<- filter(REF_DATA, coverage==max(REF_DATA$coverage))
    value<-c(as.numeric(unlist(REF_DATA[1,1])))
    REF<-c(REF,value)
  }
  
  if (REF[5]==0){
    REF[5]<-1
  }
  
  n = length(REF)
  
  for (j in 1:n){
    nameref <- paste("Ref",j, sep="")
    dhs_data[,nameref] = REF[j]
  }
  
  #Create Names Vector for Regions
  GEO_CI<- c(dhs_data$geo_names)
  
  # Crete shell dataframes to store the geographic-specific outputs
  CI_Results_GEO_Output <- data.frame(matrix(NA, nrow = length(GEO_CI), ncol = 1))
  CI_E_Results_GEO_Output<- data.frame(matrix(NA, nrow = length(GEO_CI), ncol = 1))
  AEG_Composite_Results_GEO_Output<- data.frame(matrix(NA, nrow = length(GEO_CI), ncol = 1))
  Coverage_Results_GEO_Output<- data.frame(matrix(NA, nrow = length(GEO_CI), ncol = 1))
  
  for (i in VACCINES){
    
    #i = "HIB1"
    
    # convert vaccines to to numeric vector
    dhs_data[,i] <- as.numeric(unlist(dhs_data[,i]))
    
    #Create Vaccine-Specific Variables
    underage_name <- paste("underage_",i, sep="")
    q <- dhs_data[,underage_name]
    dhs_data<-cbind(dhs_data,q)
    #dhs_data$underage_i = dhs_data[,underage_name]
    outcome_name <- paste("",i, sep="")
    q<- dhs_data[, outcome_name]
    dhs_data<-cbind(dhs_data,q)
    #dhs_data$outcome =  dhs_data[, outcome_name]
    #q<-  dhs_data$outcome
    #dhs_data<-cbind(dhs_data,q)
    
    # Create Design Matrix
    #design <- svydesign(data = dhs_data,
    #                    ids = ~caseid,
    #                    weights = ~v005)
    
    # Subset Data to be only data where outcome data is available
    data_i <- subset(dhs_data, dhs_data[,i]==1|dhs_data[,i]==0)
    data_i<- rename(data_i, outcome = outcome_name)
    data_i<- rename(data_i, underage_i = underage_name)
    
    # Create Design Matrix
    design <- svydesign(data = data_i,
                        ids = ~caseid,
                        weights = ~v005)
    
    #Set Up Logistic Regression
    
    mean<-mean(data_i$underage_i)
    
    if (i=="BCG"|i=="FULL"|mean[1]==0) {
      logit_i <- svyglm(outcome ~ 
                          relevel(factor(v101), ref = Ref1[1]) + 
                          relevel(factor(v025), ref = Ref2[1]) + 
                          relevel(factor(v106), ref = Ref3[1]) + 
                          relevel(factor(v190), ref = Ref4[1]) + 
                          relevel(factor(b4), ref = Ref5[1]) + 
                          relevel(factor(v481), ref = Ref6[1]), 
                        design = design, family = binomial(link="logit"), data = data_i)
    } else{
      logit_i <- svyglm(outcome ~ factor(underage_i) + 
                          relevel(factor(v101), ref = Ref1[1]) + 
                          relevel(factor(v025), ref = Ref2[1]) + 
                          relevel(factor(v106), ref = Ref3[1]) + 
                          relevel(factor(v190), ref = Ref4[1]) + 
                          relevel(factor(b4), Ref5[1]) + 
                          relevel(factor(v481), ref = Ref6[1]), 
                        design = design, family = binomial(link="logit"), data = data_i)
    }
    
    summary(logit_i)
    
    # Calculate Predicted Probabilities & Store in Pred_Probs Object
    pred_probs <- data.frame(hci=predict(logit_i, data = dhs_data, type="response",na.action = na.exclude))
    
    # Computing Direct Unfairness Metric
    # Compute predicted probability holding child age (fair equity) at reference category (DTP1_underage==0)
    
    data_i_unfair <- data.frame(data_i)
    #data_i_unfair$underage_i <- replace(data_i_unfair[,underage_i],data_i_unfair[,underage_i]==1,0)
    data_i_unfair$underage_i <- replace(data_i_unfair$underage_i,data_i_unfair$underage_i==1,0)
    
    # Create New Predictions on Direct Unfairness Dataset
    pred_probs_2 <- data.frame(hci_du = predict(logit_i, data_i_unfair, type="response"))
    
    # Computing Predicted Fairness Metric    
    # Compute predicted probability holding all unfair variables at reference levels (fairness equity) at reference category levels
    
    data_i_fair <- data.frame(data_i)
    
    # Create Reference Level Variable Dataset for all of the Fair Predictors
    # Set Unfair Predictors to Reference Levels
    
    data_i_fair$v101 <- replace(data_i_fair$v101,data_i_fair$v101!=0,0)
    data_i_fair$v025 <- replace(data_i_fair$v025,data_i_fair$v025!=0,0)
    data_i_fair$v106 <- replace(data_i_fair$v106,data_i_fair$v106!=0,0)
    data_i_fair$v190 <- replace(data_i_fair$v190,data_i_fair$v190!=0,0)
    data_i_fair$b4 <- replace(data_i_fair$b4,data_i_fair$b4!=0,0)
    data_i_fair$v481 <- replace(data_i_fair$v481,data_i_fair$v481!=0,0)
    
    
    data_i_fair$v101 <- replace(data_i_fair$v101,data_i_fair$v101!=REF[1],REF[1])
    data_i_fair$v025 <- replace(data_i_fair$v025,data_i_fair$v025!=REF[2],REF[2])
    data_i_fair$v106 <- replace(data_i_fair$v106,data_i_fair$v106!=REF[3],REF[3])
    data_i_fair$v190 <- replace(data_i_fair$v190,data_i_fair$v190!=REF[4],REF[4])
    data_i_fair$b4 <- replace(data_i_fair$b4,data_i_fair$b4!=REF[5],REF[5])
    data_i_fair$v481 <- replace(data_i_fair$v481,data_i_fair$v481!=REF[6],REF[6])
    
    # Set Constant (Unfair) Predictors as Factors
    data_i_fair$v101 <- as.numeric(data_i_fair$v101)
    data_i_fair$v025 <- as.numeric(data_i_fair$v025)
    data_i_fair$v106 <- as.numeric(data_i_fair$v106)
    data_i_fair$v190 <- as.numeric(data_i_fair$v190)
    data_i_fair$b4 <- as.numeric(data_i_fair$b4)
    data_i_fair$v481 <- as.numeric(data_i_fair$v481)
    
    # Create New Predictions on Predicted Fairness Dataset
    data_i_fair$predict = summary(logit_i)$coefficients[1] + summary(logit_i)$coefficients[2]*data_i_fair$underage_i
    
    logit2prob <- function(logit){
      odds <- exp(logit)
      prob <- odds / (1 + odds)
      return(prob)
    }
    
    prob <- logit2prob(data_i_fair$predict)
    
    pred_probs_3 <- data.frame(hci_fair = prob)
    
    # Calculating the Direct Concentration Index
    direct_ci <- ci(y = data_i[,"outcome"], x = pred_probs_2$hci_du.response, type = "CI")
    CI_1 <- concentration_index(direct_ci)
    CI_1
 
    # Calculating the Erreygers Corrected Composite Concentration Index 
    CIE <- ci(y = data_i[,"outcome"], x = pred_probs_2$hci_du.response, type = "CIc")
    CI_E <- concentration_index(CIE)
    CI_E
    
    # Calculating the Composite Absolute Equity Gap
    rank<-c(pred_probs_2$hci_du.response)
    quantiles<- quantile(rank, probs = seq(0, 1, 0.20))
    comp_data<- data.frame(cbind(data_i[,"outcome"],data_i[,"v005"],rank))
    names(comp_data)[names(comp_data) == "V1"] <- paste("",i, sep="")
    names(comp_data)[names(comp_data) == "V2"] <- paste("v005","", sep="")
    AEG_Composite <- round(weighted.mean(comp_data[,i][rank>=quantiles[5]],comp_data[,"v005"][rank>=quantiles[5]]) - weighted.mean(comp_data[,i][rank>=quantiles[2]],comp_data[,"v005"][rank>=quantiles[2]]), digits=3)
    
    assign(paste("CI_1_",i, sep=""),CI_1)
    assign(paste("CI_E_",i, sep=""),CI_E)
    assign(paste("AEG_Composite_",i, sep=""),AEG_Composite)
    
    #Subset Coverage to be only among correct for age-group for Vaccine (older than schedule age)
    data_ic <- subset(data_i, data_i[,"underage_i"]==0)
    coverage<- c(weighted.mean(data_ic$outcome,data_ic$v005))
    
    #Name other CI Results
    CI_Results<- c(CI_Results, assign(paste("CI_1_",i, sep=""),CI_1))
    CI_E_Results<- c(CI_E_Results, assign(paste("CI_E_",i, sep=""),CI_E))
    AEG_Composite_Results <- c(AEG_Composite_Results, assign(paste("AEG_Composite_",i, sep=""), AEG_Composite))
    Coverage_Results<- c(Coverage_Results, assign(paste("Coverage_",i, sep=""),coverage))
    
    #Decomposition Analysis
    design <- svydesign(data = data_i,
                        ids = ~caseid,
                        weights = ~v005)
    
    logit_de <- svyglm(outcome ~ underage_i+v101+v025+v106+v190+b4+v481, 
                       design = design, family = binomial(link="logit"), data = data_i)
    
    ###Loop Decomposition###
    new_data <- new_data <- data_i %>% select("underage_i",all_of(FACT))
    decomposition <- data.frame()
    decomposition_x <-data.frame() 
    for (x in colnames(new_data)) {
      decomposition_x <- tibble(
        name_t = c(x),
        cindex_t = c(concentration_index(ci(y = new_data[[x]], x = pred_probs_2$hci_du.response, type = "CI"))),
        elastic_t = c(mean(coef(logit_de)[x]*logit_de$fitted.values*(1-logit_de$fitted.values)*new_data[[x]]/logit_de$fitted.values)),
        contribution_t = abs(as.numeric(elastic_t)*as.numeric(cindex_t)),
        percent_t = 100*abs(as.numeric(contribution_t)/CI_1),
      )
      decomposition <- rbind.data.frame(decomposition,  decomposition_x)
    }
    
    resid_c_t <- CI_1-sum(decomposition$contribution_t, na.rm=TRUE)
    
    residual <- tibble(name_t = "Residual",
                       cindex_t = "" ,
                       elastic_t = "",
                       contribution_t = resid_c_t,
                       percent_t = 100*abs(resid_c_t/CI_1))
    
    decomposition <- rbind.data.frame(decomposition,  residual)
    
    decomposition$contribution_t[is.na(decomposition$contribution_t)] <- 0
    decomposition$percent_t[is.na(decomposition$percent_t)] <- 0
    
    if (sum(decomposition$percent_t)>100){
      decomposition$percent_t<-decomposition$percent_t/(sum(decomposition$percent_t))*100
    }
 
  #Create Decomposition Pie Charts:   
    decomposition_pie <- decomposition %>% 
      mutate(end = 2 * pi * cumsum(percent_t)/sum(percent_t),
             start = lag(end, default = 0),
             middle = 0.5 * (start + end),
             hjust = ifelse(middle > pi, 1, 0),
             vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))
    
    gtitle<- (ifelse(i=="ZERO","Decomposition of Zero-Dose Inequity",
                     ifelse(i=="FULL","Decomposition of Fully Immunized for Age Equity", 
                            ifelse(i=="COMPLETE","Decomposition of Equity in Completing Vaccination Schedule",
                                   paste(paste("Decomposition of",i, sep=" "),"Coverage Equity", sep=" ")))))
    
    pie_i<- ggplot(decomposition_pie, aes(ymax=end, ymin=start, xmax=4.25, xmin=3, fill = factor(name_t))) +
      geom_rect() +
      ggtitle(gtitle) +
      theme(plot.title = element_text(size=14,face="bold",hjust = 0),
            legend.title=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y = element_blank(),
            panel.border = element_blank(),
            panel.grid=element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank()) +
      coord_polar(theta="y") +
      xlim(c(2, 4.25)) + #Turns the Pie into a Donut Chart
      scale_fill_manual(values=c("#74add1","grey80","#d73027", "#f46d43", "#fdae61","#fee090","#e0f3f8","#4575b4","#4575b4"),
                        breaks=c("underage_i","v101", "v025", "v106", "v190", "b4", "v481","Residual"),
                        labels=c(paste(paste("Underage:",round(decomposition_pie$percent_t[1],1), sep=" "),"%", sep=""),
                                 paste(paste(paste(GEO,":",sep=""),round(decomposition_pie$percent_t[2],1), sep=" "),"%", sep=""),
                                 paste(paste("Urban/Rural:",round(decomposition_pie$percent_t[3],1), sep=" "),"%", sep=""),
                                 paste(paste("Maternal Education Level:",round(decomposition_pie$percent_t[4],1), sep=" "),"%", sep=""),
                                 paste(paste("Wealth Quintile:",round(decomposition_pie$percent_t[5],1), sep=" "),"%", sep=""),
                                 paste(paste("Sex of Child:",round(decomposition_pie$percent_t[6],1), sep=" "),"%", sep=""),
                                 paste(paste("Covered by Health Insurance:",round(decomposition_pie$percent_t[7],2), sep=" "),"%", sep=""),
                                 paste(paste("Unexplained Variation:",round(decomposition_pie$percent_t[8],1), sep=" "),"%", sep="")))
    
    output_list[[paste("pie_",i, sep="")]] <- assign(paste("pie_",i, sep=""),pie_i)
    
    total <- tibble(name_t = "Total",
                    cindex_t = "" ,
                    elastic_t = "",
                    contribution_t = sum(decomposition$contribution_t, na.rm=TRUE),
                    percent_t = sum(decomposition$percent_t,na.rm=TRUE))
    
    decomposition <- rbind.data.frame(decomposition,total)
    
    output[,1] <- decomposition$name_t
    output <- cbind.data.frame(output,round(decomposition$percent_t,digits=2))
    names(output)[names(output) == "round(decomposition$percent_t, digits = 2)"] <- paste("percent_",i, sep="")
    names(output)[1] <- "factors"
    
    #State-Level CIs
    CI_Results_GEO<- c()
    CI_E_Results_GEO <- c()
    AEG_Composite_Results_GEO <- c()
    Coverage_Results_GEO<- c()
  
    #Compute Concentration Indiices By Geographic Unit 
    for (k in 1:length(GEO_CI)){
      
      # Subset Data to be only data from GEO region outcome data is available
      data_k <- subset(data_i, (data_i$v101==k))
      
      # Use only that GEO region's data for equity calculation
      pred_probs_k <- subset(pred_probs, (data_i$v101==k))
      pred_probs_2_k <- subset(pred_probs_2, (data_i$v101==k))
      pred_probs_3_k <- subset(pred_probs_3, (data_i$v101==k))
      
      # Calculating the Direct Concentration Index
      direct_ci <- ci(y = data_k[,"outcome"], x = pred_probs_2_k$hci_du.response, type = "CI")
      CI_1 <- concentration_index(direct_ci)
      CI_1
      
      # Calculating the Erreygers Corrected Concentration Index 
      CIE <- ci(y = data_k[,"outcome"], x = pred_probs_2_k$hci_du.response, type = "CIc")
      CI_E <- concentration_index(CIE)
      CI_E
      
      # Calculating the Composite Absolute Equity Gap
      rank<-c(pred_probs_2_k$hci_du.response)
      quantiles<- quantile(rank, probs = seq(0, 1, 0.20))
      comp_data_k<- data.frame(cbind(data_k[,"outcome"],data_k[,"v005"],rank))
      names(comp_data_k)[names(comp_data_k) == "V1"] <- paste("",i, sep="")
      names(comp_data_k)[names(comp_data_k) == "V2"] <- paste("v005","", sep="")
      AEG_Composite <- round(weighted.mean(comp_data_k[,i][rank>=quantiles[5]],comp_data_k[,"v005"][rank>=quantiles[5]]) - weighted.mean(comp_data_k[,i][rank>=quantiles[2]],comp_data_k[,"v005"][rank>=quantiles[2]]), digits=3)
      
      # Calculating State-Level Coverage
      data_kc <- subset(data_k, data_k[,"underage_i"]==0)
      coverage_GEO<- c(weighted.mean(data_kc$outcome,data_kc$v005))
      
      #Storing Results
      CI_Results_GEO<- c(CI_Results_GEO, round(CI_1, digits=3))
      CI_E_Results_GEO<- c(CI_E_Results_GEO, round(CI_E, digits=3))
      AEG_Composite_Results_GEO<- c(AEG_Composite_Results_GEO, round(AEG_Composite, digits=3))
      Coverage_Results_GEO<- c(Coverage_Results_GEO, round(coverage_GEO, digits=3))
  
    }
    
    #Create Equity-Efficiency Plane
    invisible(capture.output(efficiency_data <- data_ic %>% 
                               group_by(data_ic$v024) %>% 
                               summarise(coverage=weighted.mean(outcome,v005),
                                         underage = weighted.mean(underage_i,v005),
                                         urban_rural = weighted.mean(v025,v005)*100/2,
                                         maternal_edu = weighted.mean(v106,v005)*100/3,
                                         wealth = weighted.mean(v190,v005)*100/5,
                                         sex = weighted.mean(b4,v005)*100/2,
                                         insurance = weighted.mean(v481,v005))))
    
        GEO_UNIT <- c(dhs_data$v101)
        GEO_NAMES <- c(dhs_data$geo_names)
        GEO_LABEL<- paste(GEO_UNIT, GEO_NAMES, sep=" = ")
    
    efficiency <- cbind.data.frame(CI_Results_GEO, GEO_UNIT, GEO_NAMES, GEO_LABEL, efficiency_data)
    efficiency$equity_axis <- abs(ifelse(efficiency$CI_Results_GEO>=0, 1 - efficiency$CI_Results_GEO, (efficiency$CI_Results_GEO*-1) - 1))
    
    etitle<- (ifelse(i=="ZERO","Zero-Dose Equity-Prevalence Plane",
                     ifelse(i=="FULL","Fully Immunized for Age Equity-Coverage Plane", 
                            ifelse(i=="COMPLETE","Completed Vaccination Schedule Equity-Coverage Plane",
                                   paste(i,"Vaccination Equity-Coverage Plane", sep=" ")))))
    
    efficiency_list[[paste("",i, sep="")]] <- assign(paste("",i, sep=""),efficiency)
    
    color_vector<-rep("#000000", length(GEO_UNIT))
    
    if (i=="ZERO"){
      efficiency_i<- ggplot(efficiency_list[[i]],aes(x= CI_Results_GEO,y=coverage*100, colour = GEO_LABEL)) +
        geom_point(alpha = .4,shape=20, color="blue3") +
        scale_y_continuous(name="Prevalence") +
        scale_x_continuous(name = "Inequality: Composite Index") +
        theme_classic(base_size = 16) +
        geom_label_repel(label=GEO_UNIT, size=2, max.overlaps = Inf) +
        ggtitle(etitle) +
        theme(axis.title.y = element_text(angle = 0)) +
        guides(colour = guide_legend(ncol = 3)) +
        theme(legend.position="bottom") +
        labs(colour = GEO) +
        scale_colour_manual(values=color_vector)
      
    } else{
      efficiency_i<- ggplot(efficiency_list[[i]],aes(x= equity_axis,y=coverage*100, colour = GEO_LABEL)) +
        geom_point(alpha = .4,shape=20, color="blue3") +
        scale_y_continuous(name="% Coverage") +
        scale_x_continuous(name = "Equity: (1 - Composite Index)") +
        theme_classic(base_size = 16) +
        geom_label_repel(label=GEO_UNIT, size=2, max.overlaps = Inf) +
        ggtitle(etitle)+
        theme(axis.title.y = element_text(angle = 0)) + 
        guides(colour = guide_legend(ncol = 3)) +
        theme(legend.position="bottom") +
        labs(colour = GEO) +
        scale_colour_manual(values=color_vector)
      
    }
    
    
    output_list[[paste("efficiency_",i, sep="")]] <- assign(paste("efficiency_",i, sep=""),efficiency_i)
    
    # Create Output
    CI_Results_GEO_Output<- cbind.data.frame(CI_Results_GEO_Output, CI_Results_GEO)
    CI_E_Results_GEO_Output<- cbind.data.frame(CI_E_Results_GEO_Output, CI_E_Results_GEO)
    AEG_Composite_Results_GEO_Output<- cbind.data.frame(AEG_Composite_Results_GEO_Output, AEG_Composite_Results_GEO)
    Coverage_Results_GEO_Output<- cbind.data.frame(Coverage_Results_GEO_Output, Coverage_Results_GEO)
    
    names(CI_Results_GEO_Output)[names(CI_Results_GEO_Output) == "CI_Results_GEO"] <- paste("CI_GEO_",i, sep="")
    names(CI_E_Results_GEO_Output)[names(CI_E_Results_GEO_Output) == "CI_E_Results_GEO"] <- paste("CI_E_GEO_",i, sep="")
    names(AEG_Composite_Results_GEO_Output)[names(AEG_Composite_Results_GEO_Output) == "AEG_Composite_Results_GEO"] <- paste("AEG_Composite_GEO_",i, sep="")
    names(Coverage_Results_GEO_Output)[names(Coverage_Results_GEO_Output) == "Coverage_Results_GEO"] <- paste("Coverage_GEO_",i, sep="")
    
    print(paste(i,"analysis completed", sep=" "))
  }
  
  
  CI_Results_GEO_Output[,1]<- GEO_CI
  CI_E_Results_GEO_Output[,1]<- GEO_CI
  AEG_Composite_Results_GEO_Output[,1]<- GEO_CI
  Coverage_Results_GEO_Output[,1]<- GEO_CI

  names(CI_Results_GEO_Output)[1] <-GEO
  names(CI_E_Results_GEO_Output)[1]<- GEO
  names(AEG_Composite_Results_GEO_Output)[1]<- GEO
  names(Coverage_Results_GEO_Output)[1]<- GEO
  
  #Create Output of Function
  output[,1] <- c("Underage",GEO,"Rural","Maternal Education","Wealth Index","Sex of Child","Health Insurance", "Residual","Total")
  
  output_list <- c(output_list, output)
  
  equity_national <- data.frame(VACCINES, Coverage_Results, CI_Results, CI_E_Results, AEG_Composite_Results)
  
  equity_state<- data.frame(Coverage_Results_GEO_Output,CI_Results_GEO_Output,CI_E_Results_GEO_Output,AEG_Composite_Results_GEO_Output)
  
  reference<-data.frame(FACTORS,REF)
  
  output_list <- c(output_list, equity_national, Coverage_Results_GEO_Output, CI_Results_GEO_Output, CI_E_Results_GEO_Output, AEG_Composite_Results_GEO_Output, reference)
  
  return(output_list)
}

#Run the code
results <- VERSE(DATA,COUNTRY,VACCINES,SCHEDULE,FACTORS,GEO)


