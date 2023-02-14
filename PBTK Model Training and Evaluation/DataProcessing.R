
###############################################################################
# This function pre-processes kinetic data                
# Authors: Anja Lehmann, Christoph Hethey                                                       
# Date: 2023-02-13                                            
###############################################################################

GetDataProcessing <- function(dataset,deagg) {

  # set seed for reproducibility of random processes
  set.seed(10)
  
  #--- 'DE-AGGREGATION' OF SUMMARY DATA -----------------------------------------
  
  if(deagg==TRUE){
    
    dataset_summarydata <- subset(dataset, dataset$YDESCRIPTION == "mean" & !is.na(dataset$YSD))
    dataset_without_summarydata <- subset(dataset, dataset$YDESCRIPTION !="mean" | is.na(dataset$YDESCRIPTION) )
    
    for(i in seq_along(dataset_summarydata$ID)) {
      
      if(dataset_summarydata$YSD[i]==0){
        
        # replicate data with YSD = 0
        row               <- dataset_summarydata[i,]
        row_rep           <- row[rep(seq_len(nrow(row)), each = N), ]
        row_rep$YDESCRIPTION <- rep(("deaggregated"), each = N)
        
        # bind row_rep to dataset_without_summarydata
        dataset_without_summarydata <- rbind(dataset_without_summarydata,row_rep)
      }
      
      else{
      ## data ID
      #ID <- dataset_summarydata$ID[i]
      
      # sample mean E
      E  <- dataset_summarydata$YORIGINAL[i]
      
      # standard deviation of the sample 
      sd <- dataset_summarydata$YSD[i]
      
      # sample size
      N  <- dataset_summarydata$N[i]
      
      # variance of cooresponding normal distribution
      sigma_square <- log(((sd^2)/(E^2)+1))
      
      # mean of cooresponding normal distribution 
      mu        <- log(E)-sigma_square/2
      
      # sample
      y <- rnorm(N)*sqrt(sigma_square)+mu
        
      # scale
      y <- y/sd(y)*sqrt(sigma_square)
      
      # center
      y <- y-mean(y)+mu
      
      # transform
      y <- exp(y)
      
      # replace summary data in dataset by sampled data
      row               <- dataset_summarydata[i,]
      row_rep           <- row[rep(seq_len(nrow(row)), each = N), ]
      row_rep$YORIGINAL <- y
      row_rep$YDESCRIPTION <- rep(("deaggregated"), each = N)
      
      # bind row_rep to dataset_without_summarydata
      dataset_without_summarydata <- rbind(dataset_without_summarydata,row_rep)
      }
      

    }
    dataset <- dataset_without_summarydata
    
  }
  else{}
  
  #############################################################################
  #transform XORIGINAL and XORIGINALUNIT into h                                                             
  #############################################################################

  #create new colums X and XUNIT with NA and length as column XORIGINAL or XORIGINALUNIT
  dataset$X <- rep_len(NA,length(dataset$XORIGINAL))
  dataset$XUNIT <- rep_len(NA,length(dataset$XORIGINALUNIT))

  #create vector is_... that shows no of rows for which XORIGINALUNIT = "..."
  is_week     <- which((dataset$XORIGINALUNIT=="week"),TRUE)
  is_rest     <- which(!((dataset$XORIGINALUNIT=="week")))

  #for rows is_week multiplicate XORIGINAL with 24*7 and give results into column X
  dataset$X[is_week]     <- dataset$XORIGINAL[is_week]*24*7
  dataset$XUNIT[is_week] <- "h"

  #copy rest
  dataset$X[is_rest]     <- dataset$XORIGINAL[is_rest]
  dataset$XUNIT[is_rest] <- dataset$XORIGINALUNIT[is_rest]
  
  #############################################################################
  #transform YORIGINAL and YORIGINALUNIT                                       
  #############################################################################

  #create new colums Y and YUNIT with NA and length as column YORIGNAL and YORIGINALUNIT
  dataset$Y <- rep_len(NA, length(dataset$YORIGINAL))
  dataset$YUNIT <- rep_len(NA, length(dataset$YORIGINALUNIT))

  #create vectors is_... which show rownumber for which ID==...
  is_White_1977_retrorsine_rat_bile          <- which((dataset$ID=="White_1977_retrorsine_rat_bile"),TRUE)
  is_Chu_1991_retrorsine_mouse_urine         <- which((dataset$ID=="Chu_1991_retrorsine_mouse_urine"),TRUE)
  is_Chu_1991_retrorsine_rat_urine           <- which((dataset$ID=="Chu_1991_retrorsine_rat_urine"),TRUE)
  is_Yang_2017_retrorsine_mouse_plasma       <- which((dataset$ID=="Yang_2017_retrorsine_mouse_plasma"),TRUE)
  is_Yang_2017_proteinadducts_mouse_liver    <- which((dataset$ID=="Yang_2017_proteinadducts_mouse_liver"),TRUE)
  is_Yang_2017_GSHconjugates_mouse_liver     <- which((dataset$ID=="Yang_2017_GSHconjugates_mouse_liver"),TRUE)
  is_Zhu_2017_DNAadducts_mouse_liver         <- which((dataset$ID=="Zhu_2017_DNAadducts_mouse_liver"),TRUE)
  is_Yang_2018_retrorsine_mouse_plasma       <- which((dataset$ID=="Yang_2018_retrorsine_mouse_plasma"),TRUE)
  is_Yang_2018_retrorsine_mouse_liver        <- which((dataset$ID=="Yang_2018_retrorsine_mouse_liver"),TRUE)
  is_Yang_2018_GSHconjugates_mouse_liver     <- which((dataset$ID=="Yang_2018_GSHconjugates_mouse_liver"),TRUE)
  is_Zhu_2017_DNAadducts_mouse_liver_10      <- which((dataset$ID=="Zhu_2017_DNAadducts_mouse_liver_10"),TRUE)
  is_Zhu_2017_DNAadducts_mouse_liver_20      <- which((dataset$ID=="Zhu_2017_DNAadducts_mouse_liver_20"),TRUE)
  is_Zhu_2017_DNAadducts_mouse_liver_40      <- which((dataset$ID=="Zhu_2017_DNAadducts_mouse_liver_40"),TRUE)
  is_Zhu_2017_DNAadducts_mouse_liver_60      <- which((dataset$ID=="Zhu_2017_DNAadducts_mouse_liver_60"),TRUE)
  is_Yang_2018_proteinadducts_mouse_liver    <- which((dataset$ID=="Yang_2018_proteinadducts_mouse_liver"),TRUE)
  is_Li_2022_retrorsine_mouse_serum          <- which((dataset$ID=="Li_2022_retrorsine_mouse_serum"),TRUE)
  is_Li_2022b_retrorsine_mouse_serum         <- which((dataset$ID=="Li_2022b_retrorsine_mouse_serum"),TRUE)
  
  
  is_rest                     <- which(!((dataset$ID=="White_1977_retrorsine_rat_bile")|
                                         (dataset$ID=="Chu_1991_retrorsine_mouse_urine")|
                                         (dataset$ID=="Chu_1991_retrorsine_rat_urine")|
                                         (dataset$ID=="Yang_2017_retrorsine_mouse_plasma")|
                                         (dataset$ID=="Yang_2017_proteinadducts_mouse_liver")|
                                         (dataset$ID=="Yang_2017_GSHconjugates_mouse_liver")|
                                         (dataset$ID=="Zhu_2017_DNAadducts_mouse_liver")|
                                         (dataset$ID=="Yang_2018_retrorsine_mouse_plasma")|
                                         (dataset$ID=="Yang_2018_retrorsine_mouse_liver")|
                                         (dataset$ID=="Yang_2018_GSHconjugates_mouse_liver")|
                                         (dataset$ID=="Zhu_2017_DNAadducts_mouse_liver_10")|
                                         (dataset$ID=="Zhu_2017_DNAadducts_mouse_liver_20")|
                                         (dataset$ID=="Zhu_2017_DNAadducts_mouse_liver_40")|
                                         (dataset$ID=="Zhu_2017_DNAadducts_mouse_liver_60")|
                                         (dataset$ID=="Yang_2018_proteinadducts_mouse_liver")|
                                         (dataset$ID=="Li_2022_retrorsine_mouse_serum")|
                                         (dataset$ID=="Li_2022b_retrorsine_mouse_serum")
                                           
                                         ))

  # Convert RET bile, uri from cumulativepercentofdose to nmol 
  bw_rat <- 0.25  #kg       (Brown et al. 1997)
  D      <- 40    #mg/kg bw (study-specific)
  M      <- 351.4 #g/mol    (NCBI 2021)  
  dataset$Y[is_White_1977_retrorsine_rat_bile]     <- dataset$YORIGINAL[is_White_1977_retrorsine_rat_bile]*(D*bw_rat)/100/10^3/M*10^9
  dataset$YUNIT[is_White_1977_retrorsine_rat_bile] <- "nmol"

  bw_mouse <- 0.025 #kg       (Brown et al. 1997)
  D        <- 25    #mg/kg bw (study-specific)
  M        <- 351.4 #g/mol    (NCBI 2021)
  dataset$Y[is_Chu_1991_retrorsine_mouse_urine]     <- dataset$YORIGINAL[is_Chu_1991_retrorsine_mouse_urine]*(D*bw_mouse)/100/10^3/M*10^9
  dataset$YUNIT[is_Chu_1991_retrorsine_mouse_urine] <- "nmol"
  
  bw_rat <- 0.25  #kg       (Brown et al. 1997)
  D      <- 25    #mg/kg bw (study-specific)
  M      <- 351.4 #g/mol    (NCBI 2021)
  dataset$Y[is_Chu_1991_retrorsine_rat_urine]     <- dataset$YORIGINAL[is_Chu_1991_retrorsine_rat_urine]*(D*bw_rat)/100/10^3/M*10^9
  dataset$YUNIT[is_Chu_1991_retrorsine_rat_urine] <- "nmol"

  # Convert RET liv, pla from µg/mL to nmol
  Vpla <- 6.124070e-04 #L equals Vblo(1-hct) (Baskurt et al. 2007, Diehl et al. 2001)
  dataset$Y[is_Yang_2017_retrorsine_mouse_plasma] <- dataset$YORIGINAL[is_Yang_2017_retrorsine_mouse_plasma]*Vpla*1000/10^6/M*10^9
  dataset$YUNIT[is_Yang_2017_retrorsine_mouse_plasma] <- "nmol"
  
  dataset$Y[is_Yang_2018_retrorsine_mouse_plasma] <- dataset$YORIGINAL[is_Yang_2018_retrorsine_mouse_plasma]*Vpla*1000/10^6/M*10^9
  dataset$YUNIT[is_Yang_2018_retrorsine_mouse_plasma] <- "nmol"
  
  dataset$Y[is_Li_2022_retrorsine_mouse_serum] <- dataset$YORIGINAL[is_Li_2022_retrorsine_mouse_serum]*Vpla*1000/10^6/M*10^9
  dataset$YUNIT[is_Li_2022_retrorsine_mouse_serum] <- "nmol"
  
  dataset$Y[is_Li_2022b_retrorsine_mouse_serum] <- dataset$YORIGINAL[is_Li_2022b_retrorsine_mouse_serum]*Vpla*1000/10^6/M*10^9
  dataset$YUNIT[is_Li_2022b_retrorsine_mouse_serum] <- "nmol"
  
  Vliv <- 1.372500e-03 #L (Brown et al. 1997)
  DF   <- 11           #- (derived from Yang et al. 2018)
  dataset$Y[is_Yang_2018_retrorsine_mouse_liver] <- dataset$YORIGINAL[is_Yang_2018_retrorsine_mouse_liver]*Vliv*1000*DF/10^6/M*10^9
  dataset$YUNIT[is_Yang_2018_retrorsine_mouse_liver] <- "nmol"
  
  #convert DHR:PROT liv from AanalyteperAispermgprotein to nmol
  Wliv     <- 1.3725*10^3 #mg            (Brown et al. 1997)
  Cprotliv <- 23.84       #mg/100 mg liv (Verma and Chakraborty 2008)
  WIS      <- 1.25*10^-8  #g             (personal correspondance X. Yang)
  MIS      <- 284.74      #g/mol         (NCBI 2021)  
  Rprot    <- 0.787       #ratio of protein amount of supernatant 
                          #to liver homogenate (personal correspondance)
  
  dataset$Y[is_Yang_2017_proteinadducts_mouse_liver] <- dataset$YORIGINAL[is_Yang_2017_proteinadducts_mouse_liver]*Cprotliv/100*Wliv*WIS/MIS*Rprot*10^9
  dataset$YUNIT[is_Yang_2017_proteinadducts_mouse_liver] <- "nmol"

  dataset$Y[is_Yang_2018_proteinadducts_mouse_liver] <- dataset$YORIGINAL[is_Yang_2018_proteinadducts_mouse_liver]*Cprotliv/100*Wliv*WIS/MIS*Rprot*10^9
  dataset$YUNIT[is_Yang_2018_proteinadducts_mouse_liver] <- "nmol"
   
  # Convert DHR:GSH liv from AanalyteperAispermgprotein to nmol 
  Wliv     <- 1.3725*10^3 #mg            (Brown et al. 1997)
  Cprotliv <- 23.84       #mg/100 mg liv (Verma and Chakraborty 2008)
  WIS      <- 1.25*10^-8   #g             (personal correspondance X. Yang)
  MIS      <- 391.5       #g/mol         (NCBI 2021)  
  Rprot    <- 0.787       #ratio of protein amount of supernatant 
                          #to liver homogenate (personal correspondance X. Yang)
  dataset$Y[is_Yang_2017_GSHconjugates_mouse_liver] <- dataset$YORIGINAL[is_Yang_2017_GSHconjugates_mouse_liver]*Cprotliv/100*Wliv*WIS/MIS*Rprot*10^9
  dataset$YUNIT[is_Yang_2017_GSHconjugates_mouse_liver] <- "nmol"

  WIS      <- 1.6*10^-8  #g             (personal correspondance X. Yang)
  dataset$Y[is_Yang_2018_GSHconjugates_mouse_liver] <- dataset$YORIGINAL[is_Yang_2018_GSHconjugates_mouse_liver]*Cprotliv/100*Wliv*WIS/MIS*Rprot*10^9
  dataset$YUNIT[is_Yang_2018_GSHconjugates_mouse_liver] <- "nmol"
  
  #Convert DHR:DNA from adducts/10^8 nucleotides to adducts/liver (nmol)
  Na       <- 6.02214076*10^23 #1/mol                         (IUPAC 1997)
  SFgenome <- 2728222451       #base pairs/mouse haploid cell (NCBI 2021) 
  SFliv    <- 128*10^6         #cells/g liver                 (Ring et al. 2011) 
  Wliv     <- 1.3725           #g                             (Brown et al. 1997)
  #2       (diploid cell)
  #2       (base pairs)
  #10^9    (mol to nmol)
  dataset$Y[is_Zhu_2017_DNAadducts_mouse_liver]        <- dataset$YORIGINAL[is_Zhu_2017_DNAadducts_mouse_liver]*SFgenome*2*2/10^8*SFliv*Wliv/Na*10^9 
  dataset$YUNIT[is_Zhu_2017_DNAadducts_mouse_liver]    <- "nmol"
  
  dataset$Y[is_Zhu_2017_DNAadducts_mouse_liver_10]     <- dataset$YORIGINAL[is_Zhu_2017_DNAadducts_mouse_liver_10]*SFgenome*2*2/10^8*SFliv*Wliv/Na*10^9 
  dataset$YUNIT[is_Zhu_2017_DNAadducts_mouse_liver_10] <- "nmol"
  
  dataset$Y[is_Zhu_2017_DNAadducts_mouse_liver_20]     <- dataset$YORIGINAL[is_Zhu_2017_DNAadducts_mouse_liver_20]*SFgenome*2*2/10^8*SFliv*Wliv/Na*10^9 
  dataset$YUNIT[is_Zhu_2017_DNAadducts_mouse_liver_20] <- "nmol"
  
  dataset$Y[is_Zhu_2017_DNAadducts_mouse_liver_40]     <- dataset$YORIGINAL[is_Zhu_2017_DNAadducts_mouse_liver_40]*SFgenome*2*2/10^8*SFliv*Wliv/Na*10^9 
  dataset$YUNIT[is_Zhu_2017_DNAadducts_mouse_liver_40] <- "nmol"
  
  dataset$Y[is_Zhu_2017_DNAadducts_mouse_liver_60]     <- dataset$YORIGINAL[is_Zhu_2017_DNAadducts_mouse_liver_60]*SFgenome*2*2/10^8*SFliv*Wliv/Na*10^9 
  dataset$YUNIT[is_Zhu_2017_DNAadducts_mouse_liver_60] <- "nmol"

  #copy rest
  dataset$Y[is_rest]     <- dataset$YORIGINAL[is_rest]
  dataset$YUNIT[is_rest] <- dataset$YORIGINALUNIT[is_rest]


  #############################################################################
  #transform LOQORIGINAL and LODORIGINAL                                 
  #############################################################################
  
  #create new colums 
  dataset$LOQ     <- rep_len(NA,length(dataset$LOQORIGINAL))
  dataset$LOQUNIT <- rep_len(NA,length(dataset$LOQORIGINALUNIT))
                                        
  dataset$LOD     <- rep_len(NA,length(dataset$LODORIGINAL))
  dataset$LODUNIT <- rep_len(NA,length(dataset$LODORIGINALUNIT))  
  
  #convert Yang_2017_retrorsine_mouse_plasma LOQORGINAL and LOQORIGINAL from µg/mL to nmol
  Vpla <- 6.124070e-04*1000 #mL
  dataset$LOQ[is_Yang_2017_retrorsine_mouse_plasma] <- dataset$LOQORIGINAL[is_Yang_2017_retrorsine_mouse_plasma]*Vpla/10^6/351.399*10^9
  dataset$LOQUNIT[is_Yang_2017_retrorsine_mouse_plasma] <- "nmol"
  
  #convert Yang_2017_retrorsine_mouse_plasma LODORGINAL and LODORIGINAL from µg/mL to nmol
  Vpla <- 6.124070e-04*1000 #mL
  dataset$LOD[is_Yang_2017_retrorsine_mouse_plasma] <- dataset$LODORIGINAL[is_Yang_2017_retrorsine_mouse_plasma]*Vpla/10^6/351.399*10^9
  dataset$LODUNIT[is_Yang_2017_retrorsine_mouse_plasma] <- "nmol"
  
  #############################################################################
  #transform DOSEORIGINAL and DOSEORIGINALUNIT                                 
  #############################################################################

  #create new colums DOSE and DOSEUNIT with NA and length as column DOSEORIGINAL or DOSEORIGINALUNIT
  dataset$DOSE     <- rep_len(NA,length(dataset$DOSEORIGINAL))
  dataset$DOSEUNIT <- rep_len(NA,length(dataset$DOSEORIGINALUNIT))

  #create vector is_mgperkg that shows no of rows for which DOSEORIGINALUNIT = ...
  is_mgperkg   <- which((dataset$DOSEORIGINALUNIT=="mgperkg"),TRUE)
  is_umolperkg <- which((dataset$DOSEORIGINALUNIT=="umolperkg"),TRUE)
  is_rest      <- which(!((dataset$DOSEORIGINALUNIT=="mgperkg")|(dataset$DOSEORIGINALUNIT=="umolperkg")))

  #covert is_mgperkg from mg/kg bw to nmol/kg bw
  dataset$DOSE[is_mgperkg]     <- dataset$DOSEORIGINAL[is_mgperkg]/10^3/351.399*10^9
  dataset$DOSEUNIT[is_mgperkg] <- "nmolperkg"

  #covert is_umolperkg from µmol/kg to nmol/kg  
  dataset$DOSE[is_umolperkg]     <- dataset$DOSEORIGINAL[is_umolperkg]*10^3
  dataset$DOSEUNIT[is_umolperkg] <- "nmolperkg"

  #copy rest
  dataset$DOSE[is_rest]     <- dataset$DOSEORIGINAL[is_rest]
  dataset$DOSEUNIT[is_rest] <- dataset$DOSEORIGINALUNIT[is_rest]

return(dataset)

}






