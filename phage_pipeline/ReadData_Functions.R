
## This script contains functions for reading in the csv files for all populations.
## These functions are used to prepare the csv files for all downstream analysis

##libraries
library("tidyverse")
library("reshape2")

## Reading in all possible PAM+protospacer regions, together with BIM targeted positions.
PAM_positions <- read_csv("../data/PAM_positions.csv")
BIM_positions <- readxl::read_xlsx("../data/Blastn_BIM_PAM_and_sequences.xlsx", sheet = 2)

## The BIM targeted regions are found based on BIM targeted position and the sign of the strand, i.e. "+" or "-"
BIM_positions$start<-ifelse(BIM_positions$`s. start`< BIM_positions$`s. end`,
                            BIM_positions$`s. start`-6, BIM_positions$`s. end`)
BIM_positions$end<-ifelse(BIM_positions$`s. start`>BIM_positions$`s. end`,
                          BIM_positions$`s. start`+7,BIM_positions$`s. end`+1)

## All escape regions, BIM targeted or not, are 37 nt long. 
## Mutations falling in these regions are Escape mutations. 

###########################################################

## Input:   A messy dataframe from either of the two functions: plotting_variants() or reading_csv()
## Output:  A tidy dataframe with extra columns needed for further analysis
tidying<-function(data){
  data<-data %>% 
    select(-1) %>% 
    mutate_at(c("TYPE", "TIME", "ALT", "REF"), as.character)
  data[data=="R1"] <- "T1" # In R7 the time points are called R1-R4 insted of T1-T4 (typo)
  data[data=="R2"] <- "T2"
  data[data=="R3"] <- "T3"
  data[data=="R4"] <- "T4"
  
  ## Adding "escape" = T/F, column to data, if position of variant falls in escape region
  data$escape <- ifelse(sapply(data$POS, function(p) 
    any(PAM_positions$escape_start <= p & PAM_positions$escape_end >= p)),TRUE, FALSE)
  
  ## Adding column with position in escape region: 1 <= escape_pos <= 37 
  ## If variant is in multiple escape regions, it gets the lowest escape_pos
  ## Whether it is on a "+" or "-" strand, is incoorporated
  escape_pos <- c()
  for (pos in data %>% .$POS) {
    if (sum(PAM_positions$escape_start <= pos & pos <= PAM_positions$escape_end) == 0) {
      escape_pos <- append(escape_pos,NA)
    }
    if (sum(PAM_positions$escape_start <= pos & pos <= PAM_positions$escape_end) != 0) {
      best_pos <- 37
      for (row in 1:nrow(PAM_positions[which(PAM_positions$escape_start <= pos & PAM_positions$escape_end >= pos),])){
        info <- PAM_positions[which(PAM_positions$escape_start <= pos & PAM_positions$escape_end >=pos)[row],]
        escape_position <- ifelse(info$strand == "-", pos-info$escape_start+1, info$escape_end-pos+1)
        if (escape_position < best_pos) {
          best_pos <- escape_position
        }
      }
      escape_pos <- append(escape_pos, best_pos)
    }
  }
  
  data <- data %>% mutate(escape_pos = escape_pos)
  
  ## Adding BIM_targeted = T/F, column to data, if position of variant falls in BIM targeted escape region
  data$BIM_targeted <- ifelse(sapply(data$POS, function(p) 
    any(BIM_positions$start <= p & BIM_positions$end >= p)),TRUE, FALSE)
  
  ## Adding BIMname column to data
  BIMname <- c()
  for (pos in data %>% .$POS) {
    if (sum(BIM_positions$start <= pos & pos <= BIM_positions$end) == 0) {
      BIMname <- append(BIMname, NA)
    }
    if (sum(BIM_positions$start <= pos & pos <= BIM_positions$end) != 0) {
      name <- BIM_positions$`#query acc.ver`[which(BIM_positions$start <= pos & pos <= BIM_positions$end)]
      BIMname <- append(BIMname, name)
    }
  }
  data <- data %>% mutate(BIMname = BIMname)
  
  
  return(data)
}

### Function used for almost all downstream analyses

## Input:   Loops through a list of file names, e.g. c("W1","W2","W3","W4","W5","W6","W7","W8"), and fetch csv's
## Output   A combined, big, tidy dataframe with all information for all samples in the list of files
##          Columns: 
##            AO:           Observed alternative allele count 
##            DP:           Read depth
##            TYPE:         Type of variant: snp, mnp, ins, del, complex
##            TIME:         Time point, T1,T2,T3,T4
##            ALT:          Alternative allele
##            REF:          Reference allele
##            POS:          Genomic position of variant
##            FREQ:         Frequency of variant
##            escape:       Escape or not, T/F
##            escape_pos:   Escape position, 1 <= pos <= 37
##            BIM_targeted: BIM targeted or not, T/F
##            BIMname:      Name of BIM if BIM_targeted == T
##            Sample:       Name of sample, W1-W8 or R1-R8
##            VarCaller     FreeBayes (originally added when trying out multiple variant callers)

reading_csv <- function(liste_of_files, path, method){
  Whole_df=tibble()
  for (j in liste_of_files){
    df<-read_tsv(paste(path,j,"_data.csv", sep=""))
    df=tidying(df)
    df=df %>% mutate(Sample=j) %>% mutate(VarCaller=method)
    Whole_df<-rbind(Whole_df,df)
  }
  return(Whole_df)
}