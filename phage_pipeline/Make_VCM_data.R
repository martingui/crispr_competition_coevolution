source("ReadData_functions.R")
## Function that makes the input for VeryCoolMethod

## Input:     samples     A list of populations to loop through
##            experiment  "W" or "R", then it knows what T0 data to use

## Output:    mutation_time_series data
##              - id    ALT:POS, the id of the variant, a combination of Alternative allele and its position
##              - site  The position 
##              - f_i    the frequency of the variant at generation i (if 5 generations then 5 columns: f1,f2,f3,f4,f5)
##              - c_i    The coverage of the variant at generation i (same idea as f_i, but for coverage)

##            The output generated is located in the folder "Steps/VCM_data/Input/"
Make_VCM_data <- function(samples, experiment) {
  
  if (experiment == "W") {
    T0<- reading_csv("T0_W",
                     path = "../data/FreeBayes/Other/", 
                     method = "FreeBayes")
  }
  if (experiment == "R") {
    T0<- reading_csv("T0_R",
                     path = "../data/FreeBayes/Other/", 
                     method = "FreeBayes")
  }
  
  T0$TIME = "T0"

  for (s in 1:length(samples)) {
    data <- reading_csv(c(samples[s]),
                        paste("../data/FreeBayes/",experiment,"_seq/", sep = ""), 
                        method = "FreeBayes")
    
    
    data <- bind_rows(T0,data)
    
    data$id <- paste(data$ALT, data$POS, sep=":")
    
    Freq <- data %>%
      select(id, site = POS, TIME, FREQ) %>%
      spread(key = TIME, value = FREQ) %>%
      rename(f_0 = T0, f_1 = T1, f_2 = T2, f_3 = T3, f_4 = T4)
    
    Cov <- data %>%
      select(id, site = POS, TIME, DP) %>%
      spread(key = TIME, value = DP) %>%
      rename(c_0 = T0, c_1 = T1, c_2 = T2, c_3 = T3, c_4 = T4)
    
    data <- bind_cols(Freq, Cov[3:7])
    data[is.na(data)] <- 0
    
    write_tsv(data,paste("../Steps/VCM_data/Input/",samples[s],".tsv", sep = ""))
  }
  
}


## Actually make the data

R <- c("R1","R2","R3","R4","R5","R6","R7","R8")
W <- c("W1","W2","W3","W4","W5","W6","W7","W8")

Make_VCM_data(R, "R")
Make_VCM_data(W, "W")
