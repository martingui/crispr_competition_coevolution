import numpy as np
import pandas as pd
import vcf

## Defining classes
vcf_file_list = []

class ImportVCF2:
    '''
    Imports VCF experimental file:
    - Full records: selfdico
    - Index of contents: self.ctrl_index (type = list)
    '''

    def __init__(self):#, name):
        # self.name = name
        self.index = []
        self.dico = {}
        self.CHROM = {}
        self.POS = {}
        self.ID = {}
        self.REF = {}
        self.ALT = {}
        self.QUAL = {}
        self.FILTER = {}
        self.INFO = {}
        self.FORMAT = {}

    def load_vcf(self, file_name):
        vcf_file_list.append(file_name)
        vcf_reader = vcf.Reader(open(file_name, 'r'))
        for record in vcf_reader:
            # print(record)
            self.index.append( record.CHROM + ';' + str(record.POS) )
            self.dico[record.CHROM + ';' + str(record.POS)]=record
            self.CHROM[record.CHROM + ';' + str(record.POS)]=record.CHROM
            self.POS[record.CHROM + ';' + str(record.POS)]=record.POS
            self.ID[record.CHROM + ';' + str(record.POS)]=record.ID
            self.ALT[record.CHROM + ';' + str(record.POS)]=record.ALT
            self.REF[record.CHROM + ';' + str(record.POS)]=record.REF
            self.QUAL[record.CHROM + ';' + str(record.POS)]=record.QUAL
            self.FILTER[record.CHROM + ';' + str(record.POS)]=record.FILTER
            self.INFO[record.CHROM + ';' + str(record.POS)]=record.INFO
            self.FORMAT[record.CHROM + ';' + str(record.POS)]=record.FORMAT

class VCFdealer:
    '''
    Overall function of this class is to deal with VCF file objects. 

    General order of the different methods:
    
    make_summary_df -- Makes a Pandas DataFrame with the wanted information needed to 
    process the information in the VCF file even further, e.g. using R-Studio. 
    '''

    def __init__(self, vcf_object, time):
        
        self.vcf_object = vcf_object
        self.time = time
        self.summary_df = self.make_summary_df()
        
    def make_summary_df(self):
        '''
        Makes a Pandas DataFrame with the columns: 
        -- POS: The genomic position of the variant
        -- TIME: T1, T2, T3 or T4
        -- ALT: Alternative allele
        -- REF: Reference allele
        -- AO: Alternate allele observation count
        -- DP: Read depth
        -- TYPE: "snp", "mnp", ins", "del" or "complex"
        -- FREQ: The allele frequency, filter >= 0.05
        '''
        
        ALT_dict = self.vcf_object.ALT
        
        REF_dict = self.vcf_object.REF
        for (key, value) in REF_dict.items():
            REF_dict[key] = [value]

        INFO_dict = self.vcf_object.INFO
        
        Info_list = ["AO","DP","TYPE"] 

        summary_dic = {}

        for (name, dic) in INFO_dict.items():
            pos = name.rsplit(';', 1)[1]
            summary_dic[pos] = {}

            for key in dic:
                if key in Info_list:
                    summary_dic[pos][key] = dic[key]
            
            ## FreeBayes call "mnp" (multi nucleotide polymorphisms) when multiple mutations
            ## are called in same (or close) positions. These positions are dealt with by adding
            ## i extra positions, e.g. pos-i, in the position column, for all i extra mutations
            ## in one position. 
            for i in range(len(ALT_dict[name])):
                if i == 0:
                    summary_dic[pos]["TIME"] = self.time
                    summary_dic[pos]["ALT"] = ALT_dict[name][i]
                    summary_dic[pos]["AO"] = INFO_dict[name]["AO"][i]
                    summary_dic[pos]["DP"] = INFO_dict[name]["DP"]
                    summary_dic[pos]["TYPE"] = INFO_dict[name]["TYPE"][i]
                
                else:
                    summary_dic[pos+"-"+str(i)] = {}
                    summary_dic[pos+"-"+str(i)]["TIME"] = self.time
                    summary_dic[pos+"-"+str(i)]["AO"] = INFO_dict[name]["AO"][i]
                    summary_dic[pos+"-"+str(i)]["DP"] = INFO_dict[name]["DP"]
                    summary_dic[pos+"-"+str(i)]["TYPE"] = INFO_dict[name]["TYPE"][i]
                    summary_dic[pos+"-"+str(i)]["ALT"] = ALT_dict[name][i]

            for i in range(len(ALT_dict[name])):
                if i == 0:
                    summary_dic[pos]["REF"] = REF_dict[name][0]
                elif i != 0 and len(REF_dict[name]) == 1:
                    summary_dic[pos+"-"+str(i)]["REF"] = REF_dict[name][0]
                else:
                    summary_dic[pos+"-"+str(i)]["REF"] = REF_dict[name][i][0] 
        
        summary_df = pd.DataFrame.from_dict(summary_dic, orient = "index")
        summary_df["POS"]=summary_df.index
        summary_df['POS'] = summary_df['POS'].str.split('-').str[0] # remove "-" in positions
        summary_df['FREQ'] = summary_df['AO']/summary_df['DP'] # Adding column with calculated allele frequency
        summary_df = summary_df[(summary_df['FREQ']>=0.025)] # Only keeping variants with allele frequencies equal to or above 0.025 
    
        return summary_df



## Imports and runs VCFdealer on multiple VCF file, outputs a csv file.
## list_files should contain the name of the files that should be in 
## same csv file, e.g. all W1 files. 

def run_on_multiple_files(list_files, path_to_vcfs, out_path, out_name):
    all_summary_df = pd.DataFrame() 

    for file in list_files:
        time = file[2:4]
        sample_indexes_list = []
        temp = ImportVCF2()
        temp.load_vcf(path_to_vcfs+file)

        value = VCFdealer(temp, time) 
        
        all_summary_df = pd.concat([all_summary_df,value.summary_df])
    
    all_summary_df.to_csv(out_path+out_name+".csv", sep = "\t")

## Defining arguments for the function: run_on_mutilple_files()
## W files

w1 =['W1T1_S18.vcf', 'W1T2_S35.vcf', 'W1T3_S51.vcf', 'W1T4_S67.vcf']
w2 =['W2T1_S19.vcf', 'W2T2_S36.vcf', 'W2T3_S52.vcf', 'W2T4_S68.vcf']
w3 =['W3T1_S20.vcf', 'W3T2_S37.vcf', 'W3T3_S53.vcf', 'W3T4_S69.vcf']
w4 =['W4T1_S21.vcf', 'W4T2_S38.vcf', 'W4T3_S54.vcf', 'W4T4_S70.vcf']
w5 =['W5T1_S22.vcf', 'W5T2_S39.vcf', 'W5T3_S55.vcf', 'W5T4_S71.vcf']
w6 =['W6T1_S23.vcf', 'W6T2_S40.vcf', 'W6T3_S56.vcf', 'W6T4_S72.vcf']
w7 =['W7T1_S24.vcf', 'W7T2_S41.vcf', 'W7T3_S57.vcf', 'W7T4_S73.vcf']
w8 =['W8T1_S26.vcf', 'W8T2_S42.vcf', 'W8T3_S58.vcf', 'W8T4-B_S92.vcf']

path_to_vcfs = 'steps/snpcall_freebayes/W_seq/'

outpath = 'steps/csv_freebayes/W_seq/'

run_on_multiple_files(list_files = w1, 
                      path_to_vcfs = path_to_vcfs, 
                      out_path = outpath,
                      out_name = "W1_data")

run_on_multiple_files(list_files = w2, 
                      path_to_vcfs = path_to_vcfs, 
                      out_path = outpath,
                      out_name = "W2_data")

run_on_multiple_files(list_files = w3, 
                      path_to_vcfs = path_to_vcfs, 
                      out_path = outpath,
                      out_name = "W3_data")

run_on_multiple_files(list_files = w4, 
                      path_to_vcfs = path_to_vcfs, 
                      out_path = outpath,
                      out_name = "W4_data")

run_on_multiple_files(list_files = w5, 
                      path_to_vcfs = path_to_vcfs, 
                      out_path = outpath,
                      out_name = "W5_data")

run_on_multiple_files(list_files = w6, 
                      path_to_vcfs = path_to_vcfs, 
                      out_path = outpath,
                      out_name = "W6_data")

run_on_multiple_files(list_files = w7, 
                      path_to_vcfs = path_to_vcfs, 
                      out_path = outpath,
                      out_name = "W7_data")

run_on_multiple_files(list_files = w8, 
                      path_to_vcfs = path_to_vcfs, 
                      out_path = outpath,
                      out_name = "W8_data")

## Defining arguments for the function: run_on_mutilple_files()
## R files

r1 = ['R1T1_S27.vcf', 'R1T2_S43.vcf', 'R1T3_S59.vcf', 'R1T4_S75.vcf']
r2 = ['R2T1_S28.vcf', 'R2T2_S44.vcf', 'R2T3_S60.vcf', 'R2T4_S76.vcf']
r3 = ['R3T1_S29.vcf', 'R3T2_S45.vcf', 'R3T3_S61.vcf', 'R3T4_S77.vcf']
r4 = ['R4T1_S30.vcf', 'R4T2_S46.vcf', 'R4T3_S62.vcf', 'R4T4_S78.vcf']
r5 = ['R5T1_S31.vcf', 'R5T2_S47.vcf', 'R5T3_S63.vcf', 'R5T4_S79.vcf']
r6 = ['R6T1_S32.vcf', 'R6T2_S48.vcf', 'R6T3_S64.vcf', 'R6T4_S80.vcf']
r7 = ['R7R1_S33.vcf', 'R7R2_S49.vcf', 'R7R3_S65.vcf', 'R7R4_S81.vcf']
r8 = ['R8T1_S34.vcf', 'R8T2_S50.vcf', 'R8T3_S66.vcf', 'R8T4_S82.vcf']


path_to_vcfs = 'steps/snpcall_freebayes/R_seq/'

outpath = 'steps/csv_freebayes/R_seq/'

run_on_multiple_files(list_files = r1, 
                      path_to_vcfs = path_to_vcfs, 
                      out_path = outpath,
                      out_name = "R1_data")

run_on_multiple_files(list_files = r2, 
                      path_to_vcfs = path_to_vcfs, 
                      out_path = outpath,
                      out_name = "R2_data")

run_on_multiple_files(list_files = r3, 
                      path_to_vcfs = path_to_vcfs, 
                      out_path = outpath,
                      out_name = "R3_data")

run_on_multiple_files(list_files = r4, 
                      path_to_vcfs = path_to_vcfs, 
                      out_path = outpath,
                      out_name = "R4_data")

run_on_multiple_files(list_files = r5, 
                      path_to_vcfs = path_to_vcfs, 
                      out_path = outpath,
                      out_name = "R5_data")

run_on_multiple_files(list_files = r6, 
                      path_to_vcfs = path_to_vcfs, 
                      out_path = outpath,
                      out_name = "R6_data")

run_on_multiple_files(list_files = r7, 
                      path_to_vcfs = path_to_vcfs, 
                      out_path = outpath,
                      out_name = "R7_data")

run_on_multiple_files(list_files = r8, 
                      path_to_vcfs = path_to_vcfs, 
                      out_path = outpath,
                      out_name = "R8_data")

f = open("steps/csv_freebayes/csv_done.txt", "w")
f.write("csv making is done.\n")
f.close

