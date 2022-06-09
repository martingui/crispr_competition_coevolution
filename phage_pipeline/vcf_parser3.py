#! /home/enrique/envs/biopython/bin/python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import vcf


vcf_file = '/mnt/alpha_raid/work/Gandon/phages/test_freebayes/results/W5T4_S71_phag.vcf'
ctrl_file = '/mnt/alpha_raid/work/Gandon/phages/test_freebayes/results/TO-WT_S83_phag.vcf'
vcf_file_list = []



class Dog:
    '''Example class to test'''

    def __init__(self):#, name):
        # self.name = name
        self.tricks = []    # creates a new empty list for each dog

    def add_trick(self, trick):
        self.tricks.append(trick)


class ImportVCF:
    '''
    *** DEPRECATED: TO BE REMOVED ON THE NEXT VERSION ***
    *** ImportVCF2 will replace it ***
    Imports VCF experimental file:
    - Full records: self.ctrl_dico
    - Index of contents: self.ctrl_index (type = list)
    '''

    def __init__(self):#, name):
        # self.name = name
        self.vcf_index = []
        self.vcf_dico = {}
        self.vcf_CHROM = {}
        self.vcf_POS = {}
        self.vcf_ALT = {}
        self.vcf_REF = {}

    def load_vcf(self, file_name):
        vcf_file_list.append(file_name)
        vcf_reader = vcf.Reader(open(file_name, 'r'))
        for record in vcf_reader:
            # print(record)
            self.vcf_index.append( record.CHROM + ';' + str(record.POS) )
            self.vcf_dico[record.CHROM + ';' + str(record.POS)]=record
            self.vcf_CHROM[record.CHROM + ';' + str(record.POS)]=record.CHROM
            self.vcf_POS[record.CHROM + ';' + str(record.POS)]=record.POS
            self.vcf_ALT[record.CHROM + ';' + str(record.POS)]=record.ALT
            self.vcf_REF[record.CHROM + ';' + str(record.POS)]=record.REF


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


class ImportControl:
    '''
    *** THIS CLASS WILL DISAPEAR ON THE NEXT VERSION ***
    *** ImportVCF2 will replace it ***
      Imports control VCF file.
    - Full records: self.ctrl_dico (type = list)
    - Index: self.ctrl_index
    - Chromosome: self.ctrl_CHROM
    - Position: self.ctrl_POS
    - Reference allele: self.ctrl_REF
    - Alternative allele: self.ctrl_ALT
    '''
    def __init__(self):#, name):
        # self.name = name
        self.ctrl_index = []
        self.ctrl_dico = {}
        self.ctrl_CHROM = {}
        self.ctrl_POS = {}
        self.ctrl_ALT = {}
        self.ctrl_REF = {}

    def load_vcf(self, file_name):
        vcf_reader = vcf.Reader(open(file_name, 'r'))
        for record in vcf_reader:
            # print(record)
            self.ctrl_index.append( record.CHROM + ';' + str(record.POS) )
            self.ctrl_dico[record.CHROM + ';' + str(record.POS)]=record
            self.ctrl_CHROM[record.CHROM + ';' + str(record.POS)]=record.CHROM
            self.ctrl_POS[record.CHROM + ';' + str(record.POS)]=record.POS
            self.ctrl_ALT[record.CHROM + ';' + str(record.POS)]=record.ALT
            self.ctrl_REF[record.CHROM + ';' + str(record.POS)]=record.REF



class ImportControl2:
    '''
    *** THIS CLASS WILL DISAPEAR ON THE NEXT VERSION ***
    *** ImportVCF2 will replace it ***
      Imports control VCF file.
    - Full records: self.dico (type = list)
    - Index: self.index
    - Chromosome: self.CHROM
    - Position: self.POS
    - Reference allele: self.REF
    - Alternative allele: self.ALT
    '''
    def __init__(self):#, name):
        # self.name = name
        self.index = []
        self.dico = {}
        self.CHROM = {}
        self.POS = {}
        self.ALT = {}
        self.REF = {}

    def load_vcf(self, file_name):
        vcf_reader = vcf.Reader(open(file_name, 'r'))
        for record in vcf_reader:
            # print(record)
            self.index.append( record.CHROM + ';' + str(record.POS) )
            self.dico[record.CHROM + ';' + str(record.POS)]=record
            self.CHROM[record.CHROM + ';' + str(record.POS)]=record.CHROM
            self.POS[record.CHROM + ';' + str(record.POS)]=record.POS
            self.ALT[record.CHROM + ';' + str(record.POS)]=record.ALT
            self.REF[record.CHROM + ';' + str(record.POS)]=record.REF



class RemoveCtrlMutations:
    '''
    Compare VCF to ctrl.
    - Removes mutations found in control IF:
        - They are on the same Chromosomoe & Position
        - If the mutation is the same one
    Only report the mutations which differ
    '''

    def __init__(self):
        self.var = []
        self.cln_index = []
        self.cln_dico = {}
        self.cln_CHROM = {}
        self.cln_POS = {}
        self.cln_ALT = {}
        self.cln_REF = {}
        self.cln_INFO = {}

    def check_ctrl_in_exp(self, ctrl_object, exp_object):
        'Checks if *and removes* the intersection between control and experimental VCFs'

        print("\nThese are the common indexes between control and experimental VCFs")
        counter = 0
        ctrl_not_in_exp = []
        val = False

        for i in ctrl_object.index:
            if i in exp_object.index:
                counter += 1
                print("{0} in exp_object".format(i))
                val = self.compare_ALT(ctrl_object, exp_object, i)
                if val:
                    ### Create a new dico witn only the ones we find. or remove from the original
                    self.cln_index.append(i)
                    self.cln_dico[i] = exp_object.dico[i]
                    self.cln_CHROM[i] = exp_object.CHROM[i]
                    self.cln_POS[i] = exp_object.POS[i]
                    self.cln_ALT[i] = exp_object.ALT[i]
                    self.cln_REF[i] = exp_object.REF[i]
                    self.cln_INFO[i] = exp_object.INFO[i]
            else:
                ctrl_not_in_exp.append(i)
                
        print( "\n{0} occurrences of control in explerimental (len = {1})".format(
                counter
                    , len(exp_object.index
                        )
                     ) )
        print("### List of control not in experimental: {}".format(ctrl_not_in_exp))


    def compare_ALT(self, ctrl_object, exp_object, iteration):
        'Compares if ALT is the same in control and experimental. Returns True or False'
        
        # print ("control ALT = {0} ; exp ALT = {1}". format(ctrl_object.ALT[iteration], exp_objectALT[iteration]))
        # print ("control REF = {0} ; exp REF = {1}". format(ctrl_objectREF[iteration], exp_objectREF[iteration]))

        value = False
        if ctrl_object.ALT[iteration] == exp_object.ALT[iteration]:
            value = True

        return value


class PresenceAbsence:
    '''
    Creates a binary pandas dataframe.
    The columns are the experimental conditions
    The lines are the mutation.
    It'll create a yes/no matrix for a set of conditions
    The method plot_binary_matrix will make a png graphic of that matrix
    '''

    def __init__(self, list_experiments, list_headers, control_list):
        self.potatoe = 'potatoe'
        self.list_experiments = list_experiments
        self.list_headers = list_headers
        self.control_list = control_list
        # self.exp_list = []
        self.sorted_single_full_list = self.make_full_list()
        self.no_control_list = self.remove_control()
        # self.binary_dataframe = self.binary_lists(self.sorted_full_list)
        self.binary_df = self.binary_lists()



    def potatoes(self):
        print("Do potatoes potate?")


    def make_full_list(self):
        # parsing_list = self.list_sorter_and_replace_strings()
        temp_list = []
        for i in self.list_experiments:
            temp_list = temp_list + i 
        parsing_list = self.list_sorter_and_replace_strings(temp_list)
        # print("Before set sort, type = {1}:\n{0}\n".format(temp_list, type(temp_list.sort())))
        # temp_list = list(set(sorted(temp_list)))
        # temp_list.sort()
        # print("After set sort, type = {1}:\n{0}\n".format(temp_list, type(temp_list)))
        return parsing_list


    def list_sorter_and_replace_strings(self, li):
        '''
        Input is a list which contains STR elements "CHROM_POS"
        as the method list.sort() is either numeric or alphabetic
        it is required to remove the "CHROM_", convert the POS to INT
        And finally sort again the whole thing.
        I can also re-add a string instead of the "CHROM".
        In case of having multiple samples, the sample name may be more pertinent
        '''
        temp_list = []
        for i in range(0,len(li)):
            # print("{0}\t{1}\t{2}".format(i, li[i], li[i].split(';')[1] ))
            temp_list.append(int(li[i].split(';')[1]))
        temp_list2 = list(set(sorted(temp_list)))
        temp_list2.sort()
        # print("{1} elements in temp_list2: {0}".format(temp_list2, len(temp_list2)))
        for i in range(0,len(temp_list2)):
            temp_list2[i] = "NC_007019.1;"+str(temp_list2[i])
            # print("++++>", temp_list2[i])
        return temp_list2


    def remove_control(self):
        '''
        Removes mutations found in the control VCF
        Control VCF could mean (re)sequencing the control strain
        used in this experiment.
        '''
        no_control_list = []
        for i in self.sorted_single_full_list :
            if i not in self.control_list:
                no_control_list. append(i)
        print(">> {0} elements in the list, {1} elements in the control. {2} elements remaining".format(len(self.sorted_single_full_list), len(self.control_list), len(no_control_list) ))
        return no_control_list


    def binary_lists(self):
        '''
        Creates a 0/1 matrix of abscence/prensence of mutated loci.
        '''
        bin_dico = {}
        for li, head in zip(self.list_experiments, self.list_headers):
            temp_list = []
            for elem in self.no_control_list:
                if elem in li:
                    temp_list.append(1)
                else:
                    temp_list.append(0)
            bin_dico[head] = temp_list
        # print(bin_dico)
        mydf = pd.DataFrame(index = self.no_control_list, data = bin_dico)
        return mydf


    def plot_binary_matrix(self, outfile, data_name):
        '''
        Plots a binary matrix of presence/abscence of mutations generated by
        "binary_lists".
        '''

        #################
        ## MAKE FIGURE
        #################

        # print("{0} Y labels, listed here:\n{1}".format(len(self.no_control_list), self.no_control_list ))
        print("{0} Y labels, listed here".format(len(self.no_control_list) ))

        only_coordinates_ytick_labels = []
        for i in self.no_control_list:
            only_coordinates_ytick_labels.append(i.split(";")[1])

        only_coordinates_xtick_labels = []
        for i in range(0, len(self.list_headers)):
            only_coordinates_xtick_labels.append(i)
        # only_coordinates_xtick_labels.append(i+1)
        print(self.binary_df)

        ## Creates an empty figure
        fig = plt.figure(figsize=(4,35))
        fig.suptitle("Presence/Absence of mutations")   # Figure title, useful in case of multi-ple graphic figures

        ## Creates the proportions the graphic will have (left, bottom, width, height)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

        ## Title for that graphic
        ax.set_title(data_name)

        ## Set x axis label and ticks and tick labels
        ax.set_xlabel('Time')
        # ax.set_xticks([0, 1, 2, 3])
        ax.set_xticks(only_coordinates_xtick_labels)
        ax.set_xticklabels(self.list_headers)

        ## Set y axis label and tick labels
        ax.set_ylabel('Coordinate by experiment')
        ax.set_yticks(range(0,len(self.no_control_list)))
        # ax.set_yticklabels(self.no_control_list)
        ax.set_yticklabels(only_coordinates_ytick_labels)

        ## Do the graph
        plt.imshow(self.binary_df
            , cmap = 'binary'     ## Color MAP
            # , interpolation='nearest'
            )

        ## Show the graph -- for dev and debug
        # plt.show()

        ## Save file to file -- for script
        fig.savefig(outfile)
        plt.close(fig)


class FrequencyOneSNP:
    '''
    Follows the sampe principle as PresenceAbsence.
    Instead of working in based on the presence or absence of a locus,
    It'll filter for only the SNPs,
    It will only work with vcf lines containing only one haplotype.

    General order of the different methods
    potatoe -- Check if the class is working
    make_full_list -- Makes a sorted single list of mutations
    remove_control -- removes the mutations present in the control sample
    '''

    def __init__(self, list_experiments, list_headers, control_list, vcf_object, kind_of_mutation, outpath, outfile_root):
        self.potatoe = 'potatoe'
        self.list_experiments = list_experiments
        self.list_headers = list_headers
        self.control_list = control_list
        self.vcf_object = vcf_object
        self.outpath = outpath
        self.df_file_name = outfile_root
        # self.exp_list = []
        self.sorted_single_full_list = self.make_full_list()
        self.no_control_list = self.remove_control()
        self.kind_of_mutation = kind_of_mutation
        # self.binary_dataframe = self.binary_lists(self.sorted_full_list)
        self.binary_df = self.make_frequency_matrix(self.kind_of_mutation)


    def potatoes(self):
        print("Do potatoes potate?")


    def make_full_list(self):
        # parsing_list = self.list_sorter_and_replace_strings()
        temp_list = []
        for i in self.list_experiments:
            temp_list = temp_list + i 
        parsing_list = self.list_sorter_and_replace_strings(temp_list)
        # print("Before set sort, type = {1}:\n{0}\n".format(temp_list, type(temp_list.sort())))
        # temp_list = list(set(sorted(temp_list)))
        # temp_list.sort()
        # print("After set sort, type = {1}:\n{0}\n".format(temp_list, type(temp_list)))
        return parsing_list


    def list_sorter_and_replace_strings(self, li):
        '''
        Input is a list which contains STR elements "CHROM_POS"
        as the method list.sort() is either numeric or alphabetic
        it is required to remove the "CHROM_", convert the POS to INT
        And finally sort again the whole thing.
        I can also re-add a string instead of the "CHROM".
        In case of having multiple samples, the sample name may be more pertinent
        '''
        temp_list = []
        for i in range(0,len(li)):
            # print("{0}\t{1}\t{2}".format(i, li[i], li[i].split(';')[1] ))
            temp_list.append(int(li[i].split(';')[1]))
        temp_list2 = list(set(sorted(temp_list)))
        temp_list2.sort()
        # print("{1} elements in temp_list2: {0}".format(temp_list2, len(temp_list2)))
        for i in range(0,len(temp_list2)):
            temp_list2[i] = "NC_007019.1;"+str(temp_list2[i])
            # print("++++>", temp_list2[i])
        return temp_list2


    def remove_control(self):
        '''
        Removes mutations found in the control VCF
        Control VCF could mean (re)sequencing the control strain
        used in this experiment.
        '''
        no_control_list = []
        for i in self.sorted_single_full_list :
            if i not in self.control_list:
                no_control_list. append(i)
        print(">> {0} elements in the list, {1} elements in the control. {2} elements remaining".format(len(self.sorted_single_full_list), len(self.control_list), len(no_control_list) ))
        return no_control_list


    def make_frequency_matrix(self, kind_of_mutation):
        '''
        If len(TYPE) == 1 && TYPE == snp
        Creates a 0/1 matrix of abscence/prensence of mutated loci.
        '''
        binary_dico = {}
        mut_types_dico = {}
        mut_freq_dico = {}
        for li, head in zip(self.list_experiments, self.list_headers):
            ## I added some number coding for different elements:
            ## 0 and 1 go to binary matrix
            ## 2 is for mutation == SNP
            ## 3 is for mutation != SNP
            ## 4 Number of mutations > 1
            ## 5 This mutation is not present in the sample
            # Maybe re-think how the matrix is built
            temp_list_bin = []
            temp_list_mut_types = []
            temp_list_mut_freq = []

            # print(self.no_control_list)
            for elem in self.no_control_list:
                ### Start by choosing INFO['AC'] == 1; else pass
                ### This is the place to chose only the type of mutations I want
                if elem in li:      ## for each mutation on a given VCF file
                    temp_list_bin.append(1)

                    try:
                        ## Test if element (mutation) exists in this list
                        print(self.vcf_object.INFO[elem]['TYPE'])
                        # self.vcf_object.INFO[elem]['TYPE'] == True
                        if len(self.vcf_object.INFO[elem]['TYPE']) == 1:        ## Select only one mutation on that coordinate. Further cases to be explored later.
                            if self.vcf_object.INFO[elem]['TYPE'][0] == kind_of_mutation:   ## Select for SNPs 
                                temp_list_mut_types.append(2) # mutation is only an SNP in this sample
                                freq = ( self.vcf_object.INFO[elem]['AO'][0] / self.vcf_object.INFO[elem]['DP'] )    #### AF, Allele Frequency Calculation !!
                                print("===== FREQ = {0} \t --- Mutation = {1}".format(freq, elem) )
                                temp_list_mut_freq.append(freq)
                            else:
                                temp_list_mut_types.append(3) # mutation is not a snp
                                temp_list_mut_freq.append(0)
                        else:
                            temp_list_mut_types.append(4) # Too many mutations (>1)
                            temp_list_mut_freq.append(0)
                    except:
                        temp_list_mut_types.append(5) # mutation not present in this sample
                        temp_list_mut_freq.append(0)
                    #     print(0)
                    #     temp_list_mut.append(0)
                    # temp_list_mut.append(2)

                else:
                    temp_list_bin.append(0)
                    temp_list_mut_types.append(0)
                    temp_list_mut_freq.append(0)

            binary_dico[head] = temp_list_bin
            mut_types_dico[head] = temp_list_mut_types
            mut_freq_dico[head] = temp_list_mut_freq        ## Add a definition to the dictionary

        my_binary_df = pd.DataFrame(index = self.no_control_list, data = binary_dico)

        my_types_df = pd.DataFrame(index = self.no_control_list, data = mut_types_dico)
        my_types_df.to_csv(self.outpath + self.df_file_name + '_mut_types.csv', sep='\t')

        my_freq_df = pd.DataFrame(index = self.no_control_list, data = mut_freq_dico)
        my_freq_df.to_csv(self.outpath + self.df_file_name + '_mut_freq.csv', sep='\t')

        return my_freq_df


    def plot_binary_matrix(self, outfile, data_name):
        '''
        Plots a binary matrix of presence/abscence of mutations generated by
        "binary_lists".
        '''

        #################
        ## MAKE FIGURE
        #################

        # print("{0} Y labels, listed here:\n{1}".format(len(self.no_control_list), self.no_control_list ))
        print("{0} Y labels, listed here".format(len(self.no_control_list) ))

        only_coordinates_ytick_labels = []
        for i in self.no_control_list:
            only_coordinates_ytick_labels.append(i.split(";")[1])

        only_coordinates_xtick_labels = []
        for i in range(0, len(self.list_headers)):
            only_coordinates_xtick_labels.append(i)
        # only_coordinates_xtick_labels.append(i+1)
        print(self.binary_df)

        ## Creates an empty figure
        fig = plt.figure(figsize=(4,35))
        fig.suptitle("Presence/Absence of mutations")   # Figure title, useful in case of multi-ple graphic figures

        ## Creates the proportions the graphic will have (left, bottom, width, height)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

        ## Title for that graphic
        ax.set_title(data_name)

        ## Set x axis label and ticks and tick labels
        ax.set_xlabel('Time')
        # ax.set_xticks([0, 1, 2, 3])
        ax.set_xticks(only_coordinates_xtick_labels)
        ax.set_xticklabels(self.list_headers)

        ## Set y axis label and tick labels
        ax.set_ylabel('Coordinate by experiment')
        ax.set_yticks(range(0,len(self.no_control_list)))
        # ax.set_yticklabels(self.no_control_list)
        ax.set_yticklabels(only_coordinates_ytick_labels)

        ## Do the graph
        plt.imshow(self.binary_df
            , cmap = 'binary'     ## Color MAP
            # , interpolation='nearest'
            )

        ## Show the graph -- for dev and debug
        # plt.show()

        ## Save file to file -- for script
        fig.savefig(outfile)
        plt.close(fig)


    def plot_frequency_matrix(self, outfile, data_name):
        '''
        Plots a binary matrix of presence/abscence of mutations generated by
        "binary_lists".
        '''

        #################
        ## MAKE FIGURE
        #################

        # print("{0} Y labels, listed here:\n{1}".format(len(self.no_control_list), self.no_control_list ))
        print("{0} Y labels, listed here".format(len(self.no_control_list) ))

        only_coordinates_ytick_labels = []
        for i in self.no_control_list:
            only_coordinates_ytick_labels.append(i.split(";")[1])

        only_coordinates_xtick_labels = []
        for i in range(0, len(self.list_headers)):
            only_coordinates_xtick_labels.append(i)
        # only_coordinates_xtick_labels.append(i+1)
        print(self.binary_df)

        ## Creates an empty figure
        fig = plt.figure(figsize=(4,35))
        fig.suptitle("Frequency of SNPs")   # Figure title, useful in case of multi-ple graphic figures

        ## Creates the proportions the graphic will have (left, bottom, width, height)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

        ## Title for that graphic
        ax.set_title(data_name)

        ## Set x axis label and ticks and tick labels
        ax.set_xlabel('Time')
        # ax.set_xticks([0, 1, 2, 3])
        ax.set_xticks(only_coordinates_xtick_labels)
        ax.set_xticklabels(self.list_headers)

        ## Set y axis label and tick labels
        ax.set_ylabel('Coordinate by experiment')
        ax.set_yticks(range(0,len(self.no_control_list)))
        # ax.set_yticklabels(self.no_control_list)
        ax.set_yticklabels(only_coordinates_ytick_labels)

        ## Do the graph
        plt.imshow(self.binary_df
            , cmap = 'viridis'     ## Color MAP
            # , interpolation='nearest'
            )

        ## Show the graph -- for dev and debug
        # plt.show()

        ## Save file to file -- for script
        fig.savefig(outfile)
        plt.close(fig)

