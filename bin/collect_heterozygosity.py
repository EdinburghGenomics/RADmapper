#!/usr/bin/env python
'''
Created on June 6, 2016

@author: jrisse

Creates a table from a vcf file and family description containing information about the hetero-/homozygosity of markers containing the following fields
contig,sire_samples,valid_sire_samples,dam_samples,valid_dam_samples,male_samples,valid_male_samples,female_samples,valid_female_samples,sire_het,dam_het,male_het,female_het,sire_hom,dam_hom,male_hom,female_hom

Inputs: vcf
family file 4 column tab delimited:
#individual family generation sex 
Amale   4A      father  M
Afemale 4A      mother  F
A_female_1      4A      offspring       F
A_female_10     4A      offspring       F

TODO: Take input filenames from commandline parameters
'''
import sys
import logging
import os
from collections import defaultdict
from utils import utils_logging
from IO_interface import vcfIO




# borrowed from RAD_sex_specific_markers.py
def count_homo_het_geno(vcf_record, list_sample):
    support_het = support_hom = 0
    genotypes=vcf_record.get_valid_genotype_per_sample(genotype_quality_threshold=20, minimum_depth=6, sample_list=list_sample)
    for genotype in genotypes:
        if len(set(genotype.split('/')))!=1:
            support_het+=len(genotypes.get(genotype))
        else:
            support_hom+=len(genotypes.get(genotype))
    return support_hom, support_het

# borrowed from RAD_sex_specific_markers.py
def read_family_file(family_file):
    sample2attibutes = defaultdict(dict)
    family2samples = defaultdict(dict)
    with open(family_file) as open_file:
        for line in open_file:
            sp_line = line.strip().split()
            if len(sp_line) > 2:
                sample = sp_line[0]
                family = sp_line[1]
                status = sp_line[2]
                sample2attibutes[sample]['family']=family
                sample2attibutes[sample]['status']=status
                sample2attibutes[sample]['sex']=sp_line[3]
                if len(sp_line) > 4: sample2attibutes[sample]['additional']=sp_line[4:] 
                if status == 'offspring':
                    if not family2samples[family].has_key(status):
                        family2samples[family][status] = []
                    family2samples[family][status].append(sample)
                else:
                    family2samples[family][status] = sample
    return sample2attibutes, family2samples    

def get_male_female_samples(family_file):
    sample2attibutes, family2samples = read_family_file(family_file)
    # list required to test for conversion for potential sex markers
    list_male = []
    list_father = []
    list_female = []
    list_mother = []
    #initiate the lines for each samples
    for family in family2samples:
        mother = family2samples.get(family).get('mother')
        list_mother.append(mother)
        father = family2samples.get(family).get('father')
        list_father.append(father)
        for sample in family2samples.get(family).get('offspring'):
            if sample2attibutes.get(sample).get('sex') == "M":
                list_male.append(sample)
            elif sample2attibutes.get(sample).get('sex') == "F":
                list_female.append(sample)
    return list_male, list_female,list_father, list_mother

def get_male_female_heterozygosity(vcf_file, list_male,list_female,list_father,list_mother):
    file_handle = utils_logging.open_input_file(vcf_file, pipe=False)
    reader = vcfIO.VcfReader(file_handle)
    het_table = {}
    for vcf_records in reader:
        homozygous_male, heterozygous_male = count_homo_het_geno(vcf_records, list_male)
        homozygous_female, heterozygous_female = count_homo_het_geno(vcf_records, list_female)
        homozygous_father, heterozygous_father = count_homo_het_geno(vcf_records, list_father)
        homozygous_mother, heterozygous_mother = count_homo_het_geno(vcf_records, list_mother)
        name = vcf_records.get_reference()
        sys.stderr.write("Name:%s\n" % name)
        pos = vcf_records.get_position()
        het_table[(name,pos)] = (len(list_father),homozygous_father+heterozygous_father,len(list_mother),homozygous_mother+ heterozygous_mother,len(list_male), heterozygous_male+homozygous_male, len(list_female), heterozygous_female+homozygous_female, heterozygous_father,heterozygous_mother, heterozygous_male, heterozygous_female, homozygous_father,homozygous_mother,homozygous_male,homozygous_female)
    return het_table

def main():

    argparser = _prepare_argparser()
    args = argparser.parse_args()
    
    vcf_file=args.vcf_file
    family_file = args.family_file
    
    
    list_male,list_female,list_father,list_mother = get_male_female_samples(family_file)
    het_table = get_male_female_heterozygosity(vcf_file,list_male,list_female,list_father, list_mother)
    
    out_file = args.get(out_file,"/dev/stdout")
    with open(out_file, 'w') as out:
        out.write("contig\tsire_samples\tvalid_sire_samples\tdam_samples\tvalid_dam_samples\tmale_samples\tvalid_male_samples\tfemale_samples\tvalid_female_samples\tsire_het\tdam_het\tmale_het\tfemale_het\tsire_hom\tdam_hom\tmale_hom\tfemale_hom\n")
        for entry in het_table:
            sys.stderr.write('entry: %s__%s\n' %(entry[0],entry[1]))
            sys.stderr.write('%s\n' % het_table[entry][0])
            out.write('%s__%s\t%s\n' %(entry[0],entry[1],'\t'.join(str(i) for i in het_table[entry])))
        
def _prepare_argparser():
    """Prepare optparser object. New arguments will be added in this
    function first.
    """
    description = """"""

    argparser = ArgumentParser(description=description)

    argparser.add_argument("-v","--vcf_file",dest="vcf_file",type=str, required=True,
                         help="Vcf file from which the stats will be gathered.")
    argparser.add_argument("-f", "--family_file", dest="family_file", type=str, required=True
                           help="File describing sex and parentage of the individuals in the vcf file")
    argparser.add_argument("-o", "--out_file", dest="out_file", type=str,
                           help="Output file destination")
    return argparser
        
if __name__ == "__main__":
    main() 
