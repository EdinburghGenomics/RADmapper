#!/usr/bin/env python
'''
Created on Mar 9, 2011

@author: tcezard,JudithR

Converts a vcf file in combination with a family description into an input file suitable for lepmap2 in the same way as RAD_vcf_to_lepmap.py.
It does however rewrite putative sex marker in males in an XO/XX reproductive system to heterozygous markers containing an additional allele as described in the two lookup tables in the code. This to be able to treat the sex chromosome linkage group as autosome.
'''
import sys
import logging
from optparse import OptionParser
import os
from collections import Counter, defaultdict
import getpass
import time
import re

from utils import utils_logging
from IO_interface import vcfIO
from utils import atgc2iupac

""" Family file layout
#sample family status sex additional
Amale   4A      father  M
Afemale 4A      mother  F
A_female_1      4A      offspring       F
"""

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

def rewrite_male(d_gt,s_gt,gt):
    lookup = {
    ("0/1","0/0","0/0"):"0/2",
    ("0/1","1/1","0/0"):"0/2",
    ("0/2","0/0","0/0"):"0/1",
    ("0/2","2/2","0/0"):"0/1",
    ("1/2","1/1","1/1"):"1/3",
    ("1/2","2/2","1/1"):"1/3",
    ("0/1","0/0","1/1"):"1/2",
    ("0/1","1/1","1/1"):"1/2",
    ("0/2","0/0","2/2"):"1/2",
    ("0/2","2/2","2/2"):"2/1",
    ("1/2","1/1","2/2"):"2/3",
    ("1/2","2/2","2/2"):"2/3"
    }
    return lookup.get((d_gt,s_gt,gt)) or gt



def rewrite_sire(d_gt,s_gt):
    lookup = {
    ("0/1","0/0"):"0/2",
    ("0/1","1/1"):"1/2",
    ("0/1","2/2"):None,
    ("0/2","0/0"):"0/1",
    ("0/2","1/1"):None,
    ("0/2","2/2"):"2/1",
    ("1/2","0/0"):None,
    ("1/2","1/1"):"1/3",
    ("1/2","2/2"):"2/3"
    }
    return lookup.get((d_gt,s_gt)) or s_gt

    
"""1st column Id for each family
2nd column individual Id
3rd Sire id
4th Dam id
5th sex (doesn't play a role just ensure your parents are coded as male
and female) 1 for male 2 for female
6th just put 0
7th onwards markers (both alleles are coded separated by space)"""

def vcf_to_lepmap(vcf_file, output_file, family_file, genotype_quality_threshold, max_prop_missing, min_het_female, max_het_male,max_het_sire):
    file_handle = utils_logging.open_input_file(vcf_file, pipe=False)
    reader = vcfIO.VcfReader(file_handle)
    all_samples_in_file = reader.get_sample_names()
    if family_file:
        sample2attibutes, family2samples = read_family_file(family_file)
        all_samples_in_file = sample2attibutes.keys()
    else:
        return

    max_missing = int(len(all_samples_in_file) * max_prop_missing)
    all_lines = {}
    nb_lines = 0
    nb_sequence = 0
    count_more_than_one_allele = 0
    count_indel = 0
    count_too_many_missing = 0
    count_non_polymorphic = 0
    nb_missing_per_sample = Counter()
    
    # list required to test for conversion for potential sex markers
    list_male = []
    list_female = []
    list_sire = []
    list_dam = []
    #initiate the lines for each samples
    for family in family2samples:
        mother = family2samples.get(family).get('mother')
        list_dam.append(mother)
        father = family2samples.get(family).get('father')
        list_sire.append(father)
        all_lines[father] = [family, father, "0", "0", "1", "0", "1 2"]
        all_lines[mother] = [family, mother, "0", "0", "2", "0", "1 1"]
        all_lines[father].extend(sample2attibutes.get(father).get('additional',[]))
        all_lines[mother].extend(sample2attibutes.get(mother).get('additional',[]))
        for sample in family2samples.get(family).get('offspring'):
            if sample2attibutes.get(sample).get('sex') == "M":
                all_lines[sample] = [family, sample, father, mother, "1" , "0",  "1 2"]
                # append male sample to list
                list_male.append(sample)
            elif sample2attibutes.get(sample).get('sex') == "F":
                all_lines[sample] = [family, sample, father, mother, "2" , "0",  "1 1"]
                # append female sample to list
                list_female.append(sample)
            else:
                all_lines[sample] = [family, sample, father, mother, "2" , "0",  "0 0"]
            all_lines[sample].extend(sample2attibutes.get(sample).get('additional',[]))
    all_markers = []
    for vcf_records in reader:
        nb_lines += 1
        if nb_lines % 10000 == 0:
            sys.stdout.write('.')
        ref_base = vcf_records.get_reference_base()
        alt_bases = vcf_records.get_alt_bases()
        
            
        # if we pass this test we will rewrite the male offspring and the sire to be X/X instead of X/.
        # the test is borrowed from RAD_sex_specific_markers.py
        convert = False
        
        
        support_for_homozygous_male, support_for_heterozygous_male = count_homo_het_geno(vcf_records, list_male)
        support_for_homozygous_female, support_for_heterozygous_female = count_homo_het_geno(vcf_records, list_female)
        support_for_homozygous_sire, support_for_heterozygous_sire = count_homo_het_geno(vcf_records, list_sire)
        support_for_homozygous_dam, support_for_heterozygous_dam = count_homo_het_geno(vcf_records, list_dam)
        #logging.info("total male: %s\nsupport for hom male: %s\nsupport for het male: %s\nsupport for het female: %s" %(len(list_male),support_for_homozygous_male,support_for_heterozygous_male,support_for_heterozygous_female))
        # if at least 80% of the male offspring are homozygous, none are heterozygous and at least one female offspring is heterozygous we pass
        #if(support_for_homozygous_male>len(list_male)*0.8 and support_for_heterozygous_male == 0 and support_for_heterozygous_female > 1):
        #        convert = True
        
        # new criteria from R plot filtering
        if(support_for_heterozygous_female>= min_het_female and support_for_heterozygous_male<= max_het_male and support_for_heterozygous_sire<= max_het_sire):
            convert=True

        # TODO sort this out as it probably gets rid of sex specific markers by accident, may be letting those that pass the convert test through would be enough
        # No idea what lepmap does with multi-allelic sites, but this is what the manual says for the filtering step:
        # By default, genotypes are outputted in bi-allelic format with alleles 1 and 2. 
        # By providing option keepAlleles=1 the original alleles are kept 
        # (If there are a lot of markers with 3-4 alleles (e.g. microsatellites), keeping the alleles is probably wise as markers with more alleles have more information about the haplotypes).
        if len(alt_bases) > 2:
            count_more_than_one_allele += 1
            continue
        if vcf_records.is_indel():
            count_indel += 1
            continue

        nb_missing = 0
        all_chars = []
        all_codes = set()

        for sample in all_samples_in_file:
            gt = vcf_records.get_genotype(sample)
            gq = vcf_records.get_genotype_quality(sample)
            
            
             # need to get the genotypes of mother and father to see what we have to work with
            family = sample2attibutes.get(sample).get('family')
            mother = family2samples.get(family).get('mother')
            father = family2samples.get(family).get('father')
            d_gt = vcf_records.get_genotype(mother)
            s_gt = vcf_records.get_genotype(father)
            # if there are missing genotypes in the parents we won't convert ,should have been caught by the ShomDhet filter from VariantFiltration, but just to be sure
            if (d_gt == './.' or s_gt == './.'): 
                convert = False
            
            # right, here goes all the conversion code
            if (convert):
                # if the the current sample is female, no conversion is required
                s_sex = sample2attibutes.get(sample).get('sex')
                # if sample is father convert according to both parents gt
                if(sample2attibutes.get(sample).get('status') == 'father'):
                    logging.info("Rewriting %s from %s" %(sample,gt))
                    gt = rewrite_sire(d_gt,gt)
                    logging.info("to %s" %gt)
                # if the sample is male offspring and has a gt, convert according to its gt and parents gt
                elif (gt and s_sex == 'M'):
                    logging.info("Rewriting %s from %s" %(sample,gt))
                    gt = rewrite_male(d_gt,s_gt,gt)
                    logging.info("to %s" %gt)
            if gt and gq > genotype_quality_threshold:
                value1, value2 = re.split('[/|]', gt)
                # TODO need to output the original vcf line to make it clear which ones ended up in the lepmap file
                code = '%s %s'%(int(value1)+1, int(value2)+1)
            else:
                nb_missing += 1
                nb_missing_per_sample[sample] += 1
                code = '0 0'
            all_chars.append(code)
            all_codes.add(code)

        if len(all_codes) == 1:
            count_non_polymorphic += 1
            continue
        if nb_missing <= max_missing:
            nb_sequence += 1
            for i, sample in enumerate(all_samples_in_file):
                all_lines[sample].append(all_chars[i])
                all_markers.append('%s:%s'%(vcf_records.get_reference(), vcf_records.get_position()))
        else:
            count_too_many_missing += 1

    if count_more_than_one_allele:
        logging.warning("%s snps remove because they had more than 2 alleles" % (count_more_than_one_allele))
    if count_indel:
        logging.warning("%s indels removed" % (count_indel))
    if count_non_polymorphic:
        logging.warning(
            "%s snps removed because no polymorphism was found between populations" % (count_non_polymorphic))
    if count_too_many_missing:
        logging.warning("%s snps removed because >%s missing samples" % (count_too_many_missing, max_missing))
    logging.info("%s snps output in Lepmap format" % (nb_sequence))
    for sample in nb_missing_per_sample:
        logging.info("%s markers missing in %s" % (nb_missing_per_sample.get(sample), sample))

    with open(output_file, 'w') as open_output:
        for sample in all_lines:
            open_output.write('%s\n'%'\t'.join(all_lines.get(sample)))
    with open(output_file + '.markers', 'w') as open_output:
        open_output.write('\n'.join(all_markers))



def main():
    # initialize the logging
    utils_logging.init_logging()
    #Setup options
    optparser = _prepare_optparser()
    (options, args) = optparser.parse_args()
    #verify options
    arg_pass = _verifyOption(options)
    if not arg_pass:
        logging.warning(optparser.get_usage())
        logging.critical("Non valid arguments: exit")
        sys.exit(1)
    if options.debug:
        utils_logging.init_logging(logging.DEBUG)

    vcf_to_lepmap(options.input_vcf_file, options.output_file, options.family_file,
                  options.geno_qual_threshold, options.max_prop_missing, options.min_female_het, options.max_male_het, options.max_sire_het)


def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog"""
    description = """"""

    optparser = OptionParser(description=description, usage=usage, add_help_option=False)
    optparser.add_option("-h", "--help", action="help", help="show this help message and exit.")
    optparser.add_option("-i", "--input_vcf_file", dest="input_vcf_file", type="string",
                         help="Path to input vcf file where the SNPs are located. Default: %default")
    optparser.add_option("-o", "--output_file", dest="output_file", type="string",
                         help="Path to the file that reformated output. Default: %default")
    optparser.add_option("-f", "--family_file", dest="family_file", type="string", default=None,
                         help="Path to the file that contains family, parents and sex information. Default: %default")
    optparser.add_option("-g", "--geno_qual_threshold", dest="geno_qual_threshold", type="int", default=20,
                         help="The genotype quality threshold above which genotypes will be used. Default: %default")
    optparser.add_option("-x", "--max_prop_missing", dest="max_prop_missing", type="float", default=.5,
                         help="The maximum of missing samples across all populations. Default: %default")
    optparser.add_option("-m", "--min_female_het", dest="min_female_het", type="int",default=10,
                         help="The minimum number of heterozygous females to be considered a sexmarker. Default: %default")
    optparser.add_option("-n", "--max_male_het", dest="max_male_het", type="int",default=1,
                         help="The maximum number of heterozygous males to be considered a sexmarker. Default: %default")
    optparser.add_option("-s", "--max_sire_het", dest="max_sire_het", type="int",default=0,
                         help="The maximum number of heterozygous sires to be considered a sexmarker. Default: %default")                         
    optparser.add_option("--phased", dest="phased", action="store_true", default=False,
                         help="Use Phasing information (if available) to create larger markers. Default: %default")
    optparser.add_option("--debug", dest="debug", action="store_true", default=False,
                         help="Set the verbosity to debug mode. Default: %default")
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass = True

    if not options.input_vcf_file or not os.path.exists(options.input_vcf_file):
        logging.error("You must specify a valid input file.")
        arg_pass = False
    if not options.output_file or not os.path.exists(os.path.dirname(os.path.abspath(options.output_file))):
        logging.error("You must specify a valid output file.")
        arg_pass = False
    return arg_pass


if __name__ == "__main__":
    main()   

