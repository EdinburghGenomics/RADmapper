#!/usr/bin/env python
'''
Created on 22 August, 2016

@author: JudithR

Utilites to convert the markers ids of linkage groups back to the #CHROM labels from the original vcf file. Works in combination with RAD_vcf_to_lepmap.py and RAD_vcf_to_lepmap_with_sexmarker_conversion.py. At the least requires a file that keeps all #CHROM labels converted to lapmap markers in that order.

First part creates the linkage file with marker id replaced with #CHROM label, second step combines the results with the table generated with collect_heterozygosity.py

TODO take filenames from command line and split into two separate parts
'''

test_nlines=99999999
contigs = []
# list of contigs were the line number corresponds to the marker number in LG*.eva.txt
with open('all_consensus_samtools1_snps_whitelist_SNP_lepmap_uniq.markers','r') as of:
    for line in of:
        line = line.strip()
        
        contigs.append(line)

# dictionary of with one entry per LG, each containing a set of contig labels corresponding to the marker number        
sets_perLG = {}

# File per LG giving marker number and distance for males and females in cM
for i in range(1,22):
    sets_perLG[i] = set()
    with open('all_consensus_samtools1_snps_whitelist_SNP_lepmap_filt_lod10_map_ordered%s.eval.txt' %i,'r') as of:
        counter = 0
        for line in of:
            
            if counter > test_nlines:
                break
            
            line = line.strip()
            if line.startswith('#'):
                continue
            counter +=1
            cols = line.split('\t')
            info_tuple = tuple([i,contigs[int(cols[0])+1].replace(':','__'),counter,cols[0]] + cols[1:3])
            sets_perLG[i].add(info_tuple)

print(str(sets_perLG))

# Table per marker with values about valid samples/heterozygosity, etc
with open('all_consensus_samtools1_snps_whitelist_heterozygosity.txt','r') as of:
    for counter,line in enumerate(of):
        line = line.strip()
        line = line.split('\t')
        if counter == 0:
            line.extend(["LG","Marker_index", "male_distance","female_distance"])
            print('\t'.join(line))
        else:
            # get list of first value of each tuple in the set of linkage group
            def scan(marker_label):
                for v in sets_perLG.values():
                    for info_tuple in v:
                        if info_tuple[1] == marker_label:
                            return info_tuple

            #lg_found = [k for k, v in sets_perLG.items() if line[0] in [v2[0] for v2 in v]]
            #assert len(lg_found) <=1
                         
            lg_info_tuple = scan(line[0])
            if lg_info_tuple:
                line.append(str(lg_info_tuple[0]))
                line.append(str(lg_info_tuple[2]))
                line.append(lg_info_tuple[4])
                line.append(lg_info_tuple[5])
                
            else:
                line.extend(['0','NA','NA','NA'])
            print('\t'.join(line))
