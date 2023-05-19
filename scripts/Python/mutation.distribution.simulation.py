# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 16:35:59 2019

@author: rmoge
"""
#%%



###############################################################################################################
###############################################################################################################
#####NOTE: IF YOU ARE GOING TO INCLUDE tRNAs AND SO FORTH, YOU HAVE TO USE A DICT OF "genes" AND NOT "CDSs"
############################################################################################################
####################################################################################################


#%%
#############
#%%
from __future__ import division
import re, os,sys, math, operator,random, copy,collections,time; import numpy as np; import pandas as pd; import csv
from itertools import groupby; import pprint as pp
import matplotlib.pyplot as plt
#%%
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
#%%
# get all sequence records for the specified genbank file
recs = [rec for rec in SeqIO.parse(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\syn3A_genome\Synthetic.bacterium_JCVI-Syn3A.gb", "genbank")]

# print the number of sequence records that were extracted
print(len(recs))

# print annotations for each sequence record
for rec in recs:
	print(rec.annotations)
    
myrec=recs[0]

# print the CDS sequence feature summary information for each feature in each
# sequence record
for rec in recs:
    feats = [feat for feat in rec.features if feat.type == "CDS"]
    for feat in feats:
        print(feat)
        
print(repr(myrec.seq))
#%%
gb_file=r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\syn3A_genome\Synthetic.bacterium_JCVI-Syn3A.gb"
gb_record = SeqIO.read(open(gb_file,"r"), "genbank")
print("Name %s, %i features" % (gb_record.name, len(gb_record.features)))
print(repr(gb_record.seq))


#%%
print(myrec.features[2].location.start)
myrec[(myrec.features[2].location.start)]#this is how you access the 0th nt of the gb file...
#for the parallel stuff, we just need lengths
#%%
def gene_to_tag(gb_record, feature_type, qualifier) :
    answer = dict()
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type=='gene':
            if qualifier in feature.qualifiers :
                #There should only be one locus_tag per feature, but there
                #are usually several db_xref entries
                for value in feature.qualifiers[qualifier] :
                    if value in answer :
                        print("WARNING - Duplicate key")
                        answer[value]['index'].append(index)
                    else :
                        answer[value] = dict()
                        answer[value]['index'] = list()
                        answer[value]['index'].append(index)
            if qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier] :
                    answer[value]['locus_tag']=feature.qualifiers['locus_tag']    
    return answer

#test = index_genbank_features(myrec, "CDS", "gene")
GD = gene_to_tag(myrec, "locus_tag", "gene")
#%%
# I think it may make more sense to build the dict all in one loop, rather than piping an updated dict from one function to the next, though. I attempt to do that in this chunk.
#We are gonna build a badass dict of dicts. Each locus tag will be a key. Its value will be a dict. In that dict, each key will be some kind of feature (e.g., length) and its value will be the corresponding value for THAT LOCUS.
def index_genbank_features(gb_record, feature_type, qualifier, qualifier_translation, AT_rel_rate, GC_rel_rate):
    max_gene_length = max([len(i) for i in gb_record.features[1:]])
    rel_rate_MAX=max(AT_rel_rate,GC_rel_rate)
    print(str(max_gene_length) + " is max gene length")
    answer = dict()
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==feature_type :
            if qualifier in feature.qualifiers :
                #There should only be one locus_tag per feature, but there
                #are usually several db_xref entries
                for value in feature.qualifiers[qualifier] :
                    if value in answer :
                        print("WARNING - Duplicate key %s for %s features %i and %i" \
                           % (value, feature_type, answer[value], index))
                    else :
                        answer[value] = {}
                        answer[value]["index"] = index
            try:
                mygg=feature.qualifiers[qualifier_translation]#attempt to access the gene name for this locus tag
            except KeyError:
                mygg=feature.qualifiers[qualifier]#if there is no gene name, just call it by the locus tag
            if qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier] :
                    answer[value]['gene'] = mygg[0]
            if qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier] :
                    answer[value]['gene_length']=len(feature)
            if qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier] :        
                    answer[value]["gene_len_rel"]= (len(feature) / max_gene_length)
            if qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier]:
                    GCcount=0
                    for (index2,pos) in enumerate(gb_record.features[index]):
                        if gb_record[pos]=="C" or gb_record[pos]=="G":
                            GCcount+=1
                        else:
                            pass
                        #print(pos)
                    answer[value]["gene_GC_prop"]= (GCcount/len(feature))
            if qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier]:
                    answer[value]['mut_rel_rate']=((answer[value]['gene_GC_prop']*GC_rel_rate) + ((1 - answer[value]['gene_GC_prop'])*AT_rel_rate)) / rel_rate_MAX
                    
            
    return answer

#test = index_genbank_features(myrec, "CDS", "gene")
DD = index_genbank_features(myrec, "gene", "locus_tag", 'gene', AT_rel_rate = 35/0.76, GC_rel_rate=486/0.24)#at_rel_rate is the relative rate of mutation of an A:T nucleotide; gc_rel_rate mutatis mutandis
#%%
#Import mutation frequency data
def freq_list(infile=""):
    mydf = pd.read_csv(infile)
    return mydf          
            
#syn3B_freqs = freq_list(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\mutation.frequencies_syn3B.csv")######THIS FILE INCLUDES SYNONYMOUS MUTATIONS, which means that your P-values may be overly conservative
#syn3B_freqs_nosyn = freq_list(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\mutation.frequencies_no.syn_syn3B.csv")###in this file all the synonymous mutations have been deleted
syn3B_freqs_onlysyn = freq_list(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\mutation.frequencies_only.synonymous_syn3B.csv")###in this file all the NONsynonymous mutations have been deleted
#syn3B_freqs.set_index('freq',inplace=True)
syn3B_freqs_raw = tuple(syn3B_freqs_onlysyn.loc[:, 'freq'].values)
#%%
roogeld = dict()
roogeld['JCVISYN3A_0001']= DD['JCVISYN3A_0001']; roogeld['JCVISYN3A_0002'] = DD['JCVISYN3A_0002']; roogeld['JCVISYN3A_0003'] = DD['JCVISYN3A_0003']

def simulate_mutations(mutation_freq_tup, indict, reps=10000):
    counter = 0; sim_master=dict(); sim_dict=dict()
    for k in indict.keys():
        sim_master[k] = list()
    while counter < reps:
        sim_temp=dict()
        for index, m in enumerate(mutation_freq_tup):
            while True:
                w = random.choice(list(indict));r = random.random();rGC=random.random()
                #print("w = " + str(w) + ", w len rel = " + str(indict[w]['gene_len_rel']) + " and r = " + str(r))
                if indict[w]['gene_len_rel'] >= r and indict[w]['mut_rel_rate'] >= rGC:
                    #do something###########################################
                    ##################################
                    #################################################
                    try:
                        sim_temp[w] += m
                    except KeyError:
                        sim_temp[w] = 0
                        sim_temp[w] += m
                    #print("value added = " + str(m))
                    break
                else:
                    continue
        sim_dict[counter] = sim_temp
        for k,v in sim_temp.items():
            sim_master[k].append(v)
        for km in sim_master.keys():
            if len(sim_master[km]) < counter+1:
                sim_master[km].append(0)
        counter += 1
        #pp.pprint(sim_temp)
    
    return sim_master
#%%
start = time.time()
smo = simulate_mutations((syn3B_freqs_raw),DD,reps = 100000)
end = time.time()
print(end - start)
#10000 simulations took 29.5 seconds. I could thus do 100 000 simulations in about 6 minutes, or 1 000 000 in 60 minutes.
#print(smo) 
###############################################################################################################
#Next steps:
#Compare simulated values to actual values. Get p values. Then, get Q values (FDR-adjusted P-values)
##could use the GENExPOP matrix to do that.
#%%
##Accessing the simulation list of values BY GENE NAME:
##homer = smo[GD['ftsZ']['locus_tag'][0]]
#o = np.count([i>=x for i in reps]), and I let reps = 10 000 in this simulation.
#for FtsZ, I get P < 0.0001      x = 4
#for FakA, I get P = 0.0016 x = 3    #note that these were calculated including synonymous mutations, which means the P-values may be overly conservative
    #tetM: P = 0.0019       x = 2.877
    #dnaN: P < 0.0001       x = 4
    #dnaA: P = 0.1388       x = 1
    #ptsG: P = 0.0264       x = 1.969
    #atpD: P < 0.0001        x = 3.052
    #rpoC: P = 0.0011       x = 4.056
    #rpoB: P = 0.3692      x = 1
    #amiF: P = 0.0023 x = 2.596 (synonymous pmsm NOT Included)
    #pyrG: P < 0.0001     x = 4
    
#recalculated without synonymous mutns in the simulation, this time with 100 000 reps
    
#%%  
  '''  #these P-values were calculated without any synonymous mutations counted, 100 000 simulations
11*np.count_nonzero([i>=4 for i in smo[GD['ftsZ']['locus_tag'][0]]])######P = 0.00001 padj < 0.00011
11*np.count_nonzero([i>=4 for i in smo[GD['dnaN']['locus_tag'][0]]])######P < 0.00001 padj = 0.00022
11*np.count_nonzero([i>=4 for i in smo[GD['pyrG']['locus_tag'][0]]])#P < 0.00001      padj < 0.00011
11*np.count_nonzero([i>=1 for i in smo[GD['rpoB']['locus_tag'][0]]])#####P = 0.33659  padj = 1
11*np.count_nonzero([i>=4.056 for i in smo[GD['rpoC']['locus_tag'][0]]])##P = 0.00031 padj = 0.00484
11*np.count_nonzero([i>=0.456 for i in smo[GD['rpoD']['locus_tag'][0]]])##P = 0.16773 padj = 1
11*np.count_nonzero([i>=3 for i in smo[GD['fakA']['locus_tag'][0]]])######P = 0.00084 padj = 0.00748
11*np.count_nonzero([i>=1.969 for i in smo[GD['ptsG']['locus_tag'][0]]]) #P = 0.02367 padj = 0.25630 <-- So the Bonferroni corrected is no longer significant.
11*np.count_nonzero([i>=3.052 for i in smo[GD['atpD']['locus_tag'][0]]]) #P < 0.00001 padj < 0.00011
11*np.count_nonzero([i>=2.877 for i in smo[GD['tetM']['locus_tag'][0]]])##P = 0.00141 padj = 0.01463
11*np.count_nonzero([i>=2.596 for i in smo[GD['amiF']['locus_tag'][0]]])#P = 0.00121  padj = 0.01386

#TO GET BONFERRONI ADJUSTED P-VALUES, YOU MULTIPLIED THESE P-VALUES BY 11. You can also use the FDR and get Q-values, or use Dunn-Sidak, to get less conservative estimates.

#np.count_nonzero([i>=1 for i in smo[GD['dnaA']['locus_tag'][0]]])######P = 0.13024
#np.count_nonzero([i>=1.074 for i in smo[GD['amiF']['locus_tag'][0]]])#P = 0.04771

'''
#############
#%%
#Now do 100 000 simulations for syn1.0
syn1recs = [rec for rec in SeqIO.parse(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\syn1.0_genome\Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.gb", "genbank")]

# print the number of sequence records that were extracted
print(len(syn1recs))

# print annotations for each sequence record
for rec in syn1recs:
	print(rec.annotations)
    
syn1rec=syn1recs[0]

# print the CDS sequence feature summary information for each feature in each
# sequence record
for rec in syn1recs:
    feats = [feat for feat in rec.features if feat.type == "CDS"]
    for feat in feats:
        print(feat)

GD1 = gene_to_tag(syn1rec, "locus_tag", "gene")
DD1 = index_genbank_features(syn1rec, "gene", "locus_tag", 'gene', AT_rel_rate=219/0.76,GC_rel_rate=1591/0.24)

#syn1_freqs_nosyn = freq_list(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\mutation.frequencies_no.syn_syn1.0.csv")###in this file all the synonymous mutations have been deleted
syn1_freqs_onlysyn = freq_list(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\mutation.frequencies_only.synonymous_syn1.0.csv")###in this file all the NONsynonymous mutations have been deleted
#syn3B_freqs.set_index('freq',inplace=True)
#syn1_freqs_raw = tuple(syn1_freqs_nosyn.loc[:, 'freq'].values)
syn1_freqs_raw = tuple(syn1_freqs_onlysyn.loc[:, 'freq'].values)

#%%


start = time.time();print(start)
smo1 = simulate_mutations((syn1_freqs_raw),DD1,reps = 100000)
end = time.time()
print(end - start)
#930 seconds for 100 000 reps with the GC code
#%%  
    #these P-values were calculated without any synonymous mutations counted, 100 000 simulations
    #IGNORE ALL THE P VALUES---they are wrong. Only read the padj values. Divide them by 12 for the P value. Note that the padj values are from Bonferroni correction. To be less conservative, you could do Q value from FDR or Dunn-Sidak
'''
12*np.count_nonzero([i>=3.079 for i in smo1[GD1['ftsZ']['locus_tag'][0]]])######P = 0.00001 padj < 0.00012
12*np.count_nonzero([i>=2.575 for i in smo1[GD1['dnaA_1']['locus_tag'][0]]])######P < 0.00001 padj = 0.00108
12*np.count_nonzero([i>=1.528 for i in smo1[GD1['rpoA']['locus_tag'][0]]])#P < 0.00001      padj < 0.02016
12*np.count_nonzero([i>=1.046 for i in smo1[GD1['rpoB']['locus_tag'][0]]])#####P = 0.33659  padj = 0.94428
12*np.count_nonzero([i>=0.782 for i in smo1[GD1['rpoC']['locus_tag'][0]]])##P = 0.00031 padj = 1
12*np.count_nonzero([i>=1.000 for i in smo1[GD1['rpoD']['locus_tag'][0]]])##P = 0.16773 padj = .72120
12*np.count_nonzero([i>=5.183 for i in smo1[GD1['lpdA']['locus_tag'][0]]])######P = 0.00084 padj < 0.00012
12*np.count_nonzero([i>=2 for i in smo1[GD1['tnpA_1']['locus_tag'][0]]]) #P = 0.02367 padj = 0.00336
12*np.count_nonzero([i>=2 for i in smo1[GD1['tnpB_1']['locus_tag'][0]]]) #P < 0.00001 padj = 0.00552
12*np.count_nonzero([i>=1.546 for i in smo1[GD1['tetM']['locus_tag'][0]]])##P = 0.00141 padj = 0.08208
12*np.count_nonzero([i>=1.335 for i in smo1['MMSYN1_0460']])##############padj = 0.33072
12*np.count_nonzero([i>=3.606 for i in smo1['MMSYN1_0641']])#P = 0.00121  padj < 0.00012
#Dunn-Sidak
1-(1-(np.count_nonzero([i>=3.079 for i in smo1[GD1['ftsZ']['locus_tag'][0]]]))/100000)**12######P = 0.00001 padj < 0.00012
1-(1-(np.count_nonzero([i>=2.575 for i in smo1[GD1['dnaA_1']['locus_tag'][0]]]))/100000)**12######P < 0.00001 padj = 001079465560347992
1-(1-(np.count_nonzero([i>=1.528 for i in smo1[GD1['rpoA']['locus_tag'][0]]]))/100000)**12#P < 0.00001      padj < 0.019974760826477422
1-(1-(np.count_nonzero([i>=1.046 for i in smo1[GD1['rpoB']['locus_tag'][0]]]))/100000)**12#####P = 0.33659  padj = 0.6260018785275774
1-(1-(np.count_nonzero([i>=0.782 for i in smo1[GD1['rpoC']['locus_tag'][0]]]))/100000)**12##P = 0.00031 padj = 0.6260018785275774
1-(1-(np.count_nonzero([i>=1.000 for i in smo1[GD1['rpoD']['locus_tag'][0]]]))/100000)**12##P = 0.16773 padj = 0.5246868876744993
1-(1-(np.count_nonzero([i>=5.183 for i in smo1[GD1['lpdA']['locus_tag'][0]]]))/100000)**12######P = 0.00084 padj < 0.00012
1-(1-(np.count_nonzero([i>=2 for i in smo1[GD1['tnpA_1']['locus_tag'][0]]]))/100000)**12 #P = 0.02367 padj = 0.00335483042639817
1-(1-(np.count_nonzero([i>=2 for i in smo1[GD1['tnpB_1']['locus_tag'][0]]]))/100000)**12 #P < 0.00001 padj = 0.005506055791773101
1-(1-(np.count_nonzero([i>=1.546 for i in smo1[GD1['tetM']['locus_tag'][0]]]))/100000)**12##P = 0.00141 padj = 0.07906148163292737
1-(1-(np.count_nonzero([i>=1.335 for i in smo1['MMSYN1_0460']]))/100000)**12##############padj = 0.2849214088321784
1-(1-(np.count_nonzero([i>=3.606 for i in smo1['MMSYN1_0641']]))/100000)**12#P = 0.00121  padj < 0.00012
'''

#%%
#o = np.count([i>=x for i in reps])
#Let x = number of OBSERVED mutations for gene of interest.
#Let [reps] be the array of simulated mutation values for the gene of interest
#Then do this line of code for each gene (in order to do this, you need the GENExPOP matrix). Store a P-value in DD for each gene, where P = o / len(reps). We are asking, "What proportion of the time do we see AT LEAST as many mutations (as were actually observed) to occur BY CHANCE ALONE?"
#Next, you need to get Q values (FDR-adjusted P values)
##Rank the genes by P value
##Calculate Q as p_i * m / i, where i is the rank, p_i is the P value for the gene, and m is the total number of P values in your rankings

#That said--- the Q value and FDR are for GWAS where you are fishing for hits. If I just submit my candidate genes, such as ftsZ, I am NOT MAKING 500 TESTS. I am only testing the ~10 genes that I really think might be involved in adaptation.

#%%
#A next step could be to add the gene LENGTH to the dict. Later when youre worrying about dNdS, you would add the gene (DNA) SEQUENCE to the dict.OR MAYBE THE gff3 file would have an easier way to access teh DNA sdequence; worth looking at.
            
def looptest(dicter, reps=2):
    counter = 0
    while counter < reps:
        while True:
            w = random.choice(list(dicter)); r = random.random()
            print("w = " + str(w) + ", w len rel = " + str(dicter[w]['gene_len_rel']) + " and r = " + str(r))
            if dicter[w]['gene_len_rel'] >= r:
                print("sucksess")
                break
            else:
                continue
        counter += 1    
    
    
#looptest(DD)
    
    
    
#%%
def output_genes(gb_record, outfile=""):
    tempdf=pd.DataFrame()
    counter=0
    for (index,feature) in enumerate(gb_record.features):
        if feature.type == 'gene':
            if 'gene' in feature.qualifiers:
                for v in feature.qualifiers['gene']:
                    tempdf[str(feature.qualifiers['gene'][0])] = ""
            else:#~147 genes do not have assigned gene names and are known ONLY by the locus tag.
                for v in feature.qualifiers['locus_tag']:
                    tempdf[str(feature.qualifiers['locus_tag'][0])] = ""
            counter+=1
    print(counter)            

    tempdf.to_csv(index=True,path_or_buf=outfile)    
    return tempdf

gxp = output_genes(myrec,r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\gxp_headers_syn3B.csv")
#pp.pprint(gxp)
##
#%%
#Now output the simulation results
simulation_3B_100000 = pd.DataFrame.from_dict(smo)
simulation_1_100000 = pd.DataFrame.from_dict(smo1)

simulation_3B_100000.to_csv(index=True,path_or_buf=r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\simulations_syn3B_GC_only.synonymous_100000.csv")
simulation_1_100000.to_csv(index=True,path_or_buf=r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\simulations_syn1.0_GC_only.synonymous_100000.csv")
#%%
#You can import simulation results from a file.
start = time.time();print(start)
recs = [rec for rec in SeqIO.parse(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\syn3A_genome\Synthetic.bacterium_JCVI-Syn3A.gb", "genbank")]
myrec=recs[0]
syn1recs = [rec for rec in SeqIO.parse(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\syn1.0_genome\Synthetic.Mycoplasma.mycoides.JCVI-syn1.0_CP002027.1.gb", "genbank")]  
syn1rec=syn1recs[0]

GD = gene_to_tag(myrec, "locus_tag", "gene")
DD = index_genbank_features(myrec, "gene", "locus_tag", 'gene', AT_rel_rate = 35/0.76, GC_rel_rate=486/0.24)
GD1 = gene_to_tag(syn1rec, "locus_tag", "gene")
DD1 = index_genbank_features(syn1rec, "gene", "locus_tag", 'gene', AT_rel_rate=219/0.76,GC_rel_rate=1591/0.24)


#smopd = pd.read_csv(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\simulations_syn3B_100000.csv")
#smo1pd = pd.read_csv(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\simulations_syn1.0_100000.csv")
smopd = pd.read_csv(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\simulations_syn3B_GC_only.synonymous_100000.csv")
smo1pd = pd.read_csv(r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\simulations_syn1.0_GC_only.synonymous_100000.csv")


smo = pd.DataFrame.to_dict(smopd)
smo1 = pd.DataFrame.to_dict(smo1pd)
end = time.time()
print(end - start)
#2405 seconds for 100 000 reps with the GC code
#%%
#Now do the analyses on the imported simulation results
#Below: SYNONYMOUS ONLY
my3bpraw=dict()
#my3bpraw['amiF']=np.count_nonzero([i>=0.633 for i in smo[GD['amiF']['locus_tag'][0]]])#P = 0.01794  padj = 0.01463
#my3bpraw['lpdA']=np.count_nonzero([i>=0.363 for i in smo[GD['lpdA']['locus_tag'][0]]])##############padj = 0.02915
#my3bpraw['lgt']=np.count_nonzero([i>=2.079 for i in smo[GD['lgt']['locus_tag'][0]]])#P = 0.00002

#%%
mys1praw=dict()
#Below: SYNONYMOUS ONLY
mys1praw['secA']=np.count_nonzero([i>=0.509 for i in smo1[GD1['secA']['locus_tag'][0]]]) #P=0.03738
#mys1praw['MMSYN1_0086']=np.count_nonzero([i>=0.303 for i in smo1['MMSYN1_0086']])#P=.02202
mys1praw['MMSYN1_0676']=np.count_nonzero([i>=0.379 for i in smo1['MMSYN1_0676']])#P=.01345

#%%
#Now for the nonsynonymous
#Now do the analyses on the imported simulation resulsts
my3bpraw_nonsyn=dict()
my3bpraw['ftsZ']=np.count_nonzero([i>=4 for i in smo[GD['ftsZ']['locus_tag'][0]].values()])######P = 0.00001 padj < 0.00011
my3bpraw['dnaN']=np.count_nonzero([i>=4 for i in smo[GD['dnaN']['locus_tag'][0]].values()])######P < 0.00001 padj = 0.00011
my3bpraw['pyrG']=np.count_nonzero([i>=4 for i in smo[GD['pyrG']['locus_tag'][0]].values()])#P < 0.00001      padj = 0.00044
my3bpraw['rpoC']=np.count_nonzero([i>=4.056 for i in smo[GD['rpoC']['locus_tag'][0]].values()])##P = 0.00031 padj = 0.00242
my3bpraw['fakA']=np.count_nonzero([i>=3 for i in smo[GD['fakA']['locus_tag'][0]].values()])######P = 0.00084 padj = 0.00737
my3bpraw['ptsG']=np.count_nonzero([i>=1.969 for i in smo[GD['ptsG']['locus_tag'][0]].values()]) #P = 0.02367 padj = 0.25532 <-- So the Bonferroni corrected is no longer significant.
my3bpraw['atpD']=np.count_nonzero([i>=3.052 for i in smo[GD['atpD']['locus_tag'][0]].values()]) #P < 0.00001 padj < 0.00011
my3bpraw['tetM']=np.count_nonzero([i>=2.877 for i in smo[GD['tetM']['locus_tag'][0]].values()])##P = 0.00141 padj = 0.01408
my3bpraw['amiF']=np.count_nonzero([i>=2.596 for i in smo[GD['amiF']['locus_tag'][0]].values()])#P = 0.00121  padj = 0.01463
my3bpraw['JCVISYN3A_0430']=np.count_nonzero([i>=2 for i in smo['JCVISYN3A_0430'].values()])##############padj = 0.00768
my3bpraw['gyrB']=np.count_nonzero([i>=0.325 for i in smo[GD['gyrB']['locus_tag'][0]].values()])#P = 0.23633
my3bpraw['clsA']=np.count_nonzero([i>=2 for i in smo[GD['clsA']['locus_tag'][0]].values()])
my3bpraw['pgpA']=np.count_nonzero([i>=0.295 for i in smo[GD['pgpA']['locus_tag'][0]].values()])
my3bpraw['JCVISYN3A_0373']=np.count_nonzero([i>=2 for i in smo['JCVISYN3A_0373'].values()])
my3bpraw['polC']=np.count_nonzero([i>=1.378 for i in smo[GD['polC']['locus_tag'][0]].values()])
my3bpraw['relA']=np.count_nonzero([i>=1.306 for i in smo[GD['relA']['locus_tag'][0]].values()])
my3bpraw['clpB']=np.count_nonzero([i>=1.143 for i in smo[GD['clpB']['locus_tag'][0]].values()])
my3bpraw['lgt']=np.count_nonzero([i>=2.082 for i in smo[GD['lgt']['locus_tag'][0]].values()])
my3bpraw['nadK']=np.count_nonzero([i>=0.136 for i in smo[GD['nadK']['locus_tag'][0]].values()])
my3bpraw['rpsC']=np.count_nonzero([i>=1.054 for i in smo[GD['rpsC']['locus_tag'][0]].values()])
my3bpraw['dnaB']=np.count_nonzero([i>=0.214 for i in smo[GD['dnaB']['locus_tag'][0]].values()])
my3bpraw['JCVISYN3A_0691']=np.count_nonzero([i>=1.906 for i in smo['JCVISYN3A_0691'].values()])
#%%
mys1praw_nonsyn=dict()
mys1praw['dnaA_1']=np.count_nonzero([i>=2.575 for i in smo1[GD1['dnaA_1']['locus_tag'][0]].values()])
mys1praw['tnpA_1']=np.count_nonzero([i>=2 for i in smo1[GD1['tnpA_1']['locus_tag'][0]].values()])
mys1praw['tnpB_1']=np.count_nonzero([i>=2 for i in smo1[GD1['tnpB_1']['locus_tag'][0]].values()])
mys1praw['MMSYN1_0030']=np.count_nonzero([i>=1.085 for i in smo1['MMSYN1_0030'].values()])
mys1praw['MMSYN1_0032']=np.count_nonzero([i>=0.993 for i in smo1['MMSYN1_0032'].values()])
mys1praw['rpsG']=np.count_nonzero([i>=0.172 for i in smo1[GD1['rpsG']['locus_tag'][0]].values()])

mys1praw['MMSYN1_0187']=np.count_nonzero([i>=1.050 for i in smo1['MMSYN1_0187'].values()])
mys1praw['lpdA']=np.count_nonzero([i>=5.183 for i in smo1[GD1['lpdA']['locus_tag'][0]].values()])
mys1praw['MMSYN1_0253']=np.count_nonzero([i>=0.576 for i in smo1['MMSYN1_0253'].values()])
mys1praw['MMSYN1_0414']=np.count_nonzero([i>=0.244 for i in smo1['MMSYN1_0414'].values()])
mys1praw['MMSYN1_0469']=np.count_nonzero([i>=0.932 for i in smo1['MMSYN1_0469'].values()])
mys1praw['MMSYN1_0490']=np.count_nonzero([i>=0.461 for i in smo1['MMSYN1_0490'].values()])
mys1praw['ftsZ']=np.count_nonzero([i>=3.079 for i in smo1[GD1['ftsZ']['locus_tag'][0]].values()])
mys1praw['parC']=np.count_nonzero([i>=0.907 for i in smo1[GD1['parC']['locus_tag'][0]].values()])
mys1praw['MMSYN1_0641']=np.count_nonzero([i>=3.606 for i in smo1['MMSYN1_0641'].values()])
mys1praw['rpoA']=np.count_nonzero([i>=1.528 for i in smo1[GD1['rpoA']['locus_tag'][0]].values()])
mys1praw['rpoB']=np.count_nonzero([i>=1.046 for i in smo1[GD1['rpoB']['locus_tag'][0]].values()])
mys1praw['prs']=np.count_nonzero([i>=0.907 for i in smo1[GD1['prs']['locus_tag'][0]].values()])
mys1praw['MMSYN1_0892']=np.count_nonzero([i>=1.123 for i in smo1['MMSYN1_0892'].values()])
mys1praw['MMSYN1_0898']=np.count_nonzero([i>=1.184 for i in smo1['MMSYN1_0898'].values()])
mys1praw['tetM']=np.count_nonzero([i>=1.546 for i in smo1[GD1['tetM']['locus_tag'][0]].values()])
mys1praw['pyk']=np.count_nonzero([i>=1.284 for i in smo1[GD1['pyk']['locus_tag'][0]].values()])
mys1praw['MMSYN1_0339']=np.count_nonzero([i>=1.302 for i in smo1['MMSYN1_0339'].values()])
mys1praw['rpoD']=np.count_nonzero([i>=1.000 for i in smo1[GD1['rpoD']['locus_tag'][0]].values()])
mys1praw['MMSYN1_0412']=np.count_nonzero([i>=1.077 for i in smo1['MMSYN1_0412'].values()])
mys1praw['MMSYN1_0460']=np.count_nonzero([i>=1.335 for i in smo1['MMSYN1_0460'].values()])
mys1praw['MMSYN1_0751']=np.count_nonzero([i>=2 for i in smo1['MMSYN1_0751'].values()])
mys1praw['MMSYN1_0257']=np.count_nonzero([i>=0.868 for i in smo1['MMSYN1_0257'].values()])
mys1praw['MMSYN1_0312']=np.count_nonzero([i>=0.130 for i in smo1['MMSYN1_0312'].values()])
mys1praw['MMSYN1_0471']=np.count_nonzero([i>=1.467 for i in smo1['MMSYN1_0471'].values()])
mys1praw['MMSYN1_0567']=np.count_nonzero([i>=0.641 for i in smo1['MMSYN1_0567'].values()])
mys1praw['MMSYN1_0691']=np.count_nonzero([i>=0.903 for i in smo1['MMSYN1_0339'].values()])

#%%
p_3b_100000 = pd.DataFrame.from_dict(my3bpraw,orient='index')
p_s1_100000 = pd.DataFrame.from_dict(mys1praw,orient='index')

#p_3b_100000.to_csv(index=True,path_or_buf=r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\pvals_simulations_syn3B_100000.csv")
#p_s1_100000.to_csv(index=True,path_or_buf=r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\pvals_simulations_syn1.0_100000.csv")
p_3b_100000.to_csv(index=True,path_or_buf=r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\pvals_simulations_syn3B_GC_only.synonymous_100000.csv")
p_s1_100000.to_csv(index=True,path_or_buf=r"C:\Users\rmoge\OneDrive - Indiana University\Mycoplasma_OneDrive\Strains\pvals_simulations_syn1.0_GC_only.synonymous_100000.csv")


#use similar code to output the nonsynonymous tables
#%%
def compare_recs(gb_record1,gb_record2):
    answer = dict()
    bads=list()
    count=0
    for (index, feature) in enumerate(gb_record1.features) :
        if feature.type=='CDS':
            if "locus_tag" in feature.qualifiers :
                lt1 = feature.qualifiers['locus_tag']
                num1 = lt1[0][-5:]
                #print(num1)
                try:
                    psq1 = feature.qualifiers['translation']
                except KeyError:
                    pass
                #try:
                 #   print(psq)
                #except UnboundLocalError:
                    #print("")
                for (index2,feature2) in enumerate(gb_record2.features):
                    if feature2.type=='CDS':
                        if "locus_tag" in feature2.qualifiers:
                            lt2 = feature2.qualifiers['locus_tag']
                            num2 = lt2[0][-5:]
                            if num1 == num2:
                                #print('match!')
                                count+=1
                                try:
                                    psq2 = feature2.qualifiers['translation']
                                    if psq1 != psq2:
                                        print('not the same')
                                        bads+=lt1
                                    else:
                                        print('same')
                                except:
                                    pass
                                       
                                       
    return bads

#test = index_genbank_features(myrec, "CDS", "gene")
cd=compare_recs(myrec, syn1rec)
print(cd)