# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 15:42:46 2019

@author: rmoge
"""

#%%
from __future__ import division
import re, os,sys, math, operator,random, copy,collections; import numpy as np; import pandas as pd
from itertools import groupby; import pprint as pp
import matplotlib.pyplot as plt
#%%
def calc_targetsize(p=1/6,q=1/6,r=1/6,s=1/6,t=1/6,possATtoCG=0,possATtoGC=0,possATtoTA=0,possCGtoGC=0,possCGtoTA=0,possGCtoTA=0):
    u=1-p-q-r-s-t
    NorS=6*(p*possATtoCG + q*possATtoGC + r*possATtoTA + s*possCGtoGC + t*possCGtoTA + u*possGCtoTA)
    return NorS

def calc_dNdS(N=1,S=1,DN=0,DS=0):
    dN = -3/4*math.log(1-(4/3*DN/N))
    dS = -3/4*math.log(1-(4/3*DS/S))
    print("dN = " + str(dN)); print("dS = " + str(dS))
    print("dN/dS = " + str(dN/dS))
    return dN/dS

def calc_dNdS_fromtuple_indiv_pseudocount(intup):
    N=intup[0];S=intup[1];DN=intup[2];DS=intup[3]
    DS_dummy = 0#turned this feature off because I am going to 
    dN = -3/4*math.log(1-(4/3*DN/N))
    dS = -3/4*math.log(1-(4/3*DS/S))
    dS_dummy = -3/4*math.log(1-(4/3*DS_dummy/S))
    #print("dN = " + str(dN)); print("dS = " + str(dS))
    try:
        meow = dN/dS
        print(dN/dS)
    except ZeroDivisionError:
        meow = dN/dS_dummy
        print("dN/dS is undefined because dS = 0.  Adding a dummy synonymous substitution so that DS=1 gives " + str(dN/dS_dummy))
    #print("dN/dS = " + str(dN/dS))
    return meow

def calc_dNdS_fromtuple(intup):
    N=intup[0];S=intup[1];DN=intup[2];DS=intup[3]
    dN = -3/4*math.log(1-(4/3*DN/N))
    dS = -3/4*math.log(1-(4/3*DS/S))
    #print("dN = " + str(dN)); print("dS = " + str(dS))
    #print("dN/dS = " + str(dN/dS))
    return dN/dS


#calc_dNdS_fromtuple(mm9)
#%%
#mydnds = calc_dNdS(1000,1000,2,2)

#%%
def freq_list(infile=""):
    mydf = pd.read_csv(infile)
    mydf=mydf.fillna(-9)
    return mydf          


            
syn3B_N = freq_list(r"~\Documents\GitHub\MinimalCell\datafiles\3B_N_DNDS_freqs.csv")######
syn3B_S = freq_list(r"~\Documents\GitHub\MinimalCell\datafiles\3B_S_DNDS_freqs.csv")######
syn1_N = freq_list(r"~\Documents\GitHub\MinimalCell\datafiles\1.0_N_DNDS_freqs.csv")######
syn1_S = freq_list(r"~\Documents\GitHub\MinimalCell\datafiles\1.0_S_DNDS_freqs.csv")######
#values were calculated ignoring stop codons

def calc_number(indf):
    thedict=dict()
    for col in indf.iteritems():
        thedict[col[0]]=col[1]
    outcomedict=dict()
    for name,l in thedict.items():
        counter=0
        for freq in l:
            r=random.random()
            #pp.pprint(r);pp.pprint(freq)
            if freq > r:
                counter+=1
            else:
                pass
        outcomedict[name]=counter    
    return outcomedict
 #%%
  
NDICT_3B=calc_number(syn3B_N)
SDICT_3B=calc_number(syn3B_S)
NDICT_1=calc_number(syn1_N)
SDICT_1=calc_number(syn1_S)
pp.pprint(SDICT_3B);pp.pprint(SDICT_1)#check whether any of the lines have S = 0
#m13 has S = 0. Therefore, need to add a pseudocount of S = 1 for all of the samples
#%%
Nmin=calc_targetsize(p=0.013435701,q=0.015355086,r=0.038387716,s=0.009596929,t=0.654510557,possATtoCG=303273,possATtoGC=204218,possATtoTA=265192,possCGtoGC=113261,possCGtoTA=96506,possGCtoTA=98885)
Smin=calc_targetsize(p=0.013435701,q=0.015355086,r=0.038387716,s=0.009596929,t=0.654510557,possATtoCG=52292,possATtoGC=157137,possATtoTA=61057,possCGtoGC=2053,possCGtoTA=11721,possGCtoTA=3100)
Nnonmin=calc_targetsize(p=0.011049724,q=0.072375691,r=0.037569061,s=0.10441989,t=0.419337017,possATtoCG=585378,possATtoGC=397994,possATtoTA=513439,possCGtoGC=215178,possCGtoTA=180519,possGCtoTA=186790)
Snonmin=calc_targetsize(p=0.011049724,q=0.072375691,r=0.037569061,s=0.10441989,t=0.419337017,possATtoCG=98053,possATtoGC=296696,possATtoTA=113194,possCGtoGC=4631,possCGtoTA=25669,possGCtoTA=6581)

essN=calc_targetsize(p=0.011049724,q=0.072375691,r=0.037569061,s=0.10441989,t=0.419337017,possATtoCG=303273,possATtoGC=204218,possATtoTA=265192,possCGtoGC=113261,possCGtoTA=96506,possGCtoTA=98885)
essS=calc_targetsize(p=0.011049724,q=0.072375691,r=0.037569061,s=0.10441989,t=0.419337017,possATtoCG=52292,possATtoGC=157137,possATtoTA=61057,possCGtoGC=2053,possCGtoTA=11721,possGCtoTA=3100)
neN=calc_targetsize(p=0.011049724,q=0.072375691,r=0.037569061,s=0.10441989,t=0.419337017,possATtoCG=585378-303273,possATtoGC=397994-204218,possATtoTA=513439-265192,possCGtoGC=215178-113261,possCGtoTA=180519-96506,possGCtoTA=186790-98885)
neS=calc_targetsize(p=0.011049724,q=0.072375691,r=0.037569061,s=0.10441989,t=0.419337017,possATtoCG=98053-52292,possATtoGC=296696-157137,possATtoTA=113194-61057,possCGtoGC=4631-2053,possCGtoTA=25669-11721,possGCtoTA=6581-3100)

mm9 = (Nmin,Smin,NDICT_3B["9"],(1+SDICT_3B["9"]))#adding a pseudocount of 1 to all samples because some of them have 0 fixed synonymous changes
mm10 = (Nmin,Smin,NDICT_3B["10"],(1+SDICT_3B["10"]))
mm11 = (Nmin,Smin,NDICT_3B["11"],(1+SDICT_3B["11"]))
mm13 = (Nmin,Smin,NDICT_3B["13"],(1+SDICT_3B["13"]))
mm9v=calc_dNdS_fromtuple(mm9)
mm10v=calc_dNdS_fromtuple(mm10)
mm11v=calc_dNdS_fromtuple(mm11)
mm13v=calc_dNdS_fromtuple(mm13)
pp.pprint(mm9v);pp.pprint(mm10v);pp.pprint(mm11v);pp.pprint(mm13v)
#0.9301089372678766
#1.7273521191214276
#1.0629829533148558
#3.18896319067647

mm1=(Nnonmin,Snonmin,NDICT_1["1"],(1+SDICT_1["1"]))
mm3=(Nnonmin,Snonmin,NDICT_1["3"],(1+SDICT_1["3"]))
mm4=(Nnonmin,Snonmin,NDICT_1["4"],(1+SDICT_1["4"]))
mm6=(Nnonmin,Snonmin,NDICT_1["6"],(1+SDICT_1["6"]))
mm1v=calc_dNdS_fromtuple(mm1)
mm3v=calc_dNdS_fromtuple(mm3)
mm4v=calc_dNdS_fromtuple(mm4)
mm6v=calc_dNdS_fromtuple(mm6)
pp.pprint(mm1v);pp.pprint(mm3v);pp.pprint(mm4v);pp.pprint(mm6v)
#1.3974940550995192
#0.5502594533166861
#0.6550732766684936
#1.5721807194227904
#%%

#%%
####DEPRECATED CELL
#values for N and S were calculated using gdtools COUNT -b. This is not the most exact possible value because i) mycoplasma uses a modified codon table wherein UGA codes for Trp, not STOP ii) we can incorporate data from the mutation spectrum to improve our predictions of syn and nonsyn target sizes---see Yang 2006 book, pg 51
#However, this may not be necessary. As indicated by M. W. Hahn (pers. comm.) the values for N and S generated by the different methods are highly correlated. This is especially true because all of the organisms are so closely related.
mm9 = (1081335,287360,11,3)
mm10 = (1081335,287360,10,1)#the true value is 0. I've added 1 dummy synonymous difference so that I don't divide by 0
mm11 = (1081335,287360,13,2)
mm13 = (1081335,287360,10,1)#this true synonymous value is also 0
#values were calculated ignoring stop codons

calc_dNdS_fromtuple(mm9)
calc_dNdS_fromtuple(mm10)
calc_dNdS_fromtuple(mm11)
calc_dNdS_fromtuple(mm13)

mm1=(2079386,544861,5,2)
#mm2=(2079386,544861,0,0);calc_dNdS_fromtuple(mm2)
mm3=(2079386,544861,8,1)#this 1 is false. actually 0.
mm4=(2079386,544861,5,1)
mm6=(2079386,544861,10,1)

calc_dNdS_fromtuple(mm1)
calc_dNdS_fromtuple(mm3)
calc_dNdS_fromtuple(mm4)
calc_dNdS_fromtuple(mm6)
####END DEPRECATED CELL

#%%
#Addressing concerns of R3
#First, the essential partition of the genome
#Values from the DNDS_freqs_COPY.FOR.R3.xlsx file
#Pseudocount of 1 added to all S values
mm1e=(essN,essS,8,3)
mm3e=(essN,essS,12,3)
mm4e=(essN,essS,8,3)
mm6e=(essN,essS,8,1)
#now, the nonessential partition of the genome
#Values from the DNDS_freqs_COPY.FOR.R3.xlsx file
#Pseudocount of 1 added to all S values
mm1ne=(neN,neS,5,1)
mm3ne=(neN,neS,7,6)
mm4ne=(neN,neS,7,4)
mm6ne=(neN,neS,6,2)

mm1ve=calc_dNdS_fromtuple(mm1e)
mm3ve=calc_dNdS_fromtuple(mm3e)
mm4ve=calc_dNdS_fromtuple(mm4e)
mm6ve=calc_dNdS_fromtuple(mm6e)

mm1vne=calc_dNdS_fromtuple(mm1ne)
mm3vne=calc_dNdS_fromtuple(mm3ne)
mm4vne=calc_dNdS_fromtuple(mm4ne)
mm6vne=calc_dNdS_fromtuple(mm6ne)

pp.pprint(mm1ve);pp.pprint(mm3ve);pp.pprint(mm4ve);pp.pprint(mm6ve)
pp.pprint(mm1vne);pp.pprint(mm3vne);pp.pprint(mm4vne);pp.pprint(mm6vne)