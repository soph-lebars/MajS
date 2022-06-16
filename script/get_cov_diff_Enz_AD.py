#!bin/bash/python3

''''
Working on AD case study will compare iggy and MajS on these case study for the different Benchmarks and will givev out tables of comparison between them benchmark

input, argument to pass for each benchmark : 
- data 1 :Dataframe of a benchmark given by MajS they are in ../ADcasestudy/out/Df/ 
- data 2 : Iggy output of same benchmark which are in ../ADcasestudy/ in the form of Iggy_*.txt
- data 3 : Observation file of same benchmark in  ../ADcasestudy/ in the form of *.obs

Output :

- Table 1 : Difference between MajS and Iggy in term of predicted node, coverage, node in commun and different
- Table 2 : All the different node between MajS and Iggy
- Table 3 : All the enzymes predicted by Iggy and/or MajS (only Majs for Benchmark1_minus pay attention to line 220:224



'''

import csv

import os, sys


from operator import itemgetter

import numpy as np

import pandas as pd

from tabulate import tabulate


data1=pd.read_csv(sys.argv[1]) #The data frame obtained by MajS

df=data1.set_index("Name")


data2=open(sys.argv[2])


obs=open(sys.argv[3])


RNAmymethod=[]
compwithother=[]

observation=[]
for ob in obs :
	ob=ob.split("=")
	observation.append(ob[0].strip())



MajS=[]

for i in range(0,df.shape[0]):

	name=df.index[i].upper()
	
	sign=df.loc[df.index[i],"SignMaj"].replace(" ","")
	
	if "0-" in sign  :
		sign="notPlus"
	if "+0" in sign:
		sign="notMinus"
	if "+-" in sign:
		sign="CHANGE"
	MajS.append([name,sign])
	

MajS=sorted(MajS, key=itemgetter(0))

RNAiggy=[]
for line in data2:
	line=line.strip()
	line=line.split("=")
	
	name=line[0].strip()
	sign=line[1].strip()
	RNAiggy.append([name,sign])
	
	

RNAiggy=sorted(RNAiggy, key=itemgetter(0))





comparedata={}


RNAiggy_set = set(tuple(x) for x in RNAiggy)
MajS_set = set(tuple(x) for x in MajS)



Nbcommun=len(RNAiggy_set.intersection(MajS_set))

if len(RNAiggy_set) < len(MajS_set):

	Nbdiff=len(RNAiggy_set.difference(MajS_set))
	
else :


	Nbdiff=len(MajS_set.difference(RNAiggy_set))





diff=[]

for el in RNAiggy_set.difference(MajS_set):

	diff.append(el[0])



pred=[]
for el in MajS:
	pred.append(el[0])

predmymethod=len(list(set(pred) - set(observation)))

coverage= int(predmymethod)/(94-len(observation))*100



prediggy=[]
for el in RNAiggy:
	prediggy.append(el[0])

prediggytot=len(list(set(prediggy)-set(observation)))


coverageiggy= round(int(prediggytot)/(94-len(observation))*100)



content2=[["Predicted node by MajS", predmymethod, "Predicted node by iggy", prediggytot],
["Coverage of predicted node by MajS", str(coverage)+" %", "Coverage of predicted node by Iggy", str(coverageiggy)+" %"],
["Number of predicted node in common between Iggy and MajS", Nbcommun-len(observation), "Number of predicted node different between Iggy and MajS",Nbdiff]]

table=tabulate(content2, tablefmt="tsv")
print(tabulate(content2, tablefmt="fancy_grid"))



name=os.path.basename(sys.argv[1]).replace("Df_","")

text_file=open("../ADcasestudy/out/results/table_AD_MajS_"+name,"w")
text_file.write(table)
text_file.close()




diff= [each_string.lower() for each_string in diff]


dfdiff=df.loc[diff]



		
dfdiff.index = dfdiff.index.str.upper()


dfRNAiggy=pd.DataFrame(RNAiggy, columns =["Name","SignIggy"])


dfRNAiggy=dfRNAiggy.set_index("Name")

diff= [each_string.upper() for each_string in diff]

dfdiffiggy=dfRNAiggy.loc[diff]
	


dfalldiff = pd.concat([dfdiff, dfdiffiggy], axis=1)

dfalldiff=dfalldiff.sort_index(axis=0)

dfalldiff.to_csv("../ADcasestudy/out/results/df_diff_iggy_MajS_"+name,index=True)	
print(dfalldiff)


enzyme=['LDHA','SLC2A1','PFKL','HK2','HK3','ENO1','ENO3','ALDOA','GAPDH','PGK1','ENO2','HK1','PDHA2','PDHB','PDHA1'] 





enzyme= [each_string.lower() for each_string in enzyme]


dfenzyme=df.loc[enzyme]



dfenzyme.index = dfenzyme.index.str.upper()


dfRNAiggy=pd.DataFrame(RNAiggy, columns =["Name","SignIggy"])


dfRNAiggy=dfRNAiggy.set_index("Name")

enzyme= [each_string.upper() for each_string in enzyme]


#NaN = np.nan


print("YOU HAVE TO DO SOMETHING for line 220:224 if benchmark1_minus only are enzyme not observed by Iggy")
dfenzymeiggy=dfRNAiggy.loc[enzyme]  #COMMENT for Benchmark1_minus as enzyme not present in IGGY 

dfallenzyme=pd.concat([dfenzyme,dfenzymeiggy], axis =1)  #COMMENT for Benchmark1_minus as enzyme not present in IGGY 

#dfallenzyme=dfenzyme  #DECOMMENT for Benchmark1_minus as enzyme not present in IGGY 

dfallenzyme=dfallenzyme.sort_index(axis=0)

dfallenzyme.to_csv("../ADcasestudy/out/results/df_Enz_iggy_MajS_"+name,index=True)   
print(dfallenzyme)

