#!bin/bash/python3

import csv

import os, sys

import os.path
from operator import itemgetter


import pandas as pd

from tabulate import tabulate




DF=pd.read_csv(sys.argv[1]) #Dataframe produce by Majs located in *casesstudy/out/df format csv



filename=os.path.basename(sys.argv[1]).replace(".csv","").replace("Df_","")
print(filename)


Iggyout=open(sys.argv[2]) #Output produce by Iggy located in *casestudy/ format Iggy_*.txt





obs=open(sys.argv[3]) # Observation used for the benchmark located in *casestudy/ format Benchmark*.obs



df=DF.set_index("Name")

RNAmymethod=[]
compwithother=[]

observation=[]
for ob in obs :
	ob=ob.split("=")
	observation.append(ob[0].strip())

#print(observation)




#print("====================================================")


MajS=[]

for i in range(0,df.shape[0]):

	name=df.index[i].upper()
	
	sign=df.loc[df.index[i],"SignMaj"].replace(" ","")
	
	if "0-" in sign  :
		sign="notPlus"
	if "+0" in sign:
		sign="notMinus"
	MajS.append([name,sign])
	

MajS=sorted(MajS, key=itemgetter(0))

#print(MajS)
		
		
#print("====================================================")
RNAiggy=[]
for line in Iggyout:
	line=line.strip()
	line=line.split("=")
	
	name=line[0].strip()
	sign=line[1].strip()
	RNAiggy.append([name,sign])
	
	

RNAiggy=sorted(RNAiggy, key=itemgetter(0))
#print(RNAiggy)




comparedata={}


RNAiggy_set = set(tuple(x) for x in RNAiggy)
MajS_set = set(tuple(x) for x in MajS)



Nbcommun=len(RNAiggy_set.intersection(MajS_set))

if len(RNAiggy_set) < len(MajS_set):

	Nbdiff=len(RNAiggy_set.difference(MajS_set))
	
else :


	Nbdiff=len(MajS_set.difference(RNAiggy_set))



#print("=================================================")


#print("nombre de noeuds en commun : ", Nbcommun, "nombre de noeuds différents : ", Nbdiff)



diff=[]

for el in RNAiggy_set.difference(MajS_set):

	diff.append(el[0])
#print("liste des noeuds différents : ", diff)


pred=[]
for el in MajS:
	pred.append(el[0])

predmymethod=len(list(set(pred) - set(observation)))

coverage= int(predmymethod)/(81-len(observation))*100


#print("nombre de noeud prédit non observé MajS : " ,predmymethod , " coverage MajS : ", coverage )


prediggy=[]
for el in RNAiggy:
	prediggy.append(el[0])

prediggytot=len(list(set(prediggy)-set(observation)))


coverageiggy= round(int(prediggytot)/(81-len(observation))*100)

#print("nombre de noeud prédit non observé Iggy: " ,prediggytot, " coverage Iggy : ", coverageiggy )


content2=[["Predicted node by MajS", predmymethod, "Predicted node by iggy", prediggytot],
["Coverage of predicted node by MajS", str(coverage)+" %", "Coverage of predicted node by Iggy", str(coverageiggy)+" %"],
["Number of predicted node in common between Iggy and MajS", Nbcommun-len(observation), "Number of predicted node different between Iggy and MajS",Nbdiff]]

table=tabulate(content2, tablefmt="tsv")
print(tabulate(content2, tablefmt="fancy_grid"))

text_file=open("../RNAcasestudy/out/results/table_RNA_MajS_"+filename+".csv","w")
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

dfalldiff.to_csv("../RNAcasestudy/out/results/df_diff_iggy_MajS_"+filename+".csv",index=True)	
print(dfalldiff)
