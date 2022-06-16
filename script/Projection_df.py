#!bin/bash/python3

import os, sys
import pandas as pd
import numpy as np
import statistics as st
import os.path
import glob

#change to AD or RNA depending on the case study


'''




'''

file_obs=glob.glob("../RNAcasestudy/out/*.txt")


path = "../RNAcasestudy/out/Df"







def data_to_def(file):

	data=open(file)
	dat={}


	def div(a,b):
		if b!=0:
			res=a/b
		else :
			res=0
		return round(res)

	signmajo={}


	def mymean(x):

		try :
			res=round(st.mean(x))
			
		except :
			res=0
		
		
		return res
		
	def myst(x):

		try :
			res=round(st.stdev(x))
		except:
			res=0
		return res








	ct=0
	for line in data :
		if "Answer: 1" in line:
			ct=ct+1	
		
		
		if "pred" in line and ct>1:
			line=line.strip().replace("(","").replace(")","")
			line=line.split("pred")
			for el in line:
				
				if el!="":
				
					el=el.split(",")
					key=el[0]
					
					if key in dat.keys():
						
						dat[key]=dat[key]+","+el[1].strip()+","+el[2].strip()
					else:
						
						dat[key]=el[1].strip()+","+el[2].strip()
						


	df = pd.DataFrame(data=[val.split(",") for val in dat.values()], index=dat.keys())

	dfSignMaj = pd.DataFrame(columns=["Name","SignMaj","Sign1","Sign-1","Sign0"])	
	for i in range(0,df.shape[0]):
		ctplus=0
		weightplus=[]
		ctmoins=0
		weightmoins=[]
		ctzero=0
		weightzero=[]
		k=0
		for j  in range(0,df.shape[1]):
			
			if j%2==0:
				if float(df.iloc[i,j])==1:
					ctplus=ctplus+1
					weightplus.append(int(df.iloc[i,j+1]))
				elif float(df.iloc[i,j])==-1:
					ctmoins=ctmoins+1
					weightmoins.append(int(df.iloc[i,j+1]))
				elif float(df.iloc[i,j])==0:
					ctzero=ctzero+1
					weightzero.append(int(df.iloc[i,j+1]))
		
		dic_count={'+':ctplus,'0':ctzero,'-':ctmoins}
		
		signmaj=max(dic_count.values())
		
		signemajoritaire=""
		
		
		for key,val in dic_count.items():
			if val==signmaj :
				signemajoritaire=signemajoritaire+" "+str(key)
			

		
		dfSignMaj.loc[df.index[i],"Name"]=df.index[i]
		dfSignMaj.loc[df.index[i],"SignMaj"]=signemajoritaire
		dfSignMaj.loc[df.index[i],"Sign1"]=[ctplus,mymean(weightplus),myst(weightplus)]
		dfSignMaj.loc[df.index[i],"Sign-1"]=[ctmoins,mymean(weightmoins),myst(weightmoins)]
		dfSignMaj.loc[df.index[i],"Sign0"]=[ctzero,mymean(weightzero),myst(weightzero)]
			
		signmajo[df.index[i].upper()]=signemajoritaire.strip()
		

	return dfSignMaj
	

for i in range(0,len(file_obs)) :
	
	name_obs =  file_obs[i]
	print(name_obs)
	name=os.path.basename(name_obs).replace(".txt","")
	res= data_to_def(name_obs)

	file_name="Df_"+name
	filename='%s.csv' % file_name
	res.to_csv(os.path.join(path, filename),index=False)
	

