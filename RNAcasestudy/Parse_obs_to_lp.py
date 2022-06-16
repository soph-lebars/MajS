#!bin/bash/python3

import os, sys
from clyngor import ASP, solve


import glob

file_obs=glob.glob("*.obs")

print(file_obs)

path = "./instance/"
def converttoinstance(network,obs):

	net=open(network)
	obs=open(obs)
	

	Vertex=set()
	EdgeV=set()
	ObservedV=set()

	for el in net:	
		if '!'in el:
			el=el.replace("!","")
			
			signe=-1
		else :
			signe =1
		line=el.split("->")
		
		parent=line[0].lower().strip()
		
		enfant=line[1].lower().strip()
		if parent not in Vertex:
			
			Vertex.add(parent)
		if enfant not in Vertex:
			
			Vertex.add(enfant)
		edge="observedE("+parent+","+enfant+","+str(signe)+")"
		EdgeV.add(edge)
		
	for el in obs: 
		print(el)
		line=el.split("=")
		node=line[0].lower().strip()
		sign=line[1].strip()
		if "-" in sign:
			sign2="-1"
		if "+" in sign:
			sign2="1"
		if "0" in sign : 
			sign2="0"
			
		observ="observedV("+node+","+str(sign2)+",100)"
		ObservedV.add(observ)


	observedE=f'{".".join([v for v in EdgeV])}.'


	observedV=f'{".".join([v for v in ObservedV])}.'


	vertex =f'vertex({";".join([v for v in Vertex])}).'
	
	return observedE+observedV+vertex


network= "regfromdagRNA.cif" 


for i in range(0,len(file_obs)) :
	
	name_obs =  file_obs[i]
	
	res= converttoinstance(network,name_obs)
	file_name=file_obs[i].replace(".obs","")
	filename='%s.lp' % file_name
	o = open(os.path.join(path, filename), 'w')
	
	o.write(res)
	
	o.close()


