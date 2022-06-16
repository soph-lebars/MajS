# MajS

---

## Installation to reproduce study results

To use MajS :  
`git clone` in a repository created on your computer
> Install clingo : version 5.5.0 was used for this study the installation for ubuntu can be found on : (https://launchpad.net/~potassco) <br>
> Install python : version 3.8.5 was used for this study 

> The python libraries used in this study are : os,sys,numpy,pandas, statistic, os.path, glob, operator, tabulate

To compare with Iggy :
> Install Iggy : version 2.1.1 was used for this study, installation can be found at : (https://github.com/bioasp/iggy/releases)

---
## Toycasestudy



`cd BMC_MajS/toycasestudy | ls`

You will find the instance used for toynetwork study, the vertex, edges and the partial observation in **instancetoynetwork.lp**

You can run the logical rules on this instance and obtained the optimal solutions by doing :

`./get_otimal_solution_toynet.sh` 

In more details, the command line in this bash script is :

`clingo -n 0 script/logical_rules.lp toycasestudy/instancetoynetwork.lp --const k=2 --opt-mode=optN,2 --project -W none`

Here we are using clingo by given our logical rules in lp format and the instance of our problem, we fix the added artificial parent with --const k, here k=2 and we want to have only optimal solutions that have only 2 added artificial parents so we use --opt-mode=optN,2 the --project is here so that we don't have duplicate prediction (the -W none, is used to prevent warning on the computation of weight in logical rules) 


The results of this bash script is present in **out_optsolution_toynet.txt**, you can find these optimal solutions in the toy case study table in Figure 1.


---
## AD case study


# input
`cd BMC_MajS/ADcasestudy/` 

> In obs format, you can find all the partial observation of benchmarks used in AD case study,  these observation were obtained as explain in Section methods of our article, the script to get this observation can be found in the gitlab of our precedent study : (https://gitlab.univ-nantes.fr/E19D080G/comparing_iggy_prob.git)

> In cif format, you can find the IG derived from HIF-1 signaling pathway and AD dataset used for this case study, for a more detail explanation and the script please refer to Section Methods and our precedent article, link just above.


> The file Iggy_*.txt are the output of Iggy, predicted nodes for each benchmark it is possible to obtain them directly by running :

> `./Iggy_output.sh`

> The instance used for AD case are all in instance repository `BMC_MajS/ADcasestudy/instance` they were obtained by running the script :

> `python Parse_obs_to_lp.py ` that take as input all observation files *.obs and the cif file and will provide instance for all benchmarks

# output
All the script to produce the outputs are in :

`cd BMC_MajS/script/`

All ouput for AD case study are in :

`cd BMC_MajS/ADcasestudy/out` 

The file *.txt are all the optimal solutions given by the logical rules. To have them, you can run the bash script in directory script
by doing :

`./MajS.sh Benchmark*.lp`  do not forget to **change RNA to AD** if you are in AD case study.  It will take as argument one of the instance benchmark in lp format that can be found in BMC_MajS/ADcasestudy/instance/. 
 

In out directory, you can find a new directory which is Df in this one you have the projection of all the optimal solution for each benchmark in the csv format this were obtain by running the python script in script directory.

`Projection_df.py` do not forget to **change RNA to AD** it will takes all of the txt files which are optimal solutions for all AD case study benchmarks in directory out and it wil return the dataframe representing the projection of all optimal solutions.

At last but not least, in the out directory you will find the results directory  in which you can find many of the tables for AD case study that were used in our article produce the Table 3. To obtain this table, you have to run in script directory :

`python get_cov_diff_Enz_AD.py` 


will take as arguments :
<ul>
<li> data 1 :Dataframe of a benchmark given by MajS they are in BMC_MajS/ADcasestudy/out/Df/ </li>
<li> data 2 : Iggy output of same benchmark which are in BMC_MajS/ADcasestudy/ in the form of Iggy_*.txt</li>
<li> data 3 : Observation file of same benchmark in  BMC_MajS/ADcasestudy/ in the form of *.obs</li>
</ul>

Will give as output :

<ul>
<li> Table 1 : Difference between MajS and Iggy in term of predicted node, coverage, node in commun and different</li>
<li> Table 2 : All the different node between MajS and Iggy </li>
<li> Table 3 : All the enzymes predicted by Iggy and/or MajS </li>

</ul>


# plot


To reproduce Figure 4 of our article :

`cd BMC_MajS/plot`

The figure is called comp_prob_majs_enz.png it was obtained with the prediction data for each enzyme in probregnet across three perturbation of HIF1A, this data are : res_prediction_enzyme_prob.csv 
> You can found this data in our precedent article
The other data used were the prediction of each enzyme across all perturbation with MajS (obtained by running get_cov_diff_Enz_AD.py) and joining all benchmarks and removed Iggy prediction.
The command line to obtain this plot is :

`Rscript Enzprobplot.R`

## RNA case study

# input

`cd BMC_MajS/RNAcasestudy/`

> All files in this directory and script are the same than for AD case study but for RNA, there is the two file observation in obs format to perform the 2 benchmarks used in RNA case study, the network in cif format reduced for RNA dataset the output of Iggy for this two benchmarks called "Iggy_*.txt" we can obtain them by using the Iggy_output.sh (move this script in the repository with mv command). 

> The instance used by MajS are also in instance repository `BMC_MajS/RNAcasestudy/instance` and we can run the script Parse_obs_to_lp.py to obtain them.


# output

All ouput for RNA case study are in :

`cd BMC_MajS/RNAcasestudy/out` 

Once again, the *.txt files are the optimal solutions of each benchmark obtain by running in script repository :

` MajS.sh Benchmark*.lp `  **change AD to RNA** 


The csv files in Df subdirectory of out are obtain by running and represent the projection of all optimal solutions :

`Projection_df.py`**change AD to RNA** 

The csv files in results subdirectory of out are the one used to obtain Table 1 and 2 on our study. You have to run :


`get_table_cov_diff_RNA.py` 


will take as arguments :
<ul>
<li> data 1 :Dataframe of a benchmark given by MajS they are in BMC_MajS/ADcasestudy/out/Df/ </li>
<li> data 2 : Iggy output of same benchmark which are in BMC_MajS/ADcasestudy/ in the form of Iggy_*.txt</li>
<li> data 3 : Observation file of same benchmark in  BMC_MajS/ADcasestudy/ in the form of *.obs</li>
</ul>

Will give as output :

<ul>
<li> Table 1 : Difference between MajS and Iggy in term of predicted node, coverage, node in commun and different</li>
<li> Table 2 : All the different node between MajS and Iggy </li>

</ul>


