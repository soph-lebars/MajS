var1=$1

##########################change to AD or RNA depending on the case study
#### RUN this script to have the row output of MajS for a benchmark pase in argument in lp format that will be  converted to dataframe during the projection phase 

shopt -s lastpipe
var2=$(clingo -n 0 -t 10  logical_rules.lp  ../RNAcasestudy/instance/$var1 --opt-mode=opt -W none --const k=4  | grep  "Optimization" | tail -1 | awk -F': ' '{print $2}' )
#var2 is here to have the exact optimum number for this benchmark it will serve as an argument to have all optimal solutions

echo "$var2"

	nohup clingo -n 0 -t 10 logical_rules.lp ../RNAcasestudy/instance/$var1 --const k=4 --opt-mode=optN,$((var2)) -W none --project > ../RNAcasestudy/out/${var1/lp/txt}   & # Put a function in the background
	#we call our script on the benchmark that we want ex: var1= Benchmark1.lp  whith 10 cores (-t 10) and we are applying all the logical rules of our program (logical_rules.lp) we then fixed k artificial influences (--const k=) and we want to have optimal solution (--opt-mode) and we don't want dduplicate (--projection)


wait 
echo "job "$var1" done"
