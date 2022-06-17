BENCHMARK_NB=$1

hc_values=($(seq 0.01 0.01 0.1))
lc_values=($(seq 0.1 0.1 1))

if [ ! -d "experimentations/standard_deviation/" ] ; then
    mkdir "experimentations/standard_deviation/" && mkdir "experimentations/standard_deviation/data"
fi

for hc_val in "${hc_values[@]}"
do
    for lc_val in "${lc_values[@]}"
    do
        hc_val=$(echo $hc_val | tr ',' '.')
        lc_val=$(echo $lc_val | tr ',' '.')
        outdir="experimentations/standard_deviation/data/hc_${hc_val}_lc_${lc_val}"
        if [ ! -d $outdir ] ; then
            mkdir $outdir
        fi
        if [ ! -d $outdir/benchmark$BENCHMARK_NB ] ; then
            mkdir $outdir/benchmark$BENCHMARK_NB
        fi
        
        cmd="python MIComp.py --majs-file ../out/Df/Df_Benchmark$BENCHMARK_NB.csv --iggy-file ../Iggy_RNA_benchmark$BENCHMARK_NB.txt --obs-file ../RNA_benchmark$BENCHMARK_NB.obs --rna-file dataRNA.csv --export  --outdir $outdir/ -lc $lc_val -hc $hc_val > $outdir/benchmark$BENCHMARK_NB/out.log"
        echo $cmd
        eval $cmd
    done
done