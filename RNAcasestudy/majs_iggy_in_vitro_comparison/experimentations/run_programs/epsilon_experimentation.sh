BENCHMARK_NB=$1

vals=($(seq 0.001 0.001 0.1))

for val in "${vals[@]}"
do
    echo "$val"
    cmd="python MIComp.py --majs-file ../out/Df/Df_Benchmark$BENCHMARK_NB.csv --iggy-file ../Iggy_RNA_benchmark$BENCHMARK_NB.txt --obs-file ../RNA_benchmark$BENCHMARK_NB.obs --rna-file dataRNA.csv --export  -eps $val --outdir experimentations/epsilon/data/$val/"
    cmd_ok=$(echo $cmd | tr ',' '.')
    echo $cmd_ok
    eval $cmd_ok
done