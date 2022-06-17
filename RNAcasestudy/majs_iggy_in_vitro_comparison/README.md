# Comparison of MajS' and Iggy's predictions

The aim of the program `MIComp.py` is to compare discrete predictions of MajS and Iggy methods with continue values (fold change). Using Normal distribution mixture, we compute significance score of both methods to predict a fold change.


## Folder content

This folder contains all data and program for the *in vitro* comparison of Iggy and MajS (Section 2.3.5 in our article). The HUVECS benchmark composed of logFC data is present in `dataRNA.csv`. The `All-logFC_RNA.csv` file contains the complete RNA-Seq dataset used for the construction of the HUVECS benchmark. In the `results` folder, we can observe mixtures of each gene, for both benchmark 1 and 2. The `experimentations` folder contains all data of the experimentations leaded in Supplementary materials of the article.


## Requirements

See `requirements.txt`.

## Usage

```
usage: MIComp.py [-h] --outdir OUTDIR --majs-file MAJS_FILE --iggy-file IGGY_FILE --rna-file
                 RNA_FILE --obs-file OBS_FILE
                 [--high_confidence_coefficient HIGH_CONFIDENCE_COEFFICIENT]
                 [--low_confidence_coefficient LOW_CONFIDENCE_COEFFICIENT]
                 [--epsilon EPSILON] [--export]

optional arguments:
  -h, --help            show this help message and exit
  --outdir OUTDIR, -o OUTDIR
                        out directory name
  --majs-file MAJS_FILE, -mf MAJS_FILE
                        majS data file (csv)
  --iggy-file IGGY_FILE, -if IGGY_FILE
                        iggy data file (csv)
  --rna-file RNA_FILE, -rf RNA_FILE
                        rNA data file (csv)
  --obs-file OBS_FILE, -of OBS_FILE
                        observations data file (csv)
  --high_confidence_coefficient HIGH_CONFIDENCE_COEFFICIENT, -min HIGH_CONFIDENCE_COEFFICIENT
                        high confidence coefficient (float)
  --low_confidence_coefficient LOW_CONFIDENCE_COEFFICIENT, -max LOW_CONFIDENCE_COEFFICIENT
                        low confidence coefficient (float)
  --epsilon EPSILON, -eps EPSILON
                        epsilon value (float)
  --export, -e          export all predictions in a JSON file
```



## Reproduce our results

### Benchmark 1

```
python MIComp.py --majs-file ../out/Df/Df_Benchmark1.csv --iggy-file ../Iggy_RNA_benchmark1.txt --obs-file ../RNA_benchmark1.obs --rna-file dataRNA.csv --outdir results --export  
```

### Benchmark 2

```
python MIComp.py --majs-file ../out/Df/Df_Benchmark2.csv --iggy-file ../Iggy_RNA_benchmark2.txt --obs-file ../RNA_benchmark2.obs --rna-file dataRNA.csv  --outdir results --export
```

Produced results are stored in `results/` folder.

## Supplementary materials experimentations

### Epsilon parameter

Generate data:
```
bash experimentations/run_programs/epsilon_experimentation.sh 1 # For Benchmark1
bash experimentations/run_programs/epsilon_experimentation.sh 2 # For Benchmark2
```

Analyze produced data:
```
cd experimentations/epsilon
python epsilon_analyze.py -d ./data -o ./results -b 1 # For Benchmark1
python epsilon_analyze.py -d ./data -o ./results -b 2 # For Benchmark2
```

Analyze results are stored in `experimentations/epsilon/` folder.

### Low and high confidence parameters

Generate data:
```
bash experimentations/run_programs/standard_deviation_experimentation.sh 1 # For Benchmark1
bash experimentations/run_programs/standard_deviation_experimentation.sh 2 # For Benchmark2
```

Analyze produced data:
```
cd experimentations/standard_deviation
python std_analyze.py -d ./data -o results -b 1 # For Benchmark1
python std_analyze.py -d ./data -o results -b 2 # For Benchmark2
```

Analyze results are stored in `experimentations/standard_deviation/` folder.