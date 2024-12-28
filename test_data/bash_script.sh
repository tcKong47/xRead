## 1. The commands and parameters used for assembly pipelines on simulated datasets
## (1) For simulated ONT datasets with an average accuracy of 87%:
## xRead
$xRead -k 15 -l 11 -w 5 -x 5 -X 10 -a 2000 -n 3 -m 500 -b 100 -t 30 -M 24 {sim_dataset} -f {out_file}
## NextDenovo
## Note: The {genome_size} option was assigned as 4.6m, 120m, 140m and 3.1g for E. coli, A. thaliana, D. melanogaster and H. sapiens, respectively.
$NextDenovo NextDenovo.cfg
## Wtdbg2
## For nanopore data and genome size < 1G:
$wtdbg2 -x preset2 -g {genome_size} -i {sim_dataset} -t 30 -fo {out_file1}
$wtpoa_cns -t 30 -i{out_file1}.ctg.lay.gz -fo {out_file2}
## For nanopore data and genome size >= 1G:
$wtdbg2 -x preset3 -g {genome_size} -i {sim_dataset} -t 30 -fo {out_file1}
$wtpoa_cns -t 30 -i{out_file1}.ctg.lay.gz -fo {out_file2}
## Flye
$flye --nano-raw {sim_dataset} --out-dir {out_dir} --threads 30
## Shasta
$shasta --input {sim_dataset} --config Nanopore-OldGuppy-Seq2020 --assemblyDirectory {out_dir} --threads 30

## (2) For simulated ONT datasets with an average accuracy of 94%:
## xRead
$xRead -k 15 -l 11 -w 13 -x 5 -X 10 -a 1500 -n 3 -m 800 -b 300 -t 30 -M 24 {sim_dataset} -f {out_file}
## NextDenovo
## The 'read_type' in config file needs to be changed to: read_type = ont
$NextDenovo NextDenovo.cfg
## Wtdbg2
## For nanopore data and genome size < 1G:
$wtdbg2 -x preset2 -g {genome_size} -i {sim_dataset} -t 30 -fo {out_file1}
$wtpoa_cns -t 30 -i{out_file1}.ctg.lay.gz -fo {out_file2}
## For nanopore data and genome size >= 1G:
$wtdbg2 -x preset3 -g {genome_size} -i {sim_dataset} -t 30 -fo {out_file1}
$wtpoa_cns -t 30 -i{out_file1}.ctg.lay.gz -fo {out_file2}
## Flye
$flye --nano-hq {sim_dataset} --out-dir {out_dir} --threads 30
## Shasta
$shasta --input {sim_dataset} --config Nanopore-May2022 --assemblyDirectory {out_dir} --threads 30

## (3) For simulated PacBio HiFi datasets with average accuracy > 99.5%:
## xRead 
$xRead -k 15 -l 11 -w 17 -x 5 -X 10 -a 1000 -n 3 -m 1000 -b 600 -t 30 -M 24 {sim_dataset} -f {out_file}
## NextDenovo
## The 'read_type' in config file needs to be changed to: read_type = hifi
$NextDenovo NextDenovo.cfg
## Wtdbg2
$wtdbg2 -x preset4 -g {genome_size} -i {sim_dataset} -t 30 -fo {out_file1}
$wtpoa_cns -t 30 -i{out_file1}.ctg.lay.gz -fo {out_file2}
## Flye
$flye --pacbio-hifi {sim_dataset} --out-dir {out_dir} --threads 30
## Shasta
$shasta --input {sim_dataset} --config HiFi-Oct2021 --assemblyDirectory {out_dir} --threads 30

## 2. The commands and parameters used for assembly pipelines on real datasets
## (1) For ONT datasets in fast base-calling mode:
## xRead
$xRead -k 15 -l 11 -w 5 -x 5 -X 10 -a 2000 -n 3 -m 500 -b 100 -t 30 -M 24 {real_dataset} -f {out_file}
## NextDenovo
## The 'read_type' in config file needs to be changed to: read_type = ont
$NextDenovo NextDenovo.cfg
## Wtdbg2
## For nanopore data and genome size < 1G:
$wtdbg2 -x preset2 -g {genome_size} -i {real_dataset} -t 30 -fo {out_file1}
$wtpoa_cns -t 30 -i{out_file1}.ctg.lay.gz -fo {out_file2}
## For nanopore data and genome size >= 1G:
$wtdbg2 -x preset3 -g {genome_size} -i {real_dataset} -t 30 -fo {out_file1}
$wtpoa_cns -t 30 -i{out_file1}.ctg.lay.gz -fo {out_file2}
## Flye
$flye --nano-raw {real_dataset} --out-dir {out_dir} --threads 30
## Shasta
$shasta --input {real_dataset} --config Nanopore-OldGuppy-Seq2020 --assemblyDirectory {out_dir} --threads 30

## (2) For human PacBio HiFi datasets:
## xRead
$xRead -k 15 -l 11 -w 17 -x 5 -X 10 -a 1000 -n 3 -m 1000 -b 600 -t 30 -M 24 {real_dataset} -f {out_file}
## NextDenovo
## The 'read_type' in config file needs to be changed to: read_type = hifi
$NextDenovo NextDenovo.cfg
## Wtdbg2
$wtdbg2 -x preset4 -g {genome_size} -i {real_dataset} -t 30 -fo {out_file1}
$wtpoa_cns -t 30 -i{out_file1}.ctg.lay.gz -fo {out_file2}
## Flye
$flye --pacbio-corr {real_dataset} --out-dir {out_dir} --threads 30
## Shasta
$shasta --input {real_dataset} --config HiFi-Oct2021 --assemblyDirectory {out_dir} --threads 30
