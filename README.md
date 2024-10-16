# psi-co-mAFiA

Here we provide a brief walkthrough to run mAFiA, using the example of chromosome X.

## 0. Preliminary
- The following steps are tested with **python 3.9**.
  Get code and activate virtual environment, e.g.:
```
git clone git@github.com:ADHDrian/psi-co-mAFiA.git
cd psi-co-mAFiA
python3 -m venv mafia-venv
source mafia-venv/bin/activate
```
If you pip version is <21.3, then upgrade it to a newer version:
```
python3 -m pip install --upgrade pip
```
Install package
```
pip install -e .
```
${mafia} will be your repo path.
- Download data from [here](https://zenodo.org/record/8321727)
    - The folder "data" contains a subset of input data on chrX:
        - fast5_chrX: dRNA-Seq raw data from HEK293 WT mRNA
- Assume that data is unzipped to ${data} respectively. Your output directory is ${output}


## 1. Basecalling
The basecalling script is adapted from the [RODAN](https://github.com/biodlab/RODAN) repository.
```
fast5dir="${data}/fast5_chrX"
backbone="${mafia}/models/RODAN_HEK293_IVT.torch"

basecall \
--fast5dir ${fast5dir} \
--model ${backbone} \
--batchsize 4096 \
--outdir ${output}
```
On a reasonably modern GPU machine, this step should take less than 30 mins. Unless otherwise specified, basecalling results will be written to "${output}/rodan.fasta"

## 2. Alignment
Align basecalling results to a reference genome (eg, GRCh38_102.fasta). Filter, sort, and index BAM file.
```
minimap2 --secondary=no -ax splice -uf -k14 -t 36 --cs ${ref} ${basecall} \
| samtools view -bST ${ref} -q50 - \
| samtools sort - > ${bam}

samtools index ${bam}
```

## 3. mAFiA read-level prediction
After the standard procedures, we can now predict modification probabilities of single nucleotides on each read.
```
classifiers="${mafia}/models/psi-co-mAFiA"

process_reads \
--bam_file ${bam} \
--fast5_dir ${fast5_dir} \
--sites ${sites} \
--backbone_model_path ${backbone} \
--classifier_model_dir ${classifiers} \
--num_jobs 4 \
--batchsize 128 \
--out_dir ${output}
```
Here ${sites} is a bed file specifying the genome / transcriptome coordinates where predictions should be performed. To exhaustively generate all the possible sites on a reference, we provide the script [WIP].
