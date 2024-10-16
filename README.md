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
Here ${sites} is a bed file specifying the genome / transcriptome coordinates where predictions should be performed. It should match the reference that was used for the alignment in step 2. To exhaustively generate all the possible sites from a reference, we provide the script [WIP].

The argument "num_jobs" is the number of parallel processes to run on the GPU. If you get an out-of-memory error, try to reduce the job number.

Unless otherwise specified, the output bam file is called "mAFiA.reads.bam" in ${outdir}.

## 4. mAFiA site-level prediction
Finally, we can aggregate the single-read probabilities to predict site-level stoichiometry.
```
pileup \
--bam_file ${outdir}/mAFiA.reads.bam \
--sites ${sites} \
--min_coverage 20 \
--out_dir ${output} \
--num_jobs 36
```
Here the input bam file should be the tagged bam generated in step 3. ""min_coverage" is the coverage cut-off for sites that should be considered. "num_jobs" are the number of cpu threads available on your computer. The output should be a bed file similar to ${sites}, but with the additional columns "coverage", "modRatio", "confidence".
