# Ψ-co-mAFiA

Here we provide a brief walkthrough of the software using a toy example.

## 0. Preliminary
Note: The following steps are tested with **python 3.9**.

Get code and activate virtual environment:
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
`${mafia}` will be your repo path.
Example data is provided in `${mafia}/example_data/input`, where
- `fast5_HSPA1A` contains a single fast5 file with reads covering <em>HSPA1A</em>
- `HSPA1A.bed` annotates the possible m<sup>6</sup>A and Ψ sites in this region 


## 1. Basecalling
The basecaller is adapted from the [RODAN](https://github.com/biodlab/RODAN) repository.
```
fast5_dir="${mafia}/example_data/input/fast5_HSPA1A"
backbone="${mafia}/models/RODAN_HEK293_IVT.torch"

basecall \
--fast5_dir ${fast5_dir} \
--model ${backbone} \
--batchsize 4096 \
--out_dir ${out_dir}
```
Basecalling results will be written to `${out_dir}/rodan.fasta`

## 2. Alignment
Align basecalling results to a reference genome `${ref}` (eg, GRCh38_102.fasta). Filter, sort, and index BAM file.
```
bam="${out_dir}/aligned.bam"

minimap2 --secondary=no -ax splice -uf -k14 -t 36 --cs ${ref} ${out_dir}/rodan.fasta \
| samtools view -bST ${ref} -q50 - \
| samtools sort - > ${bam}

samtools index ${bam}
```

## 3. Read-level prediction
After the standard procedures, we can now predict modification probabilities of single nucleotides on each read.
```
classifiers="${mafia}/models/psi-co-mAFiA"
sites="${mafia}/example_data/input/HSPA1A.bed"

process_reads \
--bam_file ${bam} \
--fast5_dir ${fast5_dir} \
--sites ${sites} \
--backbone_model_path ${backbone} \
--classifier_model_dir ${classifiers} \
--num_jobs 4 \
--batchsize 128 \
--out_dir ${out_dir}
```
On a Turing GPU, this step should finish within 3 minutes.

Among the input arguments,
- `--sites` points to a bed file specifying the genome / transcriptome coordinates where predictions should be performed. It should correspond to the reference that was used for the alignment in step 2. To exhaustively generate all possible sites covered by the models, we provide the script [WIP].
- `--num_jobs` is the number of parallel processes to run on the GPU. If you get a CUDA out-of-memory error, try to reduce the job number.

Unless otherwise specified, the output bam file will be called `mAFiA.reads.bam` in `${out_dir}`. It contains the MM / ML tags that conform to the modBAM standard. m<sup>6</sup>A and Ψ are represented by the ChEBI codes `21891` and `17802` respectively.

The modifications can be visualized in genome viewers such as [IGV](https://igv.org/). However, as of version 2.17.4, the viewer does not support RNA modifications and will treat all such tags as "Other". To display them with user-defined colors, one needs to relabel them as DNA modifications instead. For example,
```
samtools view -h ${out_dir}/mAFiA.reads.bam | sed -e 's/N+21891/N+a/g' | sed -e 's/N-21891/N-a/g' | samtools view -o ${out_dir}/relabelled.mAFiA.reads.bam
samtools index ${out_dir}/relabelled.mAFiA.reads.bam
```
will recast all m<sup>6</sup>A tags as 6mA. The resulting visualization would look like this:

![alt text](https://github.com/ADHDrian/psi-co-mAFiA/blob/main/example_data/igv_snapshot.png)

## 4. Site-level prediction
Finally, we can aggregate the single-read probabilities to predict site-level stoichiometry.
```
pileup \
--bam_file ${out_dir}/mAFiA.reads.bam \
--sites ${sites} \
--min_coverage 20 \
--out_dir ${out_dir} \
--num_jobs 36
```
Here the input bam file should be the modBAM generated in step 3. `--min_coverage` is the coverage cut-off for sites that should be considered for calculation. `--num_jobs` is the number of CPU threads available on your computer. The output is written to `${out_dir}/mAFiA.sites.bed`. It should be a bed file similar to `${sites}`, but with the additional columns
- `coverage`: Number of reads covering the specific site
- `modRatio`: Ratio of single nucleotides at the site with modification probability ≥0.5
- `confidence`: Ratio of single nucleotides at the site with modification probability ≥0.75 or <0.25.

For high-precision measurements, we recommend a minimum coverage of 50 and confidence of 80%.

The expected outputs of the entire workflow can be found in `${mafia}/example_data/output`.
