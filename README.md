# Ψ-co-mAFiA

Here we provide a brief walkthrough of the software using a toy example.

Note: The following steps are tested with **python 3.9**.

## 0. Preliminary
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
- `--sites` points to a bed file specifying the genome / transcriptome coordinates where predictions should be performed. It should correspond to the reference that was used for the alignment in step 2. To exhaustively generate all possible sites covered by the models, we provide the command `generate_sites` (see [step 5](https://github.com/ADHDrian/psi-co-mAFiA/tree/main#5-generate-annotation-from-reference) below).
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

## 5. Generate annotation from reference
To generate your own bed files as input to steps 3 and 4, do
```
generate_sites \
--ref_file ${ref} \
--out_dir ${annot_dir}
```
This will produce bed files covering all contigs in the reference, for both modifications m<sup>6</sup>A and Ψ.

To limit the output to a subset of contigs, give the additional argument
`--chroms 1,2,3` for example.

To limit to only one mod, give
`--mods m6A` or `--mods psi`.

## 6. Benchmarking results
We have generated stoichiometric predictions on the HEK293 cell-line and compared the numbers with 2 chemical assays: GLORI for m<sup>6</sup>A [Liu et al, Nat Biotechnol 41, 355–366 (2023)] and PRAISE for Ψ [Zhang et al, Nat Chem Biol 19, 1185–1195 (2023)]. Each dot in the scatter plots below represents a single modified site on the human transcriptome that is predicted by both Ψ-co-mAFiA and the orthogonal method. For Ψ-co-mAFiA, there is a minimum site coverage requirement of 10 reads and minimum prediction confidence of 80%.

<img src="https://github.com/dieterich-lab/psi-co-mAFiA/blob/main/example_data/m6A_HEK293_WT_vs_GLORI_combined_conf80_cov10.png" width="300"> <img src="https://github.com/dieterich-lab/psi-co-mAFiA/blob/main/example_data/psi_HEK293_WT_vs_PRAISE_combined_conf80_cov10.png" width="300">

- For the comparison with GLORI, there are 15928 overlapping m<sup>6</sup>A sites with correlation 0.915.
- For PRAISE, there are 180 overlapping Ψ sites with correlation 0.938.

Furthermore, we have compared our results on METTL3 knock-out (KO), TRUB1 knock-down (KD), and TRUB1 over-expressed (OE) variants of HEK293 against its wildtype (WT). METTL3 is a m<sup>6</sup>A writer, while TRUB1 is a Ψ writer for the GUUCN motif.

<img src="https://github.com/dieterich-lab/psi-co-mAFiA/blob/main/example_data/m6A_METTL3_KO_vs_HEK293_WT_combined_conf80_cov10.png" width="200"> <img src="https://github.com/dieterich-lab/psi-co-mAFiA/blob/main/example_data/psi_TRUB1_KD_vs_HEK293_WT_combined_conf80_cov10_restrict_motifs_GUUCN.png" width="200"> <img src="https://github.com/dieterich-lab/psi-co-mAFiA/blob/main/example_data/psi_TRUB1_OE_vs_HEK293_WT_combined_conf80_cov10_restrict_motifs_GUUCN.png" width="200">

- In the KO / KD scenarios, the site stoichiometries are noticeably suppressed.
- For the OE case, the values are elevated.

Data sources:
- HEK293 METTL3-KO: https://www.ebi.ac.uk/ena/browser/view/PRJEB40872
- HEK293 TRUB1-KD: https://www.ebi.ac.uk/ena/browser/view/PRJEB72637
- HEK293 TRUB1-OE: RNA sample provided by the Schwartz lab (Weizmann Institute), sequenced by I. S. Naarmann-de Vries at the Dieterich lab (Heidelberg).
