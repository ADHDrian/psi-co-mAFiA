# psi-co-mAFiA

Here we provide a brief walkthrough to run mAFiA, using the example of chromosome X.

## 0. Preliminary
- The following steps are tested with **python 3.9**.
  Get code and activate virtual environment, e.g.:
```
git clone https://github.com/dieterich-lab/mAFiA.git
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
```
backbone="${mafia}/models/RODAN_HEK293_IVT.torch"
classifiers="${mafia}/models/psi-co-mAFiA"
fast5dir="${data}/fast5_chrX"
```

## 1. Basecalling
The basecalling script is adapted from the [RODAN](https://github.com/biodlab/RODAN) repository. Assume that ${mafia} is your code directory.
```
python3 ${mafia}/RODAN/basecall.py \
--fast5dir ${fast5dir} \
--model ${backbone} \
--batchsize 4096 \
--outdir ${output}
```
On a reasonably modern GPU machine, this step should take less than 30 mins.
