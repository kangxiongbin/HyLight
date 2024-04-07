# HyLight
HyLight is a strain aware de novo assembly method based on the overlap-layout-consensus (OLC) paradigm that leverages the strengths of NGS and 3rd generation sequencing to rapidly and accurately assemble highly complex metagenomic sequencing data.


<p align="center">
<img src="https://github.com/kangxiongbin/HyLight/assets/23208764/587682ca-722f-43e1-89c6-5582b8cd1f09" alt="HiStrain_workflow" width="500px"/>
</p>

The workflow of HyLight.Broadly speaking, there are three main steps. First, an overlap graph is built using long reads, then the graph is optimized into a strain-aware graph. This graph is used to assemble long read contigs at strain-resolve. Next, short reads are aligned to the long read contigs and any short reads that align to assembled regions are removed. The remaining short reads undergo strain-aware assembly to produce short read contigs. Finally, the long read contigs and short read contigs are together used to construct a contig graph for further scaffolding and extension of the contigs into final master contigs.

## Installation and dependencies
Please note that HyLight is built for linux-based systems and python3 only. HyLight relies on the following dependencies:
HyLight relies on the following dependencies:
- [bfc](https://github.com/lh3/bfc)
- [racon](https://github.com/isovic/racon)
- [fmlrc2](https://github.com/HudsonAlpha/fmlrc2)
- [miniasm](https://github.com/lh3/miniasm)
- [minimap2](https://github.com/lh3/minimap2)
- [HaploConduct](https://github.com/HaploConduct/HaploConduct)
- g++ >=5.5.0 and with boost libraries

To install HyLight, firstly, it is recommended to intall the dependencies through [Conda](https://docs.conda.io/en/latest/):
```
conda create -n HyLight
conda activate HyLight
conda install -c bioconda python=3.6 scipy pandas minimap2 bfc fmlrc2 ropebwt2 racon
```
Subsequently, pull down the code to the directory where you want to install, and compile the code:
```
git clone https://github.com/kangxiongbin/HyLight.git
cd HyLight
sh install.sh
```
## Examples

Illumina miseq and ONT reads. The out_folder must give the full path.
```
python ../script/HyLight.py -l long_reads.fq -s short_reads.fq --nsplit 100 -t 30  -o out_folder

```
The input file must be in interleaved FASTQ format. Since the final clustering step retrieves and groups reads based on their names, read names should not contain spaces. The read file should be formatted like this:
```
@S0R0/1
TATAAGTAAGGCGTTGCGAGCGGGTCGTAAAATATTTTTGATCCGT
+
EEEEEGEDJHJ3JHKJMMMLLLKNGOOLLNLOOOMJONLOOIOLMO
@S0R0/2
TTGATTATCATGCCGGAAGTGCTGCTCTTGTTCTCTGAAAGAGAAT
+
EEEGEHHHJHFJJJJBML2MMLNLLONNLNLOLJONOLNONNNMNF
```


The output is long contigs (fasta format). Long_con_polished.fa is the assembly result of long reads. final_contigs.fa is the assembly result of both long reads and short reads. Sometimes, because the assembly result of long reads has already covered most of the genomic regions, there is no final_contigs.fa. The final assembly result is long_con_polished.fa.



## Possible issues during installation (optional)

If `g++` version of the system is not satisfied, one could try this to install:
```
conda install -c conda-forge gxx_linux-64=7.3.0
# replace the /path/to/ with your own path
ln -s /path/to/miniconda3/envs/HyLight/bin/x86_64-conda-cos6-linux-gnu-g++ /path/to/miniconda3/envs/HyLight/bin/g++
ln -s /path/to/miniconda3/envs/HyLight/bin/x86_64-conda-cos6-linux-gnu-gcc /path/to/miniconda3/envs/HyLight/bin/gcc
```
If `boost` library is not installed, you could try this to install:
```
conda install -c conda-forge boost
# set envionment variables
export LD_LIBRARY_PATH=/path/to/miniconda3/envs/HyLight/lib/:$LD_LIBRARY_PATH
export CPATH=/path/to/miniconda3/envs/HyLight/include/:$CPATH
```

If compile error occurs something like `/path/to/miniconda3/envs/HyLight/x86_64-conda_cos6-linux-gnu/bin/ld: cannot find -lboost_timer `
or `cannot find -lgomp`, 
 which means it fails to link `boost` or `libgomp` library, one could try this to solve:
```
ln -s /path/to/miniconda3/envs/HyLight/lib/libboost_* /path/to/miniconda3/envs/HyLight/x86_64-conda_cos6-linux-gnu/lib/.
ln -s /path/to/miniconda3/envs/HyLight/lib/libgomp* /path/to/miniconda3/envs/HyLight/x86_64-conda_cos6-linux-gnu/lib/.
# then re-complile and install
sh install.sh
```
