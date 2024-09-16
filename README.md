# HyLight
HyLight is a strain aware de novo assembly method based on the overlap-layout-consensus (OLC) paradigm that leverages the strengths of NGS and 3rd generation sequencing to rapidly and accurately assemble highly complex metagenomic sequencing data.


<p align="center">
<img src="https://github.com/kangxiongbin/HyLight/assets/23208764/587682ca-722f-43e1-89c6-5582b8cd1f09" alt="HiStrain_workflow" width="500px"/>
</p>

The workflow of HyLight.Broadly speaking, there are three main steps. First, an overlap graph is built using long reads, then the graph is optimized into a strain-aware graph. This graph is used to assemble long read contigs at strain-resolve. Next, short reads are aligned to the long read contigs and any short reads that align to assembled regions are removed. The remaining short reads undergo strain-aware assembly to produce short read contigs. Finally, the long read contigs and short read contigs are together used to construct a contig graph for further scaffolding and extension of the contigs into final master contigs.

## Installation and dependencies
Please note that HyLight is built for linux-based systems and python3 only. HyLight relies on the following dependencies:

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

## Inputs
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

If the FASTQ files you have are separated into R1 and R2, we recommend using [fastp](https://github.com/OpenGene/fastp) to merge the two files into an interleaved FASTQ format.

```
fastp -i input_R1.fastq -I input_R2.fastq --stdout > sample_interleaved.fastq
```
## Running HyLight

We have provided a simple test dataset of Illumina MiSeq and ONT reads in the example directory to check if the software has been installed and is running correctly. Please note that the out_folder must be specified with the full path.
```
python ../script/HyLight.py -l long_reads.fq -s short_reads.fq --nsplit 100 -t 30  -o <full path of the output folder>

```

## Outputs

- `final_contigs.fa`: Final assembly result of both long and short reads. This file is generated only if further assembly beyond `long_con_polished.fa` is required.

- `long_con_polished.fa`: Assembly result of long reads after polishing. 

- `short_stageb.fa`: Assembled contigs from short reads extended using a global overlap graph.

- `all_contigs.fa`: Contains all contigs from both short and long reads, used to construct a global overlap graph for further extension to generate final_contigs.fa.


## Docker

```
docker build -t hylight .
# 1. run directly in your path with data
docker run -v $PWD:/$PWD -w $PWD hylight python /tools/HyLight/script/HyLight.py -l long_reads.fq -s short_reads.fq --nsplit 100 -t 30 -o out_folder
# 2. start an interactive docker container session and run in your path with data
docker run -it --rm -v $PWD:/wd -w /wd -v /var/run/docker.sock:/var/run/docker.sock hylight /bin/bash
conda activate hylight
python /tools/HyLight/script/HyLight.py -l long_reads.fq -s short_reads.fq --nsplit 100 -t 30 -o /wd/out_folder

```

### Parameters:

- `-s`, `--short_reads`  
  Input short reads in **interleaved FASTQ** format.

- `-l`, `--long_reads`  
  Input long reads in **FASTQ** format.

- `-o`, `--outdir`  
  Path to output directory for the result files. A full path is required, not just a filename.

- `-t`, `--threads`  
  The number of threads to use.  
  **Default:** `20`

- `--corrected`  
  When specified, it indicates that both short reads and long reads have already been corrected, so no further correction is needed before proceeding with assembly.  
  **Default:** `False`

- `--nsplit`  
  The number of split input FASTA/FASTQ files. If the size of the long reads data exceeds 5 GB or 10 GB, we recommend splitting the data into 1000 or more files to improve processing speed.  
  **Default:** `60`

- `--min_identity`  
  The minimum identity for filtering overlaps. Higher values reduce error in the retained overlaps, but may also filter out useful overlaps.  
  **Default:** `0.95`

- `--min_ovlp_len`   
  The minimum overlap length between long reads. Larger values indicate stricter filtering.  
  **Default:** `3000`

- `--size`   
  The maximum size of the cluster for short reads. This is determined by your server's performance. Since we use an overlap graph, this process is computationally intensive. Assembling too many short reads at once may slow down performance.  
  **Default:** `15000`

- `--max_tip_len`  
  The maximum length to be removed as tips.  
  **Default:** `10000`

- `--insert_size`  
  The length of the insert size for short reads.  
  **Default:** `450`

- `--average_read_len`   
  The average length of short reads.  
  **Default:** `250`

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
