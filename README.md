Juan F. Santana, Ph.D. (juan-santana@uiowa.edu), University of Iowa, Iowa City, I.A.

# TSS-finder
Transcription start sites (TSS) are identified and their corresponding total counts and average transcript lengths are calculated. TSS are sorted based on their count numbers, with priority assigned to those with higher counts. Any TSS having lower counts and located within a user defined distance (base pairs) to major TSS are filtered out.

## Requirements

### System

- Linux OS

### Data files

- BED file output from PRO-Seq or PRO-Cap experiments.

### Usage

In the working directory, save `sampleskey.csv` and `.fastq.gz` files. 
The bed file name should correspond to a name present in `sampleskey.csv`.
Run:

```
./tss-finder <bed file> <TSS region (int)> <minimum counts per TSS (int)> <yes/no for bed output for conversion to bedgraph>
```

# Parameter description #
```
bed file: <str> Path to the BED file containing reads.

TSS region: <int> Distance, in base pairs, upstream and downstream of the major TSS used to define the region for TSS annotation, ensuring that only the major TSS at the center is included, excluding any other nearby TSSs.

minimum counts per TSS: <int> The minimum number of counts or reads required for a transcription start site (TSS) to be considered during annotation and transcript length calculation.

yes/no for bed output for conversion to bedgraph: <string> Specify 'yes' or 'no' to indicate whether the program should generate a BED output file that can later be used for conversion to a bedgraph format.

```

### Output
- BED file with TSSs:

| chromosome | TSS-start | TSS-end | TSS-counts | avg. transcript lengths | Strand |
|:----------:|:---------:|:--------|------------|-------------------------|-------:|
|    chr2    |  535      |  536    |   100      |       50.23             |  +     |


- BED file with every read pertaining to indentified TSSs (optional)
| chromosome | Read-start | Read-end | ID | column | Strand |
|:----------:|:----------:|:---------|----|--------|-------:|
| chr2       |  535       |   600    |xxx |  .     |  +     |

### Illustration

![tss-finder](https://github.com/JuanFSantana/TSS-finder/assets/38702786/719dc006-2e05-4e73-96e5-9ad899420883)