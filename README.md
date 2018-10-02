# Megaviridae Amplicon Processing System (MAPS)

Megaviridae Amplicon Processing System (MAPS) is a set of scripts,
which can automatically process paired-end MEGAPRIMER-amplicon sequencing data obtained from Illumina MiSeq systems.
This system have the ability to process reads from other platforms, however, it can only process pair-end reads currently.

## Prerequisites

Listed softwares are required by MAPS.
Please install these softwares and add their directories to the ```PATH``` environment variable or specify them in ```Pipeline.sh```.
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [pplacer](http://matsen.fhcrc.org/pplacer/) 
* [MAFFT](http://mafft.cbrc.jp/alignment/software/)
* [FLASH](https://ccb.jhu.edu/software/FLASH/)
* [Python3](https://www.python.org/) with [Biopython](http://biopython.org/wiki/Biopython) package
* [cd-hit](http://weizhongli-lab.org/cd-hit/)
* [blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [QIIME](http://qiime.org/)

Path of Trimmomatic must be specify in Pipeline, as it is hardcoded in ```Pipeline.sh```. 
Path of ```MAIN_DIR```, which should contain references and scripts, will automatically refer to the parten folder of ```Pipeline```
should be spify in Pipeline.sh as well

PBS system are used in this script.
otherwise, you can use ```parallel``` or other commond to 
If you do not use PBS system, please modify the ```qsubarray``` commond in the tail of each function. 

## Usage
Run this system by switching to ```Scripts``` directory and using ```sh Pipeline.sh``` or ```./Pipeline.sh```.
* ```-maindir MAIN_DIR``` : Path of main directory, which should at least contain ```Scripts``` directory and ```References``` directory
[default: the parten folder of this script (```Pipeline.sh```)]
* ```-i SOURCES_DIR, -in SOURCES_DIR``` :  Path of input files, whose names should be ```[barcode]_R1.fastq``` and ```[barcode]_R2.fastq```
[default: ```maindir/Sources/```)]
* ```-o | -out``` : Path to store output, output filepath in the case of will be ```[barcode]_trimmed.fna``` and ```[barcode]_trimmed.fna``` in ```OUTPUTS_DIR/7_ALIGNMENT/```
[default: ```maindir/Outputs/```)]
* ```-qusb | -pbs``` : use PBS system, otherwise use GNU ```parallel```
* ```-t | -threads``` : Specify number of threads to be used, only work for GNU ```parallel```
* ```-module``` : use GNU module system, otherwise find software in ```PATH```
* ```-bystep``` : run partical pipeline ```-bystep 123```
1.A5G40 2.CUTADAPT_G40 3.MERGE 4.DEDUPLICATION 5.FAA 6.BLASTP 7.ALIGNMENT

## List of reference and script files
### Scripts
| Filename | Description |
| ---- | :--- |
|```Decode_primer.py```|generate primer sequence option for cutadapt|
|```Trim_merge_O.py```|trim primer sequence from merged reads using "outie"|
|```Merge_Rescue.sh```|trim reads and merge again|
|```Seq_convert.py```|convert sequences from fastq to fasta|
|```Translate_rename.py```|translate DNA sequences into amino acid sequences|
|```Trim_common_region.py```|trim sequences alignment into common region|
|```Pplacer_decode.py```|decode pplacer result file|
|```Pplacer_fna_id.py```|convert names of DNA sequences into names of amino acid sequences|

### References
| Filename | Description |
| ---- | :--- |
|```primer.fna```  |primer sequences|
|```mixture.txt``` |mixture used in cocktail method (opintional)|
|```PolB_homology_search.faa``` |reference PolB sequences for homology search|
|```PolB_design_prot.aln``` and ```PolB_design_primer_nucl.aln``` |reference alignment for trimming into common region|
|```Pplacer.aln```, ```Pplacer.res``` and ```Pplacer.info``` |reference tree and alignment for pplacer pipeline|

## Authors

* **Yanze Li** - *Initial work* - yanzeli@kuicr.kyoto-u.ac.jp

## References

Li, Y.; Hingamp, P.; Watai, H.; Endo, H.; Yoshida, T.; Ogata, H.	Degenerate PCR Primers to Reveal the Diversity of Giant Viruses in Coastal Waters. Viruses 2018, 10, 496.
http://www.mdpi.com/1999-4915/10/9/496

## Acknowledgments

* Special thanks to Florian who helped in improving MAPS
