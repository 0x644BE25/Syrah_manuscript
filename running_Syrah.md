## 1. Installing Syrah

Follow the [Installation Instructions](https://github.com/0x644BE25/Syrah/blob/main/readme.md#1-installation) or [Tutorial](https://github.com/0x644BE25/Syrah/blob/main/Syrah_tutorial.md) from the main Syrah GitHub. If you follow the tutorial, you will already have one of the two necessary genome references built and one of three datasets processed. **Take note of the directory where you have installed Syrah, as we will need this information later.**

## 2. Downloading data

 - [Curio test data](https://curioseekerbioinformatics.s3.us-west-1.amazonaws.com/TestDatasets/example_input_mouse_spleen_1M.tar.gz) includes the FASTQ files and bead coordinates file (if you followed the [Tutorial](https://github.com/0x644BE25/Syrah/blob/main/Syrah_tutorial.md) you already have this)
 - [Chick data FASTQs](https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR24325152) and [bead coordinates file](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7234197&format=file&file=GSM7234197_beadbarcodes.txt.gz)
 - [Planarian data FASTQs](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289299) (prior to publication this will require a Reviewer Token) and [bead coordinates file](https://github.com/0x644BE25/Syrah_manuscript/blob/main/data/planarian_bead_coordinates.txt)

## 3. Building references

**NOTE:** if you get the error **`command not found: STAR`** <details> This means that your computer doesn't know where to find the STAR aligner you installed. After installing Syrah, the paths to different components are saved in the `paths.txt` file in the directory where you installed Syrah. To tell your computer to look in those places for software, open a terminal at the folder with `paths.txt` and run the command 
```
while IFS="" read -r p || [ -n "$p" ]; do; PATH=$p:$PATH; done < ../paths.txt
```
You will need to repeat this if you close the terminal and open a new one. You do not need to do this before running `Syrah.sh`, as it will do so itself.</details>

#### Mouse genome: GRCm39

(If you followed the tutorial, you can skip this). Download and unzip the [FASTA](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz) and [GTF](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz) files. Use [STAR to generate a genome index](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf#section.2) using the FASTA and GTF. 

#### Chick genome: 

Download and unzip the [FASTA](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz) and [GTF](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf.gz) files. Use [STAR to generate a genome index](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf#section.2) using the FASTA and GTF. 

#### Planarian transcriptome: smed_20140614

Download and unzip the [FASTA](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE72389&format=file&file=GSE72389_smed_20140614.fa.gz) and [follow these instructions](https://github.com/ejrsimr/transcriptome2star) to generate a STAR-compatible FASTA and GTF file for a transcriptome. Use [STAR to generate a genome index](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf#section.2) using the STAR-compatible FASTA and GTF.

## Running Syrah

For each of the datasets, make a copy of the [manifest.txt file](https://github.com/0x644BE25/Syrah/blob/main/manifest.txt) and fill out the [manifest parameters](https://github.com/0x644BE25/Syrah/blob/main/readme.md#manifest-parameters). If you're planning to generate the manuscript's plots afterward, it will simplify things if you male at `data` directory to write all of the results to and choose `batchName`s to match the names of [the count files in this repo](https://github.com/0x644BE25/Syrah_manuscript/tree/main/data).

Once your manifest files are ready, simply pass them to Syrah with 
```
bash Syrah.sh /path/to/my_manifest.txt
```
or, if they're all in the same folder you can use a `*` wildcard to pass them all at once with 
```
bash Syrah.sh /path/to/*_manifest.txt
```

**NOTE 1:** If you're not in the Syrah directory, you will need to use the proper path to `Syrah.sh`, like 
```
bash /path/to/Syrah_installation_folder/Syrah.sh /path/to/my_manifest.txt
```

**NOTE 2:** If you would like to save a copy of the Syrah processing log, use the [`tee` command](https://www.geeksforgeeks.org/tee-command-linux-example/) with 
```
bash Syrah.sh /path/to/my_manifest.txt | tee -a syrah_output_log.txt
```
to both display *and* save the terminal output.

## Generating Plots

Ensure that all of your counts files are in a directory together named `data` and that their names match those of [the counts files in this repo](https://github.com/0x644BE25/Syrah_manuscript/tree/main/data). Don't worry about filtering to a minimum of 10 UMIs -- [the script `generate_plots_and_xlsx.R`](https://github.com/0x644BE25/Syrah_manuscript/blob/main/generate_plots_and_xlsx.R) will handle that for you. Simply open your copy of the script, modify the `setwd()` command as needed, and run.
