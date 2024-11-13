# Exercises

1. [Setting up the cloud computing](#setting-up-the-cloud-computing)
   1. [Setting up VS Code](#setting-up-vs-code)
   2. [Connecting to the remote machine with VS Code](#connecting-to-the-remote-machine-with-vs-code)
   3. [Cloning the course's GitHub repository](#cloning-the-courses-github-repository)
2. [Getting the raw data](#getting-the-raw-data)
3. [QC and trimming](#qc-and-trimming)
   1. [QC of the raw data](#qc-of-the-raw-data)
   2. [Read trimming](#read-trimming)
   3. [QC of the trimmed data](#qc-of-the-trimmed-data)
4. [Read-based taxonomic profiling](#read-based-taxonomic-profiling)
   1. [Kraken2](#kraken2)
   2. [sourmash](#sourmash)
5. [Metagenome assembly](#metagenome-assembly)
   1. [Assembly QC](#assembly-qc)
   2. [Abundance quantification of assembled contigs](#abundance-quantification-of-assembled-contigs)
   3. [Taxonomic annotation of assembled contigs](#taxonomic-annotation-of-assembled-contigs)
7. [Assembling long reads with Flye](#assembling-long-reads-with-flye)
   1. [Polishing the assembly with Illumina reads](#polishing-the-assembly-with-illumina-reads)
6. [Automatic binning with SemiBin2](#automatic-binning-with-semibin2)
7. [Quality control and taxonomic annotation of metagenome-assembled genomes (MAGs)](#quality-control-and-taxonomic-annotation-of-metagenome-assembled-genomes-mags)
8. [Targeted functional analysis of MAGs](#targeted-functional-analysis-of-mags)

## Setting up the cloud computing

We will use the [Amazon Cloud](https://aws.amazon.com/ec2/) (AWS EC2) services for most of the analyses.  
The IP address of the remote machine will change every day, so a new IP adress will be posted in Slack each morning.  
Your username - that you have received by e-mail/Slack - will be the same for the whole course.  
We will use `ssh` to connect to the remote machine.  
We encourage the use of [VS Code](https://code.visualstudio.com/Download), but you are welcome to use any IDE or terminal emulator that you are comfortable with.  

### Setting up VS Code

1. Download `VS Code` and set it up as shown [here](Lectures/course-outline-and-practical-info.pdf)  
2. Save the `.pem` file that you have received by e-mail somewhere in your computer  
3. **Linux/MacOS users only:**  
   1. Launch the `Terminal` app  
   2. `cd` to the directory where you saved the `.pem` file
   3. run `chmod 600 userXX.pem` (remember to change `userXX.pem` by the name of your own file)
4. Back to `VS Code`, go to `View -> Command Palette`  
5. Search for `ssh config`  
6. Select `Remote-SHH: Open SSH Configuration File...`  
7. In the next dialogue box, hit `Enter/Return`
8. Copy and paste the following text:

```
Host physalia
  HostName 54.245.21.143
  User user1
  IdentityFile ~/Desktop/user1.pem
```

9. In the 2nd, 3rd, and 4th lines:  
   1. HostName: change to the IP adress of the day
   2. User: change to your own username
   3. IdentityFile: change to the location and name of the `.pem` file that you have saved in your computer 
10. Save and close the `config` file

### Connecting to the remote machine with VS Code

11. Go to `View -> Command Palette`
12. Search for `ssh connect`
13. Select `Remote-SHH: Connect to Host...`
14. Select `physalia` (a new window will open)
15. If a dialogue box opens asking the server type, select `Linux`
16. If a dialogue box opens asking if you are sure, select `Continue`
17. If a terminal does not open by default, go to `Terminal -> New Terminal`

That's it, you should now be connected to the remote machine and ready to go!  
**Remember:** every day you should redo steps 4-7 and update `HostName` to match the IP adress of the day.  

### Cloning the course's GitHub repository

Once you have connected to the remote machine, you will be in your home folder (`/users/userXX`, also represented by `~` or `$HOME`).  
**Remember:** You can check where you are with the command `pwd`.  

To have access to the course's content, let's copy the GitHub repository to your `home` folder using `git clone`:

**Do this on the first day:** 

```bash
cd ~
git clone https://github.com/NikolayOskolkov/Physalia_EnvMetagenomics_2024
```

**Do this every once in a while, at least each day before starting the activities:**  

```bash
cd ~/Physalia_EnvMetagenomics_2024
git pull
```

**Note:** All exercises will be executed inside the `Physalia_EnvMetagenomics_2024` folder that you cloned inside your own `home` folder.  
So remember to `cd ~/Physalia_EnvMetagenomics_2024` every time you connect to the remote machine.  

## Getting the raw data

Copy the raw sequencing data to your own `01_DATA` folder.  
Also copy the file `SAMPLES.txt`, which will be useful for running `for loop` and etc.  

```bash

cd ~/Physalia_EnvMetagenomics_2024
mkdir 01_DATA

cp ~/Share/toy_data/*.fastq.gz 01_DATA/
cp ~/Share/toy_data/SAMPLES.txt ./
```

Let us now explore the data a little bit. First of all, we can look inside the gzipped-file without unzipping with `zcat`:

```bash
zcat 01_DATA/*R1.fastq.gz | head
```

You should see 4 lines corresponding to each read: the first line contains the read ID (each starting with @), 
the second line corresponds to the sequence of the read, the third line is the delimiter and the fourth line contains ASCII quality scores for eac sequenced nucleotide.

Let us now count the number of reads in the fastq-files:

```bash
zcat 01_DATA/*R1.fastq.gz | grep -c @
```

How many reads do we have in the fastq-files?


## QC and trimming

Now that you have copied the raw data to your working directory, let's do some quality control.  
The sequencing process is subject to several types of problems that can introduce errors and artifacts in the sequences.  
Because of this, bioinformatics analyses usually start with the quality control of raw sequences.  
He we will use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc) and [MultiQC](https://multiqc.info/) to obtain quality reports, and [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) for trimming the Illumina data, respectively.  

### QC of the raw data

Go to your `Physalia_EnvMetagenomics_2024` folder, create a folder for the QC files, and activate the `conda` environment:  

```bash
cd ~/Physalia_EnvMetagenomics_2024
mkdir 02_QC_RAW
conda activate envmetagenomics
```

And now you're ready to run the QC on the raw data:

```bash
fastqc 01_DATA/*.fastq.gz -o 02_QC_RAW -t 4
multiqc 02_QC_RAW -o 02_QC_RAW --interactive
```

After the QC is finished, copy the `MultiQC` report (`02_QC_RAW/multiqc_report.html`) to your local machine and open it with your favourite browser.  
We will go through the report together before continuing with the pre-processing.  

**NOTE:** to move files to and from local and remote machines, you can use: 
- The command-line tool [scp](https://kb.iu.edu/d/agye)  
- A file transfer software such as [FileZilla](https://filezilla-project.org)  
- The `VS Code` built-in `Explorer` tool (`View -> Explorer`) 

Below we provide an example (please note that the IP address should be changed) of copying files to you local computer via `scp` command-line tool:

```bash
scp -r -i envmeta24.pem ubuntu@54.244.63.96:/home/ubuntu/Physalia_EnvMetagenomics_2024/02_QC_RAW/* .
```

### Read trimming

Our QC reports tell us that a significant percentage of the raw sequences contain some isses such as the presence of adapters.  
Before proceeding, it is necessary to clean up/trim the raw sequences.  
Before start trimming the data, let's create a folder for the processed data and activate the `conda` environment:  

```bash
cd ~/Physalia_EnvMetagenomics_2024
mkdir 03_TRIMMED
conda activate envmetagenomics
```

For the Illumina data, we will use a `for loop` to process each of the samples one after the other:  

```bash
for sample in $(cat SAMPLES.txt); do
  cutadapt 01_DATA/${sample}_R1.fastq.gz \
           01_DATA/${sample}_R2.fastq.gz \
           -o 03_TRIMMED/${sample}_R1.fastq.gz \
           -p 03_TRIMMED/${sample}_R2.fastq.gz \
           -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
           -A CTGTCTCTTATACACATCTGACGCTGCCGACGA \
           -m 50 \
           -q 20 \
           -j 4 > 03_TRIMMED/${sample}.log
done
```

While `Cutadapt` is running: looking at the [online manual](https://cutadapt.readthedocs.io/en/stable/index.html) or running `cutadapt --help`, answer:  

- What do the `-o`, `-p`, `-a`, `-A`, `m`, `-q`, and `-j` flags mean?  
- How did we choose the values for `-m` and `-q`?  
- What is the purpose of the redirection (`> 03_TRIMMED/${sample}.log`)?  


### QC of the trimmed data

Now the data has been trimmed, it would be a good idea to run `FastQC` and `MultiQC` again.  
Modify the [commands used for the raw data](#qc-of-the-raw-data) to match the trimmed data and run the two QC softwares.  

While you wait, take a look at the `Cutadapt` logs.  
When `Cutadapt` runs, it prints lots of interesting information to the screen, which we lose once we logout of the remote machine.  
Because we used redirection (`>`) to capture the standard output (`stdout`) of `Cutadapt`, this information is now stored in a file (`03_TRIMMED/${sample}.log`).  
Take a look at the log file for one of the samples using the program `less`:  

**NOTE:** You can scroll up and down using the arrow keys on your keyboard, or move one "page" at a time using the spacebar.  
**NOTE:** To quit `less`, hit the `q` key.  
**NOTE:** If you have set it up, you can also access the files using the `Explorer` tab on `VS Code` (`View -> Explorer`).  

By looking at the `Cutadapt` log, can you answer:  
- How many read pairs we had originally?  
- How many reads contained adapters?  
- How many read pairs were removed because they were too short?  
- How many base calls were quality-trimmed?  
- Overall, what is the percentage of base pairs that were kept?  

When `FastQC` and `MultiQC` have finished, copy the `MultiQC` report to your local machine and open it with a browser.  
Compare this with the report obtained earlier for the raw data.  
Do the data look better now?  


## Bonus step: host removal

Even if you work with environmental samples, it is quite likely that human DNA is also present in your sample, in this sense it is considered as contamination. 
Therefore, to be on a safe side, it is a good practice to explicitely clean your data from it. 
If you work with host-associated microbiome, i.e. human microbiome, this is a mandatory step, please see [here](https://retractionwatch.com/2024/06/26/all-authors-agree-to-retraction-of-nature-article-linking-microbial-dna-to-cancer/)
what can happen if you do not properly clean your data from human DNA. Here, we demonstrate how to practically perfrom the host removal step.


```bash
# DO NOT RUN THIS (ALREADY DONE): download and index human reference genome
# wget http://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
# bowtie2-build --large-index hg38.fa.gz hg38.fa.gz --threads 20

cd ~/Physalia_EnvMetagenomics_2024
mkdir 04_HOST_REMOVAL

for sample in $(cat SAMPLES.txt); do
	bowtie2 --large-index -x ~/Share/Databases/hg38.fa.gz --end-to-end --threads 4 --very-sensitive \
	-1 03_TRIMMED/${sample}_R1.fastq.gz -2 03_TRIMMED/${sample}_R2.fastq.gz | samtools view -bS -h -@ 4 - \
	> 04_HOST_REMOVAL/${sample}_aligned_to_hg38.bam
	
	samtools sort 04_HOST_REMOVAL/${sample}_aligned_to_hg38.bam -@ 4 \
	> 04_HOST_REMOVAL/${sample}_aligned_to_hg38.sorted.bam
	
	samtools index 04_HOST_REMOVAL/${sample}_aligned_to_hg38.sorted.bam
	samtools view -b -f 4 04_HOST_REMOVAL/${sample}_aligned_to_hg38.sorted.bam \
	> 04_HOST_REMOVAL/${sample}_unaligned_to_hg38.bam

	samtools bam2fq 04_HOST_REMOVAL/${sample}_unaligned_to_hg38.bam | gzip \
	> 04_HOST_REMOVAL/${sample}_unaligned_to_hg38.fastq.gz
done
```
Above, we constructed a fastq-file which is free from human DNA. This was done by aligning the trimmed reads to the human reference genome and extracting unaligned reads only.

## Read-based taxonomic profiling

There are many different tools and approaches for obtaining taxonomic profiles from metagenomes.  
Here we will use a popular read-based taxonomic profiler [Kraken2](https://github.com/DerrickWood/kraken2) and [sourmash](https://sourmash.readthedocs.io/en/latest/).  
What is the basic approach that each of these tools use and how they can impact the results?  
Well, let's find out!  

First let's create a folder to store the results:  

```bash
cd ~/Physalia_EnvMetagenomics_2024
mkdir 05_TAXONOMIC_PROFILE
```

### Kraken2

And now let's run `Kraken2`. Kraken2 is a taxonomic sequence classifier that assigns taxonomic labels to DNA sequences. 
Kraken2 examines the k-mers within a query sequence and uses the information within those k-mers to query a database. 
That database maps k-mers to the lowest common ancestor (LCA) of all genomes known to contain a given k-mer.  

```bash
conda activate envmetagenomics

for sample in $(cat SAMPLES.txt); do
  kraken2 --db ~/Share/Databases/minikraken2_v2_8GB_201904_UPDATE \
	  --paired 03_TRIMMED/sample1_ILM_R1.fastq.gz 03_TRIMMED/sample1_ILM_R2.fastq.gz \
	  --output 05_TAXONOMIC_PROFILE/${sample}_sequences.kraken \
	  --report 05_TAXONOMIC_PROFILE/${sample}_kraken.output \
	  --report-minimizer-data --use-names --threads 4
done
```

Now that we have got our hands into some tables describing the abundance of the different taxa in our metagenome, it is time to make sense of the data.  
One way to do this is making summaries, plots, statistical tests, etc, as you would normally do for any kind of species distribution data.  
Here you are free to use whichever tool you are most familiar with (but we all know that there is only one co`R`rect tool for this).  

The idea here is to: 
- Learn what are the main (most abundant) taxa in our samples  
- Learn about potential differences in community composition between the samples  
- Learn what fraction of the community we were actually able to identify at, let's say, the genus level  
- Compare the taxonomic profiles obtainted from Illumina and Nanopore data  

Hopefully you will be able to learn a bit about these metagenomic datasets.  
And realise that there is so much that still remains unknown...  

We recommend to use [R / Rstudio](https://posit.co/download/rstudio-desktop/) for visualization of microbial abundances in your sample. For example, one can use [Pavian](https://github.com/fbreitwieser/pavian) tool:

```R
# explore abundance in Pavian https://github.com/fbreitwieser/pavian
if (!require(remotes)) { install.packages("remotes") }
remotes::install_github("fbreitwieser/pavian")
pavian::runApp(port=5000)
```

Another popular way of visualization of your data is [Krona](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-385) which can be installed and used from https://github.com/marbl/Krona.


### sourmash

There are many different appraoches for taxonomic profiling of metagenomes, each of them with their own up- and downsides.  
Let's now try a `sourmash`. Sourmash is k-mer-based, similar to Kraken2, Bracken, and Centrifuge, but uses k-mers from across the entire dataset, rather than individual reads, to find best-match
genomes. In this way, it is able to leverage longer-range information present in a dataset, though
not across reads themselves.

```bash
conda activate envmetagenomics

for sample in $(cat SAMPLES.txt); do
  sourmash sketch dna 03_TRIMMED/${sample}_R?.fastq.gz \
                      -p k=31,scaled=1000,abund \
                      -o 05_TAXONOMIC_PROFILE/${sample}.sig.zip \
                      --merge ${sample}

  sourmash gather 05_TAXONOMIC_PROFILE/${sample}.sig.zip \
                  ~/Share/Databases/gtdb-rs207.genomic-reps.dna.k31.zip \
                  -k 31 --threshold-bp 10 \
                  -o 05_TAXONOMIC_PROFILE/${sample}.gather.csv
done

# Gather results
sourmash tax metagenome -g 05_TAXONOMIC_PROFILE/*.gather.csv \
                        -t ~/Share/Databases/gtdb-rs207.taxonomy.with-strain.csv.gz \
                        --output-dir 05_TAXONOMIC_PROFILE \
                        --output-base sourmash.phylum \
                        --output-format lineage_summary \
                        --rank phylum

sourmash tax metagenome -g 05_TAXONOMIC_PROFILE/*.gather.csv \
                        -t ~/Share/Databases/gtdb-rs207.taxonomy.with-strain.csv.gz \
                        --output-dir 05_TAXONOMIC_PROFILE \
                        --output-base sourmash.genus \
                        --output-format lineage_summary \
                        --rank genus
```

Now analyse the results from `sourmash` in `R` or other data analysis tool of your preference.  
Are there differences between the taxonomic profiles obtained by the two different tools?  

## Metagenome assembly

Now it's time to move forward to metagenome assembly. For the assembly we will use [Megahit](https://github.com/voutcn/megahit) which is is an ultra-fast and memory-efficient NGS assembler. 
It is optimized for metagenomes, but also works well on generic single genome assembly (small or mammalian size) and single-cell assembly.

Before you start the assembly, have a look at the [Megahit usage examples](https://github.com/voutcn/megahit?tab=readme-ov-file#usage), and the [Megahit publication](https://academic.oup.com/bioinformatics/article/31/10/1674/177884).  

__What options do we need?__  
We have only given the output directory in the script below; modify it as necessary and run `megahit`:  

```bash 
cd ~/Physalia_EnvMetagenomics_2024
conda activate envmetagenomics

for sample in $(cat SAMPLES.txt); do
megahit -r 04_HOST_REMOVAL/${sample}_unaligned_to_hg38.fastq.gz --out-dir 06_ASSEMBLY --min-contig-len 100 -t 4
done
```

### Assembly QC

Now we are going to explore the assembled contigs. First, we will run QC for the assembled contigs and compute N50 value.  

```bash
mkdir 07_ASSEMBLY_QC

conda activate envmetagenomics

# check assembly quality statistics with callN50 JavaScript script 
# that requires the k8 JavaScript shell (or node) to be installed
# download callN50: wget https://raw.githubusercontent.com/lh3/calN50/master/calN50.js

# if you need to install k8 please run
# wget -O- https://github.com/attractivechaos/k8/releases/download/v1.2/k8-1.2.tar.bz2 | tar -jxf -
# however, k8 is already installed for you in ~/Share/k8-1.2/ 

~/Share/k8-1.2/./k8-x86_64-Linux ~/Share/calN50.js 06_ASSEMBLY/final.contigs.fa > 07_ASSEMBLY_QC/assemstats.txt
```

N50 has a complex meaning. It is some sort of "average" (or representative) contig length but not exactly.

If we have contigs with length: 2,3,4,5,6,7,8,9,10 then the total assembled length is 2+3+4+5+6+7+8+9+10=54, then the largest contigs of length 10+9+8=27 make half of assembled length, therfore N50=8 and L50=3, 
i.e. 8 is the length of the smallest contig which in the sum with larger contigs make 50% of total assembled length, and 3 is the number of contigs of lengths greater or equal 8 which together make 50% of assembled length.


In this tutorial, you should see that we assembled 11 contigs of total length 1465 bp and N50=136 bp which is the "average" / "typical" or "median" contig length, and L50=5 contigs with length greater or equal than 136 bp make 50% of total assembled length.

### Abundance quantification of assembled contigs

Now, when we have assembled contigs, we might wonder what organisms they correspond to and how abundant these organisms are in our samples.

To quantify abundance of each assembled contig, let us now align the trimmed reads back to assembled contigs. We will do it with `Bowtie2` aligner, and we will first have to build the index for the assembled contigs.

```bash
bowtie2-build --large-index 06_ASSEMBLY/final.contigs.fa 06_ASSEMBLY/final.contigs.fa --threads 4

bowtie2 --large-index -x 06_ASSEMBLY/final.contigs.fa --end-to-end --threads 4 --very-sensitive \
-1 03_TRIMMED/${sample}_R1.fastq.gz -2 03_TRIMMED/${sample}_R2.fastq.gz | samtools view -bS -h -q 1 -@ 4 - \
> 07_ASSEMBLY_QC/aligned_to_assembled_contigs.bam

samtools view 07_ASSEMBLY_QC/aligned_to_assembled_contigs.bam | cut -f3 > 07_ASSEMBLY_QC/contig_count.txt
```

Above, we generated a bam-alignment where it is recorded to what contig each read is aligned. Then we used samtools to extract a list of contigs corresponding to each aligned read.
Now let us order the assembled contigs by their abundance, we will use R for this purpose:

```R
df<-scan("07_ASSEMBLY_QC/contig_count.txt",what="character")

head(sort(table(df),TRUE))

write.table(data.frame(sort(table(df),TRUE)),file="07_ASSEMBLY_QC/abund_contigs.txt",
col.names=FALSE,row.names=TRUE,quote=FALSE,sep="\t")
```
Finally, let us display top-abundant contigs:


```bash
head 07_ASSEMBLY_QC/abund_contigs.txt
```

Can you name a few most abundant contigs?

### Taxonomic annotation of assembled contigs

Now, we will figure out what organisms with available taxonomic annotation correspond to the assembled contigs. We will use Kraken2 for assigning taxa to assembled contigs:

```bash
kraken2 --db ~/Share/Databases/minikraken2_v2_8GB_201904_UPDATE --threads 4 \
--output 07_ASSEMBLY_QC/sequences.kraken_contigs --use-names \
--report 07_ASSEMBLY_QC/kraken.output_contigs 06_ASSEMBLY/final.contigs.fa
```

Please explore the taxonomic annotation of the assembled contigs and compare it with the read-based taxonomic profiling results.


## Assembling long reads with Flye


For the assembly of our Nanopore data we will use [Flye](https://github.com/fenderglass/Flye).
`Flye` is a long-read de novo assembler that can also handle metagenomic data.

Since the process will take several minutes, copy and paste the below and read the details whiles it is running.

```bash 
conda activate envmetagenomics
cp -p ~/Share/urban_soil_toy_long_reads/ONT_preprocessed.fq.gz 01_DATA/

flye \
   --nano-raw 01_DATA/ONT_preprocessed.fq.gz \
   --meta \
   -t 4 \
   --out-dir 08_ASSEMBLY_ONT
```

The options we are using are
- `--nano-raw`: the input file contains Nanopore reads
- `--meta`: metagenomic mode
- `-t 4`: number of threads (4 — use as many threads as you have cores)
- `--out-dir`: the output directory

For more information, have a look at the [Flye manual](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md), especially the parts about Nanopore data and metagenome assembly.

In real data, this would take several hours to run (even days), so we are running it on sampled dataset (0.5% of the original) and it will still take some time.

### Polishing the assembly with Illumina reads

There are several ways to deal with hybrid data (Illumina and Nanopore).
Here, we are demonstrating one of the possible ways: assemble Nanopore data with Flye (as we did above) and then polish the assembly with Illumina data.
We consider this a good approach, but handling hybrid data is an evolving area and different groups might have different preferences.

For this, we need to first align the Illumina reads to the Nanopore assembly and then polish the assembly with the aligned reads.
Since we have already done this before, we will just present the code here:

```bash
conda activate envmetagenomics

cp -p ~/Share/urban_soil_toy_long_reads/ILM_selected.pair.1.fq.gz 01_DATA/
cp -p ~/Share/urban_soil_toy_long_reads/ILM_selected.pair.2.fq.gz 01_DATA/

bowtie2-build \
    --large-index 08_ASSEMBLY_ONT/assembly.fasta 08_ASSEMBLY_ONT/assembly.fasta --threads 4

bowtie2 \
    --large-index \
    -x 08_ASSEMBLY_ONT/assembly.fasta \
    --end-to-end --threads 4 --very-sensitive \
    -1 01_DATA/ILM_selected.pair.1.fq.gz -2 01_DATA/ILM_selected.pair.2.fq.gz \
    > 08_ASSEMBLY_ONT/aligned_to_unpolished_assembly.sam
```
Now that we have obtained a SAM file, we can use `Polypolish` to polish the assembly:

```bash
conda activate envmetagenomics

polypolish polish \
        08_ASSEMBLY_ONT/assembly.fasta \
        08_ASSEMBLY_ONT/aligned_to_unpolished_assembly.sam \
        | gzip > 08_ASSEMBLY_ONT/polished.fasta.gz
```

Note that polypolish **only works with** SAM files, not BAM files.

Now, the new version of the assembly is in the file `08_ASSEMBLY_ONT/polished.fasta.gz`.
Going forward, you would use this file for further analyses.
The assembly QC could be done in the same way as for the Illumina data.


However, we are going to use a prepared assembly for the next steps (otherwise they would produce no results), so let's copy it:
```bash
cp -p ~/Share/urban_soil_toy_long_reads/polished_assembly_pre.fa 08_ASSEMBLY_ONT/
```

## Automatic binning with SemiBin2

SemiBin2 is one of the several automatic binning algorithms published.
If you want learn more, there is a [there is a publication](https://academic.oup.com/bioinformatics/article/39/Supplement_1/i21/7210480) (also, the [manuscript about the original SemiBin1](https://www.nature.com/articles/s41467-022-29843-y)). More practical reading can be found from the documentation: [https://semibin.readthedocs.io/en/latest/](https://semibin.readthedocs.io/en/latest/).  

SemiBin2 uses self-supervised learning and has some pre-trained models, which makes the computation faster. It requires as input the results from mapping reads back to assembly (sorted and indexed bam files, which we already have) and the assembly (which we also have).

```bash
mkdir 09_SEMIBIN
```

### Mapping reads back to the assembly

This is the same operation we have been doing before.

```bash

bowtie2-build \
    --large-index 08_ASSEMBLY_ONT/polished_assembly_pre.fa 08_ASSEMBLY_ONT/polished_assembly_pre.fa \
    --threads 4
bowtie2 \
    --large-index \
    -x 08_ASSEMBLY_ONT/polished_assembly_pre.fa \
    --end-to-end --threads 4 --very-sensitive \
    -1 01_DATA/ILM_selected.pair.1.fq.gz -2 01_DATA/ILM_selected.pair.2.fq.gz \
    | samtools view -bS -h -@ 4 - \
    | samtools sort -@ 4 - \
    > 09_SEMIBIN/aligned_to_assembly.sorted.bam
```

The one important step here is that we must sort the BAM file, which is done with `samtools sort`.

This should take a few minutes to run.

### Binning with SemiBin2

Once we have (1) the assembly and (2) the sorted BAM file, we can run SemiBin2.

For complex reasons, we need to switch to a different conda environment for this step (and several of the other tools in this section).

```bash
conda activate SemiBin

SemiBin2 single_easy_bin \
        --input-fasta 08_ASSEMBLY_ONT/polished_assembly_pre.fa \
        --input-bam 09_SEMIBIN/aligned_to_assembly.sorted.bam \
        --sequencing-type=long_read \
        -o 09_SEMIBIN/SemiBin2_out \
        --environment soil \
        --threads 4
```

To break down the command above:

- `--input-fasta`: the assembly file
- `--input-bam`: the sorted BAM file
- `--sequencing-type`: the type of sequencing data used (`long_read` in this case)
- `-o`: the output directory
- `--environment`: the environment type (in this case, `soil`)
- `--threads`: the number of threads to use


At the end of the run, you will have a set of bins in the `09_SEMIBIN/SemiBin2_out/output_bins` directory.
This is a toy example where almost all contigs get binned.
This is not the case in real data (but real data would take too long to run the mapping step).

## Quality control of metagenome-assembled genomes (MAGs)

Now we have obtained some bins that we think could represent genomes present in our samples. Next steps are QC and taxonomic annotation of our genomes.

We will use a program called [checkM2](https://github.com/chklovski/CheckM2) to get more precise estimates of the completeness and redundancy of these genomes and another one, called [GUNC](https://doi.org/10.1186/s13059-021-02393-0) to get an alternative check of contamination.

### checkM2

Normally, to use checkM2, you need to download the database first (by running `checkm2 database --download`).
In this case, we have already downloaded it for you, so you can just activate the environment and run the command.

```bash
conda activate checkm2

checkm2 predict \
      --input 09_SEMIBIN/SemiBin2_out/output_bins \
      -x fa.gz \
      --output-directory 09_SEMIBIN/CheckM2 \
      --database_path ~/Share/Databases/CheckM2_database/uniref100.KO.1.dmnd \
      --threads 4
```

Since this takes a while, let's start the process first and look at the details while it's running!

Command line break down (should be mostly self-explanatory, except `-x`):
- `--input`: input **directory**
- `-x`: the **extension** of the bins
- `--threads`: number of threads to use
- `--output-directory`: where to store the outputs
- `--database_path`: normally you would not need this

CheckM2 will estimate two important metrics:

- `completeness`: how much of the genome is present in the MAG
- `contamination`: how much of the MAG is not from the original genome

Once everything has run, we can inspect the file `09_SEMIBIN/CheckM2/quality_report.tsv` for these results.

Often, we will want to use them to filter the outputs as well.

### GUNC

_This takes too long to run, so we will include it here, but you can run it at your leisure later._

Like checkM2, GUNC relies on a prebuilt database to run; and, like above, we have predownloaded it for you.

```bash
conda activate gunc
mkdir 09_SEMIBIN/GUNC_out
gunc run \
    -d 09_SEMIBIN/SemiBin2_out/output_bins \
    --file_suffix .fa.gz \
    --db_file ~/Share/Databases/GUNC/gunc_db_progenomes2.1.dmnd \
    --out_dir 09_SEMIBIN/GUNC_out
```

The file `09_SEMIBIN/GUNC_out/GUNC.progenomes_2.1.maxCSS_level.tsv` will contain the results.
The most important column is `pass.GUNC` which is a boolean.
If `False`, this is a sign that that genome is probably chimeric.

You can see the outputs in `~/Share/expected_outputs/09_SEMIBIN/GUNC_out/GUNC.progenomes_2.1.maxCSS_level.tsv`

## Taxonomic annotation of metagenome-assembled genomes (MAGs)

For taxonomic annotation we will use Genome Taxonomy Database ([GTDB](https://gtdb.ecogenomic.org/)) and a tool called [GTDB-Tk](https://ecogenomics.github.io/GTDBTk/installing/index.html#installing-gtdbtk-reference-data).

### GTDB-tk

By now, it should not surprise you that GTBD-tk uses a database that we have predownloaded, but that you'd normally be expected to download.

```bash
conda activate gtdbtk

gunzip 09_SEMIBIN/SemiBin2_out/output_bins/*.fa.gz

gtdbtk classify_wf \
      --genome_dir 09_SEMIBIN/SemiBin2_out/output_bins \
      --out_dir 09_SEMIBIN/GTDB \
      -x fa \
      --cpus 4 \
      --skip_ani_screen
```

### Quality control recap

1. We have used checkM2 to estimate `completeness` and `contamination`
2. We have used GUNC to check for chimerism (basically another form of contamination)
3. We have used GTDB-Tk to annotate the genomes taxonomically

## Functional analysis of MAGs

There are several approaches we can take to annotate our MAGs and find, for example, which kind of metabolic pathways they encode.
At this we point, we stop doing metagenomics and start doing genomics, so in one way we have reached the end of our workshop.

```bash
mkdir -p 10_MAG_ANNOTATIONS/MAGs
cp -pir 09_SEMIBIN/SemiBin2_out/output_bins/SemiBin_* 10_MAG_ANNOTATIONS/MAGs
```


### Broad-scale annotation

To get a broad-scale annotation of our MAGs, we can predict genes and annotate them with a database such as [eggNOG](http://eggnog5.embl.de/) or a specific database like [CARD](https://card.mcmaster.ca/) (for antibiotic resistance genes).

Thus, we will illustrate how to do this with a single genome, but you can run this on all of your genomes.

### Predicting genes with Prodigal

You can predict genes with Prodigal by running the following command:

```bash
conda activate envmetagenomics
mkdir -p 10_MAG_ANNOTATIONS/Prodigal
prodigal -i 10_MAG_ANNOTATIONS/MAGs/SemiBin_0.fa \
    -a 10_MAG_ANNOTATIONS/Prodigal/SemiBin_0.prots.faa \
    -d 10_MAG_ANNOTATIONS/Prodigal/SemiBin_0.prots.fna \
    -o 10_MAG_ANNOTATIONS/Prodigal/SemiBin_0.prots.out
```

This should be pretty fast and produce three files:

1. `SemiBin2_0.prots.faa`: the predicted protein sequences
2. `SemiBin2_0.prots.fna`: the nucleotide sequences corresponding to the proteins
3. `SemiBin2_0.prots.out`: contains a bit more information about the predictions

In fact, almost all of the tools that we used here internally relied on Prodigal to predict genes!

### Looking for antibiotic resistance genes

One of the things we might be interested in is whether our genomes contain antibiotic resistance genes.

For this, we will use the [Resistance Gene Identifier (RGI)](https://card.mcmaster.ca/analyze/rgi) from the CARD database.

```bash
conda activate rgi
mkdir -p 10_MAG_ANNOTATIONS/RGI
rgi main \
    --input_sequence 10_MAG_ANNOTATIONS/MAGs/SemiBin_0.fa \
    --output_file 10_MAG_ANNOTATIONS/RGI/SemiBin_0.rgi \
    --clean
```

Most of arguments are self-explanatory, but the `--clean` flag tells `rgi` to remove intermediate files.

It is very important to be aware that the results of this analysis are not definitive!
They are a starting point for further investigation.

Results will be in the file `10_MAG_ANNOTATIONS/RGI/SemiBin2_0.rgi.txt` as a table.

### Running on all genomes

You can loop over all of your genomes to predict genes and annotate them with RGI.

```bash
conda activate rgi
for f in 10_MAG_ANNOTATIONS/MAGs/SemiBin*.fa ; do
    rgi main \
        --input_sequence ${f} \
        --output_file 10_MAG_ANNOTATIONS/RGI/$(basename ${f}).rgi \
        --clean
done
```

This will produce a set of files with the antibiotic resistance genes for each genome.
You should expect to get some hits, but—as mentioned above—these are not definitive results!

### Annotating with eggnog-mapper

_This step will take too long to run, so we will include it here, but you can run it at your leisure later._

Finally, we can annotate the proteins with eggNOG.
This will give us a broad-scale annotation of the proteins, which can be used to infer the functions of the genes.

As usual, this requires a database, which we have predownloaded for you.
Otherwise, you would use the `download_eggnog_data.py` script to download the database for you (it is quite large, so make sure you have the disk space).

```bash
conda activate eggnog-mapper

mkdir 10_MAG_ANNOTATIONS/SemiBin2_0.eggnog
emapper.py \
    --itype genome \
    -i 10_MAG_ANNOTATIONS/MAGs/SemiBin_0.fa \
    --output_dir 10_MAG_ANNOTATIONS/SemiBin_0.emapper \
    --data_dir ~/Share/Databases/emapper_db \
    -o SemiBin_0
```

This will produce wide-ranging results, which can be used to infer the functions of the genes.
Again, these are not definitive results, but they are a starting point for further investigation.

