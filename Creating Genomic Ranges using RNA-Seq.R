This file uses varies R packages in order to load in data from a standard Airways package, create BAM files, read coverage and create a plot with the corresponding data.
```{r}
#load necessary libraries 
library(airway)
library(ggplot2)
library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(BiocParallel)
```
```{r}
#Specify pathway in the Airway package 
dir <- system.file("extdata",package="airway",mustWork=TRUE)
#Verify that files saved properly 
list.files(dir)
#create file path for csv table
csvfile <- file.path(dir,"sample_table.csv")
#read in the data from the csv
(sampleTable <- read.csv(csvfile,row.names=1))
#create table for only bam
filenames <- file.path(dir, paste0(sampleTable$Run, "_subset.bam"))
#verify that saved properly
file.exists(filenames)
#Object that saves all data from the bam files 
bamfiles <- BamFileList(filenames, yieldSize=2000000)
#sequence information about the first bam file
seqinfo(bamfiles[1])
#load in reference genome
gtffile <- "/usr/lib64/R/library/airway/extdata/Homo_sapiens.GRCh37.75_subset.gtf"
#object that saves all transcript information from previously save gtffile
(txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))
#create genomic ranges object for only exons
(ebg <- exonsBy(txdb, by="gene"))
# A genomic ranges object for the first exon in the first gene
exonExample <- ebg[[1]][1]
```

```{r}
#Generic way to create a genomic ranges object:
gr <- GRanges(seqnames=chr,ranges=IRanges(start=range,end=range),strand=str)
# seqnames is a vector of chromosome names, start is vector of chromosomal start coordinates, end is a vector of end coordinates, and strand is a vector of "+", "-", or "*")
chrom <- paste0("chr",c(1:6))
starts <- c(4000001:4000006)
ends <- c(4000101:4000106)
strands <- rep("+",6)
grTest <- GRanges(seqnames=chrom,ranges=IRanges(start=starts,end=ends),strand=strands)
```

```{r}
#Create a genomic ranges object where each row lists the 
#coordinates of an individual base in exonExample
#the starting base for this exon is 11086580 and the end is 11087705
#set both of these coordinates
chromosome <- seqnames(exonExample) #or do <- rep(as.character(seqnames(exonExample)))
exon <- paste0(chromosome)
starts <- c(start(exonExample):end(exonExample))
ends <- starts
strands <- rep(strand(exonExample), length(starts))
grTest2 <- GRanges(seqnames=exon,ranges=IRanges(start=starts,end=ends),strand=strands)
```

```{r}
#Find the number of reads in bamfiles[1] that overlap each row in your genomic ranges object
# Use countOverlaps() or summarizeOverlaps() from GenomicRanges
#summarize = granges object
#countoverlaps = bam file
variablealn <- readGAlignments(bamfiles[[1]])
count <- countOverlaps(grTest2, variablealn, type="any", ignore.strand=FALSE)
```

```{r}
#Histogram that plots the number of reads in a single bam file across exonExample 
#genomic ranges to data frame(i.e as.dataframe(grTest2))
library(ggplot2)
overlapsdf <- as.data.frame(count)
#start and end as the same since looking at single bases
x <- starts
grhistogram <- ggplot(overlapsdf, aes(x,count)) + geom_bar(stat = "Identity") + xlab('GR coordinate positions') + ylab('BAM file coverage')
#if working on server
#pdf("~/Grhistogram.pdf")
#grhistogram
#dev.off()
```
