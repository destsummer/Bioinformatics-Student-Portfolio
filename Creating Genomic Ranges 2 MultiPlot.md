Similar to GR 1 in my repository this creates a similar genomic ranges object that looks at treated vs untreated and then plots based off reference genome.
```
#load necessary libraries 
library(airway)
library(ggplot2)
library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(BiocParallel)
```

```
#Specify pathway in the Airway package
dir <- system.file("extdata",package="airway",mustWork=TRUE)
#Verify that files saved properly 
list.files(dir)
#create file path for csv table
csvfile <- file.path(dir,"sample_table.csv")
(sampleTable <- read.csv(csvfile))
#divide treated from untreated in the csv
#Useing the "dex" column in sampleTable to identify the bam files for the treated and untreated samples.
sampletrt <- sampleTable[sampleTable$dex == "trt",]
sampleuntrt <- sampleTable[sampleTable$dex == "untrt",]
#make into bam file
filenamestrt <- file.path(dir, paste0(sampletrt$Run, "_subset.bam"))
filenamesuntrt <- file.path(dir, paste0(sampleuntrt$Run, "_subset.bam"))
#ensure that files were created properly, should return true
file.exists(filenamestrt)
file.exists(filenamesuntrt)
#run bam file in chunks based off of yield size
bamfilestrt <- BamFileList(filenamestrt, yieldSize=2000000)
bamfilesuntrt <- BamFileList(filenamesuntrt, yieldSize=2000000)
#list out the seqinfo for the bam files we just created, treated and untreated 
seqinfo(bamfilestrt[1])
seqinfo(bamfilesuntrt[1])
#load file for reference genome
gtffile <- "/usr/lib64/R/library/airway/extdata/Homo_sapiens.GRCh37.75_subset.gtf"
#take the gtf file and make into a txdb object
(txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))
#filter txdb transcripts by gene and save as variable 
(tbg <- transcriptsBy(txdb, by="gene"))
```

```
#Create a GenomicRanges object with all of the reads from all of the treated samples and another one with the reads from all of the untreated samples. 
#The samples are paired end.
treated <- readGAlignmentPairs(bamfilestrt[[1]])
untreated <- readGAlignmentPairs(filenamesuntrt[[1]])
#combine all the reads together into one gr object
comboGR <- c(treated,untreated)
```

```
# Create or load a GenomicRanges object with the coordinates of the genes(this is already done above)
# The transcripts work better than exons when visualizing the gene models
gtffile <- "/usr/lib64/R/library/airway/extdata/Homo_sapiens.GRCh37.75_subset.gtf"
(txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))
(tbg <- transcriptsBy(txdb, by="gene"))
```

```
# Create a gr for your gene of interest.
# Create a genomic ranges object where each row lists the coordinates of an individual base in the gene. 
#I am using the first gene ENSG000000009724 and first transcript ENST00000400897 to make a gr for 
#setting the start as start and end, allows to divide the coordinates up by bases
Example <- tbg[[1]][1]
chromosome <- seqnames(Example) 
trans <- paste0(chromosome)
starts <- c(start(Example):end(Example))
ends <- starts
strands <- rep(strand(Example), length(starts))
gr <- GRanges(seqnames=trans,ranges=IRanges(start=starts,end=ends),strand=strands)
```

```
# Find the number of reads from treated and untreated samples that overlap each row in your genomic ranges object.
#countOverlaps() or summarizeOverlaps() from GenomicRanges
variablealn <- readGAlignments(bamfilestrt[[1]])
counttrt <- countOverlaps(gr, variablealn, type="any", ignore.strand=FALSE)
variablealn2 <- readGAlignments(bamfilesuntrt[[1]])
countuntrt <- countOverlaps(gr, variablealn2, type="any", ignore.strand=FALSE)
```

```
# Plot the normalized number of reads in each group at each base across the gene. 
#take the overlaps table and make into a data frame
#use the starts of the coordinates to be the x axis
#the overlap read that was produced using the above code is then used as the y axis
#make two graphs for treated vs untreated
library(ggplot2)
#treated graph
overlapsdftrt <- as.data.frame(counttrt)
x <- starts
ggplot(overlapsdftrt, aes(x,counttrt)) + geom_bar(stat = "Identity") + xlab('Coordinates of Bases') + ylab('Treated Counts/1M reads')

#untreated graph
overlapsdfuntrt <- as.data.frame(countuntrt)
x <- starts
ggplot(overlapsdfuntrt, aes(x,countuntrt)) + geom_bar(stat = "Identity") + xlab('Coordinates of Bases') + ylab('Untreated Counts/1M reads')
```

```
# Plot the transcript models.
#gr object we created earlier determines the coordinates we want to use
#txdb pulls the transcripts from those coordinates
ggplot() + geom_alignment(txdb, which = gr) + xlab('Coordinates of Bases') + ylab('transcript isoforms')
```

```
# Combine the plots into a single pdf 
#assign each of the plots a variable to make it easier for calling in the multiplot function

p1 <- ggplot(overlapsdftrt, aes(x,counttrt)) + geom_bar(stat = "Identity") + xlab('Coordinates of Bases') + ylab('Treated Counts/1M reads')

p2 <- ggplot(overlapsdfuntrt, aes(x,countuntrt)) + geom_bar(stat = "Identity") + xlab('Coordinates of Bases') + ylab('Untreated Counts/1M reads')

p3 <- ggplot() + geom_alignment(txdb, which = gr) + xlab('Coordinates of Bases') + ylab('transcript isoforms')  

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#if on server: can change name of pdf 
pdf("~/multiplots.pdf")
#each graph is put into the display
#cols determine the layout of the display. 1 column means the graphs will be stacked on each other
multiplot(p1,p2,p3,cols=1)
dev.off()
dev.off()
```
