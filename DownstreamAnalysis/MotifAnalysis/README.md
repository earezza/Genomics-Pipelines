## Motif Analysis  
___  
### Software Required  
<a href="https://bedtools.readthedocs.io/en/latest/content/tools/slop.html">Bedtools slop</a> to extend peak regions surrounding peak summits.  
<a href="https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html">Bedtools getfasta</a> to map genomic sequences to the peak regions.  
<a href="https://meme-suite.org/meme/doc/meme-chip.html?man_type=web">MEME-ChIP from MEME Suite</a> to perform the motif analysis.  

### Reference Files  
<a href="https://meme-suite.org/meme/db/motifs">Motif database files</a>  
<a href="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/">UCSC hg38 genome sequence files</a>  
<a href="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/">UCSC hg19 genome sequence files</a>  
<a href="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/">UCSC mm10 genome sequence files</a>  

### Example Usage  
Extend summit range for peaks  
>bedtools slop -i summits.bed -g Reference_Files/mm10.chrom.sizes -b 500 > summits_slop500.bed

Get sequence for extended summits (use masked genome sequences)  
>bedtools getfasta -fi Reference_Files/mm10.fa.masked -bed summits_slop500.bed -fo summits_slop500.fasta

Run MEME-ChIP  
>meme-chip summits_slop500.fasta -oc MEME-OUTPUT/ -ccut 500 -meme-nmotifs 15 -meme-mod anr -minw 15 -maxw 25 -db Reference_Files/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme
