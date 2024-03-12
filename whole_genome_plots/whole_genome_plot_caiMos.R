library(stats)
library(BiocGenerics)
library(IRanges)
library(GenomeInfoDb)
library(karyoploteR)


############## kiyCaiMos genome 2022 ######################
seqIds <- c(
  "GWHBJBF00000007", "GWHBJBF00000011", "GWHBJBF00000014", "GWHBJBF00000023", "GWHBJBF00000024", "GWHBJBF00000029", "GWHBJBF00000032", "GWHBJBF00000038", "GWHBJBF00000046", "GWHBJBF00000055", "GWHBJBF00000061", "GWHBJBF00000075", "GWHBJBF00000080", "GWHBJBF00000096", "GWHBJBF00000100", "GWHBJBF00000120", "GWHBJBF00000131", "GWHBJBF00000235", "GWHBJBF00000276"
)
seqLen <- c(
  39155555, 21947970, 20669000, 7946471, 7757260, 6238420, 2768413, 192975, 106670, 62571, 53168, 43395, 38258, 32835, 31503, 25000, 23560, 3492, 1165
)
isCircular <- c(
  FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
)


SeqInfo <- Seqinfo(seqnames=seqIds, seqlengths=seqLen, isCircular=isCircular, genome="cr")

KpNoLabels <- plotKaryotype(SeqInfo, labels=NULL);

kpAddChromosomeNames(KpNoLabels, chr.names=c(
  "chr_7", "chr_10", "chr_13", "chr_22", "chr_23", "chr_28", "chr_31", "chr_37",  "scaf_46",  "scaf_55",  "scaf_61",  "scaf_75",  "scaf_80",  "scaf_96",  "scaf_100",  "scaf_120",  "scaf_131",  "scaf_235",  "scaf_276", " ", " "
)
); 
kpAddMainTitle(KpNoLabels, main="full_genome_coverage_caisMos_plsamid_and_genome"); 
kpAddBaseNumbers(KpNoLabels); 


plot <- kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/1_vacs_insert_sites/2_integration_site_fromChr13ToInsert/4_kizCaiMos1-0_ref/3_mapping/3_reads_on_plasmid_and_genome/readsCrTrimmedL600_mappedOn_j42Plasmid_and_genome.sorted.bam", col="darkorange");
kpAxis(plot, ymax=plot$latest.plot$computed.values$max.density)


###### only small ref seq ###### 

seqIds <- c(
  "GWHBJBF00000038", "GWHBJBF00000046", "GWHBJBF00000055", "GWHBJBF00000061", "GWHBJBF00000075", "GWHBJBF00000080", "GWHBJBF00000096", "GWHBJBF00000100", "GWHBJBF00000120", "GWHBJBF00000131", "GWHBJBF00000235", "GWHBJBF00000276"
)
seqLen <- c(
  1929750, 1066700, 625710, 531680, 433950, 382580, 328350, 315030, 250000, 235600, 34920, 11650
)
isCircular <- c(
  FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
)


SeqInfo <- Seqinfo(seqnames=seqIds, seqlengths=seqLen, isCircular=isCircular, genome="cr")

KpNoLabels <- plotKaryotype(SeqInfo, labels=NULL);

kpAddChromosomeNames(KpNoLabels, chr.names=c(
  "chr_37",  "scaf_46",  "scaf_55",  "scaf_61",  "scaf_75",  "scaf_80",  "scaf_96",  "scaf_100",  "scaf_120",  "scaf_131",  "scaf_235",  "scaf_276"
)
); 
kpAddMainTitle(KpNoLabels, main="full_genome_coverage_caisMos_plsamid_and_genome_short_refSeq"); 
kpAddBaseNumbers(KpNoLabels); 


kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/1_vacs_insert_sites/2_integration_site_fromChr13ToInsert/4_kizCaiMos1-0_ref/3_mapping/3_reads_on_plasmid_and_genome/readsCrTrimmedL600_mappedOn_j42Plasmid_and_genome.sorted.bam", col="violet");

####### only chr with high coverage 

seqIds <- c(
   "GWHBJBF00000024"
)
seqLen <- c(
  7757260
)
isCircular <- c(
  FALSE
)

SeqInfo <- Seqinfo(seqnames=seqIds, seqlengths=seqLen, isCircular=isCircular, genome="cr")

KpNoLabels <- plotKaryotype(SeqInfo, labels=NULL);

kpAddChromosomeNames(KpNoLabels,  chr.names="chr 23")

kpAddMainTitle(KpNoLabels, main="full_genome_coverage_caisMos_plsamid_and_chr23"); 
kpAddBaseNumbers(KpNoLabels); 


kp <- kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/1_vacs_insert_sites/2_integration_site_fromChr13ToInsert/4_kizCaiMos1-0_ref/3_mapping/3_reads_on_plasmid_and_genome/readsCrTrimmedL600_mappedOn_j42Plasmid_and_genome.sorted.bam", col="darkgreen");
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density)

kpCoverage <- kpPlotBAMCoverage(KpNoLabels, data="/home/michele/git/4_in_house_seq/1_vacs_insert_sites/2_integration_site_fromChr13ToInsert/4_kizCaiMos1-0_ref/3_mapping/3_reads_on_plasmid_and_genome/readsCrTrimmedL600_mappedOn_j42Plasmid_and_genome.sorted.bam", col="cornflowerblue", max.valid.region.size = 7757260)
kpAxis(kpCoverage, ymax=kpCoverage$latest.plot$computed.values$max.coverage)

