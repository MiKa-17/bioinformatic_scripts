library(stats)
library(BiocGenerics)
library(IRanges)
library(GenomeInfoDb)
library(karyoploteR)
 
################################################################### sample 1 ak
###### all reads ak
seqIds <- c(
  "NC_048595.1", "NC_048596.1", "NC_048597.1", "NC_048598.1", "NC_048599.1", "NC_048600.1", "NC_048601.1", "NC_048602.1", "NC_048603.1", "NC_048604.1", "NW_023276806.1", "NW_023276807.1", "NW_023276808.1", "NW_023276850.1", "NW_023276852.1", "NW_023276874.1", "NW_023276875.1", "NW_023276883.1", "NW_023276886.1", "NW_023276919.1", "NW_023276921.1", "NW_023276997.1", "NW_023277019.1", "NW_023277053.1", "NW_023277108.1", "NW_023277152.1", "NW_023277421.1", "NW_023277426.1", "NW_023277432.1"
)
seqLen <- c(
  461620117, 282827514, 231097868, 193770019, 155611870, 134359064, 99554469, 28505831, 32558357, 127255434, 275698159, 274391693, 10418059, 80947, 80724, 71633, 447013, 70062, 438294, 9941687, 59453, 278916, 231238, 215655, 187415, 25626, 682437, 111588, 595479
)
isCircular <- c(
  FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
)

SeqInfo <- Seqinfo(seqnames=seqIds, seqlengths=seqLen, isCircular=isCircular, genome="cho")

KpNoLabels <- plotKaryotype(SeqInfo, labels=NULL);
kpAxis(KpNoLabels, labels=c(
  "chr_2", " ", " ", "chr_3", " ", " ", "chr_4", " ", " ", "chr_5", " ", " ", "chr_6", " ", " ", "chr_7", " ", " ", "chr_8", " ", " ", "chr_9", " ", " ", "chr_10", " ", " ", "chr_X", " ", " ", "chr1_0", " ", " ", "chr1_1", " ", " ", "unScaf_1", " ", " ", "unScaf_137", " ", " ", "unScaf_139", " ", " ", "unScaf_159", " ", " ", "unScaf_16", " ", " ", "unScaf_167", " ", " ", "unScaf_17", " ", " ", "unScaf_2", " ", " ", "unScaf_200", " ", " ", "unScaf_27", " ", " ", "unScaf_29", " ", " ", "unScaf_32", " ", " ", "unScaf_37", " ", " ", "unScaf_409", " ", " ", "unScaf_8", " ", " ", "unScaf_84", " ", " ", "unScaf_9", " ", " "
)); 
kpAddBaseNumbers(KpNoLabels); 
###### sample 1 ak reads mapped 
kpAddMainTitle(KpNoLabels, main="covergae_cho_genome_sample1_ak_all_reads"); 
kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/4_pcld_insert_sites_transposase/1_transposase_72h_ak/3_mapping/1_mapping_to_plasmid_or_genome/reads_cho_transposase_72h_ak_trimmed_cho_plasmids.sorted.bam", col="orange");

############# reads plasmid and genome sample1 ak
seqIds <- c(
  "NC_048595.1", "NC_048596.1", "NC_048597.1", "NC_048598.1", "NC_048599.1", "NC_048600.1", "NC_048601.1", "NC_048602.1", "NC_048603.1", "NC_048604.1", "NW_023276806.1", "NW_023276807.1"
)
seqLen <- c(
  461620117, 282827514, 231097868, 193770019, 155611870, 134359064, 99554469, 28505831, 32558357, 127255434, 275698159, 274391693
)
isCircular <- c(
  FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
)

SeqInfo <- Seqinfo(seqnames=seqIds, seqlengths=seqLen, isCircular=isCircular, genome="cho")

KpNoLabels <- plotKaryotype(SeqInfo, labels=NULL);
kpAxis(KpNoLabels, labels=c(c(
  "chr_2", " ", " ", "chr_3", " ", " ", "chr_4", " ", " ", "chr_5", " ", " ", "chr_6", " ", " ", "chr_7", " ", " ", "chr_8", " ", " ", "chr_9", " ", " ", "chr_10", " ", " ", "chr_X", " ", " ", "chr1_0", " ", " ", "chr1_1", " ", " "
)));


kpAddBaseNumbers(KpNoLabels); 
kpAddMainTitle(KpNoLabels, main="covergae_cho_genome_sample1_ak_plasmid_and_genome"); 

kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/4_pcld_insert_sites_transposase/1_transposase_72h_ak/3_mapping/3_mapping_to_plasmid_and_genome/reads_cho_transposase_72h_sample1_ak_trimmed_on_plasmids_and_genome.sorted.bam", col="violet");


############# reads plasmid and genome sample1 ak insert seq TTAA
seqIds <- c(
  "NC_048595.1", "NC_048596.1", "NC_048597.1", "NC_048598.1", "NC_048599.1", "NC_048600.1", "NC_048601.1", "NC_048602.1", "NC_048604.1", "NW_023276806.1", "NW_023276807.1", "NW_023276909.1", "NW_023277286.1"
)
seqLen <- c(
  461620117, 282827514, 231097868, 193770019, 155611870, 134359064, 99554469, 28505831, 127255434, 275698159, 274391693, 62245, 154394
)
isCircular <- c(
  FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
)

SeqInfo <- Seqinfo(seqnames=seqIds, seqlengths=seqLen, isCircular=isCircular, genome="cho")

KpNoLabels <- plotKaryotype(SeqInfo, labels=NULL);
kpAxis(KpNoLabels, labels=c(
  "chr_2", " ", " ", "chr_3", " ", " ", "chr_4", " ", " ", "chr_5", " ", " ", "chr_6", " ", " ", "chr_7", " ", " ", "chr_8", " ", " ", "chr_9", " ", " ", "chr_X", " ", " ", "chr1_0", " ", " ", "chr1_1", " ", " ", "unScaf_190", " ", " ", "unScaf_53", " ", " "
));


kpAddBaseNumbers(KpNoLabels); 
kpAddMainTitle(KpNoLabels, main="covergae_cho_genome_sample1_ak_plasmid_and_genome_insert_TTAA"); 

kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/4_pcld_insert_sites_transposase/1_transposase_72h_ak/3_mapping/4_mapping_to_plasmid_genome_insert_region/2_mapping_to_plasmid_insert_region_TTAA_and_genome/reads_cho_transposase_72h_sample1_ak_trimmed_on_plasmids_insertRegionTTAA_and_genome.sorted.bam", col="red");


################################################################### sample 2 akT

###### all reads akT
seqIds <- c(
  "NC_048595.1", "NC_048596.1", "NC_048597.1", "NC_048598.1", "NC_048599.1", "NC_048600.1", "NC_048601.1", "NC_048602.1", "NC_048603.1", "NC_048604.1", "NW_023276806.1", "NW_023276807.1", "NW_023276808.1", "NW_023276850.1", "NW_023276875.1", "NW_023276883.1", "NW_023276919.1", "NW_023276997.1", "NW_023277053.1", "NW_023277108.1", "NW_023277152.1", "NW_023277426.1"
)
seqLen <- c(
  461620117, 282827514, 231097868, 193770019, 155611870, 134359064, 99554469, 28505831, 32558357, 127255434, 275698159, 274391693, 10418059, 80947, 447013, 70062, 9941687, 278916, 215655, 187415, 25626, 111588
)
isCircular <- c(
  FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
)

SeqInfo <- Seqinfo(seqnames=seqIds, seqlengths=seqLen, isCircular=isCircular, genome="cho")

KpNoLabels <- plotKaryotype(SeqInfo, labels=NULL);

kpAxis(KpNoLabels, labels=c(
  "chr_2", " ", " ", "chr_3", " ", " ", "chr_4", " ", " ", "chr_5", " ", " ", "chr_6", " ", " ", "chr_7", " ", " ", "chr_8", " ", " ", "chr_9", " ", " ", "chr_10", " ", " ", "chr_X", " ", " ", "chr1_0", " ", " ", "chr1_1", " ", " ", "unScaf_1", " ", " ", "unScaf_137", " ", " ", "unScaf_16", " ", " ", "unScaf_167", " ", " ", "unScaf_2", " ", " ", "unScaf_27", " ", " ", "unScaf_32", " ", " ", "unScaf_37", " ", " ", "unScaf_409", " ", " ", "unScaf_84", " ", " "
)); 

kpAddBaseNumbers(KpNoLabels); 

kpAddMainTitle(KpNoLabels, main="covergae_cho_genome_sample2_akT_all_reads"); 

kpPlotBAMDensity(KpNoLabels, normalize=TRUE, data="/home/michele/git/4_in_house_seq/4_pcld_insert_sites_transposase/2_transposase_72h_akT/3_mapping/1_mapping_on_plasmid_or_genome/reads_cho_transposase_72h_sample2_akT_trimmed_cho_plasmids.sorted.bam", col="darkorange");

###### reads mapped on plasmid AND genome akT
seqIds <- c(
  "NC_048595.1", "NC_048596.1", "NC_048597.1", "NC_048598.1", "NC_048599.1", "NC_048600.1", "NC_048601.1", "NC_048602.1", "NC_048603.1", "NC_048604.1", "NW_023276806.1"
)
seqLen <- c(
  461620117, 282827514, 231097868, 193770019, 155611870, 134359064, 99554469, 28505831, 32558357, 127255434, 275698159
)
isCircular <- c(
  FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
)


SeqInfo <- Seqinfo(seqnames=seqIds, seqlengths=seqLen, isCircular=isCircular, genome="cho")

KpNoLabels <- plotKaryotype(SeqInfo, labels=NULL);

kpAxis(KpNoLabels, labels=c(
  "chr_2", " ", " ", "chr_3", " ", " ", "chr_4", " ", " ", "chr_5", " ", " ", "chr_6", " ", " ", "chr_7", " ", " ", "chr_8", " ", " ", "chr_9", " ", " ", "chr_10", " ", " ", "chr_X", " ", " ", "chr1_0", " ", " "
)); 

kpAddBaseNumbers(KpNoLabels); 
kpAddMainTitle(KpNoLabels, main="covergae_cho_genome_sample2_akT_plasmid_and_genome"); 

kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/4_pcld_insert_sites_transposase/2_transposase_72h_akT/3_mapping/3_mapping_on_plasmid_and_genome/reads_cho_transposase_72h_sample2_akT_trimmed_on_plasmids_and_genome.sorted.bam", col="darkviolet");


############# reads plasmid and genome sample2 akT insert seq TTAA
seqIds <- c(
  "NC_048595.1", "NC_048596.1", "NC_048597.1", "NC_048598.1", "NC_048599.1", "NC_048600.1", "NC_048601.1", "NC_048602.1", "NC_048603.1", "NC_048604.1", "NW_023276806.1", "NW_023276807.1", "NW_023276808.1", "NW_023277286.1", "NW_023277334.1", "NW_023277422.1"
)
seqLen <- c(
  461620117, 282827514, 231097868, 193770019, 155611870, 134359064, 99554469, 28505831, 32558357, 127255434, 275698159, 274391693, 10418059, 154394, 13320, 116830
)
isCircular <- c(
  FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
)

SeqInfo <- Seqinfo(seqnames=seqIds, seqlengths=seqLen, isCircular=isCircular, genome="cho")

KpNoLabels <- plotKaryotype(SeqInfo, labels=NULL);
kpAxis(KpNoLabels, labels=c(
  "chr_2", " ", " ", "chr_3", " ", " ", "chr_4", " ", " ", "chr_5", " ", " ", "chr_6", " ", " ", "chr_7", " ", " ", "chr_8", " ", " ", "chr_9", " ", " ", "chr_10", " ", " ", "chr_X", " ", " ", "chr1_0", " ", " ", "chr1_1", " ", " ", "unScaf_1", " ", " ", "unScaf_53", " ", " ", "unScaf_573", " ", " ", "unScaf_80", " ", " "
));


kpAddBaseNumbers(KpNoLabels); 
kpAddMainTitle(KpNoLabels, main="covergae_cho_genome_sample2_akT_plasmid_and_genome_insert_TTAA"); 

kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/4_pcld_insert_sites_transposase/2_transposase_72h_akT/3_mapping/4_mapping_on_plasmid_and_genome_insert_region_ttaa/2_mapping_on_plasmid_insert_region_ttaa_and_genome/reads_cho_transposase_72h_sample2_akT_trimmed_on_plasmids_insertRegionTTAA_and_genome.sorted.bam", col="blue");

################################################# sample 1 & 2 plots for comparison


seqIds <- c(
  "NC_048595.1", "NC_048596.1", "NC_048597.1", "NC_048598.1", "NC_048599.1", "NC_048600.1", "NC_048601.1", "NC_048602.1", "NC_048603.1", "NC_048604.1", "NW_023276806.1", "NW_023276807.1", "NW_023276808.1", "NW_023277286.1", "NW_023277334.1", "NW_023277422.1"
)
seqLen <- c(
  461620117, 282827514, 231097868, 193770019, 155611870, 134359064, 99554469, 28505831, 32558357, 127255434, 275698159, 274391693, 10418059, 154394, 13320, 116830
)
isCircular <- c(
  FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
)

SeqInfo <- Seqinfo(seqnames=seqIds, seqlengths=seqLen, isCircular=isCircular, genome="cho")

KpNoLabels <- plotKaryotype(SeqInfo, labels=NULL);

kpAddChromosomeNames(KpNoLabels, chr.names=c(
  "chr_2", "chr_3", "chr_4", "chr_5", "chr_6", "chr_7", "chr_8", "chr_9", "chr_10", "chr_X", "chr1_0", "chr1_1", "unScaf_1", "unScaf_53", "unScaf_573", "unScaf_80", " ", " "
)
); 

kpAddMainTitle(KpNoLabels, main="density_cho_genome_sample1_ak&sample2_akT_plasmid_and_genome_insert_TTAA"); 
kpAddBaseNumbers(KpNoLabels); 

dataset1 = "/home/michele/git/4_in_house_seq/4_pcld_insert_sites_transposase/2_transposase_72h_akT/3_mapping/4_mapping_on_plasmid_and_genome_insert_region_ttaa/2_mapping_on_plasmid_insert_region_ttaa_and_genome/reads_cho_transposase_72h_sample2_akT_trimmed_on_plasmids_insertRegionTTAA_and_genome.sorted.bam"

dataset2 = "/home/michele/git/4_in_house_seq/4_pcld_insert_sites_transposase/1_transposase_72h_ak/3_mapping/4_mapping_to_plasmid_genome_insert_region/2_mapping_to_plasmid_insert_region_TTAA_and_genome/reads_cho_transposase_72h_sample1_ak_trimmed_on_plasmids_insertRegionTTAA_and_genome.sorted.bam"

kpPlotAkT <- kpPlotBAMDensity(KpNoLabels, data=dataset1, col= rgb(0, 0, 1, alpha = 0.5));
kpPlotAk <- kpPlotBAMDensity(KpNoLabels, data=dataset2, col=rgb(1, 0, 0, alpha = 0.5));

kpAxis(kpPlotAk, ymax=kpPlotAk$latest.plot$computed.values$max.density)

legend(x = "bottomright", fill = c("blue", "red"), legend = c("smaple1_akT", "sample2_ak"))

### coverage plot single chr 




seqIds <- c(
  "NC_048595.1", "NC_048596.1", "NC_048597.1", "NC_048598.1", "NC_048599.1", "NC_048600.1", "NC_048601.1", "NC_048602.1", "NC_048603.1", "NC_048604.1", "NW_023276806.1", "NW_023276807.1", "NW_023276808.1", "NW_023277286.1", "NW_023277334.1", "NW_023277422.1"
)
seqLen <- c(
  461620117, 282827514, 231097868, 193770019, 155611870, 134359064, 99554469, 28505831, 32558357, 127255434, 275698159, 274391693, 10418059, 154394, 13320, 116830
)
isCircular <- c(
  FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
)

# Chr 4
SeqInfo_chr4 <- Seqinfo(seqnames="NC_048597.1", seqlengths=282827514, isCircular=FALSE, genome="cho")

KpNoLabels <- plotKaryotype(SeqInfo_chr3, labels=NULL);

kpAddChromosomeNames(KpNoLabels, chr.names="chr_3"); 


kpAddMainTitle(KpNoLabels, main="covergae_cho_chr3_sample1_ak&sample2_akT_plasmid_and_genome_insert_TTAA"); 
kpAddBaseNumbers(KpNoLabels); 

dataset1 = "/home/michele/git/4_in_house_seq/4_pcld_insert_sites_transposase/2_transposase_72h_akT/3_mapping/4_mapping_on_plasmid_and_genome_insert_region_ttaa/2_mapping_on_plasmid_insert_region_ttaa_and_genome/reads_cho_transposase_72h_sample2_akT_trimmed_on_plasmids_insertRegionTTAA_and_genome.sorted.bam"

dataset2 = "/home/michele/git/4_in_house_seq/4_pcld_insert_sites_transposase/1_transposase_72h_ak/3_mapping/4_mapping_to_plasmid_genome_insert_region/2_mapping_to_plasmid_insert_region_TTAA_and_genome/reads_cho_transposase_72h_sample1_ak_trimmed_on_plasmids_insertRegionTTAA_and_genome.sorted.bam"



kpCoverageAkT <- kpPlotBAMCoverage(KpNoLabels, data=dataset1, col="cornflowerblue", max.valid.region.size = 282827514)
kpCoverageAk <- kpPlotBAMCoverage(KpNoLabels, data=dataset1, col="darkred", max.valid.region.size = 282827514)
kpAxis(kpCoverageAkT, ymax=kpCoverage$latest.plot$computed.values$max.coverage)

legend(x = "bottomright", fill = c("cornflowerblue", "darkred"), legend = c("smaple1_akT", "sample2_ak"))

