library(stats)
library(BiocGenerics)
library(IRanges)
library(GenomeInfoDb)
library(karyoploteR)
library(rtracklayer)

seqIds <- c(
  "NC_000001.11", "NC_000002.12", "NC_000003.12", "NC_000004.12", "NC_000005.10", "NC_000006.12", "NC_000007.14", "NC_000008.11", "NC_000009.12", "NC_000010.11", "NC_000011.10", "NC_000012.12", "NC_000013.11", "NC_000014.9", "NC_000015.10", "NC_000016.10", "NC_000017.11", "NC_000018.10", "NC_000019.10", "NC_000020.11", "NC_000021.9", "NC_000022.11", "NC_000023.11", "NC_000024.10", "NC_012920.1"
)

seqLen <- c(
  248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415, 16569
)

isCircular <- c(FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE
                )

SeqInfo <- Seqinfo(seqnames=seqIds, seqlengths=seqLen, isCircular=isCircular, genome="human_genome")

KpNoLabels <- plotKaryotype(SeqInfo, labels=NULL);

kpAddChromosomeNames(KpNoLabels, chr.names=c(
  "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "mitrochondion")
); 


kpAddBaseNumbers(KpNoLabels); 



################################ 9c10 clone ########################################
#### new run, new basecalled with gag primer 
# all reads
kpAddMainTitle(KpNoLabels, main="full_genome_coverage_human_9c10clone_trimmed_gag_all"); 
kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/2_vvgt_insert_sites/2_repetition_run_0823/new_basecalled/4_9c10clone_gag_primer/3_mapping/1_all_reads_trimmed/reads_trimmed_9c10clone_gag_human_genome_plasmids.sorted.bam", col="violet");
# reads on plasmid AND genome 
kpAddMainTitle(KpNoLabels, main="full_genome_coverage_human_9c10clone_trimmed_gag_plasmidAndGenome"); 
kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/2_vvgt_insert_sites/2_repetition_run_0823/new_basecalled/4_9c10clone_gag_primer/3_mapping/3_reads_on_plasmid_genome/reads_trimmed_9c10clone_on_pasmid_pLL37_kym_and_genome.sorted.bam", col="darkviolet");


### new run, new basecalled with ltr primer 
# all reads
kpAddMainTitle(KpNoLabels, main="full_genome_coverage_human_9c10clone_trimmed_ltr_all"); 
kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/2_vvgt_insert_sites/2_repetition_run_0823/new_basecalled/2_9c10clone_ltr_primer/3_mapping/1_all_reads_trimmed/reads_trimmed_9c10clone_ltr_human_genome_plasmid_pll37_kym.sorted.bam", col="darkgreen");
# reads on plasmid AND genome 
kpAddMainTitle(KpNoLabels, main="full_genome_coverage_human_9c10clone_trimmed_ltr_plasmidAndGenome"); 
kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/2_vvgt_insert_sites/2_repetition_run_0823/new_basecalled/2_9c10clone_ltr_primer/3_mapping/3_reads_on_plasmid_genome/reads_trimmed_9c10clone_on_pasmid_pll37_and_genome.sorted.bam", col="#669966");

#### new run, new basecalled with combine ltr and gag primer 
# reads on plasmid AND genome 
kpAddMainTitle(KpNoLabels, main="full_genome_coverage_human_9c10clone_trimmed_ltrAndGag_plasmidAndGenome"); 
kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/2_vvgt_insert_sites/2_repetition_run_0823/new_basecalled/2_9c10clone_ltr_primer/3_mapping/3_reads_on_plasmid_and_genome/reads_trimmed_9c10clone_on_pasmid_pll37_and_genome.sorted.bam", col="#669966");
kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/2_vvgt_insert_sites/2_repetition_run_0823/new_basecalled/4_9c10clone_gag_primer/3_mapping/3_reads_on_plasmid_and_genome/reads_trimmed_9c10clone_on_pasmid_pLL37_kym_and_genome.sorted.bam", col="darkviolet");
legend(x = "bottomright", fill = c("#669966", "darkviolet"), legend = c("ltr_primer", "gag_primer"))

######################################################################################

################################ pLVTHM ########################################
#### new run, new basecalled with gag primer 
# all reads
kpAddMainTitle(KpNoLabels, main="full_genome_coverage_human_pLVTHM_trimmed_gag_all"); 
kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/2_vvgt_insert_sites/2_repetition_run_0823/new_basecalled/3_lvthm_gag_primer/3_mapping/1_reads_on_plasmid_or_genome/reads_trimmed_plvthm_gag_human_genome_plasmid.sorted.bam", col="violet");
# reads on plasmid AND genome 
kpAddMainTitle(KpNoLabels, main="full_genome_coverage_human_pLVTHM_trimmed_gag_plasmidAndGenome"); 
kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/2_vvgt_insert_sites/2_repetition_run_0823/new_basecalled/3_lvthm_gag_primer/3_mapping/3_reads_on_plasmid_and_genome/reads_trimmed_on_plvthm_and_genome.sorted.bam", col="darkviolet");


### new run, new basecalled with ltr primer 
# all reads
kpAddMainTitle(KpNoLabels, main="full_genome_coverage_human_pLVTHM_trimmed_ltr_all"); 
kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/2_vvgt_insert_sites/2_repetition_run_0823/new_basecalled/1_lvthm_ltr_primer/3_mapping/1_reads_on_plasmid_or_genome/reads_trimmed_plvthm_ltr_human_genome_plasmid.sorted.bam", col="darkgreen");
# reads on plasmid AND genome 
kpAddMainTitle(KpNoLabels, main="full_genome_coverage_human_pLVTHM_trimmed_ltr_plasmidAndGenome"); 
kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/2_vvgt_insert_sites/2_repetition_run_0823/new_basecalled/1_lvthm_ltr_primer/3_mapping/3_reads_on_plasmid_and_genome/reads_trimmed_on_plvthm_and_genome.sorted.bam", col="#669966");

#### new run, new basecalled with combine ltr and gag primer 
# reads on plasmid AND genome 
kpAddMainTitle(KpNoLabels, main="full_genome_coverage_human_pLVTHM_trimmed_ltrAndGag_plasmidAndGenome"); 
kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/2_vvgt_insert_sites/2_repetition_run_0823/new_basecalled/1_lvthm_ltr_primer/3_mapping/3_reads_on_plasmid_and_genome/reads_trimmed_on_plvthm_and_genome.sorted.bam", col="#669966");
kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/2_vvgt_insert_sites/2_repetition_run_0823/new_basecalled/3_lvthm_gag_primer/3_mapping/3_reads_on_plasmid_and_genome/reads_trimmed_on_plvthm_and_genome.sorted.bam", col="darkviolet");
legend(x = "bottomright", fill = c("#669966", "darkviolet"), legend = c("ltr_primer", "gag_primer"))


######################################################################################

################### Lentivirus integration
### new run, new basecalled with ltr primer 
# all reads
allMappedReads <- "/home/michele/1_bioinformatic_projects/4_in_house_seq/7_vvgt_lentiviral_integration/3_mapping/1_reads_on_plasmid_or_genome/LvthmInsertReads_trimmed_l750_humanGenome_pLVTHM.sorted.bam"
kpAddMainTitle(KpNoLabels, main="full_genome_coverage_human_pLVTHM_lentiviral_all"); 
allReadsPlot <- kpPlotBAMDensity(KpNoLabels, data=allMappedReads, col="darkorange");
kpAxis(allReadsPlot, ymax=allReadsPlot$latest.plot$computed.values$max.density)

# reads on plasmid AND genome 
breakpointReads <- "/home/michele/1_bioinformatic_projects/4_in_house_seq/7_vvgt_lentiviral_integration/3_mapping/3_reads_on_plasmid_and_genome/LvthmInsertReads_trimmed_l750_on_plasmid_and_genome.sorted.bam"
kpAddMainTitle(KpNoLabels, main="full_genome_coverage_human_pLVTHM_lentiviral_breakpoints"); 
breakponitsPlot <- kpPlotBAMDensity(breakponitsPlot, data=breakpointReads, col="darkviolet");
kpAxis(breakponitsPlot, ymax=breakponitsPlot$latest.plot$computed.values$max.density)

