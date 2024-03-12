library(stats)
library(BiocGenerics)
library(IRanges)
library(GenomeInfoDb)
library(karyoploteR)
library(rtracklayer)

seqIds <- c(
            "NT_113891.3",
            "NT_167244.2",
            "NT_167245.2",
            "NT_167246.2",
            "NT_167247.2",
            "NT_167248.2",
            "NT_167249.2",
            "NW_025791765.1")

seqLen <- c(
            4795265,
            4672374,
            4604811,
            4677643,
            4827813,
            4606388,
            4929269,
            955087
           )

isCircular <- c(FALSE,
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

kpAxis(KpNoLabels, labels=c("chr6_scaffold_loci2"," "," ",
                            "chr6_scaffold_loci1"," "," ",
                            "chr6_scaffold_loci3"," "," ",
                            "chr6_scaffold_loci4"," "," ",
                            "chr6_scaffold_loci5"," "," ",
                            "chr6_scaffold_loci6"," "," ",
                            "chr6_scaffold_loci7"," "," ",
                            "chr2_patch"," "," ")); 

kpAddBaseNumbers(KpNoLabels, tick.dist=500000); 



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


