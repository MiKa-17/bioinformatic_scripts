
library(grid)
library(VennDiagram)
library(gplots)

listForDiagram4 <- list('$'(find_adapter_index_long_blast_index_1_ids_uniq, V1), '$'(find_adapter_index_long_blast_index_2_ids_uniq, V1), '$'(find_adapter_index_long_blast_index_1_reverse_ids_uniq, V1), '$'(find_adapter_index_long_blast_index_2_reverse_ids_uniq, V1))
categoryNames <- c("adapter_index_1", "adapter_index_2", "adapter_index_1_complement", "adapter_index_2_complement")
diagramOutputFile <-  '/home/michele/git/6_in_house_seq/2_vvgt_insert_sites/adapter_index_in_read.png'


transperentRed <- rgb(220, 20, 60, max=220, alpha=125, name="red50")
transperentPink <- rgb(255, 20, 147, max=255, alpha=125, name="pink50")

transperentBlue <- rgb(0, 0, 139, max=139, alpha=125, name="blue50")
transperentTeal <- rgb(0, 128, 128, max=128, alpha=125, name="teal50")


venn.diagram(x = listForDiagram4, category.names = categoryNames, filename = diagramOutputFile, col=c("#00008B", "#DC143C", "#008080", "deeppink"), fill=c(trasnparentBlue, transperentRed, transperentTeal, transperentPink) )

#### venn table 
listForDiagramTable <- list(index_1='$'(find_adapter_index_long_blast_index_1_ids_uniq, V1), index_2='$'(find_adapter_index_long_blast_index_2_ids_uniq, V1),index_1_complement='$'(find_adapter_index_long_blast_index_1_reverse_ids_uniq, V1), index_2_complement='$'(find_adapter_index_long_blast_index_2_reverse_ids_uniq, V1))
tableOutputFile_index1 <-  '/home/michele/git/6_in_house_seq/2_vvgt_insert_sites/6_find_adapter_index/reads_with_index_1_ids.csv'
tableOutputFile_index1_complement <-  '/home/michele/git/6_in_house_seq/2_vvgt_insert_sites/6_find_adapter_index/reads_with_index_1_complement_ids.csv'
tableOutputFile_index2 <-  '/home/michele/git/6_in_house_seq/2_vvgt_insert_sites/6_find_adapter_index/reads_with_index_2_ids.csv'
tableOutputFile_index2_complement <-  '/home/michele/git/6_in_house_seq/2_vvgt_insert_sites/6_find_adapter_index/reads_with_index_2_complement_ids.csv'

v.table <- venn(listForDiagramTable)
reads_index1 <- attr(v.table,"intersections")["index_1"]
reads_index1_complement <- attr(v.table,"intersections")["index_1_complement"]
reads_index2 <- attr(v.table,"intersections")["index_2"]
reads_index2_complement <- attr(v.table,"intersections")["index_2_complement"]

write.csv(reads_index1, file = tableOutputFile_index1)
write.csv(reads_index1_complement, file = tableOutputFile_index1_complement)

write.csv(reads_index2, file = tableOutputFile_index2)
write.csv(reads_index2_complement, file = tableOutputFile_index2_complement)


