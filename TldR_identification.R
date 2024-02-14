
library(tidyverse)
library(Biostrings)
library(DECIPHER)
library(biomartr)
library(ggtree)
# library(ape)

###############################################################
# dTnpB Identification

# NOTE: This file is not meant to be run as a script. There are steps that require the running of command line applications
  # NOTE (cont.): outside of R. These steps are annotated/explained as comments

setwd("/Users/tannerwiegand/Documents/Projects/TnpB/dTnpB/Code_for_Github/dTnpB_identification")
###############################################################

rm(list = ls())

############################################
# Determine which homologs in TnpB tree are dTnpB: Initial look
############################################

# Import trimmed alignment
algn <- readAAStringSet("inputs/ALL_tnpB_sequences_250aa_GB_50_clustered_mafft_alignment_ENISI_4_iterations_trimmed_g05.afa")

# Extract columns corresponding to each catalytic RuvC residue
df <- data.frame(id = names(algn), 
                 D191 = unname(as.character(narrow(algn, start = 321, end = 321))),
                 E278 = unname(as.character(narrow(algn, start = 434, end = 434))),
                 D361 = unname(as.character(narrow(algn, start = 548, end = 548))),
                 stringsAsFactors = FALSE)

# Add column that indicates whether or not a full catalytic triad is present (conservative substitutions allowed)
sim <- c("D", "E")
df$full.triad <- FALSE ; df$sum.triad <- 0
for(i in 1:nrow(df)) { 
  if(sum(df$D191[i] %in% sim,  df$E278[i] %in% sim, df$D361[i] %in% sim) == 3) {df$full.triad[i] <- TRUE}
  df$sum.triad[i] <- sum(df$D191[i] %in% sim,  df$E278[i] %in% sim, df$D361[i] %in% sim) }

# Import TnpB tree
tree <- ape::read.tree("inputs/TnpB.tre")

# Make list for OTU grouping (i.e., colored branches)
triad <- list()
triad[["Incomplete"]] <- df$id[df$full.triad == FALSE]
triad[["Complete"]] <- df$id[df$full.triad == TRUE]

# Visualize tree with colored branches for catalytic triad completeness (RuvC)
p <- ggtree(tree, layout = "circular")
groupOTU(p, triad, "Catalytic_triad") + aes(color = Catalytic_triad) + theme(legend.position = "right") #+ geom_tiplab2()

# Reorder triad info to be consistent with branch order
x <- data.frame(id = get_taxa_name(p))
triad.df <- x %>% left_join(df, by = "id")

# Reorder MSA to be consistent with branch order
y <- data.frame(id = names(algn), seq = unname(as.character(algn)), stringsAsFactors = FALSE) 
x <- x %>% left_join(y, by = "id")
triad.algn <- AAStringSet(x$seq)
names(triad.algn) <- x$id

# Export
  write_tsv(triad.df, "outputs/Catalytic_triad_completeness.firstlook.tsv")
  writeXStringSet(triad.algn, "outputs/Catalytic_triad_completeness.firstlook.orderedMSA.trimmed.afa")


# Manually inspected list of sequences/triads and tree
  # Can't tell if a second acidic residue that is misaligned in the MSA corresponds to the middle residue of the triad

# Run AlphaFold for representatives on AWS (Manually made FASTA files for each input)
  # python3 /home/ubuntu/work/alphafold_program/docker/run_docker.py \
  # --fasta_paths=ADL43404.faa,BCS83717.faa,WP_011231306.faa,WP_148882482.faa,WP_206406690.faa \
  # --data_dir=/home/ubuntu/work/alphafold_databases/ \
  # --max_template_date=2024-01-01 \
  # --output_dir=/home/ubuntu/work/alphafold/ \
  # --log_dir=/home/ubuntu/work/alphafold/ \
  # --model_preset=monomer
  
# Okay. So ALL of the structural predictions imply that the major groups on the tree I was thinking might have incomplete
  # triads actually have a compensatory residue that completes the triad, even though they're often not aligned in the MSA

# Moved files re: this first attempt at identifying catalytically inactive TnpBs to a new directory (including rank0 AFold models)
  # outputs/first_look/AFold/
  # outputs/first_look/Catalytic_triad_completeness.firstlook.tsv
  # outputs/first_look/Catalytic_triad_completeness.firstlook.orderedMSA.trimmed.afa  



############################################
# Determine which homologs in the TnpB tree are dTnpB: Strategy 1 [REMOVE sequences with gaps at expected DED positions]
############################################

# Analyze TnpB sequences with MULTIPLE mutated catalytic residues

# Import trimmed alignment
algn <- readAAStringSet("inputs/ALL_tnpB_sequences_250aa_GB_50_clustered_mafft_alignment_ENISI_4_iterations_trimmed_g05.afa")

# Extract columns corresponding to each catalytic RuvC residue: Use WP_010887311.1 (ISDra2 TnpB) as reference
df <- data.frame(id = names(algn), 
                 D191 = unname(as.character(narrow(algn, start = 321, end = 321))),
                 E278 = unname(as.character(narrow(algn, start = 434, end = 434))),
                 D361 = unname(as.character(narrow(algn, start = 548, end = 548))),
                 stringsAsFactors = FALSE)

# Add column that indicates whether or not a full catalytic triad is present (conservative substitutions allowed)
sim <- c("D", "E")
df$full.triad <- FALSE ; df$sum.triad <- 0
for(i in 1:nrow(df)) { 
  if(sum(df$D191[i] %in% sim,  df$E278[i] %in% sim, df$D361[i] %in% sim) == 3) {df$full.triad[i] <- TRUE}
  df$sum.triad[i] <- sum(df$D191[i] %in% sim,  df$E278[i] %in% sim, df$D361[i] %in% sim) }

# Remove sequences that have gaps at expected positions of catalytic resides
df.all <- df
df <- df %>% filter(!D191 == "-") %>% filter(!E278 == "-") %>% filter(!D361 == "-")

# Pull sequences with one or less catalytic residues conserved (at least as I can tell in the MSA used to build the TnpB tree) 
these <- df %>% filter(sum.triad <= 1)  

# Add sequences that have a mutation in either of the first two catalytic acidic residues in the triad (D191, E278)
# these <- rbind(these, df %>% filter(sum.triad == 2, (!D191 %in% sim | !D361 %in% sim)))

# Import TnpB tree
tree <- ape::read.tree("inputs/TnpB.tre")

# Make list for OTU grouping (i.e., colored branches)
triad <- list()
triad[["Incomplete"]] <- these$id
triad[["Complete"]] <- df.all$id[!df.all$id %in% these$id]

# Visualize tree with colored branches for catalytic triad completeness (RuvC)
df.lab <- cbind(df %>% select(label = id), df %>% select(id))
p <- ggtree(tree, layout = "circular") %<+% df.lab
groupOTU(p, triad, "Catalytic_triad") + aes(color = Catalytic_triad) + theme(legend.position = "right") +
  geom_tiplab2(label = id)

# Change values of cataltic triad: TRUE/FALSE
df.all$full.triad <- TRUE
df.all$full.triad[df.all$id %in% these$id] <- FALSE

# Reorder triad info to be consistent with branch order
x <- data.frame(id = get_taxa_name(p))
triad.df <- x %>% left_join(df.all, by = "id")

# Add position in tree to dataframe
triad.df <- triad.df %>% mutate(tree.position = 1:nrow(triad.df))

# Import cluster info
clus <- read_tsv("inputs/ALL_tnpB_sequences_250aa_55_clustered.df.tsv")

# Add cluster info to each row in the triad dataframe
triad.df <- triad.df %>% mutate(cluster.size = NA, mean.identity = NA, mean.seq.cov = NA, mean.seq.length = NA)
for(i in 1:nrow(triad.df)){
  
  # Extract only the cluster info for the representative
  x <- clus %>% filter(Group_Rep == clus$Group_Rep[clus$id == triad.df$id[i]])
  
  # Convert cluster identity and coverage from character vectors to numeric
  x$clstr_iden <- as.numeric(gsub("%", "", x$clstr_iden))
  x$clstr_cov <- as.numeric(gsub("%", "", x$clstr_cov))
  
  # Fill in each field
  triad.df$cluster.size[i] <- nrow(x)
  triad.df$mean.identity[i] <- paste0(round(mean(x$clstr_iden),2), "% [", min(x$clstr_iden), "-", max(x$clstr_iden), "%]")
  triad.df$mean.seq.cov[i] <- paste0(round(mean(x$clstr_cov),2), "% [", min(x$clstr_cov), "-", max(x$clstr_cov), "%]")
  triad.df$mean.seq.length[i] <- paste0(round(mean(x$length),2), " [", min(x$length), "-", max(x$length), "]")
  
}

# Export
write_tsv(triad.df, "outputs/Catalytic_triad_completeness.info.tsv")


############################################
# Pull genomes for each dTnpB that has more than a single representative per cluster: Strategy 1
############################################

# Start clean
rm(list = ls())

# Make new directory
dir.create("outputs/neighborhoods")

# Import cluster info
clus <- read_tsv("inputs/ALL_tnpB_sequences_250aa_55_clustered.df.tsv")

# Import triad info
df <- read_tsv("outputs/Catalytic_triad_completeness.info.tsv")

# Narrow cluster info to include only clusters with incomplete triads
x <- df$id[df$full.triad == FALSE] # & df$cluster.size > 1] # turn this on to exclude singletons
clus <- clus %>% filter(Group_Rep %in% clus$Group_Rep[clus$id %in% x])

# Export info on clusters with incomplete triads
write_tsv(clus, "outputs/neighborhoods/dtnpb.cluster.info.tsv")

# Export protein accessions
write_tsv(clus %>% select(id), "outputs/neighborhoods/dtnpb.acc.txt", col_names = FALSE)

# Use Batch-Entrez (online) to pull IPG info for each sequence
  # outputs/dtnpb.acc.txt
  # 40 records were obsolete

# Import IPG info
ipg <- read_tsv("outputs/neighborhoods/dtnpb.ipg.tsv")

# Export genomic accessions
write_tsv(ipg %>% select(`Nucleotide Accession`) %>% filter(!`Nucleotide Accession` == "N/A"),
          "outputs/neighborhoods/dtnpb.genomic.acc.txt", col_names = FALSE)

# Use Batch-Entrez to pull IPG info for each sequence
  # outputs/neighborhoods/dtnpb.genomic.fna
  # 343 records were obsolete/irretrievable

# Import genomic loci
genomes <- readDNAStringSet("outputs/neighborhoods/dtnpb.genomic.fna")

# Pull flanks of dTnpB
flanks <- DNAStringSet()
flank.size <- 20000
df.flank <- tibble()
for(i in 1:length(genomes)){
  
  # Isolate one genome
  x <- genomes[i]
  
  # Get annotations for that genome
  annot <- ipg %>% filter(`Nucleotide Accession` == word(names(x)))
  
  # Some genomes have 2 copies
  for(v in 1:nrow(annot)){
    
    # # Trim genome to only include dTnpB + max length of sequence (if 3kbp is longer than end) 
    if(annot$Stop[v] + flank.size > width(x)){ flank.end <- width(x) } else { flank.end <- annot$Stop[v] + flank.size }
    if(annot$Start[v] - flank.size < 1){ flank.start <- 1 } else { flank.start <- annot$Start[v] - flank.size }
    
    # Trim genome to only include dTnpB + 3kbp flanks
    this.flank <- narrow(x, start = flank.start, end = flank.end)
    
    # Add coordinates to FASTA header
    names(this.flank) <- paste0(word(names(this.flank)), "_", flank.start, "-", end = flank.end)
    
    # Flip sequence, if the dTnpB is on the "-" strand
    if(annot$Strand[v] == "-") { this.flank <- reverseComplement(this.flank) }
    
    # Add flank to new FASTA
    flanks <- c(flanks, this.flank) 
    
    # Add flank info to dataframe
    df.flank <- rbind(df.flank, tibble(genome = names(x), flank.id = paste0(word(names(this.flank))),
                                       length = width(x), flank.start = flank.start, flank.end = flank.end, 
                                       dtnpb.start = annot$Start[v], dtnpb.end = annot$Stop[v]))
    
  }
  
  # Start fresh, to raise errors if something goes wrong with the next genome
  rm(x, annot, this.flank, flank.start, flank.end)
  
}

# Determine flank sizes
df.flank <- df.flank %>% mutate(flank.size = flank.end - flank.start + 1) %>% distinct()

# Remove duplicates
flanks <- unique(flanks)
df.flank <- df.flank %>% filter(flank.id %in% names(flanks))

# Export flanks and info
writeXStringSet(flanks, "outputs/neighborhoods/dtnp.genomic.20kbp_flanks.fna")
write_tsv(df.flank, "outputs/neighborhoods/dtnp.genomic.20kbp_flanks.info.tsv")


############################################
# Annotate flanks with Eggnogg and visualize: Strategy 1
############################################  

# Start clean
rm(list = ls())

# Make directory for results
dir.create("outputs/neighborhoods/20kbp_flanks_Eggnog", showWarnings = FALSE)

# Run Eggnog mapper (command line)
# emapper.py --cpu 28 -o out \
# --output_dir 20kbp_flanks_Eggnog \
# --override -m diamond --dmnd_ignore_warnings \
# -i dtnp.genomic.20kbp_flanks.fna \
# --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 \
# --itype genome --genepred prodigal --tax_scope auto --target_orthologs all \
# --go_evidence non-electronic --pfam_realign none --report_orthologs \
# --decorate_gff yes \
# --excel

# Import GFF file of eggnog output
egg <- read_gff("outputs/neighborhoods/20kbp_flanks_Eggnog/out.emapper.decorated.gff")

# Import cluster info
clus <- read_tsv("outputs/neighborhoods/dtnpb.cluster.info.tsv")

# Import flank info
df.flank <- read_tsv("outputs/neighborhoods/dtnp.genomic.20kbp_flanks.info.tsv")

# Import flank fasta
seq.flank <- readDNAStringSet("outputs/neighborhoods/dtnp.genomic.20kbp_flanks.fna")

# Import IPG info
ipg <- read_tsv("outputs/neighborhoods/dtnpb.ipg.tsv")
colnames(ipg)[colnames(ipg) == "Nucleotide Accession"] <- "Genome"

# Create directory for blast comparisons of flanks
dir.create("outputs/neighborhoods/flank_blast_comparisons")

# Create directory for cluster FASTA files
dir.create("outputs/neighborhoods/flank_fasta")

# Create directory for cluster FASTA annotation files
dir.create("outputs/neighborhoods/flank_annotations")

# Initiate function to export flanks and annotations (GFF files) for each cluster
export.cluster.annotations <- function(rep.id){
  
  # Get other cluster member ids
  clus.ids <- clus$id[clus$Group_Rep == rep.id]
  
  # Get IPG info, removing genomes that couldn't be retrieved
  clus.ipg <- ipg %>% filter(Protein %in% clus.ids, Genome %in% word(df.flank$flank.id, sep = "\\_"))
  
  # Get flank info
  clus.flank <- df.flank %>% filter(word(flank.id, sep = "\\_") %in% clus.ipg$Genome) 
  
  # Get flank fasta
  clus.seq.flank <- seq.flank[names(seq.flank) %in% clus.flank$flank.id]
  
  # Get Eggnog annotation info
  clus.egg <- egg %>% filter(seqid %in% clus.flank$flank.id)
  
  # If there is more than one sequence in the cluster, BLAST flanks to find regions of homology 
  if(nrow(clus.flank) > 1) {
    
    # Check to see if a BLAST folder already exists for this cluster
    if(!rep.id %in% list.files("~/Documents/Projects/TnpB/dTnpB/outputs/neighborhoods/flank_blast_comparisons")){
      
      # Create new subdirectory for BLAST files
      dir.create(paste0("outputs/neighborhoods/flank_blast_comparisons/", rep.id))
      
      # Set up variables you'll need for BLAST search
      colz <- c("qseqid", "sseqid", "qlen", "qstart", "qend", "sstart", "send", "sstrand",
                "evalue", "length", "pident", "nident", "mismatch", "gaps")
      df <- data.frame(stringsAsFactors = FALSE)
      
      # Iterate through each sequence
      for(i in 1:length(clus.seq.flank)){
        
        # Pull sequences to assign query/subject
        qry <- clus.seq.flank[i]
        subj <- clus.seq.flank[-i]
        
        # Make temp variables for the names of query/subject
        qry.name <- paste0(names(clus.seq.flank[i]),".query.fna")
        subj.name <- paste0(names(clus.seq.flank[i]),".subject.fna")
        
        # Export temporary FASTA
        writeXStringSet(qry, paste0("outputs/neighborhoods/flank_blast_comparisons/", rep.id,"/", qry.name))
        writeXStringSet(subj, paste0("outputs/neighborhoods/flank_blast_comparisons/", rep.id,"/", subj.name))
        
        # Build BLAST db
        system(paste0("makeblastdb -in outputs/neighborhoods/flank_blast_comparisons/", rep.id,"/", subj.name, " -dbtype nucl"))
        
        # Run BLAST search
        system(paste0("blastn -query outputs/neighborhoods/flank_blast_comparisons/", rep.id,"/", qry.name,
                      " -db outputs/neighborhoods/flank_blast_comparisons/", rep.id,"/", subj.name,
                      " -outfmt '6 qseqid sseqid qlen qstart qend sstart send sstrand evalue length pident nident mismatch gaps' > ",
                      "outputs/neighborhoods/flank_blast_comparisons/", rep.id,"/", names(clus.seq.flank[i]),".out.txt"))
        
        # Import BLAST results
        blast <- read_tsv(paste0("outputs/neighborhoods/flank_blast_comparisons/", rep.id,"/", names(clus.seq.flank[i]),".out.txt"),
                          col_names = colz)
        
        # Add blast results to dataframe
        df <- rbind(df, blast)
        
      } 
      
      # Export all BLAST results
      write_tsv(df, paste0("outputs/neighborhoods/flank_blast_comparisons/", rep.id,".out.tsv"))
      
    } }
  
  # Change annotation type to be COG category
  clus.egg$type <- word(word(clus.egg$attribute, start = 2, sep = "em_COG_cat="), sep = "\\;")
  clus.egg$type[is.na(clus.egg$type)] <- "CDS"
  
  # Add dTnpB as new annotations to Eggnog output
  clus.egg <- rbind(clus.egg, tibble(seqid = paste0(clus.flank$flank.id),
                                     source = rep("IPG", nrow(clus.flank)),
                                     type = rep("dTnpB", nrow(clus.flank)),
                                     start = clus.flank$dtnpb.start - clus.flank$flank.start,
                                     end = clus.flank$dtnpb.end - clus.flank$flank.start,
                                     score = rep(".", nrow(clus.flank)),
                                     strand = rep("+", nrow(clus.flank)),
                                     phase = rep(".", nrow(clus.flank)),
                                     attribute = rep("dTnpB_from_chances_tree", nrow(clus.flank)),
  ))
  
  # Export FASTA of sequences
  writeXStringSet(clus.seq.flank, paste0("outputs/neighborhoods/flank_fasta/", rep.id, ".flanks.fna"))
  
  # Export GFF annotation files
  write_tsv(clus.egg, paste0("outputs/neighborhoods/flank_annotations/", rep.id, ".flanks.annot.gff"), col_names = FALSE)
  
}

# Make list of representatives in tree
x <- unique(clus$Group_Rep)
export.cluster.annotations(rep.id = x[1])

for(i in x) {export.cluster.annotations(rep.id = i)}



############################################
# Determine which homologs in the TnpB tree are dTnpB: Strategy 2 [KEEP sequences with gaps at expected DED positions]
############################################

# Start clean
rm(list = ls())

# Try restricting attention to TnpB sequences with MULTIPLE mutated catalytic residues: INCLUDING GAPS

# Import trimmed alignment
algn <- readAAStringSet("inputs/ALL_tnpB_sequences_250aa_GB_50_clustered_mafft_alignment_ENISI_4_iterations_trimmed_g05.afa")

# Extract columns corresponding to each catalytic RuvC residue
df <- data.frame(id = names(algn), 
                 D191 = unname(as.character(narrow(algn, start = 321, end = 321))),
                 E278 = unname(as.character(narrow(algn, start = 434, end = 434))),
                 D361 = unname(as.character(narrow(algn, start = 548, end = 548))),
                 stringsAsFactors = FALSE)

# Add column that indicates whether or not a full catalytic triad is present (conservative substitutions allowed)
sim <- c("D", "E")
df$full.triad <- FALSE ; df$sum.triad <- 0
for(i in 1:nrow(df)) { 
  if(sum(df$D191[i] %in% sim,  df$E278[i] %in% sim, df$D361[i] %in% sim) == 3) {df$full.triad[i] <- TRUE}
  df$sum.triad[i] <- sum(df$D191[i] %in% sim,  df$E278[i] %in% sim, df$D361[i] %in% sim) }
df.all <- df # Make backup

# Pull sequences with one or less catalytic residues conserved (at least as I can tell in the MSA used to build the TnpB tree) 
these <- df %>% filter(sum.triad <= 1)  

# Add sequences that have a mutation in either of the first two catalytic acidic residues in the triad (D191, E278)
# these <- rbind(these, df %>% filter(sum.triad == 2, (!D191 %in% sim | !D361 %in% sim)))

# Import TnpB tree
tree <- ape::read.tree("inputs/TnpB.tre")

# Import list of sequences that I previously looked at
not.these <- read_tsv("outputs/Catalytic_triad_completeness.info.tsv")
not.these <- not.these %>% filter(full.triad == FALSE)
these <- these %>% filter(!id %in% not.these$id)

# Make list for OTU grouping (i.e., colored branches)
triad <- list()
triad[["Incomplete, already analyzed"]] <- not.these$id
triad[["Incomplete"]] <- these$id
triad[["Complete"]] <- df.all$id[!df.all$id %in% these$id & !df.all$id %in% not.these$id]

# Visualize tree with colored branches for catalytic triad completeness (RuvC)
df.lab <- cbind(df %>% select(label = id), df %>% select(id))
p <- ggtree(tree, layout = "circular") %<+% df.lab
groupOTU(p, triad, "Catalytic_triad") + aes(color = Catalytic_triad) + theme(legend.position = "right") +
  geom_tiplab2(label = id, align = TRUE)

# Change values of catalytic triad: TRUE/FALSE
df.all$full.triad <- TRUE
df.all$full.triad[df.all$id %in% these$id] <- FALSE

# Reorder triad info to be consistent with branch order
x <- data.frame(id = get_taxa_name(p))
triad.df <- x %>% left_join(df.all, by = "id")

# Add position in tree to dataframe
triad.df <- triad.df %>% mutate(tree.position = 1:nrow(triad.df))

# Import cluster info
clus <- read_tsv("inputs/ALL_tnpB_sequences_250aa_55_clustered.df.tsv")

# Add cluster info to each row in the triad dataframe
triad.df <- triad.df %>% mutate(cluster.size = NA, mean.identity = NA, mean.seq.cov = NA, mean.seq.length = NA)
for(i in 1:nrow(triad.df)){
  
  # Extract only the cluster info for the representative
  x <- clus %>% filter(Group_Rep == clus$Group_Rep[clus$id == triad.df$id[i]])
  
  # Convert cluster identity and coverage from character vectors to numeric
  x$clstr_iden <- as.numeric(gsub("%", "", x$clstr_iden))
  x$clstr_cov <- as.numeric(gsub("%", "", x$clstr_cov))
  
  # Fill in each field
  triad.df$cluster.size[i] <- nrow(x)
  triad.df$mean.identity[i] <- paste0(round(mean(x$clstr_iden),2), "% [", min(x$clstr_iden), "-", max(x$clstr_iden), "%]")
  triad.df$mean.seq.cov[i] <- paste0(round(mean(x$clstr_cov),2), "% [", min(x$clstr_cov), "-", max(x$clstr_cov), "%]")
  triad.df$mean.seq.length[i] <- paste0(round(mean(x$length),2), " [", min(x$length), "-", max(x$length), "]")
  
}

# Export
write_tsv(triad.df, "outputs/Catalytic_triad_completeness.info2.tsv")


############################################
# Pull genomes for each dTnpB that has more than a single representative per cluster: Strategy 2
############################################

# Start clean
rm(list = ls())

# Make new directory
dir.create("outputs/neighborhoods2/", showWarnings = FALSE)

# Import cluster info
clus <- read_tsv("inputs/ALL_tnpB_sequences_250aa_55_clustered.df.tsv")

# Import triad info
df <- read_tsv("outputs/Catalytic_triad_completeness.info2.tsv")

# Narrow cluster info to include only clusters with incomplete triads
x <- df$id[df$full.triad == FALSE] # & df$cluster.size > 1] # turn this on to exclude singletons
y <- clus$Group_Rep[clus$id %in% x]
clus <- clus %>% filter(Group_Rep %in% y)

# Export info on clusters with incomplete triads: break into groups of 800 sequences, so that Batch-Entrez can handle it
write_tsv(clus, "outputs/neighborhoods2/dtnpb.cluster.info.tsv")

# Export protein accessions
x <- 1 ; y <- 1
while(x < nrow(clus)) { write_tsv(clus[x:(x+799),1], paste0("outputs/neighborhoods2/dtnpb.cluster.acc.batch",y,".txt"), col_names = FALSE)
  x <- x +800 ; y <- y +1 }

# Batch-Entrez has not been working well. Let's try E-utilities instead to pull IPG info for each sequence
# epost -db protein -input dtnpb.cluster.acc.batch1.txt | efetch -format ipg -mode text > dtnpb.cluster.batch1.ipg.txt
# epost -db protein -input dtnpb.cluster.acc.batch2.txt | efetch -format ipg -mode text > dtnpb.cluster.batch2.ipg.txt
# epost -db protein -input dtnpb.cluster.acc.batch3.txt | efetch -format ipg -mode text > dtnpb.cluster.batch3.ipg.txt
# epost -db protein -input dtnpb.cluster.acc.batch4.txt | efetch -format ipg -mode text > dtnpb.cluster.batch4.ipg.txt

# Import IPG info
ipg <- rbind(read_tsv("outputs/neighborhoods2/dtnpb.cluster.batch1.ipg.txt"), read_tsv("outputs/neighborhoods2/dtnpb.cluster.batch2.ipg.txt"),
             read_tsv("outputs/neighborhoods2/dtnpb.cluster.batch3.ipg.txt"), read_tsv("outputs/neighborhoods2/dtnpb.cluster.batch4.ipg.txt"))

# Insane amount of sequences! See if they're all coming from a well-sequenced, pathogenic genus/species
table(word(ipg$Organism)) # Mycobacterium and other Myco-s.... interesting. potentially.Checked it out: cat triad is intact
length(unique(ipg$Organism))
length(unique(ipg$Id)) # 2,022 unique IDs, so we lost 711 in the pull

#### Okay, let's try to break this up.
# Manually checked the top hits from Myco. There are a tonnnnn of the same cat triad intact. Let's get rid of them
x <- clus$id[clus$Group_Rep == "COW23552.1"]
y <- unique(ipg$Id[ipg$Protein %in% x])

ipg <- ipg %>% filter(!Id %in% y)
clus <- clus %>% filter(!Group_Rep == "COW23552.1")

# Look at the top hits with those removed
x <- data.frame(table(ipg$Id))

# Function to get representative in the tree for a given gene id
check.id <- function(ipg.id) {
  y <- unique(ipg$Protein[ipg$Id == paste0(ipg.id) & ipg$Protein %in% clus$id])
  z <- clus$Group_Rep[clus$id == y[1]]
  message(paste0(clus$id[clus$Group_Rep == z & clus$id %in% df$id])) }

# First: pull sequences that are in the RefSeq
#ref <- ipg %>% filter(Source == "RefSeq")

# Remove RefSeq sequences that are obviously truncated (dTnpB coordinate starts at 1)
#ref <- ref %>% filter(!Start == 1) 

# Check how many ids aren't covered by RefSeq sequences
#length(unique(ipg$Id)) - length(unique(ref$Id)) # 556, we'll try adding them back in

# Pick one-ish representative for each RefSeq sequence
#ref.sum <- ref %>% group_by(Id) %>% filter(Start == max(Start)) %>% ungroup()
#ref <- ref %>% filter(Protein %in% ref.sum$Protein & Assembly %in% ref.sum$Assembly)

# Pick one representative for each of the IPG ids not in a RefSeq genome
#no.ref <- ipg %>% filter(!Source == "RefSeq", !Id %in% ref$Id)
#no.ref.sum <- no.ref %>% group_by(Id) %>% filter(Start == max(Start)) %>% ungroup()
#no.ref <- no.ref %>% filter(Protein %in% no.ref.sum$Protein & Assembly %in% no.ref.sum$Assembly)

# Concatenate accessions for output
#ipg.narrowed <- rbind(ref, no.ref)
ipg.narrowed <- ipg

# Check how many clusters are covered by the resulting list of proteins
x <- clus$id[clus$id %in% ipg.narrowed$Protein]
length(unique(clus$Group_Rep[clus$id %in% x])) # 84 total clusters, missing about 30 [maybe somethine to look at later]

# Keep only unique accessions
ipg.acc <- ipg %>% select(`Nucleotide Accession`) %>% distinct

# Export genomic accessions in batches of 500
write_tsv(ipg.acc[1:500,] %>% select(`Nucleotide Accession`) %>% filter(!`Nucleotide Accession` == "N/A"),
          "outputs/neighborhoods2/dtnpb.genomic.acc.batch1.txt", col_names = FALSE)
write_tsv(ipg.acc[501:1000,] %>% select(`Nucleotide Accession`) %>% filter(!`Nucleotide Accession` == "N/A"),
          "outputs/neighborhoods2/dtnpb.genomic.acc.batch2.txt", col_names = FALSE)
write_tsv(ipg.acc[1001:1500,] %>% select(`Nucleotide Accession`) %>% filter(!`Nucleotide Accession` == "N/A"),
          "outputs/neighborhoods2/dtnpb.genomic.acc.batch3.txt", col_names = FALSE)
write_tsv(ipg.acc[1501:2000,] %>% select(`Nucleotide Accession`) %>% filter(!`Nucleotide Accession` == "N/A"),
          "outputs/neighborhoods2/dtnpb.genomic.acc.batch4.txt", col_names = FALSE)
write_tsv(ipg.acc[2001:2500,] %>% select(`Nucleotide Accession`) %>% filter(!`Nucleotide Accession` == "N/A"),
          "outputs/neighborhoods2/dtnpb.genomic.acc.batch5.txt", col_names = FALSE)
write_tsv(ipg.acc[2501:3000,] %>% select(`Nucleotide Accession`) %>% filter(!`Nucleotide Accession` == "N/A"),
          "outputs/neighborhoods2/dtnpb.genomic.acc.batch6.txt", col_names = FALSE)
write_tsv(ipg.acc[3001:3500,] %>% select(`Nucleotide Accession`) %>% filter(!`Nucleotide Accession` == "N/A"),
          "outputs/neighborhoods2/dtnpb.genomic.acc.batch7.txt", col_names = FALSE)
write_tsv(ipg.acc[3501:3985,] %>% select(`Nucleotide Accession`) %>% filter(!`Nucleotide Accession` == "N/A"),
          "outputs/neighborhoods2/dtnpb.genomic.acc.batch8.txt", col_names = FALSE)


# Use Batch-Entrez to pull IPG info for each sequence
# outputs/neighborhoods2/dtnpb.genomic1.fna to dtnpb.genomic8.fna

# Iterate through genomes to pull flanks
flanks <- DNAStringSet()
flank.size <- 20000
df.flank <- tibble()
for(z in 1:8){
  
  # Import genomic loci
  genomes <- readDNAStringSet(paste0("outputs/neighborhoods2/dtnpb.genomic",z,".fna"))
  
  # Pull flanks of dTnpB
  
  for(i in 1:length(genomes)){
    
    # Isolate one genome
    x <- genomes[i]
    
    # Get annotations for that genome
    annot <- ipg.narrowed %>% filter(`Nucleotide Accession` == word(names(x)))
    annot$Start <- as.numeric(annot$Start) ; annot$Stop <- as.numeric(annot$Stop)
    
    # Some genomes have 2 copies
    for(v in 1:nrow(annot)){
      
      # # Trim genome to only include dTnpB + max length of sequence (if 3kbp is longer than end) 
      if(annot$Stop[v] + flank.size > width(x)){ flank.end <- width(x) } else { flank.end <- annot$Stop[v] + flank.size }
      if(annot$Start[v] - flank.size < 1){ flank.start <- 1 } else { flank.start <- annot$Start[v] - flank.size }
      
      # Trim genome to only include dTnpB + 3kbp flanks
      this.flank <- narrow(x, start = flank.start, end = flank.end)
      
      # Add coordinates to FASTA header
      names(this.flank) <- paste0(word(names(this.flank)), "_", flank.start, "-", end = flank.end)
      
      # Flip sequence, if the dTnpB is on the "-" strand
      if(annot$Strand[v] == "-") { this.flank <- reverseComplement(this.flank) }
      
      # Add flank to new FASTA
      flanks <- c(flanks, this.flank) 
      
      # Add flank info to dataframe
      df.flank <- rbind(df.flank, tibble(genome = names(x), flank.id = paste0(word(names(this.flank))),
                                         length = width(x), flank.start = flank.start, flank.end = flank.end, 
                                         dtnpb.start = annot$Start[v], dtnpb.end = annot$Stop[v]))
      
    }
    
    # Start fresh, to raise errors if something goes wrong with the next genome
    rm(x, annot, this.flank, flank.start, flank.end)
    
  }
}
# Determine flank sizes
df.flank <- df.flank %>% mutate(flank.size = flank.end - flank.start + 1) %>% distinct()

# Remove duplicates
flanks <- unique(flanks)
df.flank <- df.flank %>% filter(flank.id %in% names(flanks))

# Export flanks and info
writeXStringSet(flanks, "outputs/neighborhoods2/dtnpb.genomic.20kbp_flanks.fna")
write_tsv(df.flank, "outputs/neighborhoods2/dtnpb.genomic.20kbp_flanks.info.tsv")


############################################
# Annotate flanks with Eggnogg and visualize: Strategy 2
############################################  

# Start clean
rm(list = ls())

# Run Eggnog mapper
# emapper.py --cpu 28 -o out \
# --output_dir 20kbp_flanks_Eggnog3 \
# --override -m diamond --dmnd_ignore_warnings \
# -i dtnpb.genomic.20kbp_flanks.fna \
# --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 \
# --itype genome --genepred prodigal --tax_scope auto --target_orthologs all \
# --go_evidence non-electronic --pfam_realign none --report_orthologs \
# --decorate_gff yes \
# --excel

# Import GFF file of eggnog output
egg <- read_gff("outputs/neighborhoods2/20kbp_flanks_Eggnog3/out.emapper.decorated.gff")

# Import cluster info
clus <- read_tsv("outputs/neighborhoods2/dtnpb.cluster.info.tsv")

# Import flank info
df.flank <- read_tsv("outputs/neighborhoods2/dtnpb.genomic.20kbp_flanks.info.tsv")

# Import flank fasta
seq.flank <- readDNAStringSet("outputs/neighborhoods2/dtnpb.genomic.20kbp_flanks.fna")

# Import IPG info
ipg <- rbind(read_tsv("outputs/neighborhoods2/dtnpb.cluster.batch1.ipg.txt"), read_tsv("outputs/neighborhoods2/dtnpb.cluster.batch2.ipg.txt"),
             read_tsv("outputs/neighborhoods2/dtnpb.cluster.batch3.ipg.txt"), read_tsv("outputs/neighborhoods2/dtnpb.cluster.batch4.ipg.txt"))

colnames(ipg)[colnames(ipg) == "Nucleotide Accession"] <- "Genome"

# Create directory for blast comparisons of flanks
dir.create("outputs/neighborhoods2/flank_blast_comparisons")

# Create directory for cluster FASTA files
dir.create("outputs/neighborhoods2/flank_fasta")

# Create directory for cluster FASTA annotation files
dir.create("outputs/neighborhoods2/flank_annotations")

# Remove sequences from ipg output that aren't included in the genomic flanks
ipg <- ipg %>% filter(Genome %in% word(df.flank$genome))

# Initiate function to export flanks and annotations (GFF files) for each cluster
export.cluster.annotations <- function(rep.id){
  
  # Get other cluster member ids
  clus.ids <- clus$id[clus$Group_Rep == rep.id]
  
  # Get IPG info, removing genomes that couldn't be retrieved
  clus.ipg <- ipg %>% filter(Protein %in% clus.ids, Genome %in% word(df.flank$flank.id, sep = "\\_"))
  
  # Get flank info
  clus.flank <- df.flank %>% filter(word(flank.id, sep = "\\_") %in% clus.ipg$Genome) 
  
  # Get flank fasta
  clus.seq.flank <- seq.flank[names(seq.flank) %in% clus.flank$flank.id]
  
  # Get Eggnog annotation info
  clus.egg <- egg %>% filter(seqid %in% clus.flank$flank.id)
  
  # If there is more than one sequence in the cluster, BLAST flanks to find regions of homology 
  if(nrow(clus.flank) > 1) {
    
    # Check to see if a BLAST folder already exists for this cluster
    if(!rep.id %in% list.files("~/Documents/Projects/TnpB/dTnpB/outputs/neighborhoods2/flank_blast_comparisons")){
      
      # Create new subdirectory for BLAST files
      dir.create(paste0("outputs/neighborhoods2/flank_blast_comparisons/", rep.id))
      
      # Set up variables you'll need for BLAST search
      colz <- c("qseqid", "sseqid", "qlen", "qstart", "qend", "sstart", "send", "sstrand",
                "evalue", "length", "pident", "nident", "mismatch", "gaps")
      df <- data.frame(stringsAsFactors = FALSE)
      
      # Iterate through each sequence
      for(i in 1:length(clus.seq.flank)){
        
        # Pull sequences to assign query/subject
        qry <- clus.seq.flank[i]
        subj <- clus.seq.flank[-i]
        
        # Make temp variables for the names of query/subject
        qry.name <- paste0(names(clus.seq.flank[i]),".query.fna")
        subj.name <- paste0(names(clus.seq.flank[i]),".subject.fna")
        
        # Export temporary FASTA
        writeXStringSet(qry, paste0("outputs/neighborhoods2/flank_blast_comparisons/", rep.id,"/", qry.name))
        writeXStringSet(subj, paste0("outputs/neighborhoods2/flank_blast_comparisons/", rep.id,"/", subj.name))
        
        # Build BLAST db
        system(paste0("makeblastdb -in outputs/neighborhoods2/flank_blast_comparisons/", rep.id,"/", subj.name, " -dbtype nucl"))
        
        # Run BLAST search
        system(paste0("blastn -query outputs/neighborhoods2/flank_blast_comparisons/", rep.id,"/", qry.name,
                      " -db outputs/neighborhoods2/flank_blast_comparisons/", rep.id,"/", subj.name,
                      " -outfmt '6 qseqid sseqid qlen qstart qend sstart send sstrand evalue length pident nident mismatch gaps' > ",
                      "outputs/neighborhoods2/flank_blast_comparisons/", rep.id,"/", names(clus.seq.flank[i]),".out.txt"))
        
        # Import BLAST results
        blast <- read_tsv(paste0("outputs/neighborhoods2/flank_blast_comparisons/", rep.id,"/", names(clus.seq.flank[i]),".out.txt"),
                          col_names = colz)
        
        # Add blast results to dataframe
        df <- rbind(df, blast)
        
      } 
      
      # Export all BLAST results
      write_tsv(df, paste0("outputs/neighborhoods2/flank_blast_comparisons/", rep.id,".out.tsv"))
      
    } }
  
  # Change annotation type to be COG category
  clus.egg$type <- word(word(clus.egg$attribute, start = 2, sep = "em_COG_cat="), sep = "\\;")
  clus.egg$type[is.na(clus.egg$type)] <- "CDS"
  
  # Add dTnpB as new annotations to Eggnog output
  clus.egg <- rbind(clus.egg, tibble(seqid = paste0(clus.flank$flank.id),
                                     source = rep("IPG", nrow(clus.flank)),
                                     type = rep("dTnpB", nrow(clus.flank)),
                                     start = clus.flank$dtnpb.start - clus.flank$flank.start,
                                     end = clus.flank$dtnpb.end - clus.flank$flank.start,
                                     score = rep(".", nrow(clus.flank)),
                                     strand = rep("+", nrow(clus.flank)),
                                     phase = rep(".", nrow(clus.flank)),
                                     attribute = rep("dTnpB_from_chances_tree", nrow(clus.flank)),
  ))
  
  # Export FASTA of sequences
  writeXStringSet(clus.seq.flank, paste0("outputs/neighborhoods2/flank_fasta/", rep.id, ".flanks.fna"))
  
  # Export GFF annotation files
  write_tsv(clus.egg, paste0("outputs/neighborhoods2/flank_annotations/", rep.id, ".flanks.annot.gff"), col_names = FALSE)
  
}

# Make list of representatives in tree
x <- unique(clus$Group_Rep)
#export.cluster.annotations(rep.id = x[1]) # Test

for(i in x) {export.cluster.annotations(rep.id = i)}


#####

##################################################
# Function to pull homologs of sequences in TnpB tree
##################################################

library(rentrez)
library(DECIPHER)
rm(list = ls())

dir.create("outputs/neighborhoods3", showWarnings = FALSE)
dir.create("outputs/neighborhoods3/automatic_annotations", showWarnings = FALSE) # This is where the results will be saved
dir.create("outputs/neighborhoods3/automatic_annotations/cluster_alignments", showWarnings = FALSE)

# Import Catalytic triad info from search with strategy three (<1 conserved residue, counting gaps)
triad3 <- read_tsv("outputs/Catalytic_triad_completeness.info2.tsv")

# Import Catalytic triad info from search with strategy three (<1 conserved residue, NOT counting gaps)
triad2 <- read_tsv("outputs/Catalytic_triad_completeness.info.tsv")

# Narrow to only sequences without an intact triad
triad.these <- rbind(triad2 %>% filter(full.triad == FALSE), triad3 %>% filter(full.triad == FALSE)) %>% distinct()

## Pull sequences for all cluster members of representatives with incomplete triads
# Import cluster info 
clus <- read_tsv("inputs/ALL_tnpB_sequences_250aa_55_clustered.df.tsv")

# Determine which group representative each cat dead sequence in the tree is associated with
clus.rep <- clus$Group_Rep[clus$id %in% triad.these$id]

# Fetch all ids in cluster with a cat dead representative
clus <- clus %>% filter(Group_Rep %in% clus.rep)

# Make directory for batch entrez files: accessions
dir.create("outputs/neighborhoods3/cluster_member_accessions", showWarnings = FALSE)
dir.create("outputs/neighborhoods3/cluster_member_fasta", showWarnings = FALSE)

# Export accessions [Run only once]
write_tsv((clus %>% select(id)), "outputs/neighborhoods3/cluster_member_accessions/cluster_members.acc.txt", col_names = FALSE)
 
# Try pulling from NR database
# blastdbcmd -entry_batch outputs/neighborhoods3/cluster_member_accessions/cluster_members.acc.txt \ 
# -db /Volumes/SSD/Databases/ncbi_nr/nr -outfmt %f > outputs/neighborhoods3/cluster_member_fasta/cluster_members.faa"))

# Set working directory to where results will be saved
setwd("outputs/neighborhoods3/automatic_annotations")

# Initiate function
pull.cluster.members <- function(protein.accession){
  
  # Import cluster info 
  clus <- read_tsv("../../../inputs/ALL_tnpB_sequences_250aa_55_clustered.df.tsv", show_col_types = FALSE)

  # Narrow to only the cluster of the query
  x <- clus$Group_Rep[clus$id == protein.accession]
  clus <- clus %>% filter(Group_Rep == x)
  
  ### [New approach]
  # Import FASTA of cluster members  [Manually changed two names that were too long]
  cluster.members <- readAAStringSet("../cluster_member_fasta/cluster_members.faa") 
  
  # Keep only accession as sequence header
  names(cluster.members) <- word(names(cluster.members))
  
  # Narrow to only include members in the cluster being considered
  cluster.members <- cluster.members[names(cluster.members) %in% clus$id]
  
  if(length(cluster.members) > 0) {
  
  # Print status message
  message(paste0("Sequences retrieved for ", protein.accession, " cluster..."))
  
  #### #### #### #### #### ###
  # Verify cluster is dTnpB
  #### #### #### #### #### ###
  
  # Import reference sequences
  ref <- c(readAAStringSet("../../../inputs/refs/ISDra2_TnpB.faa"),
           readAAStringSet("../../../inputs/refs/7C7L_Cas12f.faa"))
  
  # Align cluster members to reference to make sure that they are actually catalytically dead
  cluster.algn <- AlignSeqs(c(ref, cluster.members), verbose = FALSE)
  
  # Save alignment
  writeXStringSet(cluster.algn, paste0("cluster_alignments/", protein.accession, ".cluster.afa"))
  
  # Build dataframes of positions for reference sequences
  ref.tnpb.df <- tibble(res = str_split(unname(as.character(cluster.algn[1])), "")[[1]], pos.algn = 1:width(cluster.algn[1])) %>% 
    filter(!res == "-") %>% mutate(pos = 1:width(ref[1]))
  ref.cas12f.df <- tibble(res = str_split(unname(as.character(cluster.algn[2])), "")[[1]], pos.algn = 1:width(cluster.algn[2])) %>% 
    filter(!res == "-") %>% mutate(pos = 1:width(ref[2]))
  
  # Store catalytic residues
  cat.res <- c("D", "E")
  
  ## Build dataframe of residues in cluster members that correspond to catalytic residues in each reference sequence
  # Check ISDra2-like catalytic triad
  cluster.tnpb.df <- tibble(id = names(cluster.algn[3:length(cluster.algn)]), 
                            res1 = unname(as.character(narrow(cluster.algn[3:length(cluster.algn)],
                                                              start = ref.tnpb.df$pos.algn[ref.tnpb.df$pos == 191], width = 1))) %in% cat.res,
                            res2 = unname(as.character(narrow(cluster.algn[3:length(cluster.algn)],
                                                              start = ref.tnpb.df$pos.algn[ref.tnpb.df$pos == 278], width = 1))) %in% cat.res,
                            res3 = unname(as.character(narrow(cluster.algn[3:length(cluster.algn)],
                                                              start = ref.tnpb.df$pos.algn[ref.tnpb.df$pos == 361], width = 1))) %in% cat.res) %>% 
    mutate(sum.triad = res1 + res2 + res3)
  
  # Check Cas12f-like catalytic triad
  cluster.cas12f.df <- tibble(id = names(cluster.algn[3:length(cluster.algn)]), 
                              res1 = unname(as.character(narrow(cluster.algn[3:length(cluster.algn)],
                                                                start = ref.cas12f.df$pos.algn[ref.cas12f.df$pos == 326], width = 1))) %in% cat.res,
                              res2 = unname(as.character(narrow(cluster.algn[3:length(cluster.algn)],
                                                                start = ref.cas12f.df$pos.algn[ref.cas12f.df$pos == 422], width = 1))) %in% cat.res,
                              res3 = unname(as.character(narrow(cluster.algn[3:length(cluster.algn)],
                                                                start = ref.cas12f.df$pos.algn[ref.cas12f.df$pos == 510], width = 1))) %in% cat.res) %>% 
    mutate(sum.triad = res1 + res2 + res3)
  
  # Cross-reference each of the dataframes with catalytic residue info
  which.best.ref <- tibble(ref.id = c("ISDra2", "Cas12f"),
                           max.triad = c(sum(cluster.tnpb.df$sum.triad),
                                         sum(cluster.cas12f.df$sum.triad)))
  
  # Choose best reference sequence to use to assess catalytic residue conservation
  ref.df <- ref.tnpb.df ; ref.df$ref <- "ISDra2" # Automatically use ISDra2 as reference, unless another ref scores higher
  if(length(which.best.ref$ref.id[which.best.ref$max.triad == max(which.best.ref$max.triad)]) == 1) {
    if(which.best.ref$ref.id[which.best.ref$max.triad == max(which.best.ref$max.triad)] == "Cas12f") { ref.df <- ref.cas12f.df
    ref.df$ref <- "Cas12f" }}
  
  # Reference chosen, keep triad info for that reference
  if(ref.df$ref[1] == "ISDra2") { cluster.df <- cluster.tnpb.df ; cluster.df$ref <- "ISDra2" }
  if(ref.df$ref[1] == "Cas12f") { cluster.df <- cluster.cas12f.df ; cluster.df$ref <- "Cas12f" }

  # Determine if sequences are actually dTnpB
  if(sum(cluster.df$sum.triad > 1) > 0) { 
    
    # Print message if any sequences are likely catalytically active TnpB
    message(paste0(sum(cluster.df$sum.triad > 1), " of ", nrow(clus)," sequences in same cluster as ", protein.accession,
                   " have more catalytic residues than expected."))
    message("")
    message(paste0("These sequences will be saved to 'unlikely_dTnpB.faa' and their info will be saved to 'unlikely_dTnpB.tsv'"))
    
    # Save sequences to FASTA
    writeXStringSet(cluster.members[names(cluster.members) %in% cluster.df$id[cluster.df$sum.triad > 1]],
                    "unlikely_dTnpB.faa", append = TRUE)
    
    # Remove these sequences from FASTA so they are not processed further
    cluster.members <- cluster.members[!names(cluster.members) %in% cluster.df$id[cluster.df$sum.triad > 1]]
    
    # Save sequences to TSV file
    write_tsv(cluster.df %>% filter(sum.triad > 1), "unlikely_dTnpB.tsv", append = TRUE)
    
    # Remove these sequences from cluster member info so they are not processed further
    clus <- clus %>% filter(!id %in% cluster.df$id[cluster.df$sum.triad > 1])
    
    # Print message if that was the whole cluster
    if(length(cluster.members) == 0) { message(paste0("All cluster members for ", protein.accession, " were removed.")) }
    
    
  } 
  
  if(length(cluster.members) > 0) {
    
    # Print status message
    message(paste0("Catalytic triad checking complete for ", protein.accession, " cluster..."))
    
    # Add column to output indicating whether or not a cluster member was fetchable from NR
    clus <- clus %>% mutate(in.NR = FALSE)
    clus$in.NR[clus$id %in% names(cluster.members)] <- TRUE
    
    
    # Store sequences for BLAST
    writeXStringSet(cluster.members, "Likely_dTnpB.faa", append = TRUE)
    
    # Store sequence info
    write_tsv(clus, "Likely_dTnpB.cluster.info.tsv", append = TRUE)
    
  } 
  }
  message("")
  message("-----------------------------------------")
  message("")
}

# Pull members for each cluster representative
for(i in 1:nrow(triad.these)) { pull.cluster.members(protein.accession = paste0(triad.these$id[i])) }


##################################################
# BLAST dTnpB proteins to try to find homologs that were missed during the Jackhmmer search
##################################################

##### BLAST sequences #####

rm(list = ls())  # Start clean

# Make directory for results
dir.create("blast_results/",  showWarnings = FALSE)

# BLAST sequences likely dTnpB proteins against NR database
# blastp -query likely_dTnpB.faa -evalue 1e-50 -qcov_hsp_perc 80 -max_hsps 50 -num_threads 30 \
# -outfmt '6 qseqid sseqid qlen qstart qend sstart send sstrand evalue length pident nident mismatch gaps qcovs staxid' \
# -max_target_seqs 50 -db ~/databases/ncbi_nr/nr > blast_results/likely_dTnpB.blastout.v1.txt

# Import List of likely candidates
colz <- c("id", "clstr", "clstr_size", "length", "clstr_rep", "clstr_iden", "clstr_cov", "Group_Rep", "in.NR")
cand <- read_tsv("Likely_dTnpB.cluster.info.tsv", col_names = colz) %>% distinct()

# See how many we were unable to pull from NR database
table(cand$in.NR)

# Try exporting list of accessions and pulling again
write_tsv(cand %>% filter(in.NR ==FALSE) %>% select(id), "blast_results/likely_dTnpB_not_in_NR.acc.txt", col_names = FALSE)

# Try pulling again, since I can submit some entries and still get the primary sequence in non-batch entries from NR
  # blastdbcmd -entry_batch likely_dTnpB_not_in_NR.acc.txt -db /Volumes/SSD/Databases/ncbi_nr/nr -outfmt %f > likely_dTnpB_pull2.faa

# Import candidate proteins
prot <- c(readAAStringSet("Likely_dTnpB.faa"),
          readAAStringSet("likely_dTnpB_pull2.faa"))

# Fix names of second group of sequences added
for(i in 1:length(prot)) { 
  x <- gsub(">","",str_split(names(prot)[i], " ")[[1]])
  if(length(x[x %in% cand$id]) > 0) {
  names(prot)[i] <- x[x %in% cand$id] } else { prot <- prot[-i]} }

# Remove one bad one that made it through
prot <- prot[-grep("WP_007305335.1", names(prot))]

### Determine tree representative for each cluster to assess how many critical proteins are still missing
# Import cluster info
clus.all <- read_tsv("../../../inputs/ALL_tnpB_sequences_250aa_55_clustered.df.tsv")

# Import list of sequences in tree
in.tree.all <- tibble(id = 
    word(names(readAAStringSet("../../first_look/Catalytic_triad_completeness.firstlook.orderedMSA.trimmed.afa"))))

# Add Group_Rep column to list of sequences in tree
in.tree.all <- in.tree.all %>% left_join(clus.all, by = "id")

# Add tree representative column to cluster list
clus.all <- clus.all %>% left_join(in.tree.all %>% select(Group_Rep, tree_rep = id), by = "Group_Rep") %>%
  select(id, clstr, clstr.size = clstr_size, length, clstr.rep = clstr_rep, clstr.iden = clstr_iden, 
         clstr.cov = clstr_cov, group.rep = Group_Rep, tree.rep = tree_rep)
    # Added 4,718 rows to accommodate clusters with more than one representative in the tree

# Save cluster info with tree representatives
write_tsv(clus.all, "cluster.info.tsv")

# Add tree representative column to likely dTnpB candidate list
cand <- cand %>% left_join(clus.all %>% select(id, tree.rep), by = "id") %>% 
  select(id, clstr, clstr.size = clstr_size, length, clstr.rep = clstr_rep, clstr.iden = clstr_iden, 
         clstr.cov = clstr_cov, group.rep = Group_Rep, in.NR, tree.rep)
    
# Update NR column
cand$in.NR[cand$id %in% names(prot)] <- TRUE

# Check to see how many proteins that are still missing are tree representatives
table(cand$id[cand$in.NR == FALSE] == cand$tree.rep[cand$in.NR == FALSE])

# Pull missing tree representative sequences from untrimmed alignment
msa <- readAAStringSet("../../first_look/Catalytic_triad_completeness.firstlook.orderedMSA.trimmed.afa")
x <- msa[names(msa) %in% unique(cand$id[cand$in.NR == FALSE])] ; rm(msa)
y <- AAStringSet(gsub("-","", unname(as.character(x))))
names(y) <- names(x)

# Add the five tree representative sequences to the growing protein list
prot <- c(prot, y)

# Add the tree representatives to the protein FASTA that will be used for another round of BLAST
# writeXStringSet(y, "likely_dTnpB_pull2.faa", append = TRUE) ## Run only once

# Update NR column
cand$in.NR[cand$id %in% names(prot)] <- TRUE

# Determine how many missing sequences are the only one from their cluster: it shouldn't be any
table(!cand$clstr[cand$in.NR == FALSE] %in% cand$clstr[cand$in.NR == TRUE])
  # None: good, so every cluster is represented. We'll pick more up in the BLAST search step and shouldn't worry too much about this

# Determine how many we're still missing
length(unique(cand$id[cand$in.NR == FALSE])) # n = 247

# Export updated likely dTnpB candidate info
write_tsv(cand, "likely_dtnpb_candidates.intreeorcluster.tsv")

# Export final list of all queries that were used in the two blast searches
writeXStringSet(prot, "blast_results/likely_dtnpb_candidates.final_blast_queries.faa")

# BLAST sequences of likely dTnpB proteins against NR database on Titan: batch #2
# blastp -query likely_dTnpB_pull2.faa -evalue 1e-50 -qcov_hsp_perc 80 -max_hsps 50 -num_threads 30 \
# -outfmt '6 qseqid sseqid qlen qstart qend sstart send sstrand evalue length pident nident mismatch gaps qcovs staxid' \
# -max_target_seqs 50 -db ~/databases/ncbi_nr/nr > blast_results/likely_dTnpB.blastout.v2.txt


####################
# Import BLAST info
####################

# Start clean
rm(list = ls())

# Import likely dTnpB candidate info
cand <- read_tsv("likely_dtnpb_candidates.intreeorcluster.tsv")

# Import BLAST results
colz <- c("qseqid", "sseqid", "qlen", "qstart", "qend", "sstart", "send", "sstrand",
          "evalue", "length", "pident", "nident", "mismatch", "gaps", "qcov", "staxid")
blast <- rbind(read_tsv("blast_results/likely_dTnpB.blastout.v1.txt", col_names = colz),
               read_tsv("blast_results/likely_dTnpB.blastout.v2.txt", col_names = colz))

# Add column with subject accession only
blast$sacc <- word(blast$sseqid, start = 2, sep = "\\|")

# Import likely dTnpB proteins used as queries
prot <- readAAStringSet("blast_results/likely_dtnpb_candidates.final_blast_queries.faa")

# Determine how many queries had matching BLAST hits
table(names(prot) %in% blast$qseqid)

# Import initial clustering info
clus.all <- read_tsv("cluster.info.tsv")

# Check how many BLAST hits were in initial clusters
table(blast$sacc %in% unique(clus.all$id))

# Make list of unique accessions to pull IPG info
blast.unique <- blast %>% select(sacc) %>% distinct()

# Export unique accessions for IPG info [uncomment this section to run once]
# write_tsv(blast.unique[1:1000,], "blast_results/likely_dtnpb_candidates.blast.acc1.txt", col_names = FALSE)
# write_tsv(blast.unique[1001:2000,], "blast_results/likely_dtnpb_candidates.blast.acc2.txt", col_names = FALSE)
# write_tsv(blast.unique[2001:3000,], "blast_results/likely_dtnpb_candidates.blast.acc3.txt", col_names = FALSE)
# write_tsv(blast.unique[3001:4000,], "blast_results/likely_dtnpb_candidates.blast.acc4.txt", col_names = FALSE)
# write_tsv(blast.unique[4001:5000,], "blast_results/likely_dtnpb_candidates.blast.acc5.txt", col_names = FALSE)
# write_tsv(blast.unique[5001:6000,], "blast_results/likely_dtnpb_candidates.blast.acc6.txt", col_names = FALSE)
# write_tsv(blast.unique[6001:7000,], "blast_results/likely_dtnpb_candidates.blast.acc7.txt", col_names = FALSE)
# write_tsv(blast.unique[7001:8000,], "blast_results/likely_dtnpb_candidates.blast.acc8.txt", col_names = FALSE)
# write_tsv(blast.unique[8001:8895,], "blast_results/likely_dtnpb_candidates.blast.acc9.txt", col_names = FALSE)
# 
# # Export all accessions to one file to pull protein sequences from NR
# write_tsv(blast.unique, "v4/blast_results/likely_dtnpb_candidates.blast.acc.all.txt", col_names = FALSE)

# Fetch sequences
# blastdbcmd -entry_batch likely_dtnpb_candidates.blast.acc.all.txt -db /Volumes/SSD/Databases/ncbi_nr/nr -outfmt %f > 
  # likely_dtnpb_candidates.blast.faa

# Fetch IPG info with Batch-entrez for each accession list
  # blast_results/likely_dtnpb_candidates.blast.ipg1.txt  to  likely_dtnpb_candidates.blast.ipg9.txt

####################
# Check DED deterioration
####################

# Import proteins isolated in BLAST search
prot.blast <- readAAStringSet("blast_results/likely_dtnpb_candidates.blast.faa")

# Fix names
for(i in 1:length(prot.blast)) { 
  if(word(names(prot.blast)[i]) %in% blast$sacc) { names(prot.blast)[i] <- word(names(prot.blast)[i]) }
  if(!word(names(prot.blast)[i]) %in% blast$sacc) { x <- gsub(">","",str_split(names(prot.blast)[i], " ")[[1]])
  names(prot.blast)[i] <- x[x %in% blast$sacc] } }

# Combine BLAST hits with cluster members + tree representatives
prot.all <- unique(c(prot, prot.blast))

# Remove BLAST entries that I couldn't retrieve the sequence for
blast <- blast %>% filter(sacc %in% names(prot.all))

# Remove BLAST entries that we can't cross-reference with cluster members
blast <- blast %>% filter(qseqid %in% clus.all$id)

# Add tree representative to blast output
blast <- blast %>% left_join(clus.all %>% select(qseqid = id, tree.rep), by = "qseqid")

# Make list of unique tree representatives
unique.reps <- unique(blast$tree.rep)

# Make directory for outputs
dir.create("likely_dtnpb_w_blast/", showWarnings = FALSE)
dir.create("unlikely_dtnpb_w_blast/", showWarnings = FALSE)
dir.create("cluster_alignments_w_blast/", showWarnings = FALSE)


# Iterate through each tree representative
for(i in 1:length(unique.reps)){
  
  # Determine which other proteins are cluster members with the tree representative
  y <- unique(c(unique.reps[i], clus.all$id[clus.all$tree.rep %in% unique.reps[i]]))
  
  # Narrow blast results info to only the candidate being considered
  blast.these <- blast %>% filter(qseqid %in% y)
  
  # Get protein sequences of cluster members + tree representative + blast hits
  prot.these <- prot.all[names(prot.all) %in% c(y, blast.these$sacc)]
  
  # Import untrimmed alignment
  msa <- readAAStringSet("../../first_look/Catalytic_triad_completeness.firstlook.orderedMSA.trimmed.afa")
  
  # Determine which sequences are the closest in tree
  x <- grep(unique.reps[i], names(msa))
  x <- c( (x-5):(x-1), (x+1):(x+5) )
  x <- x[x > 1 & x < 2646]
  neighbors <- msa[x]
  x <- AAStringSet(gsub("-", "", unname(as.character(neighbors))))
  names(x) <- paste0(names(neighbors), " | neighbor ")
  neighbors <- x
  
  # Import reference sequences
  ref <- c(readAAStringSet("../../../inputs/refs/ISDra2_TnpB.faa"),
           readAAStringSet("../../../inputs/refs/7C7L_Cas12f.faa"))
  
  # Align sequences to closest neighbors in tree + references
  msa <- AlignSeqs(c(ref, neighbors, prot.these))
  
  # Build dataframes of positions for reference sequences
  ref.tnpb.df <- tibble(res = str_split(unname(as.character(msa[1])), "")[[1]], pos.algn = 1:width(msa[1])) %>% 
    filter(!res == "-") %>% mutate(pos = 1:width(ref[1]))
  ref.cas12f.df <- tibble(res = str_split(unname(as.character(msa[2])), "")[[1]], pos.algn = 1:width(msa[2])) %>% 
    filter(!res == "-") %>% mutate(pos = 1:width(ref[2]))
  
  # Store catalytic residues
  cat.res <- c("D", "E")
  
  ## Build dataframe of residues in cluster members that correspond to catalytic residues in each reference sequence
  # Check ISDra2-like catalytic triad
  cluster.tnpb.df <- tibble(id = names(msa[3:length(msa)]), 
                            res1 = unname(as.character(narrow(msa[3:length(msa)],
                                                              start = ref.tnpb.df$pos.algn[ref.tnpb.df$pos == 191], width = 1))) %in% cat.res,
                            res2 = unname(as.character(narrow(msa[3:length(msa)],
                                                              start = ref.tnpb.df$pos.algn[ref.tnpb.df$pos == 278], width = 1))) %in% cat.res,
                            res3 = unname(as.character(narrow(msa[3:length(msa)],
                                                              start = ref.tnpb.df$pos.algn[ref.tnpb.df$pos == 361], width = 1))) %in% cat.res) %>% 
    mutate(sum.triad = res1 + res2 + res3)
  
  # Check Cas12f-like catalytic triad
  cluster.cas12f.df <- tibble(id = names(msa[3:length(msa)]), 
                              res1 = unname(as.character(narrow(msa[3:length(msa)],
                                                                start = ref.cas12f.df$pos.algn[ref.cas12f.df$pos == 326], width = 1))) %in% cat.res,
                              res2 = unname(as.character(narrow(msa[3:length(msa)],
                                                                start = ref.cas12f.df$pos.algn[ref.cas12f.df$pos == 422], width = 1))) %in% cat.res,
                              res3 = unname(as.character(narrow(msa[3:length(msa)],
                                                                start = ref.cas12f.df$pos.algn[ref.cas12f.df$pos == 510], width = 1))) %in% cat.res) %>% 
    mutate(sum.triad = res1 + res2 + res3)
  
  # Cross-reference each of the dataframes with catalytic residue info
  which.best.ref <- tibble(ref.id = c("ISDra2", "Cas12f"),
                           max.triad = c(sum(cluster.tnpb.df$sum.triad),
                                         sum(cluster.cas12f.df$sum.triad)))
  
  # Choose best reference sequence to use to assess catalytic residue conservation
  ref.df <- ref.tnpb.df ; ref.df$ref <- "ISDra2" # Automatically use ISDra2 as reference, unless another ref scores higher
  if(length(which.best.ref$ref.id[which.best.ref$max.triad == max(which.best.ref$max.triad)]) == 1) {
    if(which.best.ref$ref.id[which.best.ref$max.triad == max(which.best.ref$max.triad)] == "Cas12f") { ref.df <- ref.cas12f.df
    ref.df$ref <- "Cas12f" }}
  
  # Reference chosen, keep triad info for that reference
  if(ref.df$ref[1] == "ISDra2") { cluster.df <- cluster.tnpb.df ; cluster.df$ref <- "ISDra2" }
  if(ref.df$ref[1] == "Cas12f") { cluster.df <- cluster.cas12f.df ; cluster.df$ref <- "Cas12f" }
  
  # Determine if sequences are actually dTnpB
  if(sum(cluster.df$sum.triad[cluster.df$id %in% names(prot.these)] > 1) > 0) { 
    
    # Write out unlikely sequences
    writeXStringSet(prot.these[names(prot.these) %in% cluster.df$id[cluster.df$sum.triad > 1]],
                    paste0("unlikely_dtnpb_w_blast/",unique.reps[i],"_related.faa"))
    
    # Remove these sequences from FASTA so they are not processed further
    prot.these <- prot.these[!names(prot.these) %in% cluster.df$id[cluster.df$sum.triad > 1]]
    
    # Save sequence info to TSV file
    write_tsv(cluster.df %>% filter(),
              paste0("unlikely_dtnpb_w_blast/",unique.reps[i],"_related.tsv"))
    
    # Remove these sequences from cluster member info so they are not processed further
    cluster.df <- cluster.df %>% filter(sum.triad <= 1)
    
  } 
  
  if(nrow(cluster.df) > 0) {
    
    # Save remaining sequences
    writeXStringSet(prot.these, paste0("likely_dtnpb_w_blast/",unique.reps[i],"_related.faa"))
    
    # Save sequence info
    write_tsv(cluster.df, paste0("likely_dtnpb_w_blast/",unique.reps[i],"_related.tsv"))
    
  }
  
  # Change names in alignment to show if they were rejected or not
  names(msa)[!names(msa) %in% names(prot.these) & !names(msa) %in% c(names(ref), names(neighbors))] <- 
    paste0(names(msa)[!names(msa) %in% names(prot.these) & !names(msa) %in% c(names(ref), names(neighbors))], " - removed*")
  
  # Save alignment
  writeXStringSet(msa, paste0("cluster_alignments_w_blast/",unique.reps[i],"_related.afa"))


}



####################
# Get genomes of likely dTnpB candidates
####################

library(XML)

# Start clean
rm(list = ls())

# Import BLAST results
colz <- c("qseqid", "sseqid", "qlen", "qstart", "qend", "sstart", "send", "sstrand",
          "evalue", "length", "pident", "nident", "mismatch", "gaps", "qcov", "staxid")
blast <- rbind(read_tsv("blast_results/likely_dTnpB.blastout.v1.txt", col_names = colz),
               read_tsv("blast_results/likely_dTnpB.blastout.v2.txt", col_names = colz))

# Add column with subject accession only
blast$sacc <- word(blast$sseqid, start = 2, sep = "\\|")

# Import initial clustering info
clus.all <- read_tsv("cluster.info.tsv")

# Import IPG info
system("ls blast_results > temp.filelist.txt")
filelist <- read_tsv("temp.filelist.txt", col_names = "file")
filelist <- filelist$file ; system("rm temp.filelist.txt")
filelist <- filelist[grep("likely_dtnpb_candidates.blast.ipg", filelist)]
ipg <- tibble()
for(i in filelist) { ipg <- rbind(ipg, read_tsv(paste0("blast_results/", i))) }

# Narrow down list by combining start-end coordinates of dTnpB with organism name
temp <- tibble(temp.id = paste0(ipg$Start,"_",ipg$Stop,"_", ipg$Organism), acc = ipg$`Nucleotide Accession`)
temp <- temp %>% group_by(temp.id) %>% summarize(these = dplyr::first(acc))
ipg <- ipg %>% filter(`Nucleotide Accession` %in% temp$these)

### Determine the length of each genomic sequence encoding dTnpB
# Make new directory
dir.create("genomic_loci/xml_files", recursive = TRUE, showWarnings = FALSE)

# Make list of accessions
temp <- tibble(acc = unique(ipg$`Nucleotide Accession`))

# Export accessions in batches
# write_tsv(temp[1:2000,], "genomic_loci/xml_files/ipg.acc.1.txt", col_names = FALSE)
# write_tsv(temp[3000:5000,], "genomic_loci/xml_files/ipg.acc.2.txt", col_names = FALSE)
# write_tsv(temp[c(2001:2999,5001:7499),], "genomic_loci/xml_files/ipg.acc.3.txt", col_names = FALSE) # oops fill in missed range
# write_tsv(temp[7500:10000,], "genomic_loci/xml_files/ipg.acc.4.txt", col_names = FALSE)
# write_tsv(temp[10001:15000,], "genomic_loci/xml_files/ipg.acc.5.txt", col_names = FALSE)
# write_tsv(temp[15001:20000,], "genomic_loci/xml_files/ipg.acc.6.txt", col_names = FALSE)
# write_tsv(temp[20001:25000,], "genomic_loci/xml_files/ipg.acc.7.txt", col_names = FALSE)
# write_tsv(temp[25001:30000,], "genomic_loci/xml_files/ipg.acc.8.txt", col_names = FALSE)
# write_tsv(temp[30001:34492,], "genomic_loci/xml_files/ipg.acc.9.txt", col_names = FALSE)

# Fetch genomic sequence info for each list of accessions as XML file with Eutils:
  # epost -db nucleotide -input ipg.acc.1.txt | efetch -format docsum > ipg.genomic.info1.xml
  # epost -db nucleotide -input ipg.acc.2.txt | efetch -format docsum > ipg.genomic.info2.xml
  # epost -db nucleotide -input ipg.acc.3.txt | efetch -format docsum > ipg.genomic.info3.xml
  # epost -db nucleotide -input ipg.acc.4.txt | efetch -format docsum > ipg.genomic.info4.xml
  # epost -db nucleotide -input ipg.acc.5.txt | efetch -format docsum > ipg.genomic.info5.xml
  # epost -db nucleotide -input ipg.acc.6.txt | efetch -format docsum > ipg.genomic.info6.xml
  # epost -db nucleotide -input ipg.acc.7.txt | efetch -format docsum > ipg.genomic.info7.xml
  # epost -db nucleotide -input ipg.acc.8.txt | efetch -format docsum > ipg.genomic.info8.xml
  # epost -db nucleotide -input ipg.acc.9.txt | efetch -format docsum > ipg.genomic.info9.xml

# Import XML files
genomic <- tibble()
system("ls genomic_loci/xml_files > temp.filelist.txt")
filelist <- read_tsv("temp.filelist.txt", col_names = "file")
filelist <- filelist$file ; system("rm temp.filelist.txt")
filelist <- filelist[grep("xml", filelist)]
for(i in filelist){
x <- read_tsv(paste0("genomic_loci/xml_files/",i), col_names = "this")
gen <- tibble(acc = x$this[grep("AccessionVersion", x$this)], len = x$this[grep('Stat type="Length" count', x$this)])
gen$acc <- word(word(gen$acc, start = 2, sep = "\\>"), sep = "\\<")
gen$len <- as.numeric(str_sub(word(gen$len, start = 2, sep = 'count="'), start = 1, end = -4))
genomic <- rbind(genomic, gen)
rm(gen, x)
}

# Remove IPG info that we couldn't retrieve a genome for
ipg <- ipg %>% filter(`Nucleotide Accession` %in% genomic$acc)

# Add Tree representative to blast results
blast <- blast %>% left_join(clus.all %>% select(qseqid = id, tree.rep), by = "qseqid")

# Add Tree representative to IPG output
ipg <- ipg %>% left_join(blast %>% select(Protein = sacc, tree.rep), by = "Protein")

# Add genomic length info to IPG output
colnames(ipg)[3] <- "gen.acc"
ipg <- ipg %>% left_join(genomic %>% select(gen.acc = acc, genome.length = len), by = "gen.acc")

# Add flank.start and flank.end to IPG output
ipg$flank.start <- ipg$Start - 20000
ipg$flank.end <- ipg$Stop + 20000
ipg$flank.start[ipg$flank.start < 1] <- 1
ipg$flank.end[ipg$flank.end > ipg$genome.length] <- ipg$genome.length[ipg$flank.end > ipg$genome.length]

# Add flank length to IPG output
ipg$flank.length <- ipg$flank.end - ipg$flank.start

# Remove flanks that are less than 15kbp
ipg <- ipg %>% filter(flank.length > 15000)

# Remove IPGs that can't be associated back to a tree representative
ipg.all <- ipg
ipg <- ipg.all %>% filter(!is.na(tree.rep))

# Export genomic accessions
temp <- tibble(gen.acc = unique(ipg$gen.acc))
write_tsv(temp[1:1500,], "genomic_loci/ipg.genacc1.txt", col_names = FALSE)
write_tsv(temp[1501:3000,], "genomic_loci/ipg.genacc2.txt", col_names = FALSE)
write_tsv(temp[3001:4500,], "genomic_loci/ipg.genacc3.txt", col_names = FALSE)
write_tsv(temp[4501:6000,], "genomic_loci/ipg.genacc4.txt", col_names = FALSE)
write_tsv(temp[6001:7500,], "genomic_loci/ipg.genacc5.txt", col_names = FALSE)
write_tsv(temp[7501:7933,], "genomic_loci/ipg.genacc6.txt", col_names = FALSE)

# Export narrowed IPG
write_tsv(ipg, "blast_results/likely_dtnpb_candidates.blast.ipg.filtered.tsv")

# Pulled genomes with Batch Entrez
  # genomic_loci/ipg.genomic1.fna
  # genomic_loci/ipg.genomic2.fna
  # genomic_loci/ipg.genomic3.fna
  # genomic_loci/ipg.genomic4.fna
  # genomic_loci/ipg.genomic5.fna
  # genomic_loci/ipg.genomic6.fna

########################
# Pull genomic flanks of BLAST dTnpB proteins
########################

# Start clean
rm(list = ls())

# Get list of genomic sequences
system("ls genomic_loci > temp.filelist.txt")
filelist <- read_tsv("temp.filelist.txt", col_names = "file")
filelist <- filelist$file ; system("rm temp.filelist.txt")
filelist <- filelist[grep("fna", filelist)]

# Import IPG info
ipg <- read_tsv("blast_results/likely_dtnpb_candidates.blast.ipg.filtered.tsv")

# Iterate through genomes to pull flanks
flank.size <- 20000
df.flank <- tibble()
for(z in filelist){
  
  flanks <- DNAStringSet()
  
  # Import genomic loci
  genomes <- readDNAStringSet(paste0("genomic_loci/", z))
  
  # Pull flanks of dTnpB
  for(i in 1:length(genomes)){
    
    # Isolate one genome
    x <- genomes[i]
    
    # Get annotations for that genome
    annot <- ipg %>% filter(gen.acc == word(names(x)))
    annot$Start <- as.numeric(annot$Start) ; annot$Stop <- as.numeric(annot$Stop)
    
    # Some genomes have 2 copies
    for(v in 1:nrow(annot)){
      
      # # Trim genome to only include dTnpB + max length of sequence (if 3kbp is longer than end) 
      flank.start <- annot$flank.start[v]
      flank.end <- annot$flank.end[v]
      
      # Trim genome to only include dTnpB + 3kbp flanks
      this.flank <- narrow(x, start = flank.start, end = flank.end)
      
      # Add coordinates to FASTA header
      names(this.flank) <- paste0(word(names(this.flank)), "_", flank.start, "-", end = flank.end)
      
      # Flip sequence, if the dTnpB is on the "-" strand
      if(annot$Strand[v] == "-") { this.flank <- reverseComplement(this.flank) 
      names(this.flank) <- paste0(names(this.flank), "_rc") }
      
      # Add column for relative start and stop coordinates of dTnpB
      temp.df <- annot[v,] %>% mutate(flank.id = names(this.flank), rel.start = NA, rel.end = NA)
      
      # Add relative coordinates for dTnpB encoded on forward strand
      temp.df$rel.start[temp.df$Strand == "+"] <-
        temp.df$Start[temp.df$Strand == "+"] - temp.df$flank.start[temp.df$Strand == "+"] + 1
      temp.df$rel.end[temp.df$Strand == "+"] <- 
        temp.df$Stop[temp.df$Strand == "+"] - temp.df$flank.start[temp.df$Strand == "+"] + 1
      
      # Add relative coordinates for dTnpB encoded on reverse strand
      temp.df$rel.start[temp.df$Strand == "-"] <-
        temp.df$flank.end[temp.df$Strand == "-"] - temp.df$Stop[temp.df$Strand == "-"] + 1
        
      temp.df$rel.end[temp.df$Strand == "-"] <-
        temp.df$flank.end[temp.df$Strand == "-"] - temp.df$Start[temp.df$Strand == "-"] + 1
        
      # Add flank to new FASTA
      flanks <- c(flanks, this.flank) 
      
      # Add flank info to dataframe
      df.flank <- rbind(df.flank, temp.df)
      
    }
    
    # Start fresh, to raise errors if something goes wrong with the next genome
    rm(x, annot, this.flank, flank.start, flank.end)
    
  }
  
  # Export flanks
  writeXStringSet(flanks, "dTnpB_BLAST_genomic.flanks.fna", append = TRUE)
  
  # Start fresh with next batch of genomes
  rm(flanks, genomes)
  
  # Print status message
  message(paste0("Exported flanks for genomes batch #", z))
  
}

# Remove duplicate flank entries
df.flank <- df.flank %>% distinct()

# Import flanks
flanks <- readDNAStringSet("dTnpB_BLAST_genomic.flanks.fna")

# Remove duplicates
flanks <- unique(flanks)
df.flank.unique <- df.flank %>% filter(flank.id %in% names(flanks))

# Export flanks and info
writeXStringSet(flanks, "dTnpB_BLAST_genomic.flanks.unique.fna")
write_tsv(df.flank.unique, "dTnpB_BLAST_genomic.flanks.unique.tsv")
write_tsv(df.flank, "dTnpB_BLAST_genomic.flanks.tsv")


########################
# Annotate flanks with Eggnog
########################

rm(list = ls())

# Run eggnog mapper on Titan
  # emapper.py --cpu 28 -o out \
  #   --output_dir flanks_annotation \
  #   --override -m diamond --dmnd_ignore_warnings \
  #   -i dTnpB_BLAST_genomic.flanks.unique.fna \
  #   --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 \
  #   --itype genome --genepred prodigal --tax_scope auto --target_orthologs all \
  #   --go_evidence non-electronic --pfam_realign none --report_orthologs \
  #   --decorate_gff yes \
  #   --excel

# Results manually added to new directory
  # outputs/neighborhoods3/automatic_annotations/Eggnog_annotations
    # Note: deleted out.orthologs, since it's a large file that I don't need

########################
# Cluster flanks at 95% nt identity
########################

# Cluster flanks
  # mmseqs easy-cluster dTnpB_BLAST_genomic.flanks.unique.fna dTnpB_BLAST_genomic_flanks_mmseqs tmp --min-seq-id 0.95 --cov-mode 0

# Import clustering info via names of representatives flank sequences (i.e., MMseqs output)
cluster <- readDNAStringSet("dTnpB_BLAST_genomic_flanks_mmseqs_rep_seq.fasta")

# Fix names by removing whitespace
names(cluster) <- gsub(" ", "", names(cluster))

# Import likely dTnpB info
likely.all <- read_tsv("dTnpB_BLAST_genomic.flanks.unique.tsv")

# Remove any spaces that got introduced into flank info
likely.all$flank.id <- gsub(" ", "", likely.all$flank.id)

# Import flank fasta
seq.flank.all <- readDNAStringSet("dTnpB_BLAST_genomic.flanks.unique.fna")

# Keep only clustered flank representatives
likely <- likely.all %>% filter(flank.id %in% names(cluster))

# Export info for clustered flanks
write_tsv(likely, "dTnpB_BLAST_genomic.flanks.unique.clustered.tsv")
writeXStringSet(cluster, "dTnpB_BLAST_genomic.flanks.unique.clustered.fna")


########################
# Format annotations
########################

# Start clean
rm(list = ls())

#### Import pertinent info
# Import GFF file of eggnog output
egg <- read_gff("Eggnog_annotations/out.emapper.decorated.gff")

# Import flank info
likely <- read_tsv("dTnpB_BLAST_genomic.flanks.unique.clustered.tsv")

# Import flank sequences
seq.flank <- readDNAStringSet("dTnpB_BLAST_genomic.flanks.unique.clustered.fna")

###### Annotate each dTnpB + Flank sequence

# Create directory for cluster FASTA files
dir.create("../flank_fasta", showWarnings = FALSE)

# Create directory for cluster FASTA annotation files
dir.create("../flank_annotations", showWarnings = FALSE)

# Initiate function to export flanks and annotations (GFF files) for each cluster
export.cluster.annotations <- function(rep.id){
  
  # First pull ids of all cluster members
  clus <- likely %>% filter(tree.rep %in% rep.id)
  
  # Pull flanks for cluster members
  clus.seq <- seq.flank[names(seq.flank) %in% clus$flank.id]
  
  # Get Eggnog annotation info
  clus.egg <- egg %>% filter(seqid %in% c(clus$flank.id))
  
  # Names got messed up for some annotations (_rc was dropped)
  these <- clus$flank.id[!clus$flank.id %in% clus.egg$seqid]
  if(length(these) > 0){
    
    possible <- egg %>% mutate(seqid = paste0(seqid, "_rc")) %>% filter(seqid %in% these)
    clus.egg <- rbind(clus.egg, possible)
    
  }
  
  # Change annotation type to be COG category
  clus.egg$type <- word(word(clus.egg$attribute, start = 2, sep = "em_COG_cat="), sep = "\\;")
  clus.egg$type[is.na(clus.egg$type)] <- "CDS"
  clus.egg$type[clus.egg$type == "None"] <- "Other"
  
  # Add dTnpB as new annotations to Eggnog output
  clus.egg <- rbind(clus.egg, tibble(seqid = clus$flank.id,
                                     source = rep("IPG", nrow(clus)),
                                     type = rep("dTnpB", nrow(clus)),
                                     start = clus$rel.start,
                                     end = clus$rel.end,
                                     score = rep(".", nrow(clus)),
                                     strand = rep("+", nrow(clus)),
                                     phase = rep(".", nrow(clus)),
                                     attribute = rep("dTnpB_from_tree", nrow(clus)),
  ))
  
  # Export FASTA of sequences
  writeXStringSet(clus.seq, paste0("../flank_fasta/", rep.id, ".flanks.fna"))
  
  # Export GFF annotation files
  write_tsv(clus.egg, paste0("../flank_annotations/", rep.id, ".flanks.annot.gff"), col_names = FALSE)
  
}

# Make list of representatives in tree
x <- unique(likely$tree.rep)
for(i in x) {export.cluster.annotations(rep.id = i)}


########################
# Viewing Info
########################

# Import Flanks and Annotations into Geneious and manually look for stably associated genes

########################
# Make final dataframe to format in Excel
########################

# Start clean
rm(list = ls())

# Import flank info
likely <- read_tsv("dTnpB_BLAST_genomic.flanks.unique.clustered.tsv")

### Import all catalytic residue info
# Likely dTnpB info
system("ls likely_dtnpb_w_blast > temp.filelist.txt")
filelist <- read_tsv("temp.filelist.txt", col_names = "file")
filelist <- filelist$file ; system("rm temp.filelist.txt")
filelist <- filelist[grep("tsv", filelist)]

res <- tibble()
for(i in 1:length(filelist)) { res <- rbind(res, read_tsv(paste0("likely_dtnpb_w_blast/", filelist[i])))}

# Unlikely dTnpB info
system("ls unlikely_dtnpb_w_blast > temp.filelist.txt")
filelist <- read_tsv("temp.filelist.txt", col_names = "file")
filelist <- filelist$file ; system("rm temp.filelist.txt")
filelist <- filelist[grep("tsv", filelist)]

for(i in 1:length(filelist)) { res <- rbind(res, read_tsv(paste0("unlikely_dtnpb_w_blast/", filelist[i])))}

# Remove all "neighbor" entries
res <- res[-grep("neighbor", res$id),]

# Keep only info for likely dTnpBs
res <- res %>% filter(id %in% likely$Protein) %>% distinct()

# Filter in the other direction to make sure that we're only keeping proteins that made it through the DED check
likely <- likely %>% filter(Protein %in% res$id) %>% distinct()

# Join dataframes, dropping unneeded columns
df <- likely %>% filter(Protein %in% res$id[res$sum.triad < 2]) %>% select(-tree.rep) %>% distinct() %>% 
  left_join(res %>% select(Protein = id, res1, res2, res3, sum.triad) %>% distinct(), by = "Protein")

# One rogue entry with two different DED residue counts, being conservative, we're gonna take the higher count
  # NZ_CP016793.1_2250687-2291508 [2 looks right from a check on HHpred against ISDra2]
df <- df %>% mutate(temp = paste0(flank.id, "_", sum.triad)) %>% filter(!temp == "NZ_CP016793.1_2250687-2291508_1") %>% 
  select(-temp)

### Import protein sequences
# Likely dTnpB info
system("ls likely_dtnpb_w_blast > temp.filelist.txt")
filelist <- read_tsv("temp.filelist.txt", col_names = "file")
filelist <- filelist$file ; system("rm temp.filelist.txt")
filelist <- filelist[grep("faa", filelist)]

prot <- AAStringSet()
for(i in 1:length(filelist)) { prot <- c(prot, readAAStringSet(paste0("likely_dtnpb_w_blast/", filelist[i])))}

# Unlikely dTnpB info
system("ls unlikely_dtnpb_w_blast > temp.filelist.txt")
filelist <- read_tsv("temp.filelist.txt", col_names = "file")
filelist <- filelist$file ; system("rm temp.filelist.txt")
filelist <- filelist[grep("faa", filelist)]

for(i in 1:length(filelist)) { prot <- c(prot, readAAStringSet(paste0("unlikely_dtnpb_w_blast/", filelist[i])))}

# Remove proteins that are likely to be TnpBs and make a dataframe of sequences
prot <- prot[names(prot) %in% df$Protein]
prot.df <- tibble(Protein = names(prot), seq = unname(as.character(prot)))

# Add sequences to master dataframe
df <- df %>% left_join(prot.df %>% distinct(), by = "Protein")

# Export master dataframe of likely dTnpB proteins
write_tsv(df, "../../../dTnpB_proteins.tsv")





