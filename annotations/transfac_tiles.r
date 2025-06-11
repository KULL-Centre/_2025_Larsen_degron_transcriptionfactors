
# Read TF scores
degron = read.csv(gzfile("../score/degrons.csv.gz"))

# Full TF tile lib (without controls data)
tf = read.csv(gzfile("../library/transfac_tiles.csv.gz"))
# Remove stop codons in AA sequences
i3 = which(substr(tf$aa,30,30)=="*")
tf[i3,"aa"] = sub("*", "", tf[i3,"aa"], fixed=T)
tf[i3,"last_resi"] = tf[i3,"last_resi"] - 1

# Controls
ctrl = read.csv("../score/controls.csv")
ctrl$gene = NA
ctrl$first_resi = 1
ctrl$last_resi = nchar(ctrl$aa)
ctrl$copyof_dna = match(ctrl$dna,tf$dna)
ctrl$copyof_aa = match(ctrl$aa,tf$aa)
ctrl$lib = 0

# Add controls to library
tf = rbind(tf, ctrl)

# Mark tiles that have zinc fingers
tf$zf = ifelse( grepl("C..C............H...H", tf$aa), "T", "F" )
print(sprintf("Marked %d of %d tiles from %s proteins for containing a zinc finger (%.1f%%)",
              sum(tf$zf=="T"), nrow(tf), length(unique(tf[which(tf$zf=="T"),"gene"])), sum(tf$zf=="T")/nrow(tf)*100))
stopifnot(all( tf[which(tf$zf=="T"),"gene"] != "" ))

# Map scores to TF library
i = match(tf$dna,degron$dna)
print(sprintf("Mapping %d scores to %d library members (%.2f%%)", nrow(degron), sum(! is.na(i)), sum(! is.na(i))/nrow(tf)*100 ))
cns = c("n","degron_score", "degron_std", "abundance_score", "abundance_std")
tf[,cns] = degron[i,cns]
stopifnot(all( is.na(tf$abundance_score) == is.na(tf$degron_score) ))

# Map PAP scores to TF library
pap = read.table(gzfile("transfac_tiles_pap.txt.gz"))
colnames(pap) = c("name","aa","score_ct","score","resi_first","resi_last","resn_cent","resi_cent")
i = match(tf$aa, pap$aa)
print(sprintf("Mapping %d pap scores to %d library members (%.2f%%)", nrow(degron), sum(! is.na(i)), sum(! is.na(i))/nrow(tf)*100 ))
tf$pap_ct = pap[i,"score_ct"]
tf$pap = pap[i,"score"]

# Map LocNES scores to TF library
lns = read.csv(gzfile("transfac_tiles_locnes.csv.gz"))
i = match(tf$name, lns$name)
print(sprintf("Mapping %d LocNES scores to %d library members (%.2f%%)", nrow(lns), sum(! is.na(i)), sum(! is.na(i))/nrow(tf)*100 ))
tf$locnes = lns[i,"locnes"]
tf$locnes_score = lns[i,"locnes_score"]

# Map structural features
struc = read.csv(gzfile("transfac_struc.csv.gz"))
i = match(tf$name, struc$name)
print(sprintf("Mapping %d structure features to %d library members (%.2f%%)", nrow(struc), sum(! is.na(i)), sum(! is.na(i))/nrow(tf)*100 ))
tf[,"tile_rasa_median"] = struc[i,"tile_rasa_median"]
tf[,"tile_plddt_median"] = struc[i,"tile_plddt_median"]

write.csv(tf, row.names=F, quote=F, file="transfac_tiles.csv")


