options(width=160)

# Get full length library
tfp = read.csv(gzfile("../library/transfac_proteins.csv.gz"))
colnames(tfp) = c("ensembl","gene","DBD","TF.assessment","binding.mode","motif.status","aa")

# Remove stop codons
i = which(substr(tfp$aa,nchar(tfp$aa),nchar(tfp$aa))=="*")
tfp[i,"aa"] = sub("*", "", tfp[i,"aa"], fixed=T)
tfp$nres = nchar(tfp$aa)

# Get PAP scores for full length proteins
pap = read.table(gzfile("transfac_proteins_pap.txt.gz"))
colnames(pap) = c("gene","aa","score_ct","score","resi_first","resi_last","resn_cent","resi_cent")

# Mark proteins that have zinc fingers
tfp$zf = ifelse( grepl("C..C............H...H", tfp$aa), "T", "F" )
print(sprintf("Marked %d of %d proteins for containing at east one zinc finger (%.1f%%)",
              sum(tfp$zf=="T"), nrow(tfp), sum(tfp$zf=="T")/nrow(tfp)*100))

# Find C-degron score per protein (all CT tiles have it)
pap$nres = nchar(tfp$aa)[match(pap$gene,tfp$gene)]
pap_ct = pap[which(pap$resi_last==pap$nres),c("gene","score_ct")]
tfp$pap_ct = pap_ct[match(tfp$gene,pap_ct$gene),"score_ct"]

# Make a string of colon-seperated PAP scores per protein
agg = aggregate(pap$score, by=list(pap$gene), paste0, collapse=";")
tfp$pap = agg[match(tfp$gene,agg[,1]),2]

# Get ADPred predictions
adp = read.table(gzfile("adpred_ss_api_full_scores.txt.gz"))
colnames(adp) = c("gene","ss","scores")
i = which(adp$gene=="AC092835.1")
adp[i,"gene"] = "ZN892"
print(sprintf("Changed %d gene in adpred name from AC092835.1 to ZN892",length(i)))
i = which(adp$gene=="T")
adp[i,"gene"] = "TBXT"
print(sprintf("Changed %d gene in adpred name from T to TBXT",length(i)))
tfp$adpred = adp[match(tfp$gene,adp$gene),"scores"]

# Get structural features
struc = read.csv(gzfile("transfac_struc.csv.gz"))
i = which(struc$gene=="AC092835.1")
struc[i,"gene"] = "ZN892"
print(sprintf("Changed %d gene name in struc from AC092835.1 to ZN892",length(i)))
i = which(struc$gene=="T")
struc[i,"gene"] = "TBXT"
print(sprintf("Changed %d gene name in struc from T to TBXT",length(i)))


# collect per residue scores from tiles and make a string per protein
py2vec = function(v) { strsplit( gsub("\\[|\\]| ", "", v), ",")[[1]] }
for (target in c("rasa","plddt")) {
    all_seqs = lapply(tfp$nres, function(n){ rep("",n) })
    names(all_seqs) = tfp$gene
    # tile data frame has controls so only rows of genes that are in the list of proteins
    for (i in which(struc$gene %in% tfp$gene)) {
        all_seqs[[struc[i,"gene"]]][struc[i,"first_resi"]:struc[i,"last_resi"]] = py2vec(struc[i,sprintf("tile_%s",target)])
    }
    stopifnot(all( names(all_seqs) == tfp$gene ))
    tfp[,target] = sapply(all_seqs, paste0, collapse=";")
}

py2vec = function(v) { strsplit( gsub("\\[|\\]| |'", "", v), ",")[[1]] }

# Get uniprot id
# select first tile of each gene for processing, skip controls
is = which(struc$first_resi==1 & struc$gene != "")
stopifnot(all( struc$gene %in% c("",struc[is,"gene"]) ))
stopifnot( length(struc[is,"gene"]) == length(unique(struc[is,"gene"])) )

# this is the uniprot used to calc exposure, if matching uniprots (col uniprot_id) does not have a structure in AFDB this col is empty
unip_list = lapply(struc[is,"uniprot_id_exposure"], py2vec)
names(unip_list) = struc[is,"gene"]

# this is all uniprot entries that match amino acid sequence
unip_all_list   = lapply(struc[is,"uniprot_id"],    py2vec)

# where missing, fill in the uniprot from the list of matching entries
i = which(sapply(unip_list,length)==0)
# check that these have zero or one entry listed
stopifnot(all( sapply(unip_all_list[i], length) <= 1 ))
unip_list[i] = sapply(unip_all_list[i], "[", 1)

# Get TF first and last position in the matched uniprot, and the full length of this
unip_first_list = lapply(struc[is,"uniprot_start"], py2vec)
unip_last_list  = lapply(struc[is,"uniprot_end"],   py2vec)
unip_nres_list  = lapply(struc[is,"uniprot_n_res"], py2vec)
unip = data.frame(gene=names(unip_list), uniprot=sapply(unip_list,"[",1))
unip_all_indices = sapply(seq(nrow(unip)), function(i){ ifelse(is.na(unip[i,"uniprot"]), NA, which(unip_all_list[[i]] == unip[i,"uniprot"])) })
unip$first_resi = as.numeric( sapply(seq(nrow(unip)), function(i){ ifelse(is.na(unip_all_indices[i]), NA, unip_first_list[[i]][unip_all_indices[i]]) } ) )
unip$last_resi  = as.numeric( sapply(seq(nrow(unip)), function(i){ ifelse(is.na(unip_all_indices[i]), NA, unip_last_list[[i]][unip_all_indices[i]]) } ) )
unip$nres       = as.numeric( sapply(seq(nrow(unip)), function(i){ ifelse(is.na(unip_all_indices[i]), NA, unip_nres_list[[i]][unip_all_indices[i]]) } ) )

i = match(tfp$gene, unip$gene)
tfp$uniprot    = unip[i,"uniprot"]
tfp$unip_first = unip[i,"first_resi"]
tfp$unip_last  = unip[i,"last_resi"]
tfp$unip_nres  = unip[i,"nres"]

write.csv(tfp, row.names=F, quote=F, file="transfac_proteins.csv")

# # minimal set to publish
# cns = c("gene","ensembl","uniprot","aa","pap","rasa","plddt","adpred")
# write.csv(tfp[,cns], row.names=F, quote=F, file="transfac_proteins.csv")

