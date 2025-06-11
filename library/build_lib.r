options(width=200)

codon_nat = list()
# T          T                      C                    A                     G
codon_nat[["TTT"]] = "F"; codon_nat[["TCT"]] = "S"; codon_nat[["TAT"]] = "Y"; codon_nat[["TGT"]] = "C" # T
codon_nat[["TTC"]] = "F"; codon_nat[["TCC"]] = "S"; codon_nat[["TAC"]] = "Y"; codon_nat[["TGC"]] = "C" # C
codon_nat[["TTA"]] = "L"; codon_nat[["TCA"]] = "S"; codon_nat[["TAA"]] = "*"; codon_nat[["TGA"]] = "*" # A
codon_nat[["TTG"]] = "L"; codon_nat[["TCG"]] = "S"; codon_nat[["TAG"]] = "*"; codon_nat[["TGG"]] = "W" # G
# C
codon_nat[["CTT"]] = "L"; codon_nat[["CCT"]] = "P"; codon_nat[["CAT"]] = "H"; codon_nat[["CGT"]] = "R" # T
codon_nat[["CTC"]] = "L"; codon_nat[["CCC"]] = "P"; codon_nat[["CAC"]] = "H"; codon_nat[["CGC"]] = "R" # C
codon_nat[["CTA"]] = "L"; codon_nat[["CCA"]] = "P"; codon_nat[["CAA"]] = "Q"; codon_nat[["CGA"]] = "R" # A
codon_nat[["CTG"]] = "L"; codon_nat[["CCG"]] = "P"; codon_nat[["CAG"]] = "Q"; codon_nat[["CGG"]] = "R" # G
# A
codon_nat[["ATT"]] = "I"; codon_nat[["ACT"]] = "T"; codon_nat[["AAT"]] = "N"; codon_nat[["AGT"]] = "S" # A
codon_nat[["ATC"]] = "I"; codon_nat[["ACC"]] = "T"; codon_nat[["AAC"]] = "N"; codon_nat[["AGC"]] = "S" # C
codon_nat[["ATA"]] = "I"; codon_nat[["ACA"]] = "T"; codon_nat[["AAA"]] = "K"; codon_nat[["AGA"]] = "R" # A
codon_nat[["ATG"]] = "M"; codon_nat[["ACG"]] = "T"; codon_nat[["AAG"]] = "K"; codon_nat[["AGG"]] = "R" # G
# G
codon_nat[["GTT"]] = "V"; codon_nat[["GCT"]] = "A"; codon_nat[["GAT"]] = "D"; codon_nat[["GGT"]] = "G" # A
codon_nat[["GTC"]] = "V"; codon_nat[["GCC"]] = "A"; codon_nat[["GAC"]] = "D"; codon_nat[["GGC"]] = "G" # C
codon_nat[["GTA"]] = "V"; codon_nat[["GCA"]] = "A"; codon_nat[["GAA"]] = "E"; codon_nat[["GGA"]] = "G" # A
codon_nat[["GTG"]] = "V"; codon_nat[["GCG"]] = "A"; codon_nat[["GAG"]] = "E"; codon_nat[["GGG"]] = "G" # G

dna2aa = function(dna_vec, codon_table=codon_nat, frame_shift=0, direction=1, truncate=0) {
    # If requested, read all sequences 3' to 5'
    if (direction < 0) dna_vec = sapply(dna_vec, function(s) {paste(rev(strsplit(s,"")[[1]]), collapse="")}, USE.NAMES=FALSE)
    translate = function(dna, frame_shift, truncate) {
        aa_seq = c()
        i_last = frame_shift + 3
        npos = nchar(dna)-truncate
        while (i_last <= npos) {
            # print(substr(dna, i_last-2, i_last))
	    codon = substr(dna, i_last-2, i_last)
	    if (! codon %in% names(codon_table)) { print(sprintf("ERROR: Codon %s not in table - skip sequence",codon)); return(NA) }
            aa_seq = append(aa_seq, codon_table[[codon]])
	    i_last = i_last + 3
        }
        paste(aa_seq, collapse="")
    }
    sapply(dna_vec, translate, USE.NAMES=FALSE, frame_shift=frame_shift, truncate=truncate)
}

fix_transcript = function(d, ensg, enst, enst_cn="transcript.version", dont_overwrite=T) {
    i = which(d$ID==ensg)
    if (length(i) != 1) {
        print(sprintf("ERROR: %s maps to %d entries: %s",ensg,length(i),paste(i, collapse=" ")))
	return(NA)
    }
    if (! is.na(d[i,enst_cn]) & dont_overwrite) {
        print(sprintf("ERROR: %s already assigned transcript %s",ensg,d[i,enst_cn]))
	return(NA)
    } else if (! is.na(d[i,enst_cn])) {
        print(sprintf("Overwriting transcript of %s from %s to %s",ensg,d[i,enst_cn],enst))
    } else {
        print(sprintf("Manually assigning transcript %s to gene %s",enst,ensg))
    }
    d[i,enst_cn] = enst
    return(d)
}



# Read data from Lambert18 supplementary. Skip first line which are header categories - see raw file for additional header info
lib = read.table("Lambert18_TableS1.csv", sep=";", quote="\"", dec=".", header=T, fill=T, skip=1, stringsAsFactors=F)

# Fix described below
i = which(lib$ID=="ENSG00000262156")
stopifnot(length(i)==1)
lib[i,"ID"] = "ENSG00000128383"
print(sprintf("Changing gene ENSG00000262156 index %d to ENSG00000128383 which should be the same protein but a more cannonical entry listed in MANE Select",i))

# entries that are not ensembles genes?
i_rm = which(substr(lib$ID, 1, 4) != "ENSG")
if (length(i_rm) > 0) {
    print(sprintf("WARNING: %d entries are not Ensembl gene ID's:",length(i_rm)))
    print(lib[i_rm,1:5])
}

print("Proteins curated to be transcription factors:")
print(table(lib$X))

# read MANE transcripts
# Downloaded May 9, 2022 from ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/mane/MANE_v1.0/
gzfile = gzfile("MANE.GRCh38.v1.0.transcripts_by_gene.tsv.gz","rt")
mane = read.table(gzfile, header=T)
close(gzfile)

# use MANE Select transcripts for genes
mane$ensg = sapply(strsplit(mane[,1],".",fixed=T), "[[", 1)
lib$transcript.version = mane[match(lib$ID,mane$ensg),"MANE_Select_Ensembl_id"]
print(sprintf("Found transcripts for %d of %d ensemle genes", sum(! is.na(lib$transcript.version)), nrow(lib)))

print("Two transcripts in MANE list seems not to excist - update")
# 727  ENSG00000251369         ZNF550     C2H2 ZF Yes  ENST00000699333.1 :  update to ENST00000457177.5 Q7Z398
# 1465 ENSG00000178928          TPRX1 Homeodomain Yes  ENST00000698655.1 :  update to ENST00000535759.1 D2CFI5  (shorter uniprot Q8N7U7 identical)
lib = fix_transcript(lib, "ENSG00000251369","ENST00000457177.5", dont_overwrite=F)
lib = fix_transcript(lib, "ENSG00000178928","ENST00000535759.1", dont_overwrite=F)

i_mis = which(is.na(lib$transcript.version))
print(sprintf("Of %d missing transcripts, %d are annotated as transcription factors", length(i_mis), sum(tolower(lib[i_mis,"X"])=="yes")))

i_mis_tf = i_mis[which(tolower(lib[i_mis,"X"])=="yes")]
# === Discarded ===
#                   ID           Name         DBD
# 68   ENSG00000250709 CCDC169-SOHLH2        bHLH  readthrough splice variant
# 161  ENSG00000267281     AC023509.3        bZIP  readthrough splice variant
# 233  ENSG00000267179     AC008770.3     C2H2 ZF  novel gene - Protein predicted
# 235  ENSG00000264668     AC138696.1     C2H2 ZF  novel gene - Protein predicted
# 572  ENSG00000249459        ZNF286B     C2H2 ZF  pseudogene - Protein uncertain
# 853  ENSG00000214534        ZNF705E     C2H2 ZF  pseudogene 
# 913  ENSG00000214189         ZNF788     C2H2 ZF  pseudogene - Protein uncertain
# 959  ENSG00000228623         ZNF883     C2H2 ZF  pseudogene - P0CG24: Experimental evidence at protein level
# 1051 ENSG00000177946        CENPBD1       CENPB  pseudogene - Protein uncertain
# 1307      DUX1_HUMAN           DUX1 Homeodomain  O43812 RefSeq NM_012146 - 89% identical to DUX4 (ENSG00000260596 ENST00000565211.1 NM_001306068.3)
# 1308      DUX3_HUMAN           DUX3 Homeodomain  Q96PT4 
# 1518 ENSG00000064489   BORCS8-MEF2B    MADS box  readthrough - Protein predicted
# 2150 ENSG00000230257           NFE4     Unknown  lncRNA (no protein) - Q86UQ8 NM_001085386.1

# === Manually mapped to ensembl cannonical transcripts ===
#                   ID           Name         DBD
# 232  ENSG00000100219           XBP1        bZIP  ENST00000344347.5   P17861-2
# 298  ENSG00000118922          KLF12     C2H2 ZF  ENST00000377669.7   Q9Y4X4       In MANE v0.95 but not in v1.0
# 472  ENSG00000152926         ZNF117     C2H2 ZF  ENST00000282869.11  Q03924       Experimental evidence at transcript level
# 500  ENSG00000188629         ZNF177     C2H2 ZF  ENST00000589262.5   Q13360
# 774  ENSG00000167962         ZNF598     C2H2 ZF  ENST00000562103.2   H3BPG6       Ensemble matched to Q86UK7
# 906  ENSG00000196381         ZNF781     C2H2 ZF  ENST00000358582.9   Q8N8C0-2     Experimental evidence at transcript level, In MANE v0.95 but not in v1.0
# 937  ENSG00000167766          ZNF83     C2H2 ZF  ENST00000597597.1   P51522       Experimental evidence at transcript level
# 952  ENSG00000178917         ZNF852     C2H2 ZF  ENST00000436261.6   Q6ZMS4
# 1183 ENSG00000204060          FOXO6    Forkhead  ENST00000641094.2 A8MYZ6
# 1203 ENSG00000220201          ZGLP1        GATA  ENST00000403903.7 P0C6A0   Experimental evidence at transcript level
# 1230 ENSG00000118418          HMGN3     HMG/Sox  ENST00000344726.9 Q15651
# 2159 ENSG00000186416           NKRF     Unknown  ENST00000304449.8   O15226       Ensembl cannonical is ENST00000688521.1 with 783 aa but not matched to uniprot. Golden is isoform 2
# 2727 ENSG00000174197            MGA       T-box  ENST00000219905.13  Q8IWI9-4     In MANE v0.95 but not in v1.0

lib = fix_transcript(lib,"ENSG00000100219","ENST00000344347.5")
lib = fix_transcript(lib,"ENSG00000118922","ENST00000377669.7")
lib = fix_transcript(lib,"ENSG00000152926","ENST00000282869.11")
lib = fix_transcript(lib,"ENSG00000188629","ENST00000589262.5")
lib = fix_transcript(lib,"ENSG00000167962","ENST00000562103.2")
lib = fix_transcript(lib,"ENSG00000196381","ENST00000358582.9")
lib = fix_transcript(lib,"ENSG00000167766","ENST00000597597.1")
lib = fix_transcript(lib,"ENSG00000178917","ENST00000436261.6")
lib = fix_transcript(lib,"ENSG00000204060","ENST00000641094.2")
lib = fix_transcript(lib,"ENSG00000220201","ENST00000403903.7")
lib = fix_transcript(lib,"ENSG00000118418","ENST00000344726.9")
lib = fix_transcript(lib,"ENSG00000186416","ENST00000304449.8")
lib = fix_transcript(lib,"ENSG00000174197","ENST00000219905.13")

i_mis_ntf = i_mis[which(tolower(lib[i_mis,"X"])=="no")]
# === Discarded ===
#                   ID       Name      DBD         ensembl            uniprot
# 611  ENSG00000168122    ZNF355P  C2H2 ZF           pseudogene  Protein uncertain
# 720  ENSG00000240225    ZNF542P  C2H2 ZF           pseudogene  Protein uncertain
# 873      ZNF73_HUMAN      ZNF73  C2H2 ZF                   NA            Deleted
# 929  ENSG00000224689    ZNF812P  C2H2 ZF           pseudogene                 NA
# 988  ENSG00000267908   ZSCAN5DP  C2H2 ZF           pseudogene  Protein uncertain
# 1048 ENSG00000212643      ZRSR1  CCCH ZF           pseudogene      Not in human?
# 1154 ENSG00000204828    FOXD4L2 Forkhead      no longer in db                 NA   May map to FOXD4L4 (ENSG00000184659) but that's already in the list index 1154
# 1216 ENSG00000258724 AC105001.2  HMG/Sox        not in RefSeq                 NA   Ensembl map to PINX1 but that seems not to ba a TF in uniprot. Match in uniprot is a predicted protein
# 1897 ENSG00000204677    FAM153C  Unknown           pseudogene                 NA
# 2398 ENSG00000176700    SCAND2P  Unknown           pseudogene  Protein uncertain
# 2431 ENSG00000135502   SLC26A10  Unknown           pseudogene  evidence at transcript level
# 2453 ENSG00000214338      SOGA3  Unknown                                           Ensembl canonical: ENST00000525778.5 (ccds and uniprot match but no evidence at protein level "Protein inferred from homology")
# 2617 ENSG00000228970     UBTFL6  Unknown pseudogene Protein uncertain
# 2648 ENSG00000188707    ZBED6CL  Unknown pseudogene Protein uncertain
# 2655 ENSG00000123870    ZNF137P  Unknown pseudogene Protein uncertain
#
# === Manually mapped to ensembl cannonical transcripts ===
#                   ID       Name      DBD         ensembl            uniprot
# 1238 ENSG00000163939      PBRM1  HMG/Sox  ENST00000296302.11  Q86U86          Poly-protein with many isoforms and poor match between Ensembl, uniprot and RefSeq. Golden is different (unip isoform4)
# 1740 ENSG00000118412   CASP8AP2  Unknown  ENST00000551025.4   Q9UKL3          Ensembl linked to A0A087WTW5 which is 100% identical to Q9UKL3 which is not linked to ensembl. 
# 1980 ENSG00000196101   HLA-DRB3  Unknown  ENST00000307137.11  P79483          Alternative gene, Not a Primary Assembly Gene. Seems not to excist an non-alternative verision?
# 2001 ENSG00000211899       IGHM  Unknown  ENST00000637539.2   P01871          IG C gene
# 2045 ENSG00000198083   KRTAP9-9  Unknown  ENST00000394008.1,  Q9BYP9-3        Experimental evidence at transcript level - otherwise well annotated, could consider?!
# 2153 ENSG00000167604     NFKBID  Unknown  ENST00000606253.5   Q8NI38-2        Transcript ENST00000396901.5 is isoform 1
# 2225 ENSG00000006576      PHTF2  Unknown  ENST00000416283.6   Q8N3S3          Experimental evidence at transcript level - otherwise well annotated, could consider?!
# 2262 ENSG00000122008       POLK  Unknown  ENST00000241436.8   Q9UBT6
# 2268 ENSG00000181222     POLR2A  Unknown  ENST00000674977.2   P24928-2        Ensembl mapper til A0A6Q8PGB0 som er 100% identisk med P24928-2
# 2388 ENSG00000079102    RUNX1T1  Unknown  ENST00000436581.6   A0A0A0MSU1      Uniprot map er længere end Q06455 som er det kanoniske protein (transcript ENST00000613302.4)
# 2390 ENSG00000163602       RYBP  Unknown  ENST00000477973.5   Q8N488
# 2461 ENSG00000157216      SSBP3  Unknown  ENST00000610401.5   Q9BWW4
#
# === Manually mapped to non-alternative ensembl entry ===
#                   ID       Name      DBD         ensembl            uniprot
# 1674 ENSG00000262156   APOBEC3A  Unknown                      P31941          Alternative ENSG00000128383, Not a Primary Assembly Gene. Ensembl canonical: ENST00000623492.3, golden NA
#
# === Methods text ===
# Protein coding sequences were obtained for all transcription factor genes via the MANE Select list of transcripts [Morales22].
# For 13 ensembl genes that were not matched to RefSeq, i.e. not in the MANE Select list, we used the Ensembl canonical
# transcript if this was matched to a uniprot entry with experimental evidence at protein level or in a few cases on
# transcript level. One Ensembl gene, entry ENSG00000262156, were changed to the more canonical version ENSG00000128383
# that refers to the same protein (uniprot P31941) and is present in the MANE select list. Two non-excisting transcripts
# in the MANE select list were updated (ENST00000699333 to ENST00000457177 and ENST00000698655 to ENST00000535759).
# The remaining 11 transcription factors without a match in the MANE Select list are either depreciated entries,
# readthrough or pesudogenes with poor evidence and thus discarded. For two transcription factors without Ensembl entry,
# protein coding DNA sequences were retrieved from the European Nucleotide Archive. In total, we obtained wild-type DNA
# sequences for 1628 of 1639 transcription factors of which 1626 have unique amino acid sequences. 

lib = fix_transcript(lib,"ENSG00000163939","ENST00000296302.11")
lib = fix_transcript(lib,"ENSG00000118412","ENST00000551025.4")
lib = fix_transcript(lib,"ENSG00000196101","ENST00000307137.11")
lib = fix_transcript(lib,"ENSG00000211899","ENST00000637539.2")
lib = fix_transcript(lib,"ENSG00000198083","ENST00000394008.1")
lib = fix_transcript(lib,"ENSG00000167604","ENST00000606253.5")
lib = fix_transcript(lib,"ENSG00000006576","ENST00000416283.6")
lib = fix_transcript(lib,"ENSG00000122008","ENST00000241436.8")
lib = fix_transcript(lib,"ENSG00000181222","ENST00000674977.2")
lib = fix_transcript(lib,"ENSG00000079102","ENST00000436581.6")
# lib = fix_transcript(lib,"ENSG00000163602","ENST00000477973.5") - transcript seems broken with "NN" and many stops
lib = fix_transcript(lib,"ENSG00000157216","ENST00000610401.5")

print(sprintf("Missing transcripts now %d of %d, %d of %d curated to be TF",
              sum(is.na(lib$transcript.version)), nrow(lib), sum(is.na(lib$transcript.version) & tolower(lib$X)=="yes"), sum(tolower(lib$X)=="yes")))
# print(lib[i_mis_tf,c(1,2,3,4,38)])

# Assign coding sequence using transcript
# Downloaded May 9 2022 from ftp.ensembl.org/pub/current_fasta/homo_sapiens/cds
gzfile = gzfile("Homo_sapiens.GRCh38.cds.all.seq.gz","rt")
ens = read.table(gzfile, header=F)
close(gzfile)
colnames(ens) = c("transcript.version","dna")

# match transcript ID's without version number
lib$transcript = sapply(strsplit(lib$transcript.version,".",fixed=T), "[[", 1)
ens$transcript = sapply(strsplit(ens$transcript.version,".",fixed=T), "[[", 1)

lib$dna = ens[match(lib$transcript,ens$transcript),"dna"]
print(sprintf("Found sequences for %d of %d transcripts", sum(! is.na(lib$dna)), sum(! is.na(lib$transcript))))

i_mis_seq = which(is.na(lib$dna) & tolower(lib$X)=="yes")
print(sprintf("Of %d missing transcripts/sequences, %d are annotated as transcription factors:", sum(is.na(lib$dna)), length(i_mis_seq)))
print(lib[i_mis_seq,c(1,2,3,4,38,40)])

# >ENA|CAA04776|CAA04776.1 Homo sapiens (human) hypothetical protein
# ATGGCCCTCCTGACAGCTTTGGACGACACCCTCCCCGAGGAAGCCCAGGGACCGGGAAGG
# CGAATGATACTCCTTTCGACCCCGAGTCAAAGTGATGCCCTGCGAGCCTGCTTTGAGCGG
# AACCTGTACCCGGGCATTGCCACCAAAGAAGAGCTGGCCCAGGGCATCGACATTCCGGAG
# CCCAGGGTCCAGATTTGGTTTCAGAATGAGAGATCATGCCAGTTGAGGCAGCACCGGCGG
# CAATCTCGGCCCTGGCCCGGGAGACGTGACCCGCAAAAAGGCAGACGAAAGCGGACTGCC
# ATCACCGGATCCCAAACCGCCCTGCTCCTCCGAGCCTTTGAGAAGGATCGCTTTCCAGGC
# ATTGCTGCCAGGGAAGAGCTGGCCAGAGAGACGGGCCTCCCGGAGTCCAGGATTCAGATC
# TGGTTTCAGAATCGAAGAGCCAGGCACCGGGGACAGTCTGGCAGGGCGCCCACGCAGGCA
# AGCATCCGGTGCAATGCAGCCCCAATTGGGTGA
# >translation
# MALLTALDDTLPEEAQGPGRRMILLSTPSQSDALRACFERNLYPGIATKEELAQGIDIPEPRVQIWFQNERSCQLRQHRRQSRPWPGRRDPQKGRRKRTAITGSQTALLLRAFEKDRFPGIAAREELARETGLPESRIQIWFQNRRARHRGQSGRAPTQASIRCNAAPIG-
# >sp|O43812|DUX1_HUMAN Double homeobox protein 1 OS=Homo sapiens OX=9606 GN=DUX1 PE=1 SV=1
# MALLTALDDTLPEEAQGPGRRMILLSTPSQSDALRACFERNLYPGIATKEELAQGIDIPEPRVQIWFQNERSCQLRQHRRQSRPWPGRRDPQKGRRKRTAITGSQTALLLRAFEKDRFPGIAAREELARETGLPESRIQIWFQNRRARHRGQSGRAPTQASIRCNAAPIG

# Use coding sequence from ENA
print("Manually assign DNA sequence to DUX1_HUMAN not listed in Ensembl")
lib[which(lib$ID=="DUX1_HUMAN"),"transcript"] = "manual_dux1"
lib[which(lib$ID=="DUX1_HUMAN"),"dna"] = "ATGGCCCTCCTGACAGCTTTGGACGACACCCTCCCCGAGGAAGCCCAGGGACCGGGAAGGCGAATGATACTCCTTTCGACCCCGAGTCAAAGTGATGCCCTGCGAGCCTGCTTTGAGCGGAACCTGTACCCGGGCATTGCCACCAAAGAAGAGCTGGCCCAGGGCATCGACATTCCGGAGCCCAGGGTCCAGATTTGGTTTCAGAATGAGAGATCATGCCAGTTGAGGCAGCACCGGCGGCAATCTCGGCCCTGGCCCGGGAGACGTGACCCGCAAAAAGGCAGACGAAAGCGGACTGCCATCACCGGATCCCAAACCGCCCTGCTCCTCCGAGCCTTTGAGAAGGATCGCTTTCCAGGCATTGCTGCCAGGGAAGAGCTGGCCAGAGAGACGGGCCTCCCGGAGTCCAGGATTCAGATCTGGTTTCAGAATCGAAGAGCCAGGCACCGGGGACAGTCTGGCAGGGCGCCCACGCAGGCAAGCATCCGGTGCAATGCAGCCCCAATTGGGTGA"

# >ENA|AAL02242|AAL02242.1 Homo sapiens (human) homeobox protein DUX3
# ATGCCGGCTGAGGTGCACGGGAGCCCGCCCGCCTCTCTCTGCCCGTGTCCGTCCGTGAAA
# TTCCGGCCGGGGCTCCCTGCGATGGCCCTCCTGACAGCTTTGGACGACACCCTCCCCGAG
# GAAGCCCAGGGACCGGGAAGGCGAATGATACTCCTTTCGACCCCGAGTCAAAGCGATGCC
# CTGCGAGCCTGCTTTGAGCGGAACCTGTACCCGGGCATTGCCACCAAAGAACAGCTGGCC
# CAGGGCATCGACATTCCGGAGCCCAGGGTCCAGATTTGGTTTCAGAATGAGAGATCATGC
# CAGTTGAGGCAGCACCGGCGGCAATCTCGGCCCTGGCCCGGGAGACGCGACCCGCAAAAA
# GGCAGACGAAAGCGGACTGCCATCACCGGATCCCAAACCGCCCTGCTCCTCCGAGCCTTT
# GAGAAGGATCGCTTTCCAGGCATTCCTGCCAGGGAAGAGCTGGCCAGAGAGACGGGCCTC
# CCGGAGTCCAGGATTCAGCTCTGGTTTCAGAATCGAAGAGCCAGGCACTGGGGACAGTCT
# GGCAGGGCGCCCACGCAGGCAAGCATCCGGTGCAATGCAGCCCCAATTGGGTGA
# >translation
# MPAEVHGSPPASLCPCPSVKFRPGLPAMALLTALDDTLPEEAQGPGRRMILLSTPSQSDALRACFERNLYPGIATKEQLAQGIDIPEPRVQIWFQNERSCQLRQHRRQSRPWPGRRDPQKGRRKRTAITGSQTALLLRAFEKDRFPGIPAREELARETGLPESRIQLWFQNRRARHWGQSGRAPTQASIRCNAAPIG-
# >sp|Q96PT4|DUX3_HUMAN Putative double homeobox protein 3 OS=Homo sapiens OX=9606 GN=DUX3 PE=2 SV=1
# MPAEVHGSPPASLCPCPSVKFRPGLPAMALLTALDDTLPEEAQGPGRRMILLSTPSQSDALRACFERNLYPGIATKEQLAQGIDIPEPRVQIWFQNERSCQLRQHRRQSRPWPGRRDPQKGRRKRTAITGSQTALLLRAFEKDRFPGIPAREELARETGLPESRIQLWFQNRRARHWGQSGRAPTQASIRCNAAPIG

# Use coding sequence from ENA
print("Manually assign DNA sequence to DUX3_HUMAN not listed in Ensembl")
lib[which(lib$ID=="DUX3_HUMAN"),"transcript"] = "manual_dux3"
lib[which(lib$ID=="DUX3_HUMAN"),"dna"] = "ATGCCGGCTGAGGTGCACGGGAGCCCGCCCGCCTCTCTCTGCCCGTGTCCGTCCGTGAAATTCCGGCCGGGGCTCCCTGCGATGGCCCTCCTGACAGCTTTGGACGACACCCTCCCCGAGGAAGCCCAGGGACCGGGAAGGCGAATGATACTCCTTTCGACCCCGAGTCAAAGCGATGCCCTGCGAGCCTGCTTTGAGCGGAACCTGTACCCGGGCATTGCCACCAAAGAACAGCTGGCCCAGGGCATCGACATTCCGGAGCCCAGGGTCCAGATTTGGTTTCAGAATGAGAGATCATGCCAGTTGAGGCAGCACCGGCGGCAATCTCGGCCCTGGCCCGGGAGACGCGACCCGCAAAAAGGCAGACGAAAGCGGACTGCCATCACCGGATCCCAAACCGCCCTGCTCCTCCGAGCCTTTGAGAAGGATCGCTTTCCAGGCATTCCTGCCAGGGAAGAGCTGGCCAGAGAGACGGGCCTCCCGGAGTCCAGGATTCAGCTCTGGTTTCAGAATCGAAGAGCCAGGCACTGGGGACAGTCTGGCAGGGCGCCCACGCAGGCAAGCATCCGGTGCAATGCAGCCCCAATTGGGTGA"

print(sprintf("Missing sequences now %d of %d, %d TF", sum(is.na(lib$dna)), nrow(lib), sum(is.na(lib$dna) & tolower(lib$X)=="yes")))

# fix some gene names 
lib[which(lib$Name=="T"),"Name"] = "TBXT"
lib[which(lib$Name=="AC092835.1"),"Name"] = "ZN892"
stopifnot( length(lib$Name) == length(unique(lib$Name)) )

# check sequences
i_dna = which(! is.na(lib$dna))
# check 
stopifnot( all( sapply(lib[i_dna,"dna"], function(s){ nchar(s)%%3==0 })))
print("Letters used in DNA sequences")
print(table(unlist(strsplit(lib[i_dna,"dna"],""))))

# Translate to amino acid and truncate the tailing stop codon
lib$aa = NA

# stop codons are in the synthesized tile lib so keeping them here makes nucl and aa numbering more consistent
lib[i_dna,"aa"] = dna2aa(lib[i_dna,"dna"], truncate=0)
# lib[i_dna,"aa"] = dna2aa(lib[i_dna,"dna"], truncate=3)

print("Amino acids used in library")
print(table(unlist(strsplit(lib[i_dna,"aa"],""))))
print("Length distribution in amino acids")
print(summary(nchar(lib[i_dna,"aa"])))

# Annotate redundant amino acid sequences
i_aa_redundant = c()
if (length(lib[i_dna,"aa"]) != length(unique(lib[i_dna,"aa"]))) {
    t = table(lib[i_dna,"aa"])
    print(sprintf("There are %d fully redundant amino acid sequences where only the first will be used",length(t[t>1])))
    for (seq in names(t[t>1])) {
        i_redun = which(lib$aa == seq)
	dna_identical = c(lib[i_redun[1],"dna"] == lib[i_redun,"dna"])
	print(cbind(lib[i_redun,c("ID","Name","X","transcript.version")], dna_identical))
	i_aa_redundant = c(i_aa_redundant, i_redun[2:length(i_redun)])
    }
}
lib$aa_unique = NA
lib[i_dna,"aa_unique"] = TRUE
lib[i_aa_redundant,"aa_unique"] = FALSE


tile_seq = function(sequence, tile_len, tile_stride, make_ct_tile=T) {
    n = nchar(sequence)
    if (n < tile_len) { return(NULL) }
    n_tiles = (n-tile_len) %/% tile_stride +1
    first = seq(0,n_tiles-1)*tile_stride+1
    last = first+tile_len-1
    # if whole tiles cannot cover the sequence, make a terminal tile if requested
    if (make_ct_tile & (n_tiles-1)*tile_stride+tile_len < n) {
        first = c(first, n-tile_len+1)
	last = c(last, n)
    }
    tiles = substring(sequence, first, last)
    
    name = names(sequence)
    if (is.null(name)) {
        names(tiles) = paste(seq_along(first), first, last, sep=".")
    } else {
        names(tiles) = paste(name, seq_along(first), first, last, sep=".")
    }
    return(tiles)
}

# TF's to use
ilib_use = which(tolower(lib$X)=="yes" & ! is.na(lib$dna) & lib$aa_unique)

# Check uniqueness of sequences, ensemble id's and gene names
stopifnot( length(ilib_use) == length(unique(lib[ilib_use,"dna"])) )
stopifnot( length(ilib_use) == length(unique(lib[ilib_use,"aa"])) )
stopifnot( length(ilib_use) == length(unique(lib[ilib_use,"ID"])) )
stopifnot( length(ilib_use) == length(unique(lib[ilib_use,"Name"])) )

# make a new vector of sequences that can have names to be returned by 'lapply' and expanded by 'unlist'
tf_seq = lib[ilib_use,"dna"]
names(tf_seq) = lib[ilib_use,"transcript"]
tiles = unlist(lapply(tf_seq, tile_seq, 90, 45))

n_s_unq = length(unique(tf_seq))
n_t_unq = length(unique(tiles))
print(sprintf("Sliced %d sequences (%d or %.0f%% unique) into %d tiles (%d or %.0f%% unique)",
              length(tf_seq), n_s_unq, n_s_unq/length(tf_seq)*100,
	      length(tiles),   n_t_unq, n_t_unq/length(tiles)*100))

name_split_list = strsplit(names(tiles), ".", fixed=T)
tile_lib = data.frame(meta = names(tiles),
                      name = NA,
                      transcript = sapply(name_split_list, '[[', 1),
                      ensembl = NA,
		      gene = NA,
                      tile_number = as.numeric( sapply(name_split_list, '[[', 2) ),
                      first_resi_dna = as.numeric( sapply(name_split_list, '[[', 3) ), 
                      last_resi_dna = as.numeric( sapply(name_split_list, '[[', 4) ),
		      dna = tiles)
i_l2t = match(tile_lib$transcript, lib$transcript)
tile_lib$ensembl = lib[i_l2t,"ID"]
tile_lib$gene = lib[i_l2t,"Name"]
tile_lib$name = sprintf("%s_tile%04d", tile_lib$ensembl, tile_lib$tile_number)
rownames(tile_lib) = NULL

# Find tiles that are redundant on DNA level
tile_lib$copyof_dna = NA
if (length(tile_lib$dna) != length(unique(tile_lib$dna))) { 
    t = table(tile_lib$dna)
    print(sprintf("There are %d DNA redundant tiles",length(t[t>1])))
    for (seq in names(t[t>1])) {
        i_redun = which(tile_lib$dna == seq)
	# dna_identical = c(lib[i_redun[1],"dna"] == lib[i_redun,"dna"])
	# print(cbind(lib[i_redun,c("ID","Name","X","transcript.version")], dna_identical))
	tile_lib[i_redun[2:length(i_redun)],"copyof_dna"] = as.character(i_redun[1])
	print(tile_lib[i_redun,c("name","transcript","gene","first_resi_dna","last_resi_dna","copyof_dna")])
    }
}
print("How many copied tiles per protein:")
agg = aggregate(tile_lib$copyof_dna, by=list(tile_lib$ensembl), function(v){sum(! is.na(v))})
protein_lib = data.frame(ensembl=agg[,1], tile_copies=agg[,2])
t = table(tile_lib$ensembl)
protein_lib$n_tiles = t[match(protein_lib$ensembl, names(t))]
i = which(protein_lib$tile_copies > 0)
print(protein_lib[i[order(protein_lib[i,"tile_copies"])],])

# Make libraries of even, odd and CT fragments
print("Split library into 3")
tile_lib$lib = NA
tile_lib[which(tile_lib$tile_number %% 2 == 0),"lib"] = 1
tile_lib[which(tile_lib$tile_number %% 2 == 1),"lib"] = 2
tile_lib[which(diff(tile_lib$tile_number) != 1),"lib"] = 3
tile_lib[length(tile_lib$tile_number),"lib"] = 3  # last fragment is also CT
print(sprintf("Library size %d split into lib1_even (%d), lib2_odd (%d), and lib3_ct (%d)",
              length(tile_lib$lib), sum(tile_lib$lib==1), sum(tile_lib$lib==2), sum(tile_lib$lib==3)))
stopifnot(length(which(is.na(tile_lib$lib))) == 0)
print(table(tile_lib$lib))

i_sub1 = which(tile_lib$lib == 1)
i_sub2 = which(tile_lib$lib == 2)
i_sub3 = which(tile_lib$lib == 3)

# Adapters, nt is 5 prime and ct is 3 prime
print("Append adapters and on sequences")
adapt = list()

adapt[['all']]['nt'] = "AAGACGCGTTCTAGAGGCAGCGGAGCCACC"
adapt[['all']]['ct'] = "TAGTAACTTAAGAATTCACCGGTCTGACCT"
tile_lib$dna_adapt = paste0(adapt[['all']]['nt'], tile_lib$dna, adapt[['all']]['ct'])


# Amino acid sequences of tiles
tile_lib$first_resi = (tile_lib$first_resi_dna-1)/3 +1
tile_lib$last_resi = tile_lib$last_resi_dna/3
tile_lib$aa = dna2aa(tile_lib$dna)

# CT tiles have a stop codon which is included in this numbering
stopifnot(all( nchar(tile_lib$aa) == tile_lib$last_resi-tile_lib$first_resi+1 ))
print("tile amino acid lengths:")
print(table(nchar(tile_lib$aa)))

# Check that AA sequences from full length protein and residue numbers are identical
tile_lib_aa2 = substr(lib[match(tile_lib$transcript,lib$transcript),"aa"], tile_lib$first_resi, tile_lib$last_resi)
stopifnot(all( tile_lib_aa2 == tile_lib$aa ))

# Check that the last residue of CT tiles are the last residues of full length sequences
tile_lib$nres = nchar(lib$aa)[match(tile_lib$transcript,lib$transcript)]
i_ct = which(tile_lib$lib==3)
stopifnot(all(  tile_lib[i_ct,"last_resi"] == tile_lib[i_ct,"nres"] ))

# Find tiles that are redundant on AA level
tile_lib$copyof_aa = NA
t = table(tile_lib$aa)
for (seq in names(t[t>1])) {
    i_redun = which(tile_lib$aa == seq)
    # reference tile must also be reference of DNA redundancy 
    stopifnot(is.na( tile_lib[i_redun[1],"copyof_dna"] ))
    # All AA copies gets the same reference tile - use this because copying DNA ref is already in the copyof_dna column
    tile_lib[i_redun[2:length(i_redun)],"copyof_aa"] = as.character(i_redun[1])
    print(tile_lib[i_redun,c("name","transcript","gene","first_resi","last_resi","copyof_dna","copyof_aa")])
}
print(sprintf("Found %d AA redundant tiles of which %d are not DNA redundant",
              sum(! is.na(tile_lib$copyof_aa)), sum(is.na(tile_lib$copyof_dna) & ! is.na(tile_lib$copyof_aa)) ))


# Output 
# Dump everything as R file
print("Dump R file transfac90.rda")
save(lib, tile_lib, adapt, file="library.rda")

# Dump tiles of individual Excel files for synthesis
print("Dump an Excel-readable csv file per library transfac90_lib*.csv")
cn = c("ensembl","tile_number","lib","dna","dna_adapt")
i1 = i_sub1[which(is.na(tile_lib[i_sub1,"copyof_dna"]))]
i2 = i_sub2[which(is.na(tile_lib[i_sub2,"copyof_dna"]))]
i3 = i_sub3[which(is.na(tile_lib[i_sub3,"copyof_dna"]))]
print(sprintf("Using %d of %d (lib1), %d of %d (lib2) and %d of %d (lib3) unique tiles",
              length(i1),length(i_sub1),length(i2),length(i_sub2),length(i3),length(i_sub3)))
write.csv2(tile_lib[i1,cn], file="transfac90_lib1_even.csv")
write.csv2(tile_lib[i2,cn], file="transfac90_lib2_odd.csv")
write.csv2(tile_lib[i3,cn], file="transfac90_lib3_ct.csv")

# Dump tiles as CSV
print(sprintf("Writing all %d tiles (%d DNA unique and %d AA unique) to transfac_tiles.csv",
              nrow(tile_lib), sum(is.na(tile_lib$copyof_dna)), sum(is.na(tile_lib$copyof_aa))))
write.csv(tile_lib[,c("name","gene","first_resi","last_resi","copyof_dna","copyof_aa","lib","aa","dna")], file="transfac_tiles.csv", row.names=F, quote=F)

# Dump proteins
print(sprintf("Writing full sequences for %d of %d transcription factors to files - in total %d residues (original file)",
              length(ilib_use),sum(tolower(lib$X)=="yes"), sum(nchar(lib[ilib_use,"aa"]))))

# CSV file
write.csv(lib[ilib_use,c("ID","Name","DBD","TF.assessment","Binding.mode","Motif.status","aa")], file="transfac_proteins.csv", row.names=F, quote=F)

# FASTA file without stop codons
aa_nostop = sub("*", "", lib[ilib_use,"aa"], fixed=T)
fasta_file = file("transfac_proteins.fasta")
write(sprintf(">%s\n%s",lib[ilib_use,"Name"],aa_nostop), fasta_file)
close(fasta_file)

# seq file for easy comparison to previous version
fh = file("transfac_dna.seq", "wt")
write(sprintf("%20s  %s", lib[ilib_use,"ID"], lib[ilib_use,"dna"]), fh)
close(fh)

