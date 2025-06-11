options(width=160)

tfp = read.csv(gzfile("../annotations/transfac_proteins.csv.gz"))
tfp$nres = nchar(tfp$aa)

tft = read.csv(gzfile("../annotations/transfac_tiles.csv.gz"))
tft$resi = ceiling((tft$first_resi+tft$last_resi)/2)

runavg = function(x, n = 3) { as.vector(filter(x, rep(1 / n, n), sides = 2)) }

plot_profile = function(genename) {
    ip = which(tfp$gene == genename)
    it = which(tft$gene == genename)
    stopifnot(length(ip)==1)
    adp = as.numeric(strsplit(tfp[ip,"adpred"],";")[[1]])
    stopifnot( length(adp) == tfp[ip,"nres"] )
    plddt = as.numeric(strsplit(tfp[ip,"plddt"],";")[[1]])
    stopifnot( length(plddt) == tfp[ip,"nres"] )
    pap = as.numeric(strsplit(tfp[ip,"pap"],";")[[1]])
    stopifnot( length(pap) == tfp[ip,"nres"] )
    # h25 = as.numeric(strsplit(tfp[ip,"h25"],";")[[1]])
    # stopifnot( length(h25) == tfp[ip,"nres"] )

    plot(0,0,col=0, xlim=c(1,tfp[ip,"nres"]), ylim=c(0,2.5), xlab="Residue", ylab="Score", main=genename)
    par(bg="white")
    lines(runavg(adp,5), lwd=2, col=1)
    lines(runavg(plddt/100.0,5), lwd=2, col=4)
    # lines(runavg(h25,5), lwd=1, col=4)
    lines(runavg(-log10(pap),5), lwd=2, col=2)
    
    # points(tft[it,"resi"], -log10(tft[it,"abundance_score"]*2), pch=20, col=3)
    screen = tft[it,"abundance_score"]*2
    cap = 0.01
    abline(h = -log10(cap), lty=2,col="grey")
    screen[which(screen<cap)] = cap
    points(tft[it,"resi"], -log10(screen), pch=20, col=3)
    
    # i_papt = match(tft[it,"name"],papt$name)
    # stopifnot(all( tft[it,"aa"] == papt[i_papt,"aa"] ))
    # # points(tft[it,"resi"], -log10(papt[i_papt,"score"]), pch=20, col=2)

    x = tft[it,"pap"]
    x[which(x<cap)] = cap
    points(tft[it,"resi"], -log10(x), pch=20, col=2)
    
    legend("topleft", c("ADPred","pLDDT","-log10 PAP","-log10 PAP tile (incl. C-deg)","-log10(screen abundance x2)"),
           pch=c(NA,NA,NA,20,20), lty=c(1,1,1,NA,NA), lwd=2, col=c(1,4,2,2,3), ncol=2)
    print(tft[it,c("gene","first_resi","last_resi","aa","abundance_score","abundance_std","pap_ct","pap")])
}

# # Plot profiles for all proteins!
# dir.create("profiles")
# for (genename in tfp$gene) {
#     fn = sprintf("profiles/profile_%s.pdf",genename)
#     pdf(fn, width=7, height=3.5, pointsize=7)
#     # quartz(width=12, height=6)
#     # par(bg="white")
#     plot_profile(genename)
#     dev.off()
#     # quartz.save(fn, type="png")
# }


# numbers for paper
n = sum(! is.na(tft$degron_score))
n05 = sum(tft$degron_score > 0.5, na.rm=T)
print(sprintf("Of %d TF tiles, %d (%.2f%%) have degron scores of which %d (%.2f%%) have score > 0.5",
              nrow(tft), n, n/nrow(tft)*100, n05, n05/n*100))
ids = which((! is.na(tft$degron_score)) & (! is.na(tft$tile_rasa_median)))
stopifnot(all(! is.na(tft[ids,"tile_plddt_median"])))
print(sprintf("Of %d TF tiles, %d (%.2f%%) have degron scores and median rASA/pLDDT scores",
              nrow(tft), length(ids), length(ids)/nrow(tft)*100))
ii_deg = which(tft[ids,"degron_score"] > 0.5)
n07 = sum(tft[ids[ii_deg],"tile_rasa_median"] > 0.7 & tft[ids[ii_deg],"tile_plddt_median"] < 70)
print(sprintf("Of %d tiles with degron scores and struc. features, %d (%.2f%%) have degron score > 0.5 of which %d (%.2f%%) have rASA > 0.7 and pLDDT < 70",
              length(ids), length(ii_deg), length(ii_deg)/length(ids)*100, n07, n07/length(ii_deg)*100))
n_rasa  = sum(tft[ids[ii_deg],"tile_rasa_median"] > 0.7)
n_plddt = sum(tft[ids[ii_deg],"tile_plddt_median"] < 70)
print(sprintf("Of %d degrons, %d (%.2f%%) have rASA > 0.7 and %d (%.2f%%) have pLDDT < 70", length(ii_deg),
              n_rasa, n_rasa/length(ii_deg)*100, n_plddt, n_plddt/length(ii_deg)*100))


# Paper plot with ADpred, pLDDT and PAP
genename = "HSF1"
ip = which(tfp$gene == genename)
it = which(tft$gene == genename)
adp = as.numeric(strsplit(tfp[ip,"adpred"],";")[[1]])
plddt = as.numeric(strsplit(tfp[ip,"plddt"],";")[[1]])
pap = as.numeric(strsplit(tfp[ip,"pap"],";")[[1]])

pdf("profile_HSF1.pdf", width=7, height=3.5, pointsize=9)
plot(0,0,col=0, xlim=c(1,tfp[ip,"nres"]), ylim=c(0,1.3), xlab="Residue", ylab="Score", main=genename,
     xaxp=c(0,500,10))
rect( 15, 0, 120, 1.0, border=NA, col="lightgray")
rect(130, 0, 203, 1.0, border=NA, col="lightgray")
rect(412, 0, 420, 1.0, border=NA, col="lightgray")
lines(runavg(adp,5), lwd=2, col=1)
lines(runavg(plddt/100.0,5), lwd=2, col=4)
lines(runavg(pap,5), lwd=2, col=2)
legend("topright", c("ADPred","pLDDT","PAP"), lty=c(1,1,1), lwd=2, col=c(1,4,2), ncol=3)
dev.off()

pdf("profile_HSF1_log.pdf", width=7, height=3.5, pointsize=9)
plot(0,0,col=0, xlim=c(1,tfp[ip,"nres"]), ylim=c(0,2.0), xlab="Residue", ylab="Score", main=genename,
     xaxp=c(0,500,10))
rect( 15, 0, 120, 1.5, border=NA, col="lightgray")
rect(130, 0, 203, 1.5, border=NA, col="lightgray")
rect(412, 0, 420, 1.5, border=NA, col="lightgray")
lines(runavg(adp,5), lwd=2, col=1)
lines(runavg(plddt/100.0,5), lwd=2, col=4)
lines(runavg(-log10(pap),5), lwd=2, col=2)
legend("topright", c("ADPred","pLDDT","-log10 PAP"), lty=c(1,1,1), lwd=2, col=c(1,4,2), ncol=3)
dev.off()


# Paper plot with rASA, pLDDT and degron score
genename = "HSF1"
ip = which(tfp$gene == genename)
it = which(tft$gene == genename)
rasa =  as.numeric(strsplit(tfp[ip,"rasa"],";")[[1]])
plddt = as.numeric(strsplit(tfp[ip,"plddt"],";")[[1]])

pdf("profile_HSF1_screen.pdf", width=7, height=3.5, pointsize=9)
plot(0,0,col=0, xlim=c(1,tfp[ip,"nres"]), ylim=c(-.05,1.35), xlab="Residue", ylab="Score", main=genename,
     xaxp=c(0,500,10))
rect(198, -0.1, 217, 1.2, border=NA, col="lightgray")
lines(runavg(rasa,5), lwd=2, col=1)
lines(runavg(plddt/100.0,5), lwd=2, col=4)
lines(tft[it,"resi"], tft[it,"degron_score"], type="b", lwd=2, pch=16, cex=1.2, col=2)
legend("topright", c("rASA","pLDDT","Degron score"), lty=c(1,1,NA), pch=c(NA,NA,16), lwd=2, pt.cex=1.2, col=c(1,4,2), ncol=3)
dev.off()


genename = "ASCL1"
ip = which(tfp$gene == genename)
it = which(tft$gene == genename)
rasa =  as.numeric(strsplit(tfp[ip,"rasa"],";")[[1]])
plddt = as.numeric(strsplit(tfp[ip,"plddt"],";")[[1]])

pdf("profile_ASCL1_screen.pdf", width=7, height=3.5, pointsize=9)
plot(0,0,col=0, xlim=c(1,tfp[ip,"nres"]), ylim=c(-.05,1.35), xlab="Residue", ylab="Score", main=genename,
     xaxp=c(0,500,10))
rect(23, -0.1, 40, 1.2, border=NA, col="lightgray")
lines(runavg(rasa,5), lwd=2, col=1)
lines(runavg(plddt/100.0,5), lwd=2, col=4)
lines(tft[it,"resi"], tft[it,"degron_score"], type="b", lwd=2, pch=16, cex=1.2, col=2)
legend("topright", c("rASA","pLDDT","Degron score"), lty=c(1,1,NA), pch=c(NA,NA,16), lwd=2, pt.cex=1.2, col=c(1,4,2), ncol=3)
dev.off()


# Paper plot with pLDDT and PAP
genename = "CRX"
ip = which(tfp$gene == genename)
plddt = as.numeric(strsplit(tfp[ip,"plddt"],";")[[1]])
pap = as.numeric(strsplit(tfp[ip,"pap"],";")[[1]])

pdf("profile_CRX.pdf", width=7, height=2.5, pointsize=9)
par(mar=c(5,4,2,1)+.1)
plot(0,0,col=0, xlim=c(1,tfp[ip,"nres"]), ylim=c(0,1.3), xlab="Residue", ylab="Score", xaxp=c(0,500,10))
rect(40, 0, 95, 1.2, border=NA, col="lightgray")
lines(runavg(plddt/100.0,5), lwd=2, col=4)
lines(runavg(-log10(pap),5), lwd=2, col=1)
legend("topright", c("pLDDT","-log10 PAP"), lty=c(1,1), lwd=2, col=c(4,1), ncol=3)
dev.off()
