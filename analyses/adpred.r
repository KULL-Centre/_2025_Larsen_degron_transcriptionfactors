options(width=160)
runavg = function(x, n = 3) { as.vector(filter(x, rep(1 / n, n), sides = 2)) }

# TF library with scores
tf = read.csv(gzfile("../annotations/transfac_tiles.csv.gz"))

# add log abundance scores
tf$logabun_score = log10( tf[,"abundance_score"] )
tf$logabun_std   = 1/(log(10)*tf[,"abundance_score"]) * tf[,"abundance_std"]

apply2tiles = function(v, list_of_vec, func=mean) {
    gene = v[1]
    first_resi = as.numeric(v[2])
    last_resi = as.numeric(v[3])
    if (! gene %in% names(list_of_vec)) { return(NA) }
    if ( last_resi > length(list_of_vec[[gene]]) )  { return(NA) }
    return( func(list_of_vec[[gene]][first_resi:last_resi]) )
}

# Read ADPred predictions of full proteins based on PSIPRED features
adp_raw = read.table(gzfile("../annotations/adpred_ss_api_full_scores.txt.gz"))
colnames(adp_raw) = c("gene","ss","score")
i = which(adp_raw$gene=="AC092835.1")
adp_raw[i,"gene"] = "ZN892"
print(sprintf("Changed %d gene name from AC092835.1 to ZN892",length(i)))
i = which(adp_raw$gene=="T")
adp_raw[i,"gene"] = "TBXT"
print(sprintf("Changed %d gene name from T to TBXT",length(i)))

adp = strsplit(adp_raw$score, ";")
adp = lapply(adp, as.numeric)
names(adp) = adp_raw$gene
tf$adp_avg = apply(tf[,c("gene","first_resi","last_resi")], MARGIN=1, apply2tiles, list_of_vec=adp, func=mean)
tf$adp_max = apply(tf[,c("gene","first_resi","last_resi")], MARGIN=1, apply2tiles, list_of_vec=adp, func=max)
# when using this on CT tiles of length 29, it means that these are in practice added a CT Gly. We bravely move on
tf$adp_cen = apply(tf[,c("gene","first_resi","last_resi")], MARGIN=1, apply2tiles, list_of_vec=adp, func=function(v){ v[16] })

# Select trans-activating tiles
adpred_cut = 0.8
iavg = which(tf$adp_avg >= adpred_cut)
imax = which(tf$adp_max >= adpred_cut)
icen = which(tf$adp_cen >= adpred_cut)
print(sprintf("Tile with score > %.2f by mean %d, max %d and center %d", adpred_cut, length(iavg), length(imax), length(icen)))

h_adp_api = hist(unlist(adp), breaks=100, plot=F)

# Plot per residue scores compared to tile scores
quartz(width=5, height=5)
# par(mfrow=c(1,2), bg="white")
plot(0,0,col=0, xlim=c(0,1), ylim=c(0,2), xlab="ADPred", ylab="Density")
lines(h_adp_api$mids, h_adp_api$density, col=1, lwd=2)
# lines(h_adp_coil$mids, h_adp_coil$density, col=2, lwd=2)
ha_avg = hist(tf$adp_avg, breaks=100, plot=F)
lines(ha_avg$mids, ha_avg$density, col=3, lwd=2)
ha_max = hist(tf$adp_max, breaks=100, plot=F)
lines(ha_max$mids, ha_max$density, col=4, lwd=2)
ha_cen = hist(tf$adp_cen, breaks=100, plot=F)
lines(ha_cen$mids, ha_cen$density, col=5, lwd=2)
legend("top", c("PSIPRED all res","Coil all res","PSIPRED tile mean","PSIPRED tile max","PSIPRED tile center"), lty=1, lwd=2, col=seq(5))
quartz.save("adpred_distributions.png", type="png")


# Plot experimental score distributions of high ADPred scoring tiles
# use ADpred score of central position and look at different score distirbutions
quartz(width=12, height=4)
par(mfrow=c(1,3), mar=c(5,4,2,2)+.1, bg="white")
breaks = 60
n = 3
for (score in c("degron","abundance","logabun")) {
    cn = sprintf("%s_score",score)

    hs = hist(tf[,cn],     breaks=breaks, plot=F)
    ha = hist(tf[icen,cn], breaks=breaks, plot=F)

    plot(0,0,col=0, xlim=hs$mids[c(1,length(hs$mids))], ylim=c(0,max(ha$density,hs$density)), xlab=cn, ylab="Density")
    lines(ha$mids, runavg( ha$density, n=n), lwd=2, col=2)
    lines(hs$mids, runavg( hs$density, n=n), lwd=2, col=1)
    legend("topleft", c("All tiles",sprintf("Center >%.2f, %d tiles (%.1f%%)",adpred_cut,length(icen),length(icen)/nrow(tf)*100)), lty=1, col=c(1,2))
}
quartz.save("adpred_high_scoring.png", type="png")


# Paper plot
pdf("adpred_high_scoring.pdf", width=3.5, height=3.5, pointsize=7)
n=5
breaks=60
cn = "degron_score"
hs = hist(tf[,cn],     breaks=breaks, plot=F)
ha = hist(tf[icen,cn], breaks=breaks, plot=F)
par(mar=c(5,4,2,2)+.1)
plot(0,0,col=0, xlim=hs$mids[c(1,length(hs$mids))], ylim=c(0,max(ha$density,hs$density)), xlab=cn, ylab="Density")
lines(ha$mids, runavg( ha$density, n=n), lwd=2, col=2)
lines(hs$mids, runavg( hs$density, n=n), lwd=2, col=1)
legend("topleft", c("All tiles","TAD"), lty=1, col=c(1,2))
dev.off()


# Experimental score distributions of high ADPred scoring tiles
#   look at different scores and smoothing strategies
for (score in c("degron","abundance","logabun")) {
cn = sprintf("%s_score",score)

breaks = 50
bw = 2e-2
hs  = density(tf[,cn],     bw=bw, na.rm=T)
hha = density(tf[iavg,cn], bw=bw, na.rm=T)
hhm = density(tf[imax,cn], bw=bw, na.rm=T)
hhc = density(tf[icen,cn], bw=bw, na.rm=T)

quartz(width=12, height=4)
par(mfrow=c(1,3), bg="white")
breaks = seq(-.5, 1.5, 0.05)

plot(0,0,col=0, xlim=hs$x[c(1,length(hs$x))], ylim=c(0,max(hha$y,hs$y)), xlab=cn, ylab="Density")
lines(hha$x, hha$y, col=2)
lines(hs$x, hs$y, col=1)
legend("topleft", c("All tiles",sprintf("Avg. >%.2f, %d tiles (%.1f%%)",adpred_cut,length(iavg),length(iavg)/nrow(tf)*100)), lty=1, col=c(1,2))

plot(0,0,col=0, xlim=hs$x[c(1,length(hs$x))], ylim=c(0,max(hhm$y,hs$y)), xlab=cn, ylab="Density")
lines(hhm$x, hhm$y, col=2)
lines(hs$x, hs$y, col=1)
legend("topleft", c("All tiles",sprintf("Max >%.2f, %d tiles (%.1f%%)",adpred_cut,length(imax),length(imax)/nrow(tf)*100)), lty=1, col=c(1,2))

plot(0,0,col=0, xlim=hs$x[c(1,length(hs$x))], ylim=c(0,max(hhc$y,hs$y)), xlab=cn, ylab="Density")
lines(hhc$x, hhc$y, col=2)
lines(hs$x, hs$y, col=1)
legend("topleft", c("All tiles",sprintf("Center >%.2f, %d tiles (%.1f%%)",adpred_cut,length(icen),length(icen)/nrow(tf)*100)), lty=1, col=c(1,2))

quartz.save(sprintf("adtile_%s.png",score), type="png")

}


# Data frame of tiles containg a TAD according to ADpred - use ADpred score of central position only
cns = c("name","gene","first_resi","last_resi","aa","adp_cen","pap_ct","pap",
        "degron_score","degron_std","abundance_score","abundance_std","logabun_score","logabun_std")
tad = tf[icen,cns]

# numbers for paper
n = sum(! is.na(tad$degron_score))
n05 = sum(tad$degron_score > 0.5, na.rm=T)
n10 = sum(tad$degron_score > 1.0, na.rm=T)
print(sprintf("Of %d TAD tiles, %d (%.2f%%) have degron scores and %d (%.2f%%) score > 0.5 and %d (%.2f%%) score > 1.0",
              nrow(tad),n,n/nrow(tad)*100,n05,n05/n*100,n10,n10/n*100))


mut_tad = function(v, list_of_vec, start_from_max=T, mut_from=c("D","E"), mut_to="V") {
    gene = v[1]
    first_resi = as.numeric(v[2])
    last_resi = as.numeric(v[3])
    
    aa_vec = strsplit(v[4],"")[[1]]
    apdpred_vec = NA
    if (start_from_max) {
        # find position with max adpred score
        if (! gene %in% names(list_of_vec)) { return(NA) }
        if ( last_resi > length(list_of_vec[[gene]]) )  { return(NA) }
	apdpred_vec = list_of_vec[[gene]][first_resi:last_resi]
        i_start = which.max(apdpred_vec)
    } else {
        # always start from central position
        i_start = 16
    }
    # search for something to mutate start by going left (from central pos 16, pos 15 is more central than 17)
    i_nearest_asp = NA
    for (neighbor in seq(0,29)) {
	if ( i_start-neighbor > 0 ) {
	    if (aa_vec[i_start-neighbor] %in% mut_from) {
                i_nearest_asp = i_start-neighbor
	        break
	    }
        }
        if ( i_start+neighbor <= length(aa_vec) ) {
            if (aa_vec[i_start+neighbor] %in% mut_from) {
                i_nearest_asp = i_start+neighbor
	        break
	    }
	}
    }
    if (is.na(i_nearest_asp)) {
        print(sprintf("%s  nothing to mutate in %s",v[1],v[4]))
        return(NA)
    } else {
        # print(sprintf("%s  %s  start %d  %s%d%s  %s", v[1], paste0(apdpred_vec,collapse=";"), i_start, aa_vec[i_nearest_asp], i_nearest_asp, mut_to, v[4]))
        aa_vec[i_nearest_asp] = mut_to
        return( paste0(aa_vec, collapse="") )
    }
}
tad$aa_mut = apply(tad[,c("gene","first_resi","last_resi","aa")], MARGIN=1, mut_tad, list_of_vec=adp, start_from_max=F, mut_from=c("D","E"), mut_to="V")
print(sprintf("Failed to mutate %d of %d tiles (%.1f%%)", sum(is.na(tad$aa_mut)), nrow(tad), sum(is.na(tad$aa_mut))/nrow(tad)*100))

write.table(tad[which(! is.na(tad$aa_mut)),c("name","aa_mut")], file="tad_mut.seq", quote=F, col.names=F, row.names=F)
tad_mut_pap = read.table(gzfile("tad_mut_pap.txt.gz"))
colnames(tad_mut_pap) = c("name","aa","pap_ct","pap","first_resi","last_resi","cen_aa","cen_resi")

i = match(tad$name,tad_mut_pap$name)
tad$pap_mut = tad_mut_pap[i,"pap"]
tad$pap_ct_mut = tad_mut_pap[i,"pap_ct"]

# numbers for paper
stopifnot(all(! is.na(tad$pap)))
n02 = sum(tad$pap < 0.2)
print(sprintf("Of %d TAD tiles with PAP scores, %d (%.2f%%) has score < 0.2", nrow(tad), n02, n02/nrow(tad)*100))
stopifnot(all(! is.na(tad$pap_mut)))
n02 = sum(tad$pap_mut < 0.2)
print(sprintf("Of %d mutated TAD tiles with PAP scores, %d (%.2f%%) has score < 0.2", nrow(tad), n02, n02/nrow(tad)*100))

# Plot PAP score distributions of AD containing tiles with and without mutations
quartz(width=10, height=5)
par(mfrow=c(1,2), mar=c(5,4,2,2)+.1, bg="white")

n=3
breaks=60
hs  = hist(tf[,"pap"],     breaks=breaks, plot=F)
hp  = hist(tad[,"pap"],    breaks=breaks, plot=F)
hpm = hist(tad[,"pap_mut"],breaks=breaks, plot=F)

plot(0,0,col=0, xlim=c(-0.2,0.8), ylim=c(0,max(hpm$density,hp$densiy,hs$densiy)), xlab="PAP score", ylab="Density")
lines(hp$mids,  runavg(hp$density ,n=n), lwd=2, col=2)
lines(hpm$mids, runavg(hpm$density,n=n), lwd=2, col=3)
lines(hs$mids,  runavg(hs$density ,n=n), lwd=2, col=1)
legend("topright", c("All tiles","TAD","Mutated TAD"), lty=1, lwd=2, col=c(1,2,3))

hs  = hist(log10( tf[,"pap"]),     breaks=breaks, plot=F)
hp  = hist(log10( tad[,"pap"]),    breaks=breaks, plot=F)
hpm = hist(log10( tad[,"pap_mut"]),breaks=breaks, plot=F)

plot(0,0,col=0, xlim=c(-3,0), ylim=c(0,max(hpm$density,hp$density,hs$density)), xlab="PAP score, log10", ylab="Density")
lines(hp$mids,  runavg(hp$density ,n=n), lwd=2, col=2)
lines(hpm$mids, runavg(hpm$density,n=n), lwd=2, col=3)
lines(hs$mids,  runavg(hs$density ,n=n), lwd=2, col=1)
legend("topleft", c("All tiles","TAD","Mutated TAD"), lty=1, lwd=2, col=c(1,2,3))

quartz.save("adpred_tad_mutated.png", type="png")


# Paper plot
pdf("adpred_tad_mutated.pdf", width=3.5, height=3.5, pointsize=7)
n=5
breaks=60
hs  = hist(tf[,"pap"],     breaks=breaks, plot=F)
hp  = hist(tad[,"pap"],    breaks=breaks, plot=F)
hpm = hist(tad[,"pap_mut"],breaks=breaks, plot=F)
par(mar=c(5,4,2,2)+.1)
plot(0,0,col=0, xlim=c(-0.2,0.8), ylim=c(0,max(hpm$density,hp$densiy,hs$densiy)), xlab="PAP score", ylab="Density")
lines(hp$mids,  runavg(hp$density ,n=n), lwd=2, col=2)
lines(hpm$mids, runavg(hpm$density,n=n), lwd=2, col=3)
lines(hs$mids,  runavg(hs$density ,n=n), lwd=2, col=1)
legend("topright", c("All tiles","TAD","Mutated TAD"), lty=1, lwd=2, col=c(1,2,3))
dev.off()

