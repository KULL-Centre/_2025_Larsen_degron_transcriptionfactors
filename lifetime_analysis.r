options(width=180)

tf = read.csv("transfac_tiles_lifetime.csv")
tfp = read.csv("transfac_proteins_lifetime.csv")

# internal PAP score
tf$pap_int = tf$pap - tf$pap_ct

# Median pLDDT per protein
tfp$plddt_med = sapply(strsplit(tfp[,"plddt"],";"), function(v){ median(as.numeric(v)) })

# Capping of high half-lifes in correlations, e.g. half-life >70h is just a stable protein
cap = list(Li21_hl=24, Liu19_hl=33, Mathieson18_hl=102, Nardone23g_psi=Inf)
stopifnot(all( names(cap) %in% colnames(tfp) ))
for (cn in names(cap)) { tfp[,cn] = sapply(tfp[,cn], min, cap[[cn]]) }

hl = read.csv("turnover.csv")
cap_hl = list(li21_hl=cap[["Li21_hl"]], liu19_hl_mean=cap[["Liu19_hl"]], mat18_hl_mean=cap[["Mathieson18_hl"]], nar23_PSI.DMSO=cap[["Nardone23g_psi"]])
stopifnot(all( names(cap_hl) %in% colnames(hl) ))
for (cn in names(cap_hl)) { hl[,cn] = sapply(hl[,cn], min, cap_hl[[cn]]) }


##
## Analysis
##

# rasa_cut = 0.4
plddt_cut = 70

predict_abundance = function(df, tile_idr_cut=70) {
    stopifnot( nrow(df) > 1 )
    i = which(df$tile_plddt_median < tile_idr_cut & ! is.na(df[,"abundance_score"]))
    # i = which(df$tile_rasa_median > tile_idr_cut & ! is.na(df[,score_cn]))
    if (length(i) > 0) {
        # Find the longest streatch of IDR tiles
        runs = rle(diff(i))
        longest_seq = max(c(0,runs$length[which(runs$values==1)])) +1
    
        # last = cumsum(runs$lengths)
        # first = last - runs$lengths +1
        # i = which( ifelse(runs$values==1, runs$lengths, 0) == (tfp[ip_gene,"longest_idr"]-1) )
    
        return( c(mean(df[i,"abundance_score"]), min(df[i,"abundance_score"]),
	          mean(df[i,"abundance_nocdeg"]), min(df[i,"abundance_nocdeg"]),
		  mean(df[i,"pap_int"]), min(df[i,"pap_int"]), length(i), longest_seq) )
    } else {
        return( c(NA,NA,NA,NA,NA,NA,0,NA) )
    }
}

tfp$mean_idr_abun = NA
tfp$min_idr_abun = NA
tfp$mean_idr_abun_nc = NA
tfp$min_idr_abun_nc = NA
tfp$mean_idr_pap = NA
tfp$min_idr_pap = NA
tfp$n_tiles = NA
tfp$n_disordered_tiles = NA
tfp$longest_idr = NA
nodata_genes = c()
for (gene in tfp[which(! grepl("nan",tfp$plddt)),"gene"]) {
    i_tiles = which(tf$gene == gene)
    ip_gene = which(tfp$gene == gene)
    # tfp[ip_gene,"mean_abun"] = mean(tf[i_tiles,"abundance_score"], na.rm=T)
    # tfp[ip_gene,"min_abun"] = min(tf[i_tiles,"abundance_score"], na.rm=T)
    mds_n = predict_abundance(tf[i_tiles,])
    stopifnot(length(ip_gene) == 1)
    tfp[ip_gene,"mean_idr_score"] = mds_n[1]
    tfp[ip_gene,"min_idr_score"] = mds_n[2]
    tfp[ip_gene,"mean_idr_score_nc"] = mds_n[3]
    tfp[ip_gene,"min_idr_score_nc"] = mds_n[4]
    tfp[ip_gene,"mean_idr_pap"] = mds_n[5]
    tfp[ip_gene,"min_idr_pap"] = mds_n[6]
    if (mds_n[7] < 1) {
        print(sprintf("Gene %s has no disorded tiles with scores (of %d)", gene, length(i_tiles)))
	# print(tf[i_tiles,cns])
	nodata_genes = c(nodata_genes,gene)
    }
    # print(sprintf("== %s tiles %d ==", gene, length(i_tiles)))
    # print(mds_n)
    tfp[ip_gene,"n_disordered_tiles"] = mds_n[7]
    tfp[ip_gene,"longest_idr"] = mds_n[8]
    tfp[ip_gene,"n_tiles"] = length(i_tiles) 
}


# Single plot pearsons
calc_rp = function(pct, exp_cn, pred_cn, cor_genes_idx) {
    i = cor_genes_idx[which(tfp[cor_genes_idx,"n_disordered_tiles"]*100.0/tfp[cor_genes_idx,"n_tiles"] >= pct)]
    rp = cor(tfp[i,exp_cn], tfp[i,pred_cn], method="pearson")
    return( c(rp,length(i)) )
}

# quartz(width=12, height=4)
pdf("lifetime_tf_correlations.pdf", width=12, height=4, pointsize=12)
par(mfrow=c(1,3), mar=c(5,4,4,4)+.1)

x = seq(0,100)
pred_cns = c("mean_idr_score","min_idr_score",
             "mean_idr_score_nc","min_idr_score_nc",
	     "mean_idr_pap","min_idr_pap")
pred_col = c(1,1, 3,3, 4,4)
pred_lty = c(1,3, 1,3, 1,3)
pred_lab = c("Mean Abun.", "Min Abun.","Mean ACC","Min ACC","Mean PAP","Min PAP")
names(pred_col) = pred_cns
names(pred_lty) = pred_cns
names(pred_lab) = pred_cns
exp_cns = c("Li21_hl","Nardone23g_psi","Mathieson18_hl")
exp_lab = c("Li21 Half-life","Nardone23 PSI","Mathieson18 Half-life")
names(exp_lab) = exp_cns
has_plddt = ! grepl("nan",tfp$plddt)
for (exp_cn in exp_cns) {
    has_exp_data = ! is.na(tfp[,exp_cn])
    ip_genes = which(has_plddt & has_exp_data)
    print(sprintf("Of %d %s, %d has pLDDT assigned", sum(has_exp_data), exp_cn, length(ip_genes)))
    # plot(0,0,col=0, xlim=c(1,100), ylim=c(0,1.2), ylab="Pearson", xlab="Protein disorder threshold [% tiles]", main=exp_lab[exp_cn])
    plot(0,0,col=0, xlim=c(1,100), ylim=c(0,1.2), ylab="", xlab="", main=exp_lab[exp_cn])
    for (pred_cn in pred_cns) {
        list_rp_l = lapply(x, calc_rp, exp_cn=exp_cn, pred_cn=pred_cn, cor_genes_idx=ip_genes)
        rp = sapply(list_rp_l, "[", 1)
        len = sapply(list_rp_l, "[", 2)
        lines(x, rp, col=pred_col[pred_cn], lty=pred_lty[pred_cn])
        lines(x, len/len[1], col=2)
    }
    axis(4, at=c(0,1), labels=c(0,len[1]), col=2)
    ax_cex = 0.8
    mtext("Number of proteins", side=4, line=2, cex=ax_cex, col=2)
    mtext("Pearson", side=2, line=2.5, cex=ax_cex, col=1)
    mtext("Protein disorder threshold [%]", side=1, line=2.5, cex=ax_cex, col=1)
    legend("top", pred_lab, lty=pred_lty, col=pred_col, ncol=3)
}
dev.off()


# Scatter plots
idp_min_frac_disorder = 0.8
pred_cn = "mean_idr_score_nc"
pdf("lifetime_tf_scatter.pdf", width=12, height=4, pointsize=12)
par(mfrow=c(1,3), mar=c(5,4,4,4)+.1)
has_plddt = ! grepl("nan",tfp$plddt)
for (exp_cn in exp_cns) {
    has_exp_data = ! is.na(tfp[,exp_cn])
    ip_genes = which(has_plddt & has_exp_data)
    i_all = ip_genes[which((! is.na(tfp[ip_genes,exp_cn])) & (! is.na(tfp[ip_genes,pred_cn])))]
    rp = cor(tfp[i_all,exp_cn], tfp[i_all,pred_cn], method="pearson")
    i_idp = ip_genes[which(tfp[ip_genes,"n_disordered_tiles"]/tfp[ip_genes,"n_tiles"] >= idp_min_frac_disorder)]
    rp_idp = cor(tfp[i_idp,exp_cn], tfp[i_idp,pred_cn], method="pearson")
    main = sprintf("Pearson %.3f (%d), IDP (%.0f%%) %.3f (%d)", rp, length(i_all), idp_min_frac_disorder*100, rp_idp, length(i_idp))
    plot(tfp[i_all,exp_cn], tfp[i_all,pred_cn], # ylim=c(0.08,0.4),
         xlab=sprintf("Experimental turnover (%s)",exp_cn), ylab=sprintf("Score for IDR tiles (%s)",pred_cn), main=main)
    points(tfp[i_idp,exp_cn], tfp[i_idp,pred_cn], pch=16, col=2)
    title(main=main, sub=sprintf("IDR tile median pLDDT < %.1f",plddt_cut))
}
dev.off()

idp_min_frac_disorder = 0.8
quartz(width=5, height=5)
exp_cn = "Li21_hl"
pred_cn = "mean_idr_score_nc"
ip_genes = which((! grepl("nan",tfp$plddt)) & (! is.na(tfp[,exp_cn])))
i_all = ip_genes[which((! is.na(tfp[ip_genes,exp_cn])) & (! is.na(tfp[ip_genes,pred_cn])))]
rp = cor(tfp[i_all,exp_cn], tfp[i_all,pred_cn], method="pearson")
i_idp = ip_genes[which(tfp[ip_genes,"n_disordered_tiles"]/tfp[ip_genes,"n_tiles"] >= idp_min_frac_disorder)]
rp_idp = cor(tfp[i_idp,exp_cn], tfp[i_idp,pred_cn], method="pearson")
main = sprintf("Pearson %.3f (%d), IDP (%.0f%%) %.3f (%d)", rp, length(i_all), idp_min_frac_disorder*100, rp_idp, length(i_idp))
plot(tfp[i_all,exp_cn], tfp[i_all,pred_cn], # ylim=c(0.08,0.4),
     xlab=sprintf("Experimental turnover (%s)",exp_cn), ylab=sprintf("Score for IDR tiles (%s)",pred_cn), main=main)
points(tfp[i_idp,exp_cn], tfp[i_idp,pred_cn], pch=16, col=2)
title(main=main, sub=sprintf("IDR tile median pLDDT < %.1f",plddt_cut))


###################################################################
### Correlations involving all turnover data
###################################################################
# Could try the all-vs-all plot using Nardone PSI averaged per gene and filtered for measurements that vary a lot between cell types

# Correlations between half-lifes - much easier when data is merged
plot_cor = function(df, pch=16, cex=.2, repnames=NA, n_tot_var=NA, ...) {
    text_rp = function(x, y, cex=1, ...) {
        i = which((! is.na(x)) & (! is.na(y)))
        rp = cor(x[i],y[i],method="pearson")    #,use="complete.obs")
        rs = cor(x[i],y[i],method="spearman")    #,use="complete.obs")
        x_mid = (max(x, na.rm=T) + min(x, na.rm=T)) /2.0
        y_mid = (max(y, na.rm=T) + min(y, na.rm=T)) /2.0
        text(x_mid, y_mid, sprintf("Pearson %.2f\nSpearman %.2f\nN=%d",rp,rs,length(i)), cex=1.6, ...)
    }
    scatter = function(x,y,...) {
        points(x,y, ...)
	abline(c(0,1), col=8)
    }
    if (any(is.na(repnames))) { repnames = colnames(df) } else { stopifnot(length(repnames)==ncol(df)) }
    if (ncol(df) == 2) {
        plot(df, xlab=repnames[1], ylab=repnames[2], pch=pch, cex=cex, ...)
	abline(c(0,1), col=8)
    } else if (ncol(df) > 2) {
        sums = apply(df, MARGIN=2, sum)
	if (is.na(n_tot_var)) { n_tot_var = nrow(df) }
        # covers = apply(df, MARGIN=2, function(v) { sum(! is.na(v))*100.0/n_tot_var })
        # labels = sprintf("%s\n%.2f%% cover",repnames,covers)
        covers = apply(df, MARGIN=2, function(v) { sum(! is.na(v)) })
        labels = sprintf("%s\nN=%d",repnames,covers)
        plot(df, labels=labels, upper.panel=text_rp, lower.panel=scatter, pch=pch, cex=cex, ...)
        # plot(df, labels=labels, upper.panel=text_rp, lower.panel=scatter, pch=pch, cex=cex, ...)
	# title(xlab="Average gate", ylab="Average gate")
    }
}


pdf("lifetime_correlations.pdf", width=6, height=6, pointsize=12)
# quartz(width=6, height=6)
plot_cor(hl[,c("li21_hl","mat18_hl_mean","nar23_PSI.DMSO")], repnames=c("Li21","Mathieson18","Nardone23"), cex=.4)
dev.off()


###################################################################
### C-degron corrected scores
###################################################################

deming_coef = c(0.05738429,1.62342637)

breaks = seq(-0.09, 1.11, 0.02)
i_use = which( (! is.na(tf$pap)) & (! is.na(tf$abundance_score)) )
ha = hist(tf[i_use,"abundance_score"], breaks=breaks, plot=F)
hc = hist(tf[i_use,"abundance_nocdeg"], breaks=breaks, plot=F)
hp = hist(tf[i_use,"pap_int"], breaks=breaks, plot=F)

pdf("cdeg_corrected_distributions.pdf", width=5, height=5, pointsize=12)
par(mar=c(5,4,2,2)+.1)
plot(ha$mids, ha$density, lwd=2, type="l", xlab="Abundance score", ylab="Density", xlim=c(-0.05,0.8))
lines(hc$mids,  hc$density,  lwd=2, col=3)
lines(hp$mids,  hp$density,  lwd=2, col=4)
legend("topright", c("Original","C-degron corrected (ACC)","PAP internal"), lty=1, lwd=2, col=c(1,3,4))
dev.off()

