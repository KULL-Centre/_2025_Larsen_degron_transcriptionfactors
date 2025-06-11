options(width=160)

args = commandArgs(trailingOnly=TRUE)
if (interactive()) {
    counts_file = "counts.rda"
} else if (length(args) != 1) {
    print("")
    print("usage: Rscript scores.r  <counts.rda>")
    quit(save="no")
} else {
    counts_file = args[1]
}
print(sprintf("Args: %s",paste0(args, collapse=" ")))


print(sprintf("Loading file with raw counts %s",counts_file))
load(counts_file)

settings = list()
settings[["version"]] = 1.0
settings[["input"]] = counts_file
settings[["norms"]] = c("psi")


#####################################################################
##     Merge replicates
#####################################################################
# Technical replica merged (tm)
trep_nreads = c()
trep_rp = c()
counts_tm = list(); for (lib in names(counts)) { counts_tm[[lib]] = counts[[lib]][,name_cols] }
# par(ask=T)
for (si in which(samples$tech_rep==1 & samples$obs)) {
    # sample row of other tech_rep 
    si2 = which(samples$lib == samples[si,"lib"]  &  samples$gate == samples[si,"gate"]  &  samples$bio_rep == samples[si,"bio_rep"]  &  
               samples$facs_rep == samples[si,"facs_rep"]  &  samples$tech_rep == 2)
    stopifnot(length(si2) == 1)
    stopifnot(samples[si2,"obs"])

    lib = samples[si,"lib"]
    cn_rep = c(samples[si,"file"], samples[si2,"file"])
    cn = sprintf("%s%d_%s_fasc%d", samples[si,"lib"], samples[si,"gate"], samples[si,"bio_rep"], samples[si,"facs_rep"])
    cn = sprintf("bio%d_facs%d_g%d", samples[si,"bio_rep"], samples[si,"facs_rep"], samples[si,"gate"])
    
    counts_tm[[lib]][,cn] = apply(counts[[lib]][,cn_rep], MARGIN=1, sum)

    # Correlations of relevant library members
    # ri = which(raw_clean$name %in% library[which(library$pool==samples[si,"lib"]),"name"])
    rp = cor(counts[[lib]][,cn_rep[1]], counts[[lib]][,cn_rep[2]], method="pearson")
    rs = cor(counts[[lib]][,cn_rep[1]], counts[[lib]][,cn_rep[2]], method="spearman")
    rp_log = cor(log(counts[[lib]][,cn_rep[1]]+1), log(counts[[lib]][,cn_rep[2]]+1), method="pearson")
    trep_rp = c(trep_rp, rp)
    print(sprintf("Pearson between technical replica %s and %s of %d lib %s members %.2f (log %.2f, spearman %.2f). Counts summaries:",
                  cn_rep[1],cn_rep[2],nrow(counts[[lib]]),samples[si,"lib"],rp,rp_log,rs))
    print(summary(counts[[lib]][,cn_rep[1]]))
    print(summary(counts[[lib]][,cn_rep[2]]))

    # plot(counts[[lib]][,cn_rep[1]], counts[[lib]][,cn_rep[2]], xlab=cn_rep[1], ylab=cn_rep[2])

    # Read counts per replica
    nreads = apply(counts[[lib]][,cn_rep], MARGIN=2, sum)
    trep_nreads = c(trep_nreads, nreads)
    print(sprintf("Number of reads for technical replica: %s", paste(nreads, collapse=", ")))
}
print("Tech replica correlations")
print(summary(trep_rp))
print("Reads per tech replica")
print(summary(trep_nreads))

print("Distribution of raw read counts per tile across all libraries and replicates (technical replicates merged)")
print(summary(unlist(sapply(counts_tm,function(df){ df[,4:ncol(df)]}))))


#####################################################################
##     Calculate PSI
#####################################################################

# Function that normalize all columns to RPM
reads_per_million = function(df) {
    col_sums = apply(df, MARGIN=2, sum)
    df * matrix(10^6/col_sums, nrow=dim(df)[1], ncol=dim(df)[2], byrow=TRUE)
}

# Calculate count frequencies with psudo_counts
settings$psudo_counts = 0
rpm = list()
for (lib in names(counts_tm)) {
    rpm[[lib]] = counts_tm[[lib]][,name_cols]
    rpm[[lib]] = cbind(rpm[[lib]], reads_per_million(counts_tm[[lib]][,(length(name_cols)+1):ncol(counts_tm[[lib]])] + settings$psudo_counts))
}

# Calc PSI values per library
psi = list()
all_rps = c()
for (lib in names(rpm)) {
    psi[[lib]] = rpm[[lib]][,name_cols]
    
    # Calc PSI for all replica in library
    cns_rep = c()
    for (si in which(samples$lib==lib & samples$gate==1 & samples$tech_rep==1 & samples$obs)) {
        cn_gates = sprintf("bio%d_facs%d_g%d", samples[si,"bio_rep"], samples[si,"facs_rep"], seq(4))
	stopifnot(all(cn_gates %in% colnames(rpm[[lib]])))
	
	cn = sprintf("bio%d_facs%d_psi", samples[si,"bio_rep"], samples[si,"facs_rep"])
	psi[[lib]][,cn] = apply(rpm[[lib]][,cn_gates], MARGIN=1, function(v){vsum=sum(v); sum( v/vsum * c(1,2,3,4) )})
	
	n_ok = sum(! is.na(psi[[lib]][,cn]))
	print(sprintf("PSI calculated for %d of %d (%.2f%%) tiles of lib %s replica %s", n_ok, nrow(psi[[lib]]), n_ok/nrow(psi[[lib]])*100, lib, cn))
	cns_rep = c(cns_rep, cn)

        # RPM sum over gates
        cn = sprintf("bio%d_facs%d_rpm_sum", samples[si,"bio_rep"], samples[si,"facs_rep"])
	psi[[lib]][,cn] = apply(rpm[[lib]][,cn_gates], MARGIN=1, sum)
	
        # Absolute number of counts summed over gates
        cn = sprintf("bio%d_facs%d_rd_sum", samples[si,"bio_rep"], samples[si,"facs_rep"])
	stopifnot(all(psi[[lib]][,"dna"] == counts_tm[[lib]][,"dna"]))
        psi[[lib]][,cn] = apply(counts_tm[[lib]][,cn_gates], MARGIN=1, sum)
    }
    
    # Pearson correlations of replica PSI - these are independt of normalization
    rp_mat = cor(psi[[lib]][,cns_rep], method="pearson", use="pairwise.complete.obs")
    rps = rp_mat[upper.tri(rp_mat)]
    print(sprintf("Lib %s PSI Pearson of %d replica pairs: mean %.4f range %.4f - %.4f", lib, length(rps), mean(rps), min(rps), max(rps)))
    all_rps = c(all_rps, rps)
}
print(sprintf("All lib PSI Pearson of %d replica mean %.4f range %.4f - %.4f", length(all_rps), mean(all_rps), min(all_rps), max(all_rps)))


#####################################################################
##     Filter unreliable PSI values
#####################################################################
# Each tile measurement should have read counts (rd_sum) or read frequency (rpm_sum) above threshold to be considered reliable
settings$rd_threshold = 100

# store one vector with all measurements (tile,replica)
tr_dev = c()
tr_rd = c()

for (lib in names(psi)) {
    cns_psi = colnames(psi[[lib]])[ which(grepl("facs[0-9]_psi", colnames(psi[[lib]]))) ]
    cns_psi = cns_psi[order(cns_psi)]

    cns_rd = colnames(psi[[lib]])[ which(grepl("facs[0-9]_rd_sum", colnames(psi[[lib]]))) ]
    cns_rd = cns_rd[order(cns_rd)]

    # scale works column-wise and apply returns a column per sweep
    zscore = data.frame(t(scale(t(psi[[lib]][,cns_psi]), center=T, scale=F)))
    
    # column-wise flattening
    x = unlist(psi[[lib]][,cns_rd])
    y = unlist(zscore)
    tr_rd = c(tr_rd, x)
    tr_dev = c(tr_dev, y)
}

# # Plot PSI deviation per measurement as a function of reads
# # Firstly, calculate standard deviation of deviance in bins of read counts
# breaks = 10^seq(-1,4,0.1)
# mids = (breaks[2:length(breaks)] + breaks[1:(length(breaks)-1)])/2
# tr_bin = cut(tr_rd, breaks=breaks)
# agg = aggregate(tr_dev, by=list(tr_bin), sd, na.rm=T)

# # Assume tr_bin is ordered and that this order is preserved by aggregate
# lines(mids[match(as.character(agg[,1]),levels(tr_bin))], agg[,2]*3, col=2, lwd=2)
# abline(h=0, col="grey")
# abline(v=settings$rd_threshold, lty=2)
# legend("topright", c("Measurements", "Bin std. dev. x3", "Threshold"), pch=c(".",NA,NA), lty=c(NA,1,2), col=c(1,2,1), lwd=2)

# Check that tr_dev has an element for all possible measurements assuming that each library has 6 experimets
n_meas = sum( sapply(psi, nrow) * sapply(psi, function(df){ sum(grepl("psi", colnames(df))) }) )
stopifnot(n_meas == length(tr_dev))

# Report coverage
i_notna = which((! is.na(tr_rd)) & (! is.na(tr_dev)))
n_reliable = sum(tr_rd[i_notna] > settings$rd_threshold)
print(sprintf("Of %d measurements, %d (%.2f%%) exists and %d (%.2f%%) are above read count threshold of %.1f",
              n_meas, length(i_notna), length(i_notna)/n_meas*100, n_reliable, n_reliable/n_meas*100, settings$rd_threshold))


# Make a new psi data frame and set unreliable measurements to NA
psi_clean = psi
n_rm = 0
for (lib in names(psi_clean)) {
    cns_psi = colnames(psi_clean[[lib]])[ which(grepl("facs[0-9]_psi", colnames(psi_clean[[lib]]))) ]
    cns_psi = cns_psi[order(cns_psi)]

    # cns_rpm = colnames(psi_clean[[lib]])[ which(grepl("facs[0-9]_rpm_sum", colnames(psi_clean[[lib]]))) ]
    # cns_rpm = cns_rpm[order(cns_rpm)]
    cns_rd = colnames(psi_clean[[lib]])[ which(grepl("facs[0-9]_rd_sum", colnames(psi_clean[[lib]]))) ]
    cns_rd = cns_rd[order(cns_rd)]

    for (ic in seq_along(cns_psi)) {
        cn_rd = cns_rd[ic]
	cn_psi = cns_psi[ic]
        ip = which(psi_clean[[lib]][,cn_rd] < settings$rd_threshold)
	psi_clean[[lib]][ip,cn_psi] = NA
	n_rm = n_rm + length(ip)
	
	# report
	n_na = sum(is.na(psi_clean[[lib]][,cn_psi]))
	print(sprintf("Removed %d (%.2f%%) PSI's resulting in a total %d (%.2f%%) NA values for lib %s replica %s. Read distribution of non-NA psi values:",
	              length(ip), length(ip)/nrow(psi_clean[[lib]])*100, n_na, n_na/nrow(psi_clean[[lib]])*100, lib, cn_psi))
    }
}
print(sprintf("Removed %d measurements with read count sum < %.1f", n_rm, settings$rd_threshold))


# Look at read frequencies of PSI values that pass the RPM threshold
for (lib in names(psi)) {
    cns_psi = colnames(psi_clean[[lib]])[ which(grepl("facs[0-9]_psi", colnames(psi_clean[[lib]]))) ]
    cns_psi = cns_psi[order(cns_psi)]

    cns_rpm = colnames(psi_clean[[lib]])[ which(grepl("facs[0-9]_rpm_sum", colnames(psi_clean[[lib]]))) ]
    cns_rpm = cns_rpm[order(cns_rpm)]

    # cns_rd = colnames(psi_clean[[lib]])[ which(grepl("facs[0-9]_rd_sum", colnames(psi_clean[[lib]]))) ]
    # cns_rd = cns_rd[order(cns_rd)]

    for (ic in seq_along(cns_psi)) {
        cn_rpm = cns_rpm[ic]
	cn_psi = cns_psi[ic]
        # cn_rd = cns_rd[ic]
	i_ok = which(! is.na(psi_clean[[lib]][,cn_psi]))
	
	print(sprintf("Lib %s replica %s read distribution of %d (%.2f%%) non-NA psi values:",
	              lib, cn_psi, length(i_ok), length(i_ok)/nrow(psi_clean[[lib]])*100))
	print(summary(psi_clean[[lib]][i_ok,cn_rpm]))
	# print(summary(psi_clean[[lib]][i_ok,cn_rd]))
    }
}


#####################################################################
##     Function to plot normalizations
#####################################################################

plot_norm = function(norm_name="psi", xlim=NA, ylim=NA, breaks=NA) {
    # Plot PSI distribution per FACS replica
    quartz(width=12, height=9)
    par(bg="white")
    
    breaks_min = NA; breaks_max = NA
    for (lib in names(psi_clean)) {
        cns = colnames(psi_clean[[lib]])[ which(grepl(sprintf("facs[0-9]_%s",norm_name), colnames(psi_clean[[lib]]))) ]
        breaks_min = min(c(breaks_min,unlist(psi_clean[[lib]][,cns])), na.rm=T)
        breaks_max = max(c(breaks_max,unlist(psi_clean[[lib]][,cns])), na.rm=T)
    }
    print(sprintf("Score range %.3f to %.3f",breaks_min,breaks_max))
    if (any(is.na(breaks))) { breaks = seq(breaks_min,breaks_max,length.out=100) }
    else if (length(breaks)==1) { breaks = seq(breaks_min,breaks_max,length.out=breaks) }
    if (breaks[3]-breaks[2] > 1) { print(sprintf("WARNING: Histogram bin size is %.1f",breaks[3]-breaks[2])) }
    if (any(is.na(xlim))) { xlim=c(breaks_min,breaks_max) }
    if (any(is.na(ylim))) { ylim=c(0,2) }
    plot(0,0,col=0, xlim=xlim, ylim=ylim, xlab="Score", ylab="Density", main=norm_name)
    for (il in seq_along(psi_clean)) {
        lib = names(psi_clean)[il]
        cns = colnames(psi_clean[[lib]])[ which(grepl(sprintf("facs[0-9]_%s",norm_name), colnames(psi_clean[[lib]]))) ]
        cns = cns[order(cns)]
        print(sprintf("Library %s has %d replicates: %s",lib,length(cns),paste(cns,collapse=", ")))
        for (ic in seq_along(cns)) {
            cn = cns[ic]
            h = hist(psi_clean[[lib]][,cn], breaks=breaks, plot=F)
	    lines(h$mids, h$density, col=il, lty=((ic-1)%/%2)+1, lwd=2)
        }
    }
    legend("topleft", paste(names(psi_clean),sapply(psi_clean,nrow),"tiles"), lty=1, lwd=2, col=seq(length(psi_clean)), ncol=2)
    legend("topright", sprintf("bio%d", seq(3)), lty=seq(3), lwd=2, col=1, ncol=2)
    # quartz.save("psi_distributions.png", type="png")
}


#####################################################################
##     Normalize each sorting replica using conntrol peptides
#####################################################################
# 
settings$norms = c(settings$norms, "cnorm")

# calculate mean PSI for each control in each lib
for (lib in names(psi_clean)) {
    cns_all = colnames(psi_clean[[lib]])
    cns_psi = cns_all[which(grepl("psi",cns_all))]
    psi_mat = psi_clean[[lib]][match(ctrl[,"name"],psi_clean[[lib]][,"name"]),cns_psi]
    ctrl[,lib] = apply(psi_mat, MARGIN=1, mean)
    # sd not informative since replicates are not normalized yet
    # ctrl[,sprintf("%s_sd",lib)] = apply(psi_mat, MARGIN=1, sd)
}
has_all_measured = apply(ctrl[,names(psi_clean)], MARGIN=1, function(v){ all(! is.na(v)) })
ctrl_ids = ctrl[which(has_all_measured),"name"]
print(sprintf("cnorm control peptides measure in all replicates: %s", paste(ctrl_ids,collapse=",")))

settings$cnorm_stab = c("pparg_1","pparg_8")
settings$cnorm_deg  = c("pparg_27","pparg_30","pparg_32")
print(sprintf("cnorm stable control peptides: %s", paste(settings$cnorm_stab,collapse=",")))
print(sprintf("cnorm degron control peptides: %s", paste(settings$cnorm_deg,collapse=",")))

norm_list = list()
for (lib in names(psi_clean)) {
    # Find indices of control peptides
    is = c()	
    for (pep_name in settings$cnorm_stab) {
        is = c(is, which(grepl(pep_name, psi_clean[[lib]][,"name"])))
    }
    stopifnot(length(is) == length(unique(is)))
    id = c()
    for (pep_name in settings$cnorm_deg) {
        id = c(id, which(grepl(pep_name, psi_clean[[lib]][,"name"])))
    }
    stopifnot(length(id) == length(unique(id)))

    cns = colnames(psi_clean[[lib]])
    cns = cns[which(grepl("facs[0-9]_psi",cns))]
    cns = cns[order(cns)]
    for (cn in cns) {
	degron_psi = mean(psi_clean[[lib]][id,cn], na.rm=T)
	stabil_psi = mean(psi_clean[[lib]][is,cn], na.rm=T)

	# Min-max normalize
        cn_norm = sub("psi","cnorm",cn)
        psi_clean[[lib]][,cn_norm] = (psi_clean[[lib]][,cn] - degron_psi)/(stabil_psi - degron_psi)

	# Store normalization values
	norm_list[[paste0(lib,"_",cn)]] = c(lib,cn,degron_psi,stabil_psi)
    }
}

# Data frame with normalization values
cnorm = data.frame(name = names(norm_list),
                   lib = sapply(norm_list, "[", 1),
                   colname = sapply(norm_list, "[", 2),
                   degron_psi = sapply(norm_list, "[", 3),
                   stabil_psi = sapply(norm_list, "[", 4))


#####################################################################
##     Normalize each sorting replica by peak position
#####################################################################
settings$norms = c(settings$norms, "pnorm")

norm_list = list()
for (lib in names(psi_clean)) {
    cns = colnames(psi_clean[[lib]])
    cns = cns[which(grepl("facs[0-9]_psi",cns))]
    cns = cns[order(cns)]
    for (cn in cns) {
	# find peak positions as median of upper and lower 25% of scores
	q1_q3 = quantile(psi_clean[[lib]][,cn], na.rm=T)[c(2,4)]
	degron_psi = median(psi_clean[[lib]][which(psi_clean[[lib]][,cn] < q1_q3[1]), cn])
	stabil_psi = median(psi_clean[[lib]][which(psi_clean[[lib]][,cn] > q1_q3[2]), cn])

	# Min-max normalize
        cn_score = sub("psi","pnorm",cn)
        psi_clean[[lib]][,cn_score] = (psi_clean[[lib]][,cn] - degron_psi)/(stabil_psi - degron_psi)

	# Store normalization values
	norm_list[[paste0(lib,"_",cn)]] = c(lib,cn,degron_psi,stabil_psi)
    }
}

# Data frame with normalization values
pnorm = data.frame(name = names(norm_list),
                   lib = sapply(norm_list, "[", 1),
                   colname = sapply(norm_list, "[", 2),
                   degron_psi = sapply(norm_list, "[", 3),
                   stabil_psi = sapply(norm_list, "[", 4))


#####################################################################
##     Normalize each sorting replica by FACS distribution
#####################################################################
# per library, plot replica PSI distributions, weighted distributions, per-bin average SD and the FACS transformed distributions
whist = function(x, w, breaks=20, ...) {
    stopifnot(length(x)==length(w))
    ret = list()

    nb = length(breaks)
    if (nb == 1) {
        breaks = seq(min(x,rm.na=T), max(x,rm.na=T), length.out=breaks)
        nb = length(breaks)
    } else {
        bin_width = breaks[2]-breaks[1]
        ret$equidist = all(abs( breaks[1:(nb-1)]+bin_width-breaks[2:nb]) < 1e-6 )
        stopifnot(ret$equidist)
    }

    bins = cut(x, dig.lab=6, breaks=breaks, ...)
    agg = aggregate(w, by=list(bins), sum)
    df = data.frame(group = levels(bins), mids = (breaks[1:(nb-1)] + breaks[2:nb])/2.0)
    df$counts = agg[match(df$group,agg[,1]),2]
    df[which(is.na(df$counts)),"counts"] = 0

    ret[["breaks"]] = breaks
    ret[["mids"]] = df$mids
    ret[["counts"]] = df$counts
    ret[["density"]] = df$counts / sum(df$counts) / (breaks[2]-breaks[1])
    return(ret)
}

settings$norms = c(settings$norms, "fnorm")
settings$norms = c(settings$norms, "fwnorm")

# Map each sorting replica to a FACS distribution that represents the library
load("facs/facs.rda")      # get facs_events, facs_table, and facs_set

# select a FACS data set per replica
facs_table$rep = NA
facs_table[which(facs_table$name == "CT Sort 1 _T1"),  "rep"] = "CT_bio1_facs1"
facs_table[which(facs_table$name == "CT Sort 1_T2"),   "rep"] = "CT_bio1_facs2"
facs_table[which(facs_table$name == "CT Sort 2_T1"),   "rep"] = "CT_bio2_facs1"
facs_table[which(facs_table$name == "CT Sort 2_T2"),   "rep"] = "CT_bio2_facs2"
facs_table[which(facs_table$name == "CT Sort 3_T1"),   "rep"] = "CT_bio3_facs1"
facs_table[which(facs_table$name == "CT Sort 3_T2"),   "rep"] = "CT_bio3_facs2"
facs_table[which(facs_table$name == "Even Sort 1_T2"), "rep"] = "Even_bio1_facs1"
facs_table[which(facs_table$name == "Even Sort 2_T1"), "rep"] = "Even_bio1_facs2"
facs_table[which(facs_table$name == "Even Sort 2_T2"), "rep"] = "Even_bio2_facs1"
facs_table[which(facs_table$name == "Even Sort 3_T1"), "rep"] = "Even_bio2_facs2"
facs_table[which(facs_table$name == "Even Sort 3_T2"), "rep"] = "Even_bio3_facs1"
facs_table[which(facs_table$name == "Even sort 1 _T1"),"rep"] = "Even_bio3_facs2"
facs_table[which(facs_table$name == "Odd 1 Sort 1_T2"),"rep"] = "Odd_bio1_facs1"
facs_table[which(facs_table$name == "Odd Sort 2_T2"),  "rep"] = "Odd_bio1_facs2"
facs_table[which(facs_table$name == "Odd Sort 3_T2"),  "rep"] = "Odd_bio2_facs1"
facs_table[which(facs_table$name == "odd sort 2 _T1"), "rep"] = "Odd_bio2_facs2"
facs_table[which(facs_table$name == "odd sort 3 _T1"), "rep"] = "Odd_bio3_facs1"
facs_table[which(facs_table$name == "odd1 sort 1 _T1"),"rep"] = "Odd_bio3_facs2"

# Make histogram bins smaller at small values to capture the narrow peak
fbreaks = 2^(seq(-13, 0, length.out=500))

# 0.01 is not enough for odd3 bio2 facs1 - some of the high random numbers transform to very negative values
settings$fnorm_prior_fac = 0.02

# # Normalize all replica to the same distribution
# xf = facs_events[["Odd 1 Sort 1_T2"]][,facs_set$cn_ratio]
# xf = xf[which(xf > 2^(-13) & xf < 1.0)]
# hf = hist(xf, breaks=fbreaks, plot=F)
# hfd = (1.0-settings$fnorm_prior_fac)*hf$density + settings$fnorm_prior_fac*dnorm(hf$mids, mean=mean(xf), sd=sd(xf))
# bin_widths = hf$breaks[2:length(hf$breaks)] - hf$breaks[1:(length(hf$breaks)-1)]
# hfdd = hfd*bin_widths
# # Function that transform uniformly distributed numbers in [0,1] to FACS distribution
# facs_quant = splinefun(cumsum(hfdd)/sum(hfdd), hf$mids)

# claculate histogram breaks with extended range
calc_breaks = function(x, n, n_extra=2) {
    x_min = min(x, na.rm=T)
    x_max = max(xs,na.rm=T)
    x_db = (x_max-x_min)/n
    breaks = seq(x_min-n_extra*x_db, x_max+n_extra*x_db, length.out=n)
    return(breaks)
}

# quartz(width=8, height=6)
# rand_unif = runif(10^6, min=0.0, max=1.0)
# par(ask=T)
for (lib in names(psi_clean)) {
    # column names of replica PSI's
    cns = colnames(psi_clean[[lib]])
    cns = cns[which(grepl("facs[0-9]_psi",cns))]
    cns = cns[order(cns)]
    for (cn in cns) {
        replica = sprintf("%s_%s", lib, sub("_psi","",cn))
	fi = which(facs_table$rep==replica)
	stopifnot(length(fi)==1)
	
        xf = facs_events[[facs_table[fi,"name"]]][,facs_set$cn_ratio]
        xf = xf[which(xf > 2^(-13) & xf < 1.0)]
        hf = hist(xf, breaks=fbreaks, plot=F)

        # the gaussian prior is added to avoid regions without density in the empirical PDF that leads to constant regions in the CDF that cannot be inverted
        hfd = (1.0-settings$fnorm_prior_fac)*hf$density + settings$fnorm_prior_fac*dnorm(hf$mids, mean=mean(xf), sd=sd(xf))
        # hist returns points on the PDF but these represent different densities when the bin width varies so compensate with bin width
        bin_widths = hf$breaks[2:length(hf$breaks)] - hf$breaks[1:(length(hf$breaks)-1)]
        hfdd = hfd*bin_widths

        # Function that transform uniformly distributed numbers in [0,1] to FACS distribution
        facs_quant = splinefun(cumsum(hfdd)/sum(hfdd), hf$mids)

	# Score distribution
        xs = psi_clean[[lib]][,cn]
	# with small bins original distribution is better removed but this requires more points. breaks=1000 works well with 20,000 scores 
	settings$fnorm_hs_nbreaks = 500
	# make sure there is an empty bin in each end of the score distribution
        hs = hist(xs, breaks=calc_breaks(xs, settings$fnorm_hs_nbreaks, 2), plot=F)
        hsd = (1.0-settings$fnorm_prior_fac)*hs$density + settings$fnorm_prior_fac*dnorm(hs$mids, mean=mean(xs,na.rm=T), sd=sd(xs,na.rm=T))
	
	# Function that transform degron-score-distributed numbers to uniform distribution in [0,1]
        score_cdf = splinefun(hs$mids, cumsum(hsd)/sum(hsd))

        # Function that transforms score distributed numbers to FACS distributed numbers
        # this produces negative numbers, can I avoid that?
        score2facs = function(x){ facs_quant(score_cdf(x)) }

	cn_facs = sub("psi","fnorm",cn)
        psi_clean[[lib]][,cn_facs] = score2facs(psi_clean[[lib]][,cn])

	# Repeat with a score histogram that is weighted by the representation among cells
	cn_cells = sub("psi","rpm_sum",cn)
        hw = whist(xs, psi_clean[[lib]][,cn_cells], breaks=hs$breaks, plot=F)
        hwd = (1.0-settings$fnorm_prior_fac)*hw$density + settings$fnorm_prior_fac*dnorm(hw$mids, mean=mean(xs,na.rm=T), sd=sd(xs,na.rm=T))
        score_w_cdf = splinefun(hw$mids, cumsum(hwd)/sum(hwd))
        scorew2facs = function(x){ facs_quant(score_w_cdf(x)) }
	cn_wfacs = sub("psi","fwnorm",cn)
        psi_clean[[lib]][,cn_wfacs] = scorew2facs(psi_clean[[lib]][,cn])

	# # Plot
        # cn_rep = sub("_psi","",cn)
        # hr = hist(facs_quant(rand_unif), breaks=200, plot=F)
        # hsf = hist(psi_clean[[lib]][,cn_facs], breaks=200, plot=F)
	# hsfw = whist(psi_clean[[lib]][,cn_wfacs], psi_clean[[lib]][,cn_cells], breaks=hsf$breaks, plot=F)
        # plot(0,0,col=0, xlim=c(0.0,1.0), ylim=c(0.0,9.0), xlab="Score", ylab="Density",
	#      main=sprintf("%s %s, %d scores, %d breaks",lib,cn_rep,sum(! is.na(psi_clean[[lib]][,cn])),settings$fnorm_hs_nbreaks))
        # lines(hr$mids, hr$density, col=3)
        # lines(hf$mids, hf$density, col=2)
	# lines(hsf$mids, hsf$density, col=1)
	# lines(hsfw$mids, hsfw$density,col=4)
	# hsf_lab = sprintf("Transformed score distribution (%.2f - %.2f)", hsf$breaks[1], hsf$breaks[length(hsf$breaks)])
	# hsfw_lab = sprintf("Transf. score distribution weighted (%.2f - %.2f)", hsfw$breaks[1], hsfw$breaks[length(hsfw$breaks)])
	# legend("top", c(hsf_lab,hsfw_lab,"Target FACS distribution","Transformed unif. random numbers"),
	#        col=c(1,4,2,3), lty=1)
    }
}


# Normalize FACS distributions to control peptides
settings$norms = c(settings$norms, "fcnorm")
for (lib in names(psi_clean)) {
    # Find indices of control peptides
    is = c()	
    for (pep_name in settings$cnorm_stab) {
        is = c(is, which(grepl(pep_name, psi_clean[[lib]][,"name"])))
    }
    stopifnot(length(is) == length(unique(is)))
    id = c()
    for (pep_name in settings$cnorm_deg) {
        id = c(id, which(grepl(pep_name, psi_clean[[lib]][,"name"])))
    }
    stopifnot(length(id) == length(unique(id)))

    cns = colnames(psi_clean[[lib]])
    cns = cns[which(grepl("facs[0-9]_fnorm",cns))]
    cns = cns[order(cns)]
    for (cn in cns) {
	degron_psi = mean(psi_clean[[lib]][id,cn], na.rm=T)
	stabil_psi = mean(psi_clean[[lib]][is,cn], na.rm=T)
        cn_norm = sub("fnorm","fcnorm",cn)

        # # Min-max normalize
        # psi_clean[[lib]][,cn_norm] = (psi_clean[[lib]][,cn] - degron_psi)/(stabil_psi - degron_psi)

	# Scale using high fluorescence control. This is how the low-throughput validations are treated
        psi_clean[[lib]][,cn_norm] = psi_clean[[lib]][,cn] /stabil_psi

	# Store normalization values
	norm_list[[paste0(lib,"_",cn)]] = c(lib,cn,degron_psi,stabil_psi)
    }
}

# Data frame with normalization values
fcnorm = data.frame(name = names(norm_list),
                    lib = sapply(norm_list, "[", 1),
                    colname = sapply(norm_list, "[", 2),
                    degron_psi = sapply(norm_list, "[", 3),
                    stabil_psi = sapply(norm_list, "[", 4))


#####################################################################
##     Normalize each sorting replica to a uniform distribution
#####################################################################
## Alternativer double-sigmoidal kurve
##   dsig = function(x,a,x0){ 1.0/( (1+exp(-a*(x-x0))) * (1+exp(a*(x-1.0+x0))) ) }
##   alt: 0.5*( tanh((x-c1)/w1) - tanh((x-c2)/w2) )
##   plot(xx, dsig(xx,50,0.1), type="l")
##
settings$norms = c(settings$norms, "unorm")

settings$unorm_prior_fac = 0.01
settings$unorm_hs_nbreaks = settings$fnorm_hs_nbreaks


# quartz(width=8, height=6)
# rand_unif = runif(10^6, min=0.0, max=1.0)
# par(ask=T)
for (lib in names(psi_clean)) {
    # column names of replica PSI's
    cns = colnames(psi_clean[[lib]])
    cns = cns[which(grepl("facs[0-9]_psi",cns))]
    cns = cns[order(cns)]
    for (cn in cns) {
        xs = psi_clean[[lib]][,cn]
        hs = hist(xs, breaks=calc_breaks(xs, settings$unorm_hs_nbreaks, 2), plot=F)
	cn_cells = sub("psi","rpm_sum",cn)
        hw = whist(xs, psi_clean[[lib]][,cn_cells], breaks=hs$breaks, plot=F)
	
        hsd = (1.0-settings$unorm_prior_fac)*hs$density + settings$unorm_prior_fac*dnorm(hs$mids, mean=mean(xs,na.rm=T), sd=sd(xs,na.rm=T))
        hwd = (1.0-settings$unorm_prior_fac)*hw$density + settings$unorm_prior_fac*dnorm(hw$mids, mean=mean(xs,na.rm=T), sd=sd(xs,na.rm=T))
	
	# Function that transform degron-score-distributed numbers to uniform distribution in [0,1]
        # score_cdf = splinefun(hw$mids, cumsum(hwd)/sum(hwd))
        score_cdf = splinefun(hs$mids, cumsum(hsd)/sum(hsd))

	cn_facs = sub("psi","unorm",cn)
        psi_clean[[lib]][,cn_facs] = score_cdf(psi_clean[[lib]][,cn])
	
        # cn_rep = sub("_psi","",cn)
        # hsf = hist(psi_clean[[lib]][,cn_facs], breaks=200, plot=F)
	# hsfw = whist(psi_clean[[lib]][,cn_facs], psi_clean[[lib]][,cn_cells], breaks=hsf$breaks, plot=F)
        # plot(0,0,col=0, xlim=c(0.0,1.0), ylim=c(0.0,2.0), main=sprintf("%s %s",lib,cn_rep))
	# lines(hsf$mids, hsf$density, col=1)
	# lines(hsfw$mids, hsfw$density, col=2)
	# legend("top", c("Score distribution","Score distribution weighted"), col=c(1,2), lty=1)
    }
}


#####################################################################
##     Average PSI per library
#####################################################################
# This is mostly done to make plots, e.g. compare library-averaged scores of tiles that are present in pairs of libraries

for (lib in names(psi_clean)) {
    cns = colnames(psi_clean[[lib]])
    cns = cns[which(grepl("facs[0-9]_psi",cns))]
    for (norm in settings$norms) {
        # Calc average PSI for library
        psi_clean[[lib]][,sprintf("avg_%s",norm)] = apply(psi_clean[[lib]][,sub("psi",norm,cns)], MARGIN=1, mean, na.rm=T)
        psi_clean[[lib]][,sprintf("std_%s",norm)] = apply(psi_clean[[lib]][,sub("psi",norm,cns)], MARGIN=1, sd, na.rm=T)

        print(sprintf("Summary of SD between %s values of biological replicates for %s library",norm,lib))
        print(summary(psi_clean[[lib]][,sprintf("std_%s",norm)]))
    }
}

# Calculate cell representations (rpm_sum is representation out of 4e6 cells)
#   same genotype may be differently represented if present in different libraries
for (libname in names(psi_clean)) {
    ci = which(grepl("rpm_sum",names(psi_clean[[libname]])))
    psi_clean[[libname]][,"cells_avg"] = apply(psi_clean[[libname]][,ci], MARGIN=1, mean)
    psi_clean[[libname]][,"cells_sd"] = apply(psi_clean[[libname]][,ci], MARGIN=1, sd)
}


# Plot tiles that are in more libraries incl. controls
quartz(width=12, height=4)
par(mfrow=c(1,3), bg="white")
cn_plot = "avg_pnorm"
il1s = rep(seq(3),each=3)
il2s = rep(seq(3),times=3)
for (iil in which(il1s<il2s)) {
    lib1 = names(psi_clean)[il1s[iil]]
    lib2 = names(psi_clean)[il2s[iil]]
    
    common_tiles = intersect(psi_clean[[lib1]][,"dna"], psi_clean[[lib2]][,"dna"])
    print(sprintf("Libraries %s and %s have %d tiles in common",lib1,lib2,length(common_tiles)))

    ip1 = which(psi_clean[[lib1]][,"dna"] %in% common_tiles)
    ip2 = which(psi_clean[[lib2]][,"dna"] %in% common_tiles)

    plot(psi_clean[[lib1]][ip1,cn_plot], psi_clean[[lib2]][ip2,cn_plot], main=sprintf("%s libs %s and %s",cn_plot,lib1,lib2), xlab=paste("Score",lib1), ylab=paste("Score",lib2))
    abline(c(0,1))
}
quartz.save(sprintf("common_tiles_%s.png",cn_plot), type="png")


#####################################################################
##     Make a collected data.frame with unique DNA sequences
#####################################################################
# Merge libraries 
all_dna = unique(unname(unlist(sapply(psi_clean, "[", which(name_cols=="dna")))))
degron = data.frame(name=NA, dna = all_dna, aa=NA, lib="none")
for (lib in names(psi_clean)) {
    cns = colnames(psi_clean[[lib]])
    i = match(degron[,"dna"],psi_clean[[lib]][,"dna"])
    for (norm in settings$norms) {
        cns_norm = cns[which(grepl(sprintf("facs[0-9]_%s",norm),cns))] 
        cns_norm = cns_norm[order(cns_norm)]
        cns_deg = sprintf("%s_%s", lib, cns_norm)
        degron[,cns_deg] = psi_clean[[lib]][i,cns_norm]
    }
    # per library, store the average representation per genotype
    degron[,sprintf("%s_cells",lib)] = psi_clean[[lib]][i,"cells_avg"]

    # get names etc without overwriting previous names with NA
    ii = which(! is.na(i))
    degron[ii,"name"] = psi_clean[[lib]][i[ii],"name"]
    degron[ii,  "aa"] = psi_clean[[lib]][i[ii],  "aa"]
    degron[ii,"lib"] = paste(degron[ii,"lib"], lib, sep=",")
}
degron$lib = gsub("none,","",degron$lib)
degron = degron[order(degron$name),]

# store the largest representation of each tile
cns_cells = sprintf("%s_cells",names(psi_clean))
degron$cells_max = apply(degron[,cns_cells], MARGIN=1, max, na.rm=T)

# Average scores over replica per normalization
for (norm in settings$norms) {
    cns	= colnames(degron)
    cns = cns[which(grepl(sprintf("facs[0-9]_%s",norm),cns))]
    print(sprintf("Averaging %s scores over %d replicates:  %s",norm,length(cns),paste(cns,collapse=", ")))
    degron[,sprintf("n_%s",norm)]     = apply(degron[,cns], MARGIN=1, function(l){ sum(! is.na(l)) })
    degron[,sprintf("score_%s",norm)] = apply(degron[,cns], MARGIN=1, function(l){ mean(l, na.rm=T) })
    degron[,sprintf("std_%s",norm)]   = apply(degron[,cns], MARGIN=1, function(l){ sd(l, na.rm=T) })

    # print(sprintf("Distribution of numbers of %s measurements per tile",norm))
    # print(table(degron[,sprintf("n_%s",norm)]))

    print(sprintf("Degron %s score std. dev. summary:",norm))
    print(summary(degron[,sprintf("std_%s",norm)]))
}

# Remove DNA sequences without measurements
settings$min_measurements = 2
stopifnot(all( degron$n_pnorm == degron$n_cnorm ))
stopifnot(all( degron$n_pnorm == degron$n_fnorm ))
degron$n = degron$n_pnorm
print("Distribution of number of replicate measurements before filtering")
print(table(degron$n))
i_remove = which(degron$n < settings$min_measurements)
print(sprintf("Discarding %d tiles with fewer than %d measurements", length(i_remove), settings$min_measurements))
degron = degron[-c(i_remove),]

print(sprintf("Final coverages is %d degron scores out of %d (%.2f%%)", nrow(degron), length(all_dna), nrow(degron)/length(all_dna)*100))

# Transform the selected normalization to the selected FACS distribution
selected_norm = "pnorm"
selected_facs = "Even_bio1_facs1"
cn_score = paste0("score_",selected_norm)
cn_n = paste0("score_",selected_norm)

fi = which(facs_table$rep==selected_facs)
xf = facs_events[[facs_table[fi,"name"]]][,facs_set$cn_ratio]
# truncate distribution to avoid very flat end regions
xf = xf[which(xf > 2^(-13) & xf < 1.0)]
fbreaks = 2^(seq(-13, 0, length.out=1000))
hf = hist(xf, breaks=fbreaks, plot=F)
hfd = (1.0-settings$fnorm_prior_fac)*hf$density + settings$fnorm_prior_fac*dnorm(hf$mids, mean=mean(xf), sd=sd(xf))
bin_widths = hf$breaks[2:length(hf$breaks)] - hf$breaks[1:(length(hf$breaks)-1)]
hfdd = hfd*bin_widths
facs_quant = splinefun(cumsum(hfdd)/sum(hfdd), hf$mids, method="hyman")

# Score distribution
xs = degron[,cn_score]
hs = hist(xs, breaks=calc_breaks(xs, 1000, 2), plot=F)
hsd = (1.0-settings$fnorm_prior_fac)*hs$density + settings$fnorm_prior_fac*dnorm(hs$mids, mean=mean(xs,na.rm=T), sd=sd(xs,na.rm=T))
score_cdf = splinefun(hs$mids, cumsum(hsd)/sum(hsd), method="hyman")

# Transform
score2facs = function(x){ facs_quant(score_cdf(x)) }
score2facs_deriv = function(x){ facs_quant(score_cdf(x), deriv=1) * score_cdf(x, deriv=1) }
cn_std = paste0("std_",selected_norm)
cn_score_facs = paste0("score_",selected_norm,"_facs")
cn_std_facs = paste0("std_",selected_norm,"_facs")
degron[,cn_score_facs] = score2facs(degron[,cn_score])
degron[,cn_std_facs]   = degron[,cn_std] * score2facs_deriv(degron[,cn_score])

# plot of the final FACS transform
quartz(width=12, height=6)
par(mfcol=c(1,2))
xx=seq(-0.1,1.1,.001)
plot(xx, score2facs(xx), type="l", lwd=2, ylim=c(0,1))
lines(xx, score2facs_deriv(xx)/2, col=2)
legend("topleft", c("score2facs function","First derivative"), lty=1, col=c(1,2))
hsf = hist(degron[,cn_score_facs], breaks=100, plot=F)
plot(0,0,col=0, xlim=c(0,1), ylim=c(0,8), xlab="FACS score", ylab="Density")
lines(hf$mids, hf$density, lwd=2, col=1)
lines(hsf$mids, hsf$density, lwd=2, col=2)
legend("top", c("FACS","Scores"), lty=1, col=c(1,2))
quartz.save("facs_transform.png", type="png")

# The degron score we have used previously where a degron score of 1 mean degron
degron$deg_score = 1 - degron[,cn_score]


# Dump everything
save(degron, psi_clean, ctrl, cnorm, pnorm, settings, file="degrons.rda")

# Unique measruements
df =   degron[,c("name","dna","aa","n",   "deg_score",      cn_std,    cn_score_facs,   cn_std_facs)]
colnames(df) = c("name","dna","aa","n","degron_score","degron_std","abundance_score","abundance_std")
write.csv(df, row.names=F, quote=F, file="degrons.csv")

