options(width=160)

args = commandArgs(trailingOnly=TRUE)
if (interactive()) {
    infiles = c("data_230428/export_81835_Q1 BV421-A- , PE-Texas Red-A+.csv")
    infiles = c("data/export_03-05-2023_Even Sort 3_T1_001_Q1 BV421-A- , PE-Texas Red-A+.csv")
} else if (length(args) < 1) {
    print("")
    print("usage: Rscript structure_facs.r  <facs1.csv>  [facs2.csv  ...]")
    quit(save="no")
} else {
    # infiles = args[2:length(args)]
    infiles = args
}
# print(sprintf("Args: %s",paste0(args, collapse=" ")))


parse_facs_csv = function(filename, n_cols=20) {
    d = read.table(filename, sep=",", row.names=NULL, fill=T, col.names=sprintf("col%.02d",seq(n_cols)))

    # keep reading until an excess of columns are read (read.table(fill=T)) only looks at first 5 lines)
    while (! all(is.na(d[,n_cols]) | n_cols > 200) ) {
        n_cols = n_cols+10
        d = read.table(filename, sep=",", row.names=NULL, fill=T, col.names=sprintf("col%.02d",seq(n_cols)))
    }

    facs = list()
    facs$sample = d[which(grepl("TUBE", d$col01)), "col02"]
    facs$date = d[which(grepl("DATE", d$col01)), "col02"]

    # find the number of rows with settings assuming a header row. Look for the row from which everything in column 1 is numeric
    # print("Search end of settings lines")
    settings_rows = 0
    while (suppressWarnings( any(is.na(as.numeric(d[(settings_rows+2):nrow(d),"col01"]))) & settings_rows < 500)) {
        settings_rows = settings_rows +1
    }
    stopifnot(settings_rows < 500)

    # # print("Search number of settings cloumns")
    # settings_cols = n_cols
    # while (all(is.na(d[1:settings_rows,settings_cols])) & settings_cols > 0) {
    #     settings_cols = settings_cols -1
    # }
    # stopifnot(settings_cols > 0)
    # settings = d[1:settings_rows,1:settings_cols]
    
    facs$events = read.table(filename, sep=",", skip=settings_rows, header=T)
    return(facs)
}


##
## Settings
##
facs_set = list()
facs_set$infiles = infiles

facs_events = list()
facs_table = list()
for (infile in infiles) {
    print(sprintf("Reading %s",infile))
    facsfile = parse_facs_csv(infile)
    print(sprintf("    Found sample %s from %s with %d events and %d channels",facsfile$sample,facsfile$date,nrow(facsfile$events),ncol(facsfile$events)))
    name = facsfile$sample
    if ( name %in% names(facs_events)) {
        i = 1
	name = sprintf("%s_%02d", facsfile$sample, i)
	while (name %in% names(facs_events)) {
	    i = i+1
	    name = sprintf("%s_%02d", facsfile$sample, i)
	}
	print(sprintf("    Renaming %s to %s", facsfile$sample, name))
    }
    facs_table[[name]] = c(facsfile$sample, facsfile$date, nrow(facsfile$events), ncol(facsfile$events))
    facs_events[[name]] = facsfile$events
}

# channels we are using
facs_set$cn_gfp = "GFP.A"
facs_set$cn_cherry = "PE.Texas.Red.A"
facs_set$cn_ratio = "Derived....GFP.A.PE.Texas.Red.A"

# only store used channels
cns = c(facs_set$cn_gfp, facs_set$cn_cherry, facs_set$cn_ratio)
for (sample in names(facs_events)) {
    stopifnot(all( cns %in% colnames(facs_events[[sample]]) ))
    facs_events[[sample]] = facs_events[[sample]][,cns]
}

# reformat as table
facs_table = data.frame(name = names(facs_table),
                        sample = sapply(facs_table, "[[", 1),
                        date = sapply(facs_table, "[[", 2),
                        events = sapply(facs_table, "[[", 3))
write.csv(facs_table, quote=F, row.names=F, file="facs_table.csv")

save(facs_events, facs_table, facs_set, file="facs.rda")

# plot all distributions
hl = list()
x_min = 0; x_max = 0; y_max = 0
for (sample in names(facs_events)) {
    h = hist(facs_events[[sample]][,facs_set$cn_ratio], breaks=200, plot=F)
    x_min = min(c(x_min, h$mids))
    x_max = max(c(x_max, h$mids))
    y_max = max(c(y_max, h$density))
    hl[[sample]] = h
}

quartz(width=15, height=6)
par(mfrow=c(1,2), bg="white")

# plot(0,0,col=0, xlim=c(x_min,x_max), ylim=c(0,y_max), xlab="GFP/mCherry", ylab="Density")
plot(0,0,col=0, xlim=c(-0.1,0.7), ylim=c(0,y_max), xlab="GFP/mCherry", ylab="Density")
for (i in seq_along(facs_events)) {
    sample = names(facs_events)[i]
    lines(hl[[sample]]$mids, hl[[sample]]$density, col=i, lty=i%/%8+1)
}
legend("topright", names(hl), col=seq(i), lty=seq(i)%/%8+1, ncol=3, cex=.8)
plot(0,0,col=0, xlim=c(-3,0), ylim=c(0,1.5), xlab="GFP/mCherry, log10", ylab="Density")
for (i in seq_along(facs_events)) {
    sample = names(facs_events)[i]
    ipos = which(facs_events[[sample]][,facs_set$cn_ratio] > 0.0)
    if (length(ipos) < nrow(facs_events[[sample]])) { print(sprintf("Not plotting %d events with non-positive fluorescense for %s",nrow(facs_events[[sample]])-length(ipos),sample)) }
    h = hist(log10(facs_events[[sample]][ipos,facs_set$cn_ratio]), breaks=200, plot=F)
    lines(h$mids, h$density, col=i, lty=i%/%8+1)
}

quartz.save("all_distributions.png", type="png")

