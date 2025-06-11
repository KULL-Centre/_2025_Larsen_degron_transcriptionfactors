options(width=160, digits=4, stringsAsFactors=F)

args = commandArgs(trailingOnly=TRUE)
if (interactive()) {
    samples_file = "samples_all.csv"
    files = c("../counts/O-2-1-3-1_S19_L001_counts.txt","../counts/O-2-1-3-1_S19_L002_counts.txt","../counts/O-2-1-3-1_S19_L003_counts.txt")
} else if (length(args) < 2) {
    print("")
    print("usage: Rscript merge_and_map.r  <samples.csv>  <counts1.txt>  [counts2.txt  ...]")
    quit(save="no")
} else {
    samples_file = args[1]
    files = args[2:length(args)]
}
print(sprintf("Args: %s",paste0(args, collapse=" ")))


##
## Settings
##
settings = list()

# number of read counts required to consider DNA covered in summary table
settings$coverage_min_counts = 20

# store unmapped reads for later analysis
settings$store_unmapped = TRUE
settings$unmapped_min_reads = 5
settings$unmapped_good_lengths = c(90,72,66)


##
## Helper functions
##

dna2aa = list()
# T          T                      C                    A                     G
dna2aa[["TTT"]] = "F"; dna2aa[["TCT"]] = "S"; dna2aa[["TAT"]] = "Y"; dna2aa[["TGT"]] = "C" # T
dna2aa[["TTC"]] = "F"; dna2aa[["TCC"]] = "S"; dna2aa[["TAC"]] = "Y"; dna2aa[["TGC"]] = "C" # C
dna2aa[["TTA"]] = "L"; dna2aa[["TCA"]] = "S"; dna2aa[["TAA"]] = "*"; dna2aa[["TGA"]] = "*" # A
dna2aa[["TTG"]] = "L"; dna2aa[["TCG"]] = "S"; dna2aa[["TAG"]] = "*"; dna2aa[["TGG"]] = "W" # G
# C
dna2aa[["CTT"]] = "L"; dna2aa[["CCT"]] = "P"; dna2aa[["CAT"]] = "H"; dna2aa[["CGT"]] = "R" # T
dna2aa[["CTC"]] = "L"; dna2aa[["CCC"]] = "P"; dna2aa[["CAC"]] = "H"; dna2aa[["CGC"]] = "R" # C
dna2aa[["CTA"]] = "L"; dna2aa[["CCA"]] = "P"; dna2aa[["CAA"]] = "Q"; dna2aa[["CGA"]] = "R" # A
dna2aa[["CTG"]] = "L"; dna2aa[["CCG"]] = "P"; dna2aa[["CAG"]] = "Q"; dna2aa[["CGG"]] = "R" # G
# A
dna2aa[["ATT"]] = "I"; dna2aa[["ACT"]] = "T"; dna2aa[["AAT"]] = "N"; dna2aa[["AGT"]] = "S" # A
dna2aa[["ATC"]] = "I"; dna2aa[["ACC"]] = "T"; dna2aa[["AAC"]] = "N"; dna2aa[["AGC"]] = "S" # C
dna2aa[["ATA"]] = "I"; dna2aa[["ACA"]] = "T"; dna2aa[["AAA"]] = "K"; dna2aa[["AGA"]] = "R" # A
dna2aa[["ATG"]] = "M"; dna2aa[["ACG"]] = "T"; dna2aa[["AAG"]] = "K"; dna2aa[["AGG"]] = "R" # G
# G
dna2aa[["GTT"]] = "V"; dna2aa[["GCT"]] = "A"; dna2aa[["GAT"]] = "D"; dna2aa[["GGT"]] = "G" # A
dna2aa[["GTC"]] = "V"; dna2aa[["GCC"]] = "A"; dna2aa[["GAC"]] = "D"; dna2aa[["GGC"]] = "G" # C
dna2aa[["GTA"]] = "V"; dna2aa[["GCA"]] = "A"; dna2aa[["GAA"]] = "E"; dna2aa[["GGA"]] = "G" # A
dna2aa[["GTG"]] = "V"; dna2aa[["GCG"]] = "A"; dna2aa[["GAG"]] = "E"; dna2aa[["GGG"]] = "G" # G

translate = function(dna) {
    n = nchar(dna)
    codons = substring(dna, seq(1, n-2, by=3), seq(3, n, by=3))
    if (length(codons)*3 != n) return(NA)
    paste0(dna2aa[codons], collapse="")
}


# read a sequence file (.seq) - a text file with one sequence per line
# return a data.frame with 2 columns: name and dna
read.seqfile = function(filename) {
    df = read.table(filename)
    if (ncol(df) == 1) {
        colnames(df) = c("dna")
    } else if (ncol(df) == 2) {
        colnames(df) = c("name","dna")
    } else {
        print(sprintf("ERROR: Seq-file %s does not have 1 or 2 columns", filename))
	return(NULL)
    }
    return(df)
}

# read a csv file with sequences
# return a data.frame with 2 columns: name and dna
read.seqcsv = function(filename) {
    df = read.csv2(filename)
    if ("ensembl" %in% colnames(df) & ! "name" %in% colnames(df)) {
        colnames(df)[which(colnames(df)=="ensembl")] = "name"
    }
    if (! "name" %in% colnames(df)) {
        print(sprintf("ERROR: Cannot find name column in %s",filename))
	return(NA)
    }
    if (! "tile_number" %in% colnames(df)) {
        print(sprintf("ERROR: Cannot find tile_number column in %s",filename))
	return(NA)
    }
    df$name = sprintf("%s_tile%04d",df$name,df$tile_number)
    if (! "dna" %in% colnames(df)) {
        print(sprintf("ERROR: Cannot find dna column in %s",filename))
	return(NA)
    }
    return(df[,c("name","dna")])
}

###
###  Libraries
###
# Read library files
# These are made from which(lib$even1), which(lib$even2), etc and should be non-redundant and including controls
raw = list()
raw[["Even"]] = read.seqcsv("../library/transfac90_lib1_even.csv")
raw[["Odd"]] = read.seqcsv("../library/transfac90_lib2_odd.csv")
raw[["CT"]] = read.seqcsv("../library/transfac90_lib3_ct.csv")

# check library and translate to amino acid sequences
for (lib in names(raw)) {
    # DNA sequences should be unique per library
    stopifnot(length(raw[[lib]][,"dna"]) == length(unique(raw[[lib]][,"dna"])))
    
    # translate
    raw[[lib]][,"aa"] = sapply(raw[[lib]][,"dna"], translate)    
    print(sprintf("Check if library %s contains unique AA sequences: %s",lib,length(raw[[lib]][,"aa"]) == length(unique(raw[[lib]][,"aa"]))))
}


###
###  Control sequences
###
# autodetect control sequences from libraries
ctrl_dna = raw[[1]][,"dna"]
for (ri in 2:length(raw)) { ctrl_dna = intersect(ctrl_dna, raw[[ri]][,"dna"]) }
print(sprintf("Found %d peptides common to all %d libraries", length(ctrl_dna), length(raw)))
if (length(ctrl_dna) > 0) { print(raw[[1]][match(ctrl_dna,raw[[1]][,"dna"]),"name"]) }
# get names from first library
# ctrl = data.frame(names = raw[[1]][match(ctrl_dna,raw[[1]][,"dna"]),"name"], dna = ctrl_dna)

ctrl_fn = "controls.csv"
print(sprintf("Read controls from %s and add to all libraries", ctrl_fn))
ctrl = read.csv(ctrl_fn)
stopifnot(all( c("name","dna","aa") %in% colnames(ctrl) ))
for (libname in names(raw)) {
    # Check if controls are already in library
    i_redun = which(raw[[libname]][,"dna"] %in% ctrl[,"dna"])
    if (length(i_redun) > 0) {
        for (ii in i_redun) {
	    ictrl = which(ctrl[,"dna"] == raw[[libname]][ii,"dna"])
            print(sprintf("WARNING: DNA sequence of control peptide %s is already in lib %s as %s - using control name", ctrl[ictrl,"name"], libname, raw[[libname]][ii,"name"]))
	}
	raw[[libname]] = raw[[libname]][-c(i_redun),]
    }
    # Add all controls 
    raw[[libname]] = rbind(raw[[libname]],ctrl)
}

print(sprintf("Found %d control sequences",nrow(ctrl)))


# translate DNA to protein
for (lib in names(raw)) { raw[[lib]][,"aa"] = sapply(raw[[lib]][,"dna"], translate)  }

# check that all library files have the same number of columns
name_cols = colnames(raw[[1]])
stopifnot(all(sapply(raw,ncol) == length(name_cols)))

# Read information about samples
samples = read.csv(samples_file)


###
### Read count files
###
raw_read_counts = list()
if (settings$store_unmapped) { unmapped_raw = list() }

# for each file, map counts to library and put them in a data frame
for (file in files) {
    # check filename and extract file_id
    fl = strsplit(file, "/")[[1]]
    stopifnot(substr(fl[length(fl)], nchar(fl[length(fl)])-10, nchar(fl[length(fl)])) == "_counts.txt")
    file_id = substr(fl[length(fl)], 1, nchar(fl[length(fl)])-11)

    # determine sample index
    if ( grepl("_L00",file_id) ) {
        si = which(samples$file == substr(file_id, 1, nchar(file_id)-5))
    } else {
	si = which(samples$file == file_id)
    }
    stopifnot(length(si) == 1)

    # library of sample
    lib = samples[si,"lib"]
    stopifnot(lib %in% names(raw))
    stopifnot(! file_id %in% colnames(raw[[lib]]))

    # read counts file: col 1 should be the DNA sequence, col 1 the counts
    cf = read.table(file)
    colnames(cf) = c("dna","counts","rpm")

    # map file read counts to library
    i_mapped = match(raw[[lib]][,"dna"],cf$dna)
    raw[[lib]][,file_id] = cf[i_mapped,"counts"]
    
    # set un-observed variants to zero counts
    i_na = is.na( raw[[lib]][,file_id] )
    raw[[lib]][i_na,file_id] = 0

    # store number of unique variants and total read counts
    n_unq_dna = nrow(cf)
    n_raw_counts = sum(cf$counts)
    raw_read_counts[[file_id]] = c(n_unq_dna, n_raw_counts)

    # Report
    n_mapped_dna = sum(raw[[lib]][,file_id] > 0)
    n_mapped_counts = sum(raw[[lib]][,file_id])
    # n_filtered_counts = sum(cf$counts)
    # print(sprintf("Mapped %d dna seq covering %.2f%% of library and %d of %d length-filtered read counts (%.2f%%)",
    #               n_mapped_dna, n_mapped_dna/nrow(raw[[lib]])*100,
    # 		  n_mapped_counts, n_filtered_counts, n_mapped_counts/n_filtered_counts*100))

    print(sprintf("Mapped %d out of %d raw counts (%.2f%%) from %s",n_mapped_counts,n_raw_counts,n_mapped_counts/n_raw_counts*100,file_id))

    # store unmapped redas for later analysis
    if (settings$store_unmapped) {
        i_above_threshold = which(cf$counts >= settings$unmapped_min_reads)
        unmapped_raw[[file_id]] = cf[setdiff(i_above_threshold,i_mapped),c("dna","counts")]
    }
    
    # print(table(unlist(strsplit(cf[,"dna"],""))))
}


# merge lane counts
counts = list()
unmapped_counts = list()
lane_names = list()
samples$obs = FALSE
for (lib in names(raw)) {
    if (ncol(raw[[lib]]) <= length(name_cols)) {
        print(sprintf("No data for %s - library will not be considered further",lib))
	next
    }
    
    # init counts data frame
    counts[[lib]] = raw[[lib]][,name_cols]
    
    # column names of columns with reads
    cn = colnames(raw[[lib]])[(length(name_cols)+1):ncol(raw[[lib]])]
    
    # column names with lane info removed
    cn_nolane = unique( sapply(cn, function(s){ if (grepl("_L00",s)) {substr(s,1,nchar(s)-5)} else {s} }) )
    
    for (sample_name in cn_nolane) {
        stopifnot(sample_name %in% samples$file)
	si = which(samples$file==sample_name)
	stopifnot(! samples[si,"obs"])
	samples[si,"obs"] = TRUE
	
        sample_lane_names = cn[which(grepl(sample_name, cn))]
	counts[[lib]][,sample_name] = apply(raw[[lib]][,sample_lane_names], MARGIN=1, sum)

        # consider merging raw_read_counts[[file_id]] = c(n_unq_dna, n_raw_counts)

	# report correlations of raw[[lib]][,sample_lane_names]
	lane_names[[sample_name]] = sample_lane_names
	print(sprintf("Merging %d lanes into %s", length(sample_lane_names), sample_name))

        # merge lanes for unmapped read counts
	if (settings$store_unmapped) {
            sample_unmapped_dna = unique(unlist(sapply(unmapped_raw[sample_lane_names], "[", "dna")))

            unmapped_counts[[sample_name]] = data.frame(dna = sample_unmapped_dna)
	    for (sn in sample_lane_names) {
	        unmapped_counts[[sample_name]][,sn] = unmapped_raw[[sn]][match(sample_unmapped_dna,unmapped_raw[[sn]][,"dna"]),"counts"]
	    }
	    # print(sprintf("Merged lanes for %s resulting in %d unmapped reads from %d unique DNA sequences",
	    #               sample_name, sum(unmapped_counts[[sample_name]][,"sum"]), length(sample_unmapped_dna)))
	}
    }
}

# counts per sample, info on lanes and unmapped reads not present
save(counts, samples, name_cols, ctrl, settings, file="counts.rda")

# data on unmapped reads and lanes
if (settings$store_unmapped) {
    save(unmapped_counts, raw_read_counts, lane_names, file="counts_unmapped.rda")
}

# table of read counts
all_lib_dna = unique(unlist(sapply(raw, "[", "dna")))
counts_summary = data.frame(sample=NULL, lib=NULL, coverage = NULL, ctrl_cover=NULL,
                            all=NULL, mapped=NULL, controls=NULL, mapped_other_lib=NULL, unmapped=NULL, unmapped_bad_len=NULL)

for (sn in sort(names(unmapped_counts))) {
    si = which(sn == samples$file)
    lib = samples[si,"lib"]
    sample_lane_names = lane_names[[sn]]

    # i_ctrl = which(counts[[lib]][,"dna"] %in% ctrl$dna) 
    i_ctrl = which(counts[[lib]][,"name"] %in% ctrl$name)
    stopifnot(length(i_ctrl) == nrow(ctrl))

    mapped_other_lib = 0
    unmapped = 0
    unmapped_bad_len = 0

    if (settings$store_unmapped) {
        # unmapped DNA that mappes to other libs
        imo = which(unmapped_counts[[sn]][,"dna"] %in% all_lib_dna)
        unmapped_counts[[sn]][,"sum"] = apply(unmapped_counts[[sn]][,sample_lane_names], MARGIN=1, sum, na.rm=T)
        mapped_other_lib = sum(unmapped_counts[[sn]][imo,"sum"])

        # completely unmapped DNA
        ium = setdiff(seq(nrow(unmapped_counts[[sn]])), imo)
        completely_unmapped = sum(unmapped_counts[[sn]][ium,"sum"])
        unmapped = sum(unmapped_counts[[sn]][ium,"sum"])

        # unmapped DNA with bad length, i.e. length not in
        unmapped_counts[[sn]][,"length"] = sapply(unmapped_counts[[sn]][,"dna"], nchar)
        ibl = ium[ which(! unmapped_counts[[sn]][ium,"length"] %in% settings$unmapped_good_lengths) ]
        unmapped_bad_len = sum(unmapped_counts[[sn]][ibl,"sum"])
    }
    
    df = data.frame(sample = sn, lib = lib,
                    coverage = sum(counts[[lib]][,sn] >= settings$coverage_min_counts)/nrow(counts[[lib]])*100,
		    ctrl_cover = sum(counts[[lib]][i_ctrl,sn] >= settings$coverage_min_counts)/nrow(ctrl)*100,
                    all = sum(sapply(raw_read_counts[sample_lane_names], "[", 2)),
                    mapped = sum(counts[[lib]][,sn]),
	            controls = sum(counts[[lib]][i_ctrl,sn]),
	            mapped_other_lib = mapped_other_lib,
		    unmapped = unmapped,
		    unmapped_bad_len = unmapped_bad_len)
    counts_summary = rbind(counts_summary, df)
}
print("")
print("Read counts summary")
print(counts_summary)
write.csv(counts_summary, file="counts_summary.csv")

pct_summary = counts_summary
cns = colnames(counts_summary)[5:ncol(counts_summary)]
for (ri in seq(nrow(pct_summary))) {
    pct_summary[ri,cns] = pct_summary[ri,cns]/pct_summary[ri,"all"]*100
}
print("")
print("Read percent summary")
print(pct_summary)

