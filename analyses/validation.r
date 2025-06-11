options(width=160)
runavg = function(x, n = 3) { as.vector(filter(x, rep(1 / n, n), sides = 2)) }

# this has 82 points which is fewer than the 89 in prev plot because some tiles didn't replicate 
ltp = read.csv("low_throughput_validations_TF_processed_data.csv", sep=";", dec=",")
# rename for systematic generation of column names
colnames(ltp)[which(colnames(ltp)=="New_Standard_Error")] = "Normalized_Error"

# get tile scores
htp = read.csv(gzfile("../annotations/transfac_tiles.csv.gz"))

htp$nres = nchar(htp$aa)
htp_use = which(htp$nres==30)

i_htp2ltp = match(ltp$dna,htp$dna)
ltp$degron_score = htp[i_htp2ltp,"degron_score"]
ltp$degron_std   = htp[i_htp2ltp,"degron_std"]

ltp$abundance_score = htp[i_htp2ltp,"abundance_score"]
ltp$abundance_std   = htp[i_htp2ltp,"abundance_std"]

htp$logabun_score = log10( htp[,"abundance_score"] )
htp$logabun_std   = 1/(log(10)*htp[,"abundance_score"]) * htp[,"abundance_std"]
ltp$logabun_score = htp[i_htp2ltp,"logabun_score"]
ltp$logabun_std   = htp[i_htp2ltp,"logabun_std"]
ltp$Log_Normalized_Mean = -log10( ltp$Normalized_Mean )
ltp$Log_Normalized_Error = 1/(log(10)*ltp[,"Normalized_Mean"]) * ltp[,"Normalized_Error"]

plot_validation = function(x_cbn, y_cbn="Normalized", ltp_fac=2.5, ltp_ticks=NA, ltp_col="brown3", xlim=NA, ylim=NA, lylab="Low-throughput score", ...) {
    cn_score = sprintf("%s_score",x_cbn)
    cn_std = sprintf("%s_std",x_cbn)
    ycn_score = sprintf("%s_Mean",y_cbn)
    ycn_std = sprintf("%s_Error",y_cbn)

    x = ltp[,cn_score]
    dx = ltp[,cn_std]
    y = ltp[,ycn_score]
    dy = max(ltp[,ycn_std], 1E-5)
    
    x_all = htp[,cn_score]
    h_htp = hist(x_all, breaks=200, plot=F)
    d_htp = density(x_all, bw=1E-2, na.rm=T)
    
    if (any(is.na(xlim))) {
        xlim=c(min(htp[,cn_score],na.rm=T),max(htp[,cn_score],na.rm=T))
	print(sprintf("Setting xlim to %s",paste0(xlim,collapse=", ")))
    }
    if (any(is.na(ylim))) {
        ylim=c(0,max(c(y*ltp_fac,d_htp$y),na.rm=T))
	print(sprintf("Setting ylim to %s",paste0(ylim,collapse=", ")))
    }
    if (any(is.na(ltp_ticks))) {
        ltp_ticks = seq(ylim[1]/ltp_fac, floor(ylim[2]/ltp_fac*10)/10, length.out=4 )
	print(sprintf("Setting LTP ticks at %s",paste0(ltp_ticks,collapse=", ")))
    }
    
    rp = cor(x, y, method="pearson")
    print(sprintf("Pearson of %d low-throughput validations is %.4f (R-sq %.4f)",nrow(ltp),rp,rp^2))
    
    # quartz(width=8, height=5)
    pdf(sprintf("validation_%s.pdf",x_cbn), width=4, height=2.5, pointsize=7)
    par(mar=c(5,4,2,5)+.1, bg="white")
    plot(0,0,col=0, xlim=xlim, ylim=ylim, ylab="Density (run. avg. size 5)", ...)
    # plot(0,0,col=0, xlim=xlim, ylim=ylim, ylab="Density (kernel bandwidth 0.01)", ...)
    axis(side=4, at=ltp_ticks*ltp_fac, labels=ltp_ticks, col=ltp_col, col.axis=ltp_col)
    mtext(lylab, side=4, col=ltp_col, line=2.5)
    lines(h_htp$mids, runavg(h_htp$density,5), lwd=2, col="darkgray")  
    # lines(d_htp, lwd=2, col="darkgray")  
    points(x, y*ltp_fac, pch=16, col=ltp_col, cex=.8)
    arrows(x, y0=(y-dy)*ltp_fac, y1=(y+dy)*ltp_fac, code=3, angle=90, length=.01, col=ltp_col)
    arrows(x0=x-dx, x1=x+dx, y0=y*ltp_fac, code=3, angle=90, length=.01, col=ltp_col)

    # HTP as function of LTP because errors are (mostly) on HTP
    fit = lm(x~y, weights=1/dx^2)
    slope = 1/coef(fit)[2]
    intercep = -1.0*coef(fit)[1]/coef(fit)[2]
    abline(intercep*ltp_fac, slope*ltp_fac, lty=2, lwd=2, col="darkgray")
    x_calc = coef(fit)[1] + y*coef(fit)[2]
    chi_sq = sum( (x-x_calc)^2 / dx^2 )
    chi_sq_red = chi_sq / (nrow(ltp)-2)
    print(sprintf("Under a weighted model: HTP = %.2f + LTP x %.2f", coef(fit)[1], coef(fit)[2]))
    print(sprintf("Fit in HTP space has chi-squared %.2f and reduced chi-square: %.2f",chi_sq,chi_sq_red))
    print(sprintf("Inverse model: LTP = %.2f + HTP x %.2f", intercep, slope))

    dev.off()
    # quartz.save(sprintf("validation_%s.png",x_cbn), type="png")
}


print(sprintf("Summary of %d high-throughput degron scores:",length(htp_use)))
print(summary(htp$degron_score))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.37620  0.07279  0.38021  0.44452  0.72561  1.43094 )
plot_validation("degron", xlim=c(-0.3,1.4), ltp_fac=2.5,
                xlab="High-throughput degron score", lylab="Low-throughput score, GFP normalized")

print(sprintf("Summary of %d high-throughput abundance scores:",length(htp_use)))
print(summary(htp$abundance_score))
plot_validation("abundance", xlim=c(-.05,0.6),ltp_fac=10, ltp_ticks=c(0,.1,.2,.3,.4,.5),
                xlab="High-throughput abundance score", lylab="Low-throughput score, GFP normalized")

print(sprintf("Summary of %d high-throughput abundance scores:",length(htp_use)))
print(summary(htp$logabun_score))
plot_validation("logabun", "Log_Normalized", ltp_fac=0.4, ltp_ticks=seq(0,3,0.5),
                xlab="High-throughput log abundance score", lylab="Low-throughput -log score, GFP normalized")





