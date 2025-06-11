options(width=160)

tfp = read.csv(gzfile("../annotations/transfac_proteins.csv.gz"))

dir.create("papmaps")

# genename = "HSF1"
for (genename in tfp$gene) {

ip = which(tfp$gene == genename)
stopifnot( length(ip) == 1 )
print(sprintf("%-4d plot %6s of %5d residues",ip,genename,nchar(tfp[ip,"aa"])))

# Get delta-pap scores
dpap = read.table(sprintf("prism/prism_qcdpred_%s.txt",genename), header=T)

plot_cn = "delta_score"
plot_lab = expression(Delta*"PAP-score")
filename = sprintf("papmaps/papmap_%s",genename)

dpap$wt = substr(dpap$variant, 1, 1)
dpap$mut = substr(dpap$variant, nchar(dpap$variant), nchar(dpap$variant))
dpap$resi = as.numeric(substr(dpap$variant, 2, nchar(dpap$variant)-1))

wt = list()
wt$aa = tfp[ip,"aa"]
wt_seq = strsplit(wt$aa,"")[[1]]
nres = length(wt_seq)

stopifnot(all(wt_seq[dpap$resi] == dpap$wt))

# Order of amino acids
aa2int = seq(20)
# names(aa2int) = strsplit("ACDEFGHIKLMNPQRSTVWY","")[[1]]
names(aa2int) = strsplit("CDEKRHNQAGSTVMLIFYWP","")[[1]]

# single mutants to plot
scale = 10
i = which(dpap$mut %in% names(aa2int))
pmap = data.frame(subst = dpap[i,"variant"], resi=dpap[i,"resi"], mut=dpap[i,"mut"])
pmap[,plot_cn] = dpap[i,plot_cn] *scale

grad_range_min = -3
grad_range_max = 3
# grad_range_center = round(pmap[iwt,plot_cn])
grad_range_center = 0
ns = (grad_range_center - grad_range_min)*100
nd = (grad_range_max - grad_range_center)*100
col_not_in_lib = "white"
col_native     = "#454545"
# -----------------------------------------------------------------------------------------------
# Variant :   destabilizing           neutral               stabilizing

#     col_destab = "#079700";  col_neutral = "#ffff00";    col_stab = "#ff5a6b"
      col_destab = "#005bc7";  col_neutral = "#ffffff";    col_stab = "#c40031"
#     col_destab = "#3094ff";  col_neutral = "#ffff00";    col_stab = "#ff5a6b"
      
# Position: stable/conserved      neutral/tolerant     unstable/engineering potential
# If bith gradients pass through brown-black the difference between slightly stab and destab disapears
# white as tolerant works well because the small effects are easier to distinguis and I get a saturation gradient fom
#     tolerant to non- tolerant if both stab and destab are saturated colors
# -----------------------------------------------------------------------------------------------

col_grad = c(colorRampPalette( c(col_stab, col_neutral), space="rgb")(ns), colorRampPalette( c(col_neutral,col_destab), space="rgb")(nd))
col_breaks = seq(grad_range_min, grad_range_max, length.out=ns+nd+1)
# gene     nres  width   nres/14+2
# NKX6-2    275     24          22
# HSF1      530     40          40
height = ceiling(length(wt_seq)/12+2)

# quartz(width=5, height=14)
# tiff("heatmap_300dpi.tiff", width=8, heigh=20, units="cm", res=300, pointsize=9, compression="lzw")
# jpeg(paste0(filename,".jpg"), width=8, heigh=height, units="cm", res=300, pointsize=9, quality=90)
pdf(paste0(filename,".pdf"), width=8/2.54 , heigh=height/2.54, pointsize=9)
# layout(matrix(c(1,1,3,2), ncol=2, ), width=c(4,1), height=c(3,1))
layout(matrix(c(1,1,3,2), ncol=2, ), width=c(4,1), height=c(round(height)-5,5))
par(mar=c(2,4,2,1)+.1)

m = matrix(1, ncol=nres, nrow=20)
# image(m, xaxt="n", yaxt="n", col=col_not_in_lib, ylim=c(1+.7/nres,-.7/nres), xlim=c(-.03,1.08))
image(m, xaxt="n", yaxt="n", col=col_not_in_lib, ylim=c(1+.7/nres,-.7/nres), xlim=c(-.025,1.025))
axis(1, seq(0, 20, length=20)/20, names(aa2int), cex.axis=.7, las=1, gap.axis=0)
axis(3, seq(0, 20, length=20)/20, names(aa2int), cex.axis=.7, las=1, gap.axis=0)
first_resn = 1
res_lab = paste(wt_seq, seq(first_resn,first_resn-1+length(wt_seq)), sep="")
mask = rep_len(c(TRUE,FALSE), length.out=nres)
axis(2,         (seq(0, nres, length=nres)/nres)[mask], labels=F, tcl=-.7, lwd=.6)
axis(2,         (seq(0, nres, length=nres)/nres)[!mask], labels=F, tcl=-2.2, lwd=.6)
axis(2,         (seq(0, nres, length=nres)/nres)[mask], res_lab[mask], cex.axis=.5, las=2, gap.axis=0, tick=F)
axis(2, line=1.4, (seq(0, nres, length=nres)/nres)[!mask], res_lab[!mask], cex.axis=.5, las=2, gap.axis=0, tick=F)

# Mark native
m[] = NA
m[cbind(aa2int[wt_seq], seq(nres))] = 1
image(m, col=col_native, add=T)

# # Gradient color of scores
m[] = NA
m[cbind(aa2int[pmap$mut], pmap$resi)] = pmap[,plot_cn]
image(m, zlim=c(grad_range_min,grad_range_max), col=col_grad, breaks=col_breaks, add=T)

# quartz(width=1.6, height=5)
par(mar=c(3,1,1,3)+.1)
image(t(col_breaks), zlim=c(grad_range_min,grad_range_max), col=col_grad, breaks=col_breaks, xaxt="n", yaxt="n")
n = grad_range_max - (grad_range_min-1)
axis(4,seq(0,n, length=n)/n, seq(grad_range_min,grad_range_max)/scale, las=2)
mtext(plot_lab, 1, 1, cex=.8)
dev.off()

}
