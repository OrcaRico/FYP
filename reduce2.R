#!/bin/env Rscript

## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
## EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
## OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
## NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
## HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
## WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
## FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
## OTHER DEALINGS IN THE SOFTWARE.

peakcut = -0.82
aligcut = 1.0
symmcut = 0.87




argv = commandArgs(trailingOnly=T)
argc = length(argv)
if (argc < 2) stop("usage: ./reduce.Rscript <morph.log> <redshift> [boot.dat] [ellipses to ignore]")
logfile = argv[1]
z = as.numeric(argv[2])
boot = NA
if (argc >= 3) boot = argv[3]
igel = integer()
if (argc >= 4) igel = eval(argv[4])





back = function(x) x[length(x)]

# ellipse definitions
nonel = 2 + 1 + 1 + 1 + 4 + 4 + 1 # center x,y; csb; filth frac; centroid variance; 2x4 power ratios; avg  central sb
perel = 8 # x, y, a, b, theta, exposed frac, avcos, avsin
getx = function(i,j,f) f[i,nonel+(j-1)*perel+1]
gety = function(i,j,f) f[i,nonel+(j-1)*perel+2]
geta = function(i,j,f) f[i,nonel+(j-1)*perel+3]
getb = function(i,j,f) f[i,nonel+(j-1)*perel+4]
gettheta = function(i,j,f) f[i,nonel+(j-1)*perel+5]
getfrac = function(i,j,f) f[i,nonel+(j-1)*perel+6]
getcos = function(i,j,f) f[i,nonel+(j-1)*perel+7]
getsin = function(i,j,f) f[i,nonel+(j-1)*perel+8]
elfrac.min = 0.5
eltrig.max = 0.4
goodfrac.min = 0.75

reduce = function(f, igel) {
 colnames(f)[1:nonel] = c('center.x', 'center.y', 'csb', 'filth.frac', 'cv', paste('PR', 1:4, sep=''), 'censbav', paste('cenPR', 1:4, sep=''))
 g = f[,c(3,5:nonel)] # all non-ellipse measurements except global center, filth frac
 nel = (ncol(f) - nonel)/perel

 # cull bad ellipse fits
 for (i in 1:nrow(f)) for (j in 1:nel) {
  if (is.na(getfrac(i,j,f)) | is.na(getcos(i,j,f)) | is.na(getsin(i,j,f))) next
  if (getfrac(i,j,f)<elfrac.min | abs(getcos(i,j,f))>eltrig.max | abs(getsin(i,j,f))>eltrig.max) f[i,nonel+(j-1)*perel+1:perel] = NA
 }

 # check how many ellipses survived
 allrows = 1:nrow(f)
 ngood = rep(0, nel)
 for (j in 1:nel) {
  ngood[j] = length(which(!is.na(geta(allrows,j,f))))
  if (j %in% igel) ngood[j] = -1
}
 useel = which(ngood/nrow(f) >= goodfrac.min)

 if (length(useel) > 1) {
  
  # ellipticity
  g$ellipticity.mean = NA
  g$ellipticity.slope = NA
  for (j in allrows) {
   thing = as.numeric(1-getb(j,useel,f)/geta(j,useel,f))
   g$ellipticity.mean[j] = mean(thing, na.rm=T)
   if (length(thing)-length(which(is.na(thing))) > 1) g$ellipticity.slope[j] = as.numeric(lm(thing~useel)$coefficients[2])
  }

  # position angle
  g$pa.slope = NA
  for (i in allrows) {
   pa = as.numeric(gettheta(i,useel,f))
   j = which(!is.na(pa))
   if (length(j) > 1) {
    pa = pa[j]
    j = useel[j]
    last = pa[1]
    for (jj in 2:length(j)) {
     while (pa[jj]-last>90) pa[jj] = pa[jj] - 180
     while (last-pa[jj]>90) pa[jj] = pa[jj] + 180
    }
    if (length(pa)-length(which(is.na(pa))) > 1) g$pa.slope[i] = as.numeric(lm(pa~j)$coefficients[2])
   }
  }

  # shifts
  for (i in allrows) {

   g$shift.center[i] = sqrt((getx(i,back(useel),f)-f$center.x[i])^2 + (gety(i,back(useel),f)-f$center.y[i])^2) / (0.5*( geta(i,back(useel),f) + getb(i,back(useel),f) ))

   el1 = useel[-length(useel)]
   el2 = useel[-1]

   thing = as.numeric(sqrt((getx(i,el1,f)-getx(i,el2,f))^2 + (gety(i,el1,f)-gety(i,el2,f))^2) / (0.25*( geta(i,el1,f) + getb(i,el1,f) + geta(i,el2,f) + getb(i,el2,f) )))
   g$shift.mean[i] = mean(thing, na.rm=T)

   thing = as.numeric(sqrt((getx(i,useel,f)-f$center.x[i])^2 + (gety(i,useel,f)-f$center.y[i])^2) / (0.5*( geta(i,useel,f) + getb(i,useel,f) )))
   g$asym4.mean[i] = mean(thing, na.rm=T)

  }

 } # >1 ellipse

 centering = sqrt(var(f$center.x, na.rm=T) + var(f$center.y, na.rm=T))

 el.centering = NA
 el.centers = matrix(NA, 1, 2)
 if (length(useel) > 0) {
  el.centers = matrix(NA, back(useel), 2)
  el.centering = rep(NA, back(useel))
  for (j in useel) {
   el.centers[j,1] = median(getx(allrows,j,f), na.rm=T)
   el.centers[j,2] = median(gety(allrows,j,f), na.rm=T)
   el.centering[j] = sqrt(var(getx(allrows,j,f), na.rm=T) + var(gety(allrows,j,f), na.rm=T))
  }
 }

 el.semimajor = NA
 el.semiminor = NA
 el.rad.mean = NA
 el.rad.sd = NA
 if (length(useel) > 0) {
   el.semimajor = rep(NA, back(useel))
   el.semiminor = rep(NA, back(useel))
   el.rad.mean = rep(NA, back(useel))
   el.rad.sd = rep(NA, back(useel))
  for (j in useel) {
   el.semimajor[j] = median(geta(allrows,j,f), na.rm=T)
   el.semiminor[j] = median(getb(allrows,j,f), na.rm=T)
   el.rad.mean[j] = mean(0.5*(geta(allrows,j,f) + getb(allrows,j,f)), na.rm=T)
   el.rad.sd[j] = sd(0.5*(geta(allrows,j,f) + getb(allrows,j,f)), na.rm=T)
  }
 }

 cov = tryCatch(cov(g, use='p'), error=function(e) matrix(NA, ncol(g), ncol(g)))

 list(table=g, filth.frac=mean(f$filth.frac, na.rm=T), ellipses=useel, N.ellipse=length(useel), global.centering=centering, ellipse.centering=el.centering, ellipse.centers=data.frame(x=el.centers[,1], y=el.centers[,2]), ellipse.semimajor = el.semimajor, ellipse.semiminor=el.semiminor, ellipse.rad.mean=el.rad.mean, ellipse.rad.sd.mn=el.rad.sd, mean=colMeans(g, na.rm=T), cov=cov)
}













con = pipe(paste("awk '/SB scaling factor/{print $NF}'",logfile), open='r')
sbscal = scan(con, quiet=T)
close(con)


if (is.na(boot)) {

 con = pipe(paste('tail -n1',logfile), open='r')
 line = scan(con, quiet=T)
 close(con)
 f = data.frame(matrix(line, 1, length(line)))

} else {

 f = read.table(boot, head=F, comment='/', colClass=double(), fill=T)

}

cluster = reduce(f, igel)
cluster$source = logfile
if (!is.na(boot)) cluster$boot = boot
cluster$sbscal = sbscal



important = c(3,1,2)
zpower = 1

if (cluster$N.ellipse < 2) important = 1

for (i in 1:ncol(cluster$table)) cluster$table[which(!is.finite(cluster$table[,i])),i] = NA
cluster$redshift = z
cluster$table2 = data.frame(peak=log10(pmax(cluster$table$censbav/cluster$sbscal*(1+cluster$redshift)^zpower, 1e-4)), alig=-log10(cluster$table$shift.mean), symm=-log10(cluster$table$asym4.mean))
cluster$table3 = cluster$table2[,important]
cluster$mean3 = colMeans(cluster$table3, na.rm=T)
cluster$cov3 = tryCatch(cov(cluster$table3, use='p'), error=function(e) matrix(NA, ncol(cluster$table3), ncol(cluster$table3)))

nn = length(which(!is.na(cluster$table3$peak) & !is.na(cluster$table3$alig) & !is.na(cluster$table3$symm)))
Prel = length(which(cluster$table3$peak >= peakcut & cluster$table3$alig >= aligcut & cluster$table3$symm >= symmcut)) / nn

if (nrow(cluster$table) < 2) {
 write(paste('#', Prel), '')
} else {
 write(Prel, '')
}
write(cluster$mean3, '')
write(cluster$cov3, '', ncol=length(important))





