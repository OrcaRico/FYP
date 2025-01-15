#!/bin/bash

## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
## EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
## OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
## NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
## HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
## WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
## FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
## OTHER DEALINGS IN THE SOFTWARE.

[ $# -eq 0 ] && {
    echo "usage: morphology.bash <whatever/morph> [nboot]"
    exit
}

# limit memory usage to ~600 MB
#ulimit -v 614400
# ~1GB
ulimit -v 1048576
# ~2GB
#ulimit -v 2097152

knhcorr=~amantz/bin/K_nH_correction.py

im=img2_0.6-2.0_mask.fits.gz
bg=bgimg2_0.6-2.0_mask.fits.gz
ex=expmap2_1.30kev_mask.fits.gz
smooth=morph_smooth.fits.gz
isoph=morph_isophotes.fits.gz
filth=morph_filth.fits.gz
mask=morph_mask_imagecoords.reg
centerfile=morph_fix_center.dat
Rfile=makeplots.R
pfrac=0.05
Nisoph=5

pixel=0.984 # arcsec
Erange=0.6-2.0

dir=$(pwd)

work=${1%/}
nboot=$2

id=$work

cd $work || exit
[ -r $ex ] || ex=expmap2_1.30kev.fits.gz; ls $im $bg $ex $scal
[ -r $im ] || exit
[ -r $bg ] || exit
[ -r $ex ] || exit
[ -r $scal ] || exit

read z < redshift.txt
[ $z ] || {
    echo "Warning: couldn't find redshift of cluster."
    exit
}

read nH < nH.txt
[ $nH ] || {
    echo "Warning: couldn't find nH of cluster."
    exit
}

read kpc < kpc.txt
[ $kpc ] || {
    echo "Warning: couldn't find pixel size in kpc."
#    exit
}

read kT < kT.txt
[ $kT ] || {
    echo "Warning: couldn't find kT of cluster."
    exit
}

read betarad < beta_rmin.txt
[ $betarad ] && {
    betarad=" --beta-rmin $betarad"
}

read -a bgscal < bg_scaling.txt
[ $bgscal ] || {
    echo "Warning: couldn't find blank-sky scaling factor."
    exit
}

# h=0.7, OmegaM=0.3, OmegaL=0.7, z<2
E=$(awk -v z=$z 'BEGIN{print 1.0+0.45619*z+0.344297*z^2-0.0406396*z^3}')

K=$($knhcorr --kT $kT --nH $nH --Eobs $Erange --Ez 0.1-50.0 $z)

# SB scaling factor
norm=$(awk -v K=$K -v kT=$kT -v z=$z -v E=$E -v p=$pixel 'BEGIN{print K*kT*E^3/(1.0+z)^4*(p/0.984)^2}')
# isophote SB levels
isoph_max=$(awk -v x=$norm 'BEGIN{print 0.05*x}')
isoph_flux=$(awk -v x=$norm 'BEGIN{print 2.0e-3*x}')
# SB level for window
flux=$(awk -v x=$norm 'BEGIN{print 0.0475*x}') # just needs to be lower than isoph_max

[ $kpc ] && {
    rhoc=$(awk -v z=$z 'BEGIN{print 1.36e+11 + 1.224e+11 * z + 1.224e+11 * z^2 + 4.08e+10 * z^3}') # Msun/Mpc^3
    r500=$(awk -v kpc=$kpc -v E=$E -v rhoc=$rhoc -v kT=$kT 'function log10(x){return log(x)/log(10.0)}BEGIN{print ( 10.^((log10(kT)-0.89)/0.49 - log10(E) + 15.) / (4.18879*500*rhoc) )^0.3333 * 1e3 / kpc}')
    kpc="--kpc $kpc --Ralt $r500 "
}

if [ -e $mask ]; then mask="--mask $mask "; else mask=; fi

if [ $nboot ]; then
    verb="--quiet"
    out=morph_boot.dat
#    rm -f $out
    boot="--boot $nboot"
else
    verb="--obnoxious"
    boot="--write-smooth $smooth --write-isophotes $isoph --write-filth $filth"
    out=morph.log
    rm -f $out
    echo "SB scaling factor = $norm" | tee -a $out
fi

[ -e $centerfile ] && fixcenter="--fix-center $(cat $centerfile) "

echo $dir/morphology $im $verb $mask$fixcenter--bg-file $bg ${bgscal[@]} --expmap $ex --peaky-flux $flux --peaky-frac $pfrac $kpc--isoph-min-level $isoph_flux --isoph-max-level $isoph_max $boot --num-isoph $Nisoph$betarad

time $dir/morphology $im $verb $mask$fixcenter--bg-file $bg ${bgscal[@]} --expmap $ex --peaky-flux $flux  $kpc--isoph-min-level $isoph_flux --isoph-max-level $isoph_max $boot --num-isoph $Nisoph$betarad | tee -a $out

[ $nboot ] && exit

centeringreg=morph_centering.reg
sbdat=morph_sb.dat
ellipsereg=morph_isoph_ellipse.reg
shiftreg=morph_isoph_shifts.reg
isolev=morph_isoph_levels.dat
isow=morph_isoph_expfrac.dat
centroids=morph_centroids.reg

awk '
/Region file format/ {f=1;for(i=0;i<4;++i){print;getline}}
f==1 && (/box/||/vector/||/^point/) {print}
f==1 && /Calculated/ {exit}
' $out > $centeringreg

awk '
/SB profile has [0-9]+ annuli/ {f=1;n=$4;for(i=0;i<6;++i)getline;next}
f!="" {print;if(f==n+1)exit;++f}
' $out > $sbdat

awk '
/Region file format/ {++f}
f==2 && /Region file format/ {for(i=0;i<4;++i){print;getline}}
f==2 && /^ellipse/ {print}
' $out > $ellipsereg

awk '
/Region file format/ {++f}
f==2 && /Region file format/ {for(i=0;i<4;++i){print;getline}next}
f==2 && /vector/ {print}
' $out > $shiftreg

awk '
/Isophote flux ranges are/ {f=1}
f==1 {
 while (1==1) {
  getline
  if ($0!~/^[0-9]/) exit
  print
 }
}
' $out > $isolev

awk '/Exposed fraction of isophote/ {
 printf $7" "$12
 getline
 print " "$9" "$11
}' $out > $isow

awk '
/Region file format/ {++f}
f==3 && /Region file format/ {for(i=0;i<4;++i){print;getline}next}
f==3 && /^point/ {print}
' $out > $centroids

echo "Making plots..."
cp $dir/$Rfile . || exit
./$Rfile $id

echo "All done."
