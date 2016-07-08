#!/bin/csh

#----------DESCRIPTION----------#
#make seismic analysis of SION sequence with ETH-SED catalogue
#1st: choose map type[country/topography as background]
#2nd: choose the prefered grd file
#3rd: choose an optimal cpt file
#4th: choose catalogue to plot seismicity
#5th: choose between projection methods
#6th: choose study area
#make script executable with "chmod +x seismicity_schweiz.csh"
#ETOPO grd file: https://www.ngdc.noaa.gov/mgg/global/global.html
#various cpt files could be downloaded at http://soliton.vm.bytemark.co.uk/pub/cpt-city/views/totp-cpt.html
#html color codes: http://html-color-codes.info
#rgb color codes: http://www.rapidtables.com/web/color/RGB_Color.htm
#----------DESCRIPTION----------#

#----------CONSOLE----------#
####MAP_TYPE####
set maptype = 0 #different countries has their own background color [0 = topo / 1 = country]
set colortype = 1 #colors to be implied [0 = time series for last 5 years / 1 = depth change]

####GRD&CPT####
set grd = "/Users/timothy.lee/polybox/Shared/SAMSFAULTZ/lib/srtm_ECOS.grd"
set cpt = "gray"
set seiscpt = "rainbow"
set focalcpt = "meca"

####CALTALOGUE####
set catalogue = "/Users/timothy.lee/polybox/Shared/SAMSFAULTZ/lib/MANULOC_List_1984_2012_ECOS_short.log.latest"

####PROJECTION####
set proj = M9i

####REGION####
set east = 7.10
set west = 7.70
set south = 46.15
set north = 46.40

####HISTOGRAM####
set binsize = 1

####COMMENT####
set comment = ""
#----------CONSOLE----------#


#----------PRE_SET----------#
if ($maptype == 1) then 
	set map = "country"
else
	set map = "topo"
endif
if ($colortype == 1) then 
	set col = "depth"
else
	set col = "time"
endif
if(! -e $catalogue) then
   echo "file not found: $catalogue"
   exit
endif
set region = "$east/$west/$south/$north"
set ps = seismicity_sion_"$map"_"$cpt"_"$col"_"$seiscpt"_"$east"_"$west"_"$south"_"$north"_"$comment".ps
set pshist = seismicity_sion_hist_"$col"_"$east"_"$west"_"$south"_"$north"_"$comment".ps
####SELECT_DATA####
gmt gmtselect $catalogue -R$region > temp_eq_list
####GMT_SET####
gmt gmtset PS_MEDIA = letter
gmt gmtset PS_PAGE_ORIENTATION = landscape
gmt gmtset PS_CHAR_ENCODING ISOLatin1+
gmt gmtset FONT_ANNOT_PRIMARY 12p,AvantGarde-Book,gray30
gmt gmtset FORMAT_GEO_MAP = D
#gmtset MAP_GRID_CROSS_SIZE_PRIMARY
#gmtset MEASURE_UNIT = inch
#gmtset HEADER_FONT_SIZE = 20
#gmtset HEADER_FONT = Courier-Bold
#gmtset LABEL_FONT_SIZE = 16
#gmtset LABEL_FONT = Helvetica-Narrow-Bold
#gmtset FRAME_WIDTH = 0.05 
gmt gmtset MAP_FRAME_TYPE = plain
#----------PRE_SET----------#

#----------MAP----------#
####BASEMAP####
gmt psbasemap -R$region -J$proj -Bxa0.1f0.05g0.1+"longitude" -Bya0.1f0.05g0.1+l"latitude" -BWNes+t"SION Seismicity" -K > $ps

####GRD####
if ($maptype == 0) then
	#cpt
	#gmt grdcut $grd -Gtemp.grd -R$region
	#gmt grd2cpt temp.grd -A50 -Cgray -R$region -S-5000/5000/100 -Z > temp_grd_gray.cpt
	gmt makecpt -A+80 -Cdem4 -T-100/+3000/100 -Z > temp_dem4.cpt
	gmt makecpt -Cgray -T-5000/+5000/100 -Z > temp_gray.cpt
	gmt makecpt -Cocean -T-1000/0/250 -Z > temp_ocean.cpt
	#grd
	gmt grdgradient $grd -Gtemp_grad.grd -A315
	#gmt grdgradient temp.grd -Gtemp_gradient.grd -A0/100 -Ne1.0
	gmt grdmath temp_grad.grd 80000 DIV = temp_grad_div.grd
	if (-e "temp_$cpt.cpt") then
		gmt grdimage $grd -R -J -C"temp_$cpt.cpt" -Itemp_grad_div.grd -nb -O -K >> $ps
		#gmt grdimage temp.grd -R -J -C$cpt -E300 -Itemp_gradient.grd -Q -O -K >> $ps
	else
		gmt grdimage $grd -R -J -C"$cpt.cpt" -Itemp_grad_div.grd -nb -O -K >> $ps
	endif
	gmt grdcontour $grd -R -J -A1000+pgray30 -C1000 -L-5000/5000 -Wa1p,gray30 -V -O -K >> $ps
	gmt psscale -C$cpt -Dx-0.7c/-0.7c+jBL+w6c/0.4c+h -Bxaf+l"topography" -By+lkm -I -O -K -V  >> $ps
endif

####PSCOAST####
if ($maptype == 1) then
	set color1='#CD5C5C@50' #switzerland
	set colorgroup1='CH'
	set color2='coral@50' #german
	set colorgroup2='DE'
	set color3='240/230/140@50' #france
	set colorgroup3='FR'
	set color4='0/36/74/4@50' #itlay
	set colorgroup4='IT'
	set color5='#8DC740@50' #austria
	set colorgroup5='AT'
	set color6='250/138/255@50' #liechtenstein
	set colorgroup6='LI'
	set color0='169@50' #other countries
	gmt pscoast -R -J -G${color0} -S -Dh -B -V -O -K >> $ps
	gmt pscoast -R -J -O -K -E${colorgroup1}+g${color1} >> $ps
	gmt pscoast -R -J -O -K -E${colorgroup2}+g${color2} >> $ps
	gmt pscoast -R -J -O -K -E${colorgroup3}+g${color3} >> $ps
	gmt pscoast -R -J -O -K -E${colorgroup4}+g${color4} >> $ps
	gmt pscoast -R -J -O -K -E${colorgroup5}+g${color5} >> $ps
	gmt pscoast -R -J -O -K -E${colorgroup6}+g${color6} >> $ps
	gmt pscoast -R -J -O -K -Dh -Ir/1p,cornflowerblue -N1/1p,,- -W1p,black >> $ps
endif

####SEISMICITY####
#[time = all / magnitude = all / quality = all / focal mech. = all]
#awk '{print $1,$2}' temp_eq_list | psxy -R -J -Sp -W1p -O -K >> $ps 
awk -v ms=0.0200 '{rms=ms*1.0;if(substr($1,1,1)!="#" && substr($12,1,3)=="SED" && $4>1.0) {print $1,$2,$4*$4*rms}}' temp_eq_list | psxy -R -J$proj -Sc -W1p,black -O -K >> $ps

##[time = 2014 / magnitude = all / quality = all / focal mech. = all]
#awk -v ms=0.0200 '{rms=ms*1.0;if(substr($1,1,1)!="#" && substr($12,1,3)=="SED" && $4>1.0 && $6 >= 2010 && $6 <= 2010) {print $1,$2,$4*$4*rms}}' temp_eq_list | psxy -R -J$proj -Sc -W1p,purple -O -K >> $ps
#
##[time = 2014 / magnitude = all / quality = all / focal mech. = all]
#awk -v ms=0.0200 '{rms=ms*1.0;if(substr($1,1,1)!="#" && substr($12,1,3)=="SED" && $4>1.0 && $6 >= 2011 && $6 <= 2011) {print $1,$2,$4*$4*rms}}' temp_eq_list | psxy -R -J$proj -Sc -W1p,43/0/255 -O -K >> $ps
#
##[time = 2014 / magnitude = all / quality = all / focal mech. = all]
#awk -v ms=0.0200 '{rms=ms*1.0;if(substr($1,1,1)!="#" && substr($12,1,3)=="SED" && $4>1.0 && $6 >= 2012 && $6 <= 2012) {print $1,$2,$4*$4*rms}}' temp_eq_list | psxy -R -J$proj -Sc -W1p,blue -O -K >> $ps
#
##[time = 2014 / magnitude = all / quality = all / focal mech. = all]
#awk -v ms=0.0200 '{rms=ms*1.0;if(substr($1,1,1)!="#" && substr($12,1,3)=="SED" && $4>1.0 && $6 >= 2013 && $6 <= 2013) {print $1,$2,$4*$4*rms}}' temp_eq_list | psxy -R -J$proj -Sc -W1p,green -O -K >> $ps
#
##[time = 2014 / magnitude = all / quality = all / focal mech. = all]
#awk -v ms=0.0200 '{rms=ms*1.0;if(substr($1,1,1)!="#" && substr($12,1,3)=="SED" && $4>1.0 && $6 >= 2014 && $6 <= 2014) {print $1,$2,$4*$4*rms}}' temp_eq_list | psxy -R -J$proj -Sc -W1p,yellow -O -K >> $ps
#
##[time = 2015 / magnitude = all / quality = all / focal mech. = all]
#awk -v ms=0.0200 '{rms=ms*1.0;if(substr($1,1,1)!="#" && substr($12,1,3)=="SED" && $4>1.0 && $6 >= 2015 && $6 <= 2015) {print $1,$2,$4*$4*rms}}' temp_eq_list | psxy -R -J$proj -Sc -W1p,orange -O -K >> $ps
#
##[time = 2016 / magnitude = all / quality = all / focal mech. = all]
#awk -v ms=0.0200 '{rms=ms*1.0;if(substr($1,1,1)!="#" && substr($12,1,3)=="SED" && $4>1.0 && $6 >= 2016 && $6 <= 2016) {print $1,$2,$4*$4*rms}}' temp_eq_list | psxy -R -J$proj -Sc -W1p,red -O -K >> $ps

if ($colortype == 0) then
	#[time = all / magnitude = all / quality = all / focal mech. = all] with automatic color scale coresponding with time 
	set seismin = `gmtinfo -C -i5 /Users/timothy.lee/polybox/Shared/SAMSFAULTZ/lib/MANULOC_List_1984_2012_ECOS_short.log.latest | awk '{print $1}'`
	set seismax = `gmtinfo -C -i5 /Users/timothy.lee/polybox/Shared/SAMSFAULTZ/lib/MANULOC_List_1984_2012_ECOS_short.log.latest | awk '{print $2}'`
	set seismin = `gmt gmtmath -Q $seismax 5 SUB =`
	set seismax = `gmt gmtmath -Q $seismax 1 ADD =`
	set seisinterv = `gmt gmtmath -Q 1 12 DIV =`
	gmt makecpt -C$seiscpt -T$seismin/$seismax/$seisinterv > temp_seis_$seiscpt.cpt
	awk -v ms=0.0200 '{rms=ms*1.0;if(substr($1,1,1)!="#" && substr($12,1,3)=="SED" && $4>1.0 && $6 >= '$seismin' && $6 <= '$seismax') {print $1,$2,$6+$7/12,$4*$4*rms}}' temp_eq_list | gmt psxy -R -J$proj -Ctemp_seis_$seiscpt.cpt -Sc -W1p -O -K >> $ps
	gmt psscale -Ctemp_seis_$seiscpt.cpt -Dx6.7c/-0.7c+jBL+w6c/0.4c+h -Bxaf+l"time" -By+lyear -I -O -K -V  >> $ps
endif

if ($colortype == 1) then
	#[time = all / magnitude = all / quality = all / focal mech. = all] with automatic color scale coresponding with depth
	set seismin = `gmtinfo -C -i2 /Users/timothy.lee/polybox/Shared/SAMSFAULTZ/lib/MANULOC_List_1984_2012_ECOS_short.log.latest | awk '{print $1}'`
	set seismax = `gmtinfo -C -i2 /Users/timothy.lee/polybox/Shared/SAMSFAULTZ/lib/MANULOC_List_1984_2012_ECOS_short.log.latest | awk '{print $2}'`
	gmt makecpt -C$seiscpt -I -T-2/20/20+ -Z > temp_seis_$seiscpt.cpt
	awk -v ms=0.0200 '{rms=ms*1.0;if(substr($1,1,1)!="#" && substr($12,1,3)=="SED" && $4>1.0) {print $1,$2,$3,$4*$4*rms}}' temp_eq_list | gmt psxy -R -J$proj -Ctemp_seis_$seiscpt.cpt -Sc -W1p -O -K >> $ps
	gmt psscale -Ctemp_seis_$seiscpt.cpt -Dx6.7c/-0.7c+jBL+w6c/0.4c+h -Bxaf+l"depth" -By+lkm -I -O -K -V  >> $ps
endif

#[time = all / magnitude = all / quality = all / focal mech. = all] with automatic color scale coresponding with focal mechanism


####FOCAL_MECHANISMS####
#gmt psmeca -J -R -CP5p -Sa1.3c -Z"$CPT.cpt" -K -O >> $PS

#add stations

#add faults
#----------MAP----------#

#--------PROFILE----------#
#project -C120.613/24.017 -E120.711/24.035 -G.5 -Dd -Q > trackA
#grdtrack trackA -GTotal80.grd | gawk "{print $3, $4 }" > trackA.dat
#--------PROFILE----------#

#--------HISTOGRAM----------#
if ($colortype == 0) then
	set histxmin = `awk '{print $6}' temp_eq_list | gmt pshistogram -I -W$binsize | awk '{print $1}'`
	set histxmin = `gmt gmtmath -Q $histxmin 1 SUB =`
	set histxmax = `awk '{print $6}' temp_eq_list | gmt pshistogram -I -W$binsize | awk '{print $2}'`
	set histxmax = `gmt gmtmath -Q $histxmax 1 ADD =`
	set histymin = `awk '{print $6}' temp_eq_list | gmt pshistogram -I -W$binsize | awk '{print $3}'`
	set histymax = `awk '{print $6}' temp_eq_list | gmt pshistogram -I -W$binsize | awk '{print $4}'`
	set histymax = `gmt gmtmath -Q $histymax 500 ADD =`
	awk '{print $6}' temp_eq_list | gmt pshistogram -Bxa5f1+l"Year" -Bya500f100+l"Frequency (counts)" -BWSne+t"Histogram"+glightblue -R$histxmin/$histxmax/$histymin/$histymax -JX6i/3i -D+r -F -Gorange -L1p -N0 -W$binsize -P -K > $pshist
	set chistxmin = `awk '{print $6}' temp_eq_list | gmt pshistogram -I -Q -W$binsize | awk '{print $1}'`
	set chistxmin = `gmt gmtmath -Q $chistxmin 1 SUB =`
	set chistxmax = `awk '{print $6}' temp_eq_list | gmt pshistogram -I -Q -W$binsize | awk '{print $2}'`
	set chistxmax = `gmt gmtmath -Q $chistxmax 1 ADD =`
	set chistymin = `awk '{print $6}' temp_eq_list | gmt pshistogram -I -Q -W$binsize | awk '{print $3}'`
	set chistymax = `awk '{print $6}' temp_eq_list | gmt pshistogram -I -Q -W$binsize | awk '{print $4}'`
	set chistymax = `gmt gmtmath -Q $histymax LOG 1 ADD =`
	awk '{print $6}' temp_eq_list | gmt pshistogram -Bxa5f1+l"Year" -Bya5f1+l"Frequency (counts)" -BWSne+t"Culmulative Histogram"+glightblue -R$chistxmin/$chistxmax/$chistymin/$chistymax -JX6i/3i -D+r -F -Gorange -L1p -N0 -W$binsize -Q -Y5i -Z4 -O -K >> $pshist
	#gmt psrose fractures.d -: -A10r -S1.8in -P -Gorange -R0/1/0/360 -X2.5i -Bx0.2g0.2 -By30g30 -B+glightblue -W1p -O >> $pshist
endif
if ($colortype == 1) then
	set histxmin = `awk '{print $3}' temp_eq_list | gmt pshistogram -I -W$binsize | awk '{print $1}'`
	set histxmin = `gmt gmtmath -Q $histxmin 1 SUB =`
	set histxmax = `awk '{print $3}' temp_eq_list | gmt pshistogram -I -W$binsize | awk '{print $2}'`
	set histxmax = `gmt gmtmath -Q $histxmax 1 ADD =`
	set histymin = `awk '{print $3}' temp_eq_list | gmt pshistogram -I -W$binsize | awk '{print $3}'`
	set histymax = `awk '{print $3}' temp_eq_list | gmt pshistogram -I -W$binsize | awk '{print $4}'`
	set histymax = `gmt gmtmath -Q $histymax 500 ADD =`
	awk '{print $3}' temp_eq_list | gmt pshistogram -Bxa5f1+l"Depth" -Bya500f100+l"Frequency (counts)" -BWSne+t"Histogram"+glightblue -R$histxmin/$histxmax/$histymin/$histymax -JX6i/3i -D+r -F -Gorange -L1p -N0 -W$binsize -P -K > $pshist
	set chistxmin = `awk '{print $3}' temp_eq_list | gmt pshistogram -I -Q -W$binsize | awk '{print $1}'`
	set chistxmin = `gmt gmtmath -Q $chistxmin 1 SUB =`
	set chistxmax = `awk '{print $3}' temp_eq_list | gmt pshistogram -I -Q -W$binsize | awk '{print $2}'`
	set chistxmax = `gmt gmtmath -Q $chistxmax 1 ADD =`
	set chistymin = `awk '{print $3}' temp_eq_list | gmt pshistogram -I -Q -W$binsize | awk '{print $3}'`
	set chistymax = `awk '{print $3}' temp_eq_list | gmt pshistogram -I -Q -W$binsize | awk '{print $4}'`
	set chistymax = `gmt gmtmath -Q $histymax LOG 1 ADD =`
	awk '{print $3}' temp_eq_list | gmt pshistogram -Bxa5f1+l"Depth" -Bya5f1+l"Frequency (counts)" -BWSne+t"Culmulative Histogram"+glightblue -R$chistxmin/$chistxmax/$chistymin/$chistymax -JX6i/3i -D+r -F -Gorange -L1p -N0 -W$binsize -Q -Y5i -Z4 -O -K >> $pshist
	#gmt psrose fractures.d -: -A10r -S1.8in -P -Gorange -R0/1/0/360 -X2.5i -Bx0.2g0.2 -By30g30 -B+glightblue -W1p -O >> $pshist
endif
#--------HISTOGRAM----------#

#----------INSERT_MAP----------#
gmt pscoast -R5/13/45.5/49 -JM2i -B0 -Bwnes+gwhite -Df -N1 -W -A5000 --MAP_FRAME_TYPE=plain -X0.2i -Y4i -O -K >> $ps
set color1='#CD5C5C@50' #switzerland
set colorgroup1='CH'
set color2='coral@50' #german
set colorgroup2='DE'
set color3='240/230/140@50' #france
set colorgroup3='FR'
set color4='0/36/74/4@50' #itlay
set colorgroup4='IT'
set color5='#8DC740@50' #austria
set colorgroup5='AT'
set color6='250/138/255@50' #liechtenstein
set colorgroup6='LI'
set color0='169@50' #other countries
gmt pscoast -R -J -G${color0} -S -Dh -B -V -O -K >> $ps
gmt pscoast -R -J -O -K -E${colorgroup1}+g${color1} >> $ps
gmt pscoast -R -J -O -K -E${colorgroup2}+g${color2} >> $ps
gmt pscoast -R -J -O -K -E${colorgroup3}+g${color3} >> $ps
gmt pscoast -R -J -O -K -E${colorgroup4}+g${color4} >> $ps
gmt pscoast -R -J -O -K -E${colorgroup5}+g${color5} >> $ps
gmt pscoast -R -J -O -K -E${colorgroup6}+g${color6} >> $ps
gmt pscoast -R -J -O -K -Dh -Ir/1p,cornflowerblue -N1/1p,,- -W1p,black >> $ps
gmt psbasemap -R5/13/45.5/49 -JM2i -D$region -F+p2p,red -O >> $ps
#----------INSERT_MAP----------#

####SAVE_BOTH_OUTPUT_&_ERROR_MESSAGE####
#gmt module > output.d 2> errors.log

####CONVERT_FROM_POSTSCRIPT_TO_PDF####
#gmt psconvert -Tf $ps

echo "Figure saved as: $ps"
echo "Histogram saved as: $pshist"

rm temp*