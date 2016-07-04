#!/bin/csh

#----------DESCRIPTION----------#
#make seismic analysis of SION sequence with ETH-SED catalogue
#1st: choose map type[country/topography based]
#2nd: choose the prefered grd file
#3rd: choose an optimal cpt file
#4th: choose catalogue to plot seismicity
#5th: choose between projection methods
#6th: choose study area
#last: make script executable with "chmod +x seismicity_schweiz.csh"

#ETOPO grd file: https://www.ngdc.noaa.gov/mgg/global/global.html
#various cpt files could be downloaded at http://soliton.vm.bytemark.co.uk/pub/cpt-city/views/totp-cpt.html
#html color codes: http://html-color-codes.info
#rgb color codes: http://www.rapidtables.com/web/color/RGB_Color.htm
#----------DESCRIPTION----------#

#----------CONSOLE----------#
####MAP_TYPE####
set maptype = 0 #different countries has their own background color [0 = topo / 1 = country]
if ($maptype == 1) then 
	set map = "country"
else
	set map = "topo"
endif

####GRD&CPT####
set grd = "/Users/timothy.lee/polybox/Shared/SAMSFAULTZ/lib/euro30m.grd"
#set grd = srtm_ECOS.grd
#set grd = srtm_ECOS.grd
#set grd = Wagner_Moho__LET_CSS_Euro
#set grd = Wagner_Moho_LET_CSS_Adria
#set grd = Wagner_Moho_LET_CSS_Liguria
#set grd = ETOPO1_Bed_g_gmt4.grd
set cpt = "temp_dem4"

####CALTALOGUE####
set catalogue = "/Users/timothy.lee/polybox/Shared/SAMSFAULTZ/lib/MANULOC_List_1984_2012_ECOS_short.log.latest"

####PROJECTION####
set proj = M9i

####set_region####
set east = 7.10
set west = 7.70
set south = 46.15
set north = 46.40
set region = "$east/$west/$south/$north"
#set region = 7.10/7.70/46.15/46.40
#----------CONSOLE----------#

#----------PRE_SET----------#
####FILE_NAME####
set comment = "2015-2016"
set ps = seismicity_sion_"$map"_"$cpt"_"$east"_"$west"_"$south"_"$north"_"$comment".ps
####GMT_SET####
gmtset PS_MEDIA = letter
gmtset PS_PAGE_ORIENTATION = landscape
gmtset FONT_ANNOT_PRIMARY 12p,Times-Bold,red
#gmtset FORMAT_GEO_MAP
#gmtset MAP_GRID_CROSS_SIZE_PRIMARY
#gmtset PAPER_MEDIA A4
#gmtset MEASURE_UNIT = inch
#gmtset HEADER_FONT_SIZE = 20
#gmtset HEADER_FONT = Courier-Bold
#gmtset LABEL_FONT_SIZE = 16
#gmtset LABEL_FONT = Helvetica-Narrow-Bold
#gmtset ANNOT_FONT_SIZE = 12
#gmtset FRAME_WIDTH = 0.05 
gmt gmtset MAP_FRAME_TYPE = plain
#gmt gmtset OUTPUT_DEGREE_FORMAT +D
#gmt gmtset PLOT_DEGREE_FORMAT +DF
#----------PRE_SET----------#

#----------MAP----------#
####BASEMAP####
gmt psbasemap -R$region -J$proj -Bxa0.1f0.05g0.1+"longitude" -Bya0.1f0.05g0.1+l"latitude" -BWNes+t"SION Seismicity" -K > $ps

####GRD&CPT####
if ($maptype == 0) then
	gmt grdcut $grd -Gtemp.grd -R$region
	#CPT
	gmt grd2cpt temp.grd -A90 -Cdem4 -R$region -S-3000/3000/100 > temp_dem4.cpt
	gmt makecpt -Cgray -T-5000/+5000/100 > temp_gray.cpt
	gmt makecpt -Cocean -T-1000/0/250 -Z > temp_ocean.cpt
	gmt makecpt -Cmby.cpt -T-8000/5100/1000 -Z > temp_mby.cpt
	#GRD
	gmt grdimage temp.grd -R -J -C$cpt -O -K >> $ps
	gmt grdcontour temp.grd -R -J -C1000 -L-6000/3000 -Wa5/25 -V -O -K >> $ps
	gmt psscale -C$cpt -D5/-1.5/7/0.2h -P -O -Ba2000:meter: -I -O -K -V  >> $ps
	#gmt grdgradient
	#gmt grdmath
endif

####pscoast####
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

#add serismicity
#[time = all / magnitude = all / quality = all / focal mech. = all]
#awk '{print $1,$2}' $catalogue | psxy -R -J -Sp -W1p -O -K >> $ps 
awk -v ms=0.0150 '{rms=ms*1.0;if(substr($1,1,1)!="#" && substr($12,1,3)=="SED" && $4>1.0) {print $1,$2,$4*$4*rms}}' $catalogue | psxy -R -J$proj -Sc -W1p,128/128/128 -O -K >> $ps

#[time = 2015 / magnitude = all / quality = all / focal mech. = all]
awk -v ms=0.0150 '{rms=ms*1.0;if(substr($1,1,1)!="#" && substr($12,1,3)=="SED" && $4>1.0 && $6 >= 2015 && $6 <= 2015) {print $1,$2,$4*$4*rms}}' $catalogue | psxy -R -J$proj -Sc -W1p,255/128/0 -O -K >> $ps

#[time = 2016 / magnitude = all / quality = all / focal mech. = all]
awk -v ms=0.0150 '{rms=ms*1.0;if(substr($1,1,1)!="#" && substr($12,1,3)=="SED" && $4>1.0 && $6 >= 2016 && $6 <= 2016) {print $1,$2,$4*$4*rms}}' $catalogue | psxy -R -J$proj -Sc -W1p,255/0/0 -O -K >> $ps


#add focal mechanisms

#add stations

#add faults
#----------MAP----------#

#--------PROFILE----------#
#project -C120.613/24.017 -E120.711/24.035 -G.5 -Dd -Q > trackA
#grdtrack trackA -GTotal80.grd | gawk "{print $3, $4 }" > trackA.dat
#--------PROFILE----------#

#--------HISTOGRAM----------#

#--------HISTOGRAM----------#


####SAVE_BOTH_OUTPUT_&_ERROR_MESSAGE####
#gmt module > output.d 2> errors.log

####CONVERT_FROM_POSTSCRIPT_TO_PDF####
#gmt psconvert -Tf $ps

echo "Figure saved as: $ps"

rm temp*