#!/Users/timothy.lee/anaconda/bin/python

#exercise to map the seismicity in switzerland, data from sed-fdsn service

#---LIBRARIES---#
from obspy import read, read_inventory, read_events, UTCDateTime, Trace, Stream
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import URL_MAPPINGS
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import UnivariateSpline
import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import numpy as np
import numpy.matlib as mat
#import scipy as Sci
#import scipy.linalg


#---PARAMETERIZATION---#
plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = 12, 8
def get_marker_color(depth):
    if depth < 5.0:
        return ('ro')
    elif depth < 10.0:
        return ('mo')
    elif depth < 20.0:
        return ('yo')
    elif depth < 30.0:
        return ('go')
    elif depth < 40.0:
        return ('co')
    else:
        return ('bo')

#---DATA_FROM_FDSN---#
#get specified event
client = Client(base_url = "http://arclink.ethz.ch")
starttime = UTCDateTime("2016-01-01")
endtime = UTCDateTime()
cat = client.get_events(starttime=starttime, endtime=endtime, minmagnitude=1)#, filename="requested_events.xml"
evtnum = cat.count()
print(type(cat)) #Catalog
print(cat) #add , CatalogObject.__str__(print_all=True) to print all events
#cat.plot() #add also resouces like: projection="local"
focaltime = []
hypolon = []
hypolat = []
hypodep = []
eqmag = []
for x in range(0, evtnum):
	event = cat[x]
	focaltime.append(event.origins[0].time)
	hypolon.append(event.origins[0].longitude)
	hypolat.append(event.origins[0].latitude)
	hypodep.append(event.origins[0].depth/1000)
	eqmag.append(event.magnitudes[0].mag)

##---MAPPING_DATA---#
##global map
#global_map = Basemap(projection='robin', lat_0=7.5, lon_0=47.5, resolution='l', area_thresh=1000)#, llcrnrlon=-136.25, llcrnrlat=56, urcrnrlon=-134.25, urcrnrlat=57.75) #also try: projection='ortho'
#global_map.drawcoastlines()
#global_map.drawcountries()
#global_map.fillcontinents(color='KHAKI') #try also: global_map.bluemarble()
#global_map.drawmapboundary()
#global_map.drawmeridians(np.arange(0, 360, 30))
#global_map.drawparallels(np.arange(-90, 90, 30))
#lons = [8.547219]
#lats = [47.378578] 
#x,y = global_map(lons, lats)
#global_map.plot(x, y, marker='o', markerfacecolor='turquoise', markersize=12) #or use abbr. like: 'bo' -means blue circle and use different properties like: , color='b', linestyle='dashed'
#labels = ['SED']
#x_offsets = [10000]
#y_offsets = [5000]
#for label, xpt, ypt, x_offset, y_offset in zip(labels, x, y, x_offsets, y_offsets):
#	plt.text(xpt+x_offset, ypt+y_offset, label)
#plt.plot()
#plt.savefig('output/py_seismicity_schweiz_global.png')
##plt.show()
#
##regional map
#regional_map = Basemap(projection='merc', lat_0=7.5, lon_0=47.5, resolution='h', area_thresh=0.1, llcrnrlon=5, urcrnrlon=13, llcrnrlat=45.5, urcrnrlat=49) #also try: projection='ortho'
#regional_map.drawcoastlines()
#regional_map.drawcountries()
#regional_map.fillcontinents(color='KHAKI')
#regional_map.drawmapboundary()
#regional_map.drawmeridians(np.arange(0, 360, 30))
#regional_map.drawparallels(np.arange(-90, 90, 30))
#lons = [08.547219]
#lats = [47.378578] 
#x,y = regional_map(lons, lats)
#regional_map.plot(x, y, marker='o', markerfacecolor='turquoise', markersize=12) #or use abbr. like: 'bo' -means blue circle and use different properties like: , color='b', linestyle='dashed'
#labels = ['SED']
#x_offsets = [10000]
#y_offsets = [5000]
#for label, xpt, ypt, x_offset, y_offset in zip(labels, x, y, x_offsets, y_offsets):
#	plt.text(xpt+x_offset, ypt+y_offset, label)
#plt.plot()
#plt.savefig('output/py_seismicity_schweiz_regional.png')
##plt.show()
#
#hypocenter map
eq_map = Basemap(projection='merc', lat_0=7.5, lon_0=47.5, resolution='h', area_thresh=0.1, llcrnrlon=5, urcrnrlon=13, llcrnrlat=45.5, urcrnrlat=49) #also try: projection='ortho'
eq_map.drawcoastlines()
eq_map.drawcountries()
eq_map.fillcontinents(color='lightgray')
eq_map.drawmapboundary()
eq_map.drawmeridians(np.arange(0, 360, 30))
eq_map.drawparallels(np.arange(-90, 90, 30))
min_marker_size = 4
for lon, lat, mag, dep in zip(hypolon, hypolat, eqmag, hypodep):
    x,y = eq_map(lon, lat)
    msize = mag * min_marker_size
    marker_string = get_marker_color(dep)
    eq_map.plot(x, y, marker_string, markersize=msize)
labels = ['SED'] #label
lons = [08.547219]
lats = [47.378578]
xpts,ypts = eq_map(lons, lats)
for label, xpt, ypt in zip(labels, xpts, ypts):
	plt.text(xpt, ypt, label)
title_string = "Earthquakes of Magnitude 1.0 or Greater\n" #title
title_string += "%s through %s" % (starttime, endtime)
plt.title(title_string)
plt.plot()
plt.savefig('output/py_seismicity_schweiz_hypo.png')
#plt.show()

##---HISTOGRAM---#
##time series
#times = []
#for i in range(1, len(focaltime)):
#    dt = UTCDateTime(focaltime[i-1]) - UTCDateTime(focaltime[len(focaltime)-1])
#    dt = dt / 86400
#    times.append(dt)
#time_n, time_bins, time_patches = plt.hist(times, 20, normed=1, alpha=0.5, facecolor = 'green')
#time_mean = np.mean(times)
#time_std = np.std(times)
#time_y = mlab.normpdf(time_bins, time_mean, time_std)
#plt.plot(time_bins, time_y, 'r--')
#plt.xlabel("Magnitude 1+ event since 2016 [days]")
#plt.ylabel("count")
#plt.savefig('output/py_seismicity_schweiz_hist_time.png')
#plt.show()
#
##depth series
#depth_n, depth_bins, depth_patches = plt.hist(hypodep, normed=1, alpha=0.5, facecolor = 'green')
#depth_mean = np.mean(hypodep)
#depth_std = np.std(hypodep)
#depth_y = mlab.normpdf(depth_bins, depth_mean, depth_std)
#plt.plot(depth_bins, depth_y, 'r--')
#plt.xlabel("Depth since 2016 [km]")
#plt.ylabel("count")
#plt.savefig('output/py_seismicity_schweiz_hist_depth.png')
#plt.show()
#
##magnitude series
#mag_n, mag_bins, mag_patches = plt.hist(eqmag, normed=1, alpha=0.5, facecolor = 'green')
#mag_mean = np.mean(eqmag)
#mag_std = np.std(eqmag)
#mag_y = mlab.normpdf(mag_bins, mag_mean, mag_std)
#plt.plot(mag_bins, mag_y, 'r--')
#plt.xlabel("Mag. since 2016 [Ml]")
#plt.ylabel("count")
#plt.savefig('output/py_seismicity_schweiz_hist_mag.png')
#plt.show()
