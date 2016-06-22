#!/Users/timothy.lee/anaconda/bin/python

#exercise to map the seismicity in switzerland, data from sed-fdsn service

#---LIBRARIES---#
from obspy import read, read_inventory, read_events, UTCDateTime, Trace, Stream
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import URL_MAPPINGS
from mpl_toolkits.basemap import Basemap
import matplotlib.pylab as plt
import numpy as np
import numpy.matlib as mat
#import scipy as Sci
#import scipy.linalg


#---PARAMETERIZATION---#
plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = 12, 8
def get_marker_color(depth):
    if depth < 5.0:
        return ('bo')
    elif depth < 10.0:
        return ('co')
    elif depth < 20.0:
        return ('go')
    elif depth < 30.0:
        return ('yo')
    elif depth < 40.0:
        return ('mo')
    else:
        return ('ro')

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
#plt.show()
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
#plt.show()
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
plt.show()

#---HISTOGRAM---#
##time series
#inter_event_times = []
#for i in range(1, len(focaltime)):
#    dt = UTCDateTime(focaltime[i-1]) - UTCDateTime(focaltime[len(focaltime)-1])
#    dt = dt / 86400
#    inter_event_times.append(dt)
#plt.hist(inter_event_times)#, bins=range(0, 200, 10))
#plt.xlabel("Magnitude 1+ interevent times since 2016 [days]")
#plt.ylabel("count")
#plt.show()

#depth series
#plt.hist(hypodep)#, bins=range(0, 200, 10))
#plt.xlabel("Depth since 2016 [km]")
#plt.ylabel("count")
#plt.show()

#magnitude series
#plt.hist(eqmag)#, bins=range(0, 200, 10))
#plt.xlabel("Mag. since 2016 [Ml]")
#plt.ylabel("count")
#plt.show()
