#!/Users/timothy.lee/anaconda/bin/python

#map the seismicity near station:SION in switzerland, data from sed-fdsn service

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
    if depth < 2.0:
        return ('ro')
    elif depth < 4.0:
        return ('mo')
    elif depth < 6.0:
        return ('yo')
    elif depth < 8.0:
        return ('go')
    elif depth < 10.0:
        return ('co')
    else:
        return ('bo')

#---DATA_FROM_FDSN---#
#fetch station:SIOM information
client = Client(base_url = "http://arclink.ethz.ch", user='sed', password='sc31007arc')
starttime = UTCDateTime("2014-01-01")
endtime = UTCDateTime()
inv = client.get_stations(network="CH", station="SIOM", starttime=starttime, endtime=endtime, level="station") #choose between other level, like: "network" "channel" "response"
print(type(inv)) #Inventory
print(inv)
network = inv[0]
print(network)
station = network[0]
print(station)
#inv.plot(projection="local")
centerlat = station.latitude
centerlong = station.longitude

#get specified event
cat = client.get_events(starttime=starttime, endtime=endtime, latitude=centerlat, longitude=centerlong, maxradius=1, minmagnitude=1)#, filename="sion_events.xml")
#print len(cat)
evtnum = cat.count()
print(type(cat)) #Catalog
print(cat) #add , CatalogObject.__str__(print_all=True) to print all events
#cat.plot() #add also resouces like: projection="local"
focaltime = []
hypolon = []
hypolat = []
hypodep = []
eqmag = []
for x in range(0, evtnum-1):
	event = cat[x]
	focaltime.append(event.origins[0].time)
	hypolon.append(event.origins[0].longitude)
	hypolat.append(event.origins[0].latitude)
	hypodep.append(event.origins[0].depth/1000)
	eqmag.append(event.magnitudes[0].mag)

##---MAPPING_DATA---#
##hypocenter map
#eq_map = Basemap(projection='merc', lat_0=centerlat, lon_0=centerlong, resolution='h', area_thresh=0.1, llcrnrlon=centerlong-1, urcrnrlon=centerlong+1 , llcrnrlat=centerlat-0.5, urcrnrlat=centerlat+0.5)
#eq_map.drawcoastlines()
#eq_map.drawcountries()
#eq_map.fillcontinents(color='lightgray')
#eq_map.drawmapboundary()
#eq_map.drawmeridians(np.arange(0, 360, 30))
#eq_map.drawparallels(np.arange(-90, 90, 30))
#min_marker_size = 4
#for lon, lat, mag, dep in zip(hypolon, hypolat, eqmag, hypodep):
#    x,y = eq_map(lon, lat)
#    msize = mag * min_marker_size
#    marker_string = get_marker_color(dep)
#    eq_map.plot(x, y, marker_string, markersize=msize)
#labels = ['SIOM'] #label
#lons = [centerlong]
#lats = [centerlat]
#xpts,ypts = eq_map(lons, lats)
#for label, xpt, ypt in zip(labels, xpts, ypts):
#	plt.text(xpt, ypt, label)
#title_string = "Earthquakes of Magnitude 1.0 or Greater\n" #title
#title_string += "%s through %s" % (starttime, endtime)
#plt.title(title_string)
#plt.plot()
#plt.savefig('output/py_seismicity_sion_sequence_hypo.png')
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
