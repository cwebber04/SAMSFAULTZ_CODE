#!/Users/timothy.lee/anaconda/bin/python

#---RESOURCES---#
#python documentation: http://docs.python.org/
#python basic: https://github.com/krischer/seismo_live/blob/master/notebooks/Python%20Courses/Python_Crash_Course-with_solutions.ipynb
#obspy documentation: https://docs.obspy.org/
#obspy baisc: https://github.com/krischer/seismo_live/blob/f87def2c59794812ecb6ef9d292d7d89b9ffa515/notebooks/ObsPy
#numpy cheatsheet: http://mathesaurus.sourceforge.net/matlab-numpy.html
#ipython notebooks: https://github.com/ipython/ipython/wiki/A-gallery-of-interesting-IPython-Notebooks
#jupyter documentation: http://jupyter.readthedocs.io/en/latest/running.html#running

#---LIBRARIES---#
#install cartopy: conda install -c scitools cartopy
#install basemap: conda install basemap
#complete dictionary of colors:
#for name, hex in matplotlib.colors.cnames.items():
#    print(name, hex)
from obspy import read, read_inventory, read_events, UTCDateTime, Trace, Stream
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import URL_MAPPINGS
from obspy.geodetics import locations2degrees
from obspy.taup import TauPyModel
from mpl_toolkits.basemap import Basemap
import matplotlib.pylab as plt
import numpy as np
import numpy.matlib as mat
import scipy as Sci
import scipy.linalg

#---PARAMETERIZATION---#
plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = 12, 8

#---DATA_FROM_XML---#
inventory = read_inventory("./data/station_PFO.xml", format="STATIONXML")

#---DATA_FROM_CATALOGUE---#
catalog = read_events("*.xml")
print(type(catalog)) #Catalog
print(catalog)
catalog.plot(projection="local");

event = catalog[0]
print(type(event)) #Event
print(event)

coords = inventory.get_coordinates("II.BFO.00.BHE")
coords
origin = catalog[0].preferred_origin()
print(origin)
distance = locations2degrees(origin.latitude, origin.longitude, coords["latitude"], coords["longitude"])
print(distance)

m = TauPyModel(model="ak135")
arrivals = m.get_ray_paths(
    distance_in_degree=distance,
    source_depth_in_km=origin.depth / 1000.0)
arrivals.plot();

first_arrival = origin.time + arrivals[0].time
print(first_arrival)

print(type(event.magnitudes)) #choose between other resources, like: event.ortigins
print(type(event.magnitudes[0]))
print(event.magnitudes[0])

largest_magnitude_events = catalog.filter("magnitude >= 7.0")
print(largest_magnitude_events)
largest_magnitude_events.write("/tmp/large_events.xml", format="QUAKEML")


#---DATA_FROM_FDSN---#
#client = Client("IRIS") #choose between other services, like: "ETH" "GFZ" "USGS"
client = Client(base_url = "http://arclink.ethz.ch")
#waveform
t = UTCDateTime("2010-02-27T06:45:00.000")
st = client.get_waveforms("IU", "ANMO", "*", "LHZ", t, t + 60 * 60, attach_response=True) #add filename="/tmp/requested_waveforms.mseed") to save the waveform data
print(type(st)) #Stream
print(st)
st.plot()
#event
starttime = UTCDateTime("2002-01-01")
endtime = UTCDateTime("2002-01-02")
cat = client.get_events(starttime=starttime, endtime=endtime, minmagnitude=6, catalog="ISC") #add filename="/tmp/requested_events.xml" to save the event data, add also other resources: , minmagnitude=7, limit=5, mindepth=400
print(type(cat)) #Catalog
print(cat)
cat.plot()
#station
starttime = UTCDateTime("2002-01-01")
endtime = UTCDateTime("2002-01-02")
lat = 9.0
lon = 47.5
#inv = client.get_stations(network="IU", station="A*", starttime=starttime, endtime=endtime, level="station") #choose between other level, like: "network" "channel" "response"
inventory = client.get_stations(longitude=lon, latitude=lat, minradius=90, maxradius=100, level="station", matchtimeseries=True) #add filename="/tmp/requested_stations.xml" to save the station data
print(type(inv)) #Inventory
print(inv)
network = inventory[0]
print(network)
station = network[0]
print(station)
channel = station[0]
print(channel)
print(channel.response)
inventory.plot()
inventory.plot_response(0.001)

#st = read("./data/*, starttime = t + 10 * 60", endtime = t + 12 * 60)
#st = read("./data/*.mseed", format="mseed")
st = read("*.SAC")
print(type(st)) #Stream
print(st.traces)
for tr in st:
	evt_time =  UTCDateTime(year = st[0].stats.sac.nzyear, julday = st[0].stats.sac.nzjday, hour = st[0].stats.sac.nzhour, minute = st[0].stats.sac.nzmin) + st[0].stats.sac.nzsec + st[0].stats.sac.nzmsec*0.1
	print("CURRENT Time is: ",UTCDateTime())
	print("EVENT Original Time is: ",evt_time)
	print("EVENT is ",(UTCDateTime() - evt_time) / 86400,"DAYS AGO")
	print("For All Components: ",st)
	print("For Z Component: ",st.select(component="Z"))
	print(st[tr].stats)
	print(st[tr].data.max())
	print(st[tr].data.mean())
	print(st[tr].data)
	print(np.abs(st[tr].data)
st.select(component="Z").plot()
st.select(component="Z").spectrogram(log=True, wlen=50)


#---BASIC_PROCESSING---#
st2 = st.copy()
st2.remove_response(inventory=inv, output='DISP', pre_filt=(0.005, 0.006, 30.0, 35.0), water_level=60)
#st2.plot()
st2.trim(st2.stats.starttime + 12 * 60, st2.stats.starttime + 14 * 60)
st2.detrend("linear")
#st2.plot()
st2.taper(type="hann", max_percentage=0.05)
st2.taper(type='cosine', max_percentage=0.05)
#st2.plot()
st2.filter("lowpass", freq=0.5)
#st2.plot()
st2.write("temp.SAC",format="SAC")

#---MAPPING_DATA---#
my_map = Basemap(projection='ortho', lat_0=50, lon_0=-100, resolution='l', area_thresh=1000.0)
my_map.drawcoastlines()
my_map.drawcountries()
my_map.fillcontinents(color='coral')
my_map.drawmapboundary()
my_map.drawmeridians(np.arange(0, 360, 30))
my_map.drawparallels(np.arange(-90, 90, 30))
plt.show()


