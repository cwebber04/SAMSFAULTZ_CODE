#!/Users/timothy.lee/anaconda/bin/python

#exercise to retrieve the waveforms from sed-fdsn service

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

#---DATA_FROM_FDSN---#
#get specified event
client = Client(base_url = "http://arclink.ethz.ch")
starttime = UTCDateTime("2016-05-20")
endtime = UTCDateTime("2016-05-22")
cat = client.get_events(starttime=starttime, endtime=endtime, minmagnitude=2, limit=5, mindepth=5)
print(type(cat)) #Catalog
print(cat)
cat.plot() #add also resouces like: projection="local"

#get stations with the event
evt = cat[1]
print(type(evt))
print(evt)
origin = evt.origins[0]
otime = origin.time
print(type(origin))
t = origin.time
inv = client.get_stations(longitude=origin.longitude, latitude=origin.latitude, maxradius=0.2, starttime=t, endtime =t+100, channel="LH*", network="CH", level="station")
print(type(inv))
print(inv)
inv.plot(projection="local")

#get the waveforms
st = Stream()
for network in inv:
    for station in network:
        try:
            st += client.get_waveforms(network.code, station.code, "*", "LH*", t - 5 * 60, t + 30 * 60, attach_response=True)
        except:
            pass
#print(type(st))
#print(st)
#st.select(component="Z").plot(bgcolor="#FF5733")


#remove instrument response
print("...processing: remove instr. response")
st.remove_response(output='DISP', pre_filt=(0.005, 0.006, 30.0, 35.0), water_level=20)
print(type(st))
print(st)
st.select(component="Z").plot(bgcolor="#FF5733")

#basic processing
for tr in st:
    print(tr.id)
print("...processing: trim")
st.trim(otime, otime+10*60)
st.select(component="Z").plot(bgcolor="#FF5733")
print("...processing: remove trend")
st.detrend("linear")
st.select(component="Z").plot(bgcolor="#FF5733")
print("...processing: tapor")
st.taper(type="hann", max_percentage=0.05)
st.select(component="Z").plot(bgcolor="#FF5733")
print("...processing: filter")
st.filter("lowpass", freq=0.5)
st.select(component="Z").plot(bgcolor="#FF5733")

#final plot in both time and frequency domain
#st.plot(bgcolor="#FF5733")
#st.spectrogram(log=True, wlen=50);

#write sac filefor tr in st: 
print("writing files...")
for tr in st: 
	tr.write(tr.id + ".SAC", format="SAC")

print("program ended")