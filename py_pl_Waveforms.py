#!/Users/timothy.lee/anaconda/envs/python2.7/bin/python

import re
import math
import time
import numpy as np
import cPickle
from obspy.core import UTCDateTime
from obspy.core import read
from obspy.core.util import gps2DistAzimuth
from obspy.taup.taup import getTravelTimes
from obspy.arclink.client import Client
from obspy.fdsn import Client as Client_WS
#import pg
import os
import shutil
import glob
import string
import matplotlib.pyplot as plt
import optparse

global VERY_SMALL_DOUBLE, SMALL_DOUBLE
VERY_SMALL_DOUBLE = 1.0e-30
SMALL_DOUBLE = 1.0e-8

#This code needs 3 input files:
# 1) Location information
# 2) Files with waveform and station to plot

###########################################################################################
# Classes:
#------------------------------------------------------------------------------------------
class hypo(object):

      def __init__(self,yy,mm,dd,hh,mi,ss,timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mtype,mnobs,merr,mmeth,rms,gap,mdist,nobs,etype,lqual,chx,chy,agency,evID1,evID2,evPID,orID1,methodID,earthmodelID,author,region,fstatus):
          self.yy = int(yy)
          self.mm = int(mm)
          self.dd = int(dd)
          self.hh = int(hh)
          self.mi = int(mi)
          self.ss = float(ss)
          self.timestamp = float(timestamp)
          self.lat = float(lat)
          self.lon = float(lon)
          self.dep = float(dep)
          self.lonerr = float(lonerr)
          self.laterr = float(laterr)
          self.deperr = float(deperr)
          self.mag = float(mag)
          self.mtype = str(mtype)
          self.mnobs = int(mnobs)
          self.merr  = float(merr)
          self.mmeth = str(mmeth)
          self.rms = float(rms)
          self.gap = int(gap)
          self.mdist = float(mdist)
          self.nobs = int(nobs)
          self.etype = str(etype)
          self.lqual = str(lqual)
          self.chx   = float(chx)                 #CH-coordinates
          self.chy   = float(chy)                 #CH-coordinates
          self.agency = str(agency)
          self.evID1 = str(evID1)                 #Filename
          self.evID2 = long(evID2)                #DD-ID
          self.evPID = str(evPID)                 #publicID of event
          self.orID1 = str(orID1)                 #publicID of origin
          self.methodID = str(methodID)
          self.earthmodelID = str(earthmodelID)
          self.author = str(author)
          self.region = str(region)
          self.fstatus = str(fstatus)


      def __str__(self):
          return '%04d/%02d/%02d %02d:%02d:%06.3f %9.4f %8.4f %6.2f %4.1f %-5s %6.3f %3d %7.1f %3d %-2s %1s %3.0f %3.0f %-7s %-15s %-35s %-3s %15d %-15s %-s %-s' % (self.yy,self.mm,self.dd,self.hh,self.mi,self.ss,self.lon,self.lat,self.dep,self.mag,self.mtype,self.rms,self.gap,self.mdist,self.nobs,self.etype,self.lqual,self.chx,self.chy,self.agency,self.author,self.region,self.fstatus,self.evID2,self.methodID[0:15],self.evID1,self.orID1)

#------------------------------------------------------------------------------------------
class stationobj(object):

      def __init__(self,staorg,staali,net,loc,lat,lon,ele):
          self.staorg = str(staorg)
          self.staali = str(staali)
          self.net = str(net)
          self.loc = str(loc) 
          self.lat = float(lat)
          self.lon = float(lon)
          self.ele = float(ele)

      def __str__(self):
          return '%-7s %-7s %-7s %-7s %8.4f %9.4f %6.3f' % (self.staorg,self.staali,self.net,self.loc,self.lat,self.lon,self.ele)

#------------------------------------------------------------------------------------------
class arrivalpick(object):

      def __init__(self,yy,mm,dd,hh,mi,ss,timestamp,phasetype,errlo,errup,quality,onsetqual,polar,used,station,network,location,channel,mode,author,pweight,arrres,arrazi,arrdist,arrtakeoff,statmag,statmagwgt,statamp,statamptimestamp):
          self.yy = int(yy)
          self.mm = int(mm)
          self.dd = int(dd)
          self.hh = int(hh)
          self.mi = int(mi)
          self.ss = float(ss)
          self.timestamp = float(timestamp)
          self.phasetype = str(phasetype)
          self.errlo = float(errlo)
          self.errup = float(errup)
          self.quality = int(quality)
          self.onsetqual = str(onsetqual)
          self.polar = str(polar)
          self.used = int(used)
          self.station = str(station)
          self.network = str(network)
          self.location = str(location)
          self.channel = str(channel)
          self.mode = str(mode)
          self.author = str(author)
          self.pweight = float(pweight)
          self.arrres = float(arrres)
          self.arrazi = float(arrazi)
          self.arrdist = float(arrdist)
          self.arrtakeoff = float(arrtakeoff)
          self.statmag = float(statmag)
          self.statmagwgt = float(statmagwgt)
          self.statamp = float(statamp)
          self.statamptimestamp = float(statamptimestamp)

      def __str__(self):
          return '%-7s %-9s %-1s %-1s %1d %04d/%02d/%02d %02d:%02d:%06.3f %-3s %6.3f %6.3f %1d' % (self.station,self.phasetype,self.onsetqual,self.polar,self.used,self.yy,self.mm,self.dd,self.hh,self.mi,self.ss,self.channel,self.errlo,self.errup,self.quality)
#------------------------------------------------------------------------------------------
class stream2extract(object):

      def __init__(self,format,wavefile,netstat,componentlist,channellist,integ,filter,scale):
          self.format   = str(format)
          self.wavefile = str(wavefile)
          self.netstat  = str(netstat)
          self.componentlist = str(componentlist)
          self.channellist = str(channellist)
          self.integ       = int(integ)
          self.filter      = str(filter)
          self.scale       = float(scale)

      def __str__(self):
          return "%-10s %-10s %-15s %-10s %-s" % (self.netstat,self.componentlist,self.channellist,self.format,self.wavefile)
#------------------------------------------------------------------------------------------
class traceindex(object):

      def __init__(self,traceidx,eventidx,integ,filter,scale):
          self.traceidx   = int(traceidx)
          self.eventidx   = int(eventidx)
          self.integ      = int(integ)
          self.filter     = str(filter)
          self.scale      = float(scale)

      def __str__(self):
          return "%10d %10d" % (self.traceidx,self.eventidx)

#------------------------------------------------------------------------------------------
class pairCC(object):

      def __init__(self,id1,id2,stat,chan1,chan2,pha1,pha2,cc,lag,csf,edi,hdi,SN1,SN2,lat1,lon1,lat2,lon2,latS,lonS,epD1,epD2,rlag,dep1,dep2,resE,resH):
          self.id1 = long(id1)
          self.id2 = long(id2)
          self.stat = str(stat)
          self.chan1 = str(chan1)
          self.chan2 = str(chan2)
          self.pha1 = str(pha1)
          self.pha2 = str(pha2)
          self.cc  = float(cc)
          self.lag = float(lag)
          self.csf = int(csf)        #Cycle-skip flag 
          self.edi = float(edi)      #InterEvent Epicentral  distance
          self.hdi = float(hdi)      #InterEvent Hypocentral distance
          self.SN1 = float(SN1)      #Signal-to-noise ratio trace 1
          self.SN2 = float(SN2)      #Signal-to-noise ratio trace 2
          self.lat1 = float(lat1)    #event   latitude  1
          self.lon1 = float(lon1)    #event   longitude 1
          self.lat2 = float(lat2)    #event   latitude  2
          self.lon2 = float(lon2)    #event   longitude 2
          self.latS = float(latS)    #station latitude  2
          self.lonS = float(lonS)    #station longitude 2
          self.epD1 = float(epD1)    #epicentral distance 1
          self.epD2 = float(epD2)    #epicentral distance 2
          self.rlag = float(rlag)    #Relative lag (effective lag not considering static offset), can be used for filtering the data
          self.dep1 = float(dep1)    #Not included in the original output of pyXCorr, added by mergeCCraw.py
          self.dep2 = float(dep2)    #Not included in the original output of pyXCorr, added by mergeCCraw.py
          self.resE = float(resE)    #Not included in the original output of pyXCorr, added by mergeCCraw.py
          self.resH = float(resH)    #Not included in the original output of pyXCorr, added by mergeCCraw.py

      def __str__(self):
          return '%9d %9d %-9s %-3s %-3s %-9s %-9s %6.3f %8.3f %2d %8.1f %8.1f %10.1f %10.1f %8.4f %9.4f %8.4f %9.4f %8.4f %9.4f %8.2f %8.2f %8.3f %7.2f %7.2f %7.3f %7.3f' % (self.id1,self.id2,self.stat,self.chan1,self.chan2,self.pha1,self.pha2,self.cc,self.lag,self.csf,self.edi,self.hdi,self.SN1,self.SN2,self.lat1,self.lon1,self.lat2,self.lon2,self.latS,self.lonS,self.epD1,self.epD2,self.rlag,self.dep1,self.dep2,self.resE,self.resH)
#------------------------------------------------------------------------------------------

##########################################################################################
#Subroutines/functions:

def getnetstatloc(netstat,chanstring):

    clist = string.split(chanstring,",",-1)
    netw = ""
    stat = ""
    loca = ""

    netstatsp = string.split(netstat,".",-1)

    if(len(netstatsp) == 2):
       case = 2
       netw = netstatsp[0]
       stat = netstatsp[1]
       loca = ""
    if(len(netstatsp) == 3):
       case = 3
       netw = netstatsp[0]
       stat = netstatsp[1]
       loca = netstatsp[2]
    if(len(netstatsp) != 2) and (len(netstatsp) != 3):
       case = 1
       netw = ""
       loca = ""
       stat = netstatsp[0]


    return netw,stat,loca,clist

##########################################################################################

def getcomponents(station,trlist):
    clist = []

    #Get the channels for current station:
    #Loop over all traces in memory
    for j in range(len(trlist)):
        #Check id station codes match
        if(trlist[j].stats.station == station):
           #Check if channel already known:
           cnew = True
           for k in range(len(clist)):
               if(trlist[j].stats.channel[2:3]==clist[k]):
                  cnew = False
                  break
           if(cnew):
              clist.append(trlist[j].stats.channel[2:3])

    return clist
#------------------------------------------------------------------------------------------
def filtertrlist(trlist,tridxlist,station,component):

    filter_trlist = []
    filter_tridxlist = []
    filter_cnttrc = 0

    for i in range(len(trlist)):

        #Check if station matches:
        if(trlist[i].stats.station == station) and (trlist[i].stats.channel[2:3] == component):

           filter_cnttrc += 1
           filter_trlist.append(trlist[i])
           filter_tridxlist.append(tridxlist[i])
           filter_tridxlist[filter_cnttrc-1].traceidx = filter_cnttrc

    return filter_trlist,filter_tridxlist,filter_cnttrc
#------------------------------------------------------------------------------------------
def plotstreams(refevlist,trlist,tridxlist,cnttrc,alignflg,ytype,epimin,epimax,ofilefig,asciidump,dbgmode):

    #Open new figure:
    fig = plt.figure()
    plcnt = 0
    pevID = ""

    if(asciidump):
       asciidata = ofilefig + ".gmt.data"
       fpAD = open(asciidata,"w")

    #All in one plot:
    ax = fig.add_subplot(1,1,1)

    #Get y-range for plot:
    if(ytype == "trcnt"):
       #ymin = 1-(1.0*ascale)
       #ymax = cnttrc+(1.0*ascale)
       ymin = 1-(1.0)
       ymax = cnttrc+(1.0)

    if(ytype == "edist"):
       ymin = epimin
       ymax = epimax

    #Mark 0 time:
    ax.plot([0,0],[ymin,ymax],"gray",linewidth=1.0, linestyle='dotted')

    for i in range(cnttrc):

        #Get current event ID from filename:
        cevID = refevlist[tridxlist[i].eventidx].evID1
        cevDD = refevlist[tridxlist[i].eventidx].evID2

        #Check if station should be plotted:
        #if(len(stat2plot) > 0):
        #   noplot = True
        #   for m in range(len(stat2plot)):
        #       if(trlist[i].stats.station == stat2plot[m]):
        #          noplot = False
        #          break
        #   #Station not included -> skip trace
        #   if(noplot):
        #      continue

        #Check if Component should be ploted:
        #if(len(comp2plot) > 0):
        #   noplot = True
        #   for m in range(len(comp2plot)):
        #       if(trlist[i].stats.channel[2:3] == comp2plot[m]):
        #          noplot = False
        #          break
        #   #Component not included -> skip trace
        #   if(noplot):
        #      continue

        #Load potential picks:
        #Prior to SC3 load it from MANUPICK file, afterwards from database?, for the momment, from MANUPICK file:
        if(cevID != pevID):
           pkdir = "%04d/%02d" % (refevlist[tridxlist[i].eventidx].yy,refevlist[tridxlist[i].eventidx].mm)
           pkfil = KPevDir+"/"+pkdir+"/"+refevlist[tridxlist[i].eventidx].evID1+".MANUPICK"
           picks = loadMANUPICK(pkfil,dbgmode)

        #Now extract picks for station and sort them by time:
        statpicks = getpicks4station(trlist[i].stats.station,trlist[i].stats.channel,picks)

        #Get distance and azimuth:
        slat,slon = getcoordinates4station(trlist[i].stats.station,stationlist)
        if(slon < -400) or (slon < -400):
           print "WARNING: Station",trlist[i].stats.station,"not found in stationfile:",statlist,"--> Trace skipped"
           continue
        #         gps2DistAzimuth(latA,lonA,latB,lonB)
        #         
        #returns: Great circle distance in m, azimuth A->B in degrees, azimuth B->A in degrees
        epid=gps2DistAzimuth(slat,slon,refevlist[tridxlist[i].eventidx].lat,refevlist[tridxlist[i].eventidx].lon)

        #Demean
        trlist[i].detrend('demean')

        #Integrate if requested (before filter, used in the original version):
        #if(tridxlist[i].integ == 1):
        #   trlist[i].integrate('cumtrapz')
        #   #trlist[i].detrend('simple')

        #Filter trace:
        if(trlist[i].stats.channel[0:2] == "HH") or (trlist[i].stats.channel[0:2] == "EH") or (trlist[i].stats.channel[0:2] == "SH") or (trlist[i].stats.channel[0:2] == "HG") or (trlist[i].stats.channel[0:2] == "EL"):
           #print "--> Filter ",trlist[i].stats.channel,trlist[i].stats.station
           print trlist[i].stats.channel,trlist[i].stats.delta
           trlist[i].filter('bandpass', freqmin=1, freqmax=30, corners=2, zerophase=False)
           #To compare old telemetry SH with digital EH/BH use lowpas of about 12 Hz:
           #trlist[i].filter('bandpass', freqmin=1, freqmax=12, corners=2, zerophase=False)

        #Integrate if requested (after filter, new version):
        if(tridxlist[i].integ == 1):
           trlist[i].integrate('cumtrapz')
           #trlist[i].detrend('simple')

        #Normalize traces:
        cabst = 1E-25*-1
        for j in range(len(trlist[i].data)):
            if(abs(trlist[i].data[j]) > cabst):
               cabst = abs(trlist[i].data[j])

        #Check how traces should be aligned:
        if(alignflg != 0) and (alignflg != 10) and (alignflg != 20):
           print "WARNING: Unsupported value of parameter alignflg --> set alignflg=0"
           alignflg=0

        #Calculate synthetic traveltime in global model for all traces:
        #Problems with shallow sources:
        if(refevlist[tridxlist[i].eventidx].dep <= 1.5):
           print "--> WARNING: Depth <= 1.5 km - fix depth to 1.5 km"
           wdepth = 1.5
        else:
           wdepth = refevlist[tridxlist[i].eventidx].dep
        tt = getTravelTimes(delta=(epid[0]/1000.0)/111.1949, depth=wdepth, model='ak135')

        #Check if reduction velocity is used, if yes, force to align with respect to origin time:
        if(vred > 0.0):
           print "--> WARNING: Reduction velocity",vred,"set alignflg=0"
           alignflg = 0

        #Check how traces should be aligned:
        if(alignflg == 0):
           at = refevlist[tridxlist[i].eventidx].timestamp + tt[0]['time']
           if(vred <= 0.0):
              #Get difference between start of seismogram and origin time
              tdiff = trlist[i].stats.starttime - UTCDateTime(refevlist[tridxlist[i].eventidx].timestamp)
              str2 = "Time after origin (s)"
           else:
              #Get difference between start of seismogram and origin time
              tdiff = trlist[i].stats.starttime - UTCDateTime(refevlist[tridxlist[i].eventidx].timestamp) - ((epid[0]/1000.0)/vred)
              str2 = "Reduced Traveltime Tr = T -To - Dist/%04.1f (s)" % (vred)

        if(alignflg == 10):
           #Check if P-pick (P, P1, Pn, Pg) is available for station. If not use predicted
           if(len(statpicks) > 0) and ((statpicks[0].phasetype == "P") or (statpicks[0].phasetype == "P1") or (statpicks[0].phasetype == "Pn") or (statpicks[0].phasetype == "Pg")):
              at = statpicks[0].timestamp
           else:
              at = refevlist[tridxlist[i].eventidx].timestamp + tt[0]['time']
           tdiff = trlist[i].stats.starttime - UTCDateTime(at)
           str2 = "Time after picked first P (s)"

        if(alignflg == 20):
           #Get difference between start of seismogram and predicted first P-arrival
           at = refevlist[tridxlist[i].eventidx].timestamp + tt[0]['time']
           tdiff = trlist[i].stats.starttime - UTCDateTime(at)
           str2 = "Time after predicted first P (s) | earth model: ak135"

        #Now add possible CC-correction
        lagfound = 0
        if(corrlag) and (cevDD != corrMaE):
           for j in range(len(cclag)):
               if(cevDD == cclag[j].id1):
                  print cevDD,corrMaE," rlag:",cclag[j].rlag," CC-coef:",cclag[j].cc
                  tdiff = tdiff+cclag[j].rlag
                  lagfound = 1
                  break
               if(cevDD == cclag[j].id2):
                  print corrMaE,cevDD," rlag:",cclag[j].rlag," CC-coef:",cclag[j].cc
                  tdiff = tdiff-cclag[j].rlag
                  lagfound = 1
                  break
        if(corrlag) and (lagfound == 0) and (cevDD != corrMaE):
           print "--> WARNING No CC-dt for pair ",corrMaE,cevDD

        #Time with respect to origin time & normalization in display window:
        cabsw = 1E-25*-1
        time = np.zeros(trlist[i].stats.npts)
        for j in range(len(time)):
            time[j]=(j*(trlist[i].stats.delta))+tdiff
            #Check if time is within displayed window:
            if(time[j] >= plstwin) and (time[j] <=  plenwin) and (abs(trlist[i].data[j]) > cabsw):
               cabsw = abs(trlist[i].data[j])

        #Create some comment strings:
        str0 = "%-7s | %-s" % (trlist[i].stats.station,cevID)
        str1 = "%-s Dist: %5.1f Km, Dep: %6.1f, Mag: %4.1f" % (trlist[i].stats.channel,epid[0]/1000.0,refevlist[tridxlist[i].eventidx].dep,refevlist[tridxlist[i].eventidx].mag)
        #str2 = "Time after Origin (s) %04d/%02d/%02d %02d\:%02d\:%05.2f" % (refevlist[tridxlist[i].eventidx].yy,refevlist[tridxlist[i].eventidx].mm,refevlist[tridxlist[i].eventidx].dd,refevlist[tridxlist[i].eventidx].hh,refevlist[tridxlist[i].eventidx].mi,refevlist[tridxlist[i].eventidx].ss)

        #Add seismogramm to plot:
        plcnt += 1

        #Individual subplots for each trace:
        #ax = fig.add_subplot(cnttrc,1,plcnt)
        #ax.plot(time,trlist[i].data/cabst,'black')

        #What is the y-axis? -> number of trace:
        if(ytype == "trcnt"):
           ypos = cnttrc-plcnt+1
        if(ytype == "edist"):
           ypos = epid[0]/1000.0

        #Plot traces normalized by maximum of entire trace:
        if(normT == 1):
           if(corrlag) and (cevDD == corrMaE):
              ax.plot(time,ypos+(trlist[i].data/cabst)*ascale*tridxlist[i].scale,'red')
           else:
              ax.plot(time,ypos+(trlist[i].data/cabst)*ascale*tridxlist[i].scale,'black')
              if(asciidump):
                 fpAD.write("> > > >\n")
        #Plot traces normalized by maximum in display window:
        if(normT == 2):
           if(corrlag) and (cevDD == corrMaE):
              ax.plot(time,ypos+(trlist[i].data/cabsw)*ascale*tridxlist[i].scale,'red')
           else:
              ax.plot(time,ypos+(trlist[i].data/cabsw)*ascale*tridxlist[i].scale,'black')
              if(asciidump):
                 fpAD.write("> > > > %15.6E %15.6E %15.6E %15.6E\n" % (plstwin,plenwin,ymin,ymax))
                 for m in range(len(time)):
                     fpAD.write("%10.3f %20.11E\n" % (time[m],ypos+(trlist[i].data[m]/cabsw)*ascale*tridxlist[i].scale))
        #Backgroundcolors not correctly displayed in interactive tool, only in final eps file...
        #plt.text(plstwin, ypos+(0.25*ascale), str0, backgroundcolor='white', color='black', fontsize=8, ha='left', va='bottom')
        #plt.text(plstwin, ypos-(0.25*ascale), str1, backgroundcolor='white', color='black', fontsize=8, ha='left', va='top')
        plt.text(plstwin, ypos+(0.10), str0, fontsize=fsizetrace, ha='left', va='bottom')
        plt.text(plstwin, ypos-(0.10), str1, fontsize=fsizetrace, ha='left', va='top')

        #Add picks:
        if(len(statpicks) > 0):
           for j in range(len(statpicks)):
               #relative time:
               pt = UTCDateTime(statpicks[j].timestamp)-UTCDateTime(trlist[i].stats.starttime)+tdiff
               #P-Picks:
               if((statpicks[j].phasetype[0:1] == "P") or (statpicks[j].phasetype[0:1] == "p")) and ((pt >= plstwin) and (pt <= plenwin)):
                   ax.plot([pt,pt],[ypos-(0.25),ypos+(0.25)],"blue",linewidth=1.0)
                   plt.text(pt,ypos+(0.25),statpicks[j].phasetype,color="blue", fontweight='regular',ha='center', va='bottom', fontsize=fsizepick)
               #S-Picks:
               if((statpicks[j].phasetype[0:1] == "S") or (statpicks[j].phasetype[0:1] == "s")) and ((pt >= plstwin) and (pt <= plenwin)):
                   ax.plot([pt,pt],[ypos-(0.25),ypos+(0.25)],"red",linewidth=1.0)
                   plt.text(pt,ypos+(0.25),statpicks[j].phasetype,color="red", fontweight='regular',ha='center', va='bottom', fontsize=fsizepick)
        else:
           #No Pick avalable:
           #relative time:
           pt = UTCDateTime(at)-UTCDateTime(trlist[i].stats.starttime)+tdiff
           #P-Picks:
           if(pt >= plstwin) and (pt <= plenwin) and (vred <= 0.0):
              ax.plot([pt,pt],[ypos-(0.25),ypos+(0.25)],"blue",linewidth=1.0)
              plt.text(pt,ypos+(0.25),"Ppre",color="blue",ha='center', va='bottom', fontsize=fsizepick)

        #Store current eventID
        pevID = cevID

    #----------------------------------------------------------------
    #Show plot:
    if(cnttrc > 0):
       plt.xlim((plstwin,plenwin))
       plt.ylim((ymin,ymax))
       plt.xlabel(str2)


       #Show and save:
       if(savefig):

          if(savefor != "eps"):
             print "Outputformat",savefor,"not supported yet, use eps instead..."

          if(savefor == "eps"):
             ofilefig = ofilefig+"."+savefor
             plt.savefig(ofilefig, format=savefor, bbox_inches='tight', pad_inches=0.10)

          print ""
          print "Figure saved as:",ofilefig

       if(showfig):
          plt.show()

    else:
       print "No traces to plot..."

    if(asciidump):
       fpAD.close()
       print "Traces written to ASCII file:",asciidata

    return
#------------------------------------------------------------------------------------------
def celleb(lon,lat):

    #**********************
    #  lon
    #  SCHWEIZ. PROJEKTIONSSYSTEM  FORMELN VON H. ODERMATT
    #  TRANSFORMATION ELLIPSOID - EBENE
    #  L,B  laenge und breite in grad
    #  Y,X LANDESKOORDINATEN IN KILO-METER y= e-w; x=n-s
    #  MY  MERIDIANKONVERGENZ ( SEXAG. SEK.)
    #
    #  from nlloc/eth_custom/get_region_name_nr.c
    #*************************

    #To be consistent with orignal code:
    l = lon      #Laenge = longitude
    b = lat      #Breite = latitude

    #Initials:
    y = 999
    x = 999

    #Define TOP 8
    TOP = 8

    bb = 169028.66
    bl = 26782.5

    a = l
    a = a*3600.0 - bl
    c = b
    c = c * 3600.0 - bb
   
    d = []
    e = []
    f = []
    rw = []
    iw = []

    d.append(float( 0.0))                 #d00
    d.append(float( 0.68382546262761))    #d01
    d.append(float(-3.91798328045E-8))    #d02
    d.append(float( 1.4965410352E-15))    #d03
    d.append(float(-8.039471422E-23))     #d04
    d.append(float( 7.0021390E-30))       #d05
    d.append(float(-5.586904E-37))        #d06
    d.append(float( 4.0402E-44))          #d07
    d.append(float(-3.06E-51))            #d08

    e.append(float( 0.0))                 #e00
    e.append(float( 2.3635916074715E-2))  #e01
    e.append(float( 0.0))                 #e02
    e.append(float( 4.527219881E-17))     #e03
    e.append(float(-3.89081120E-24))      #e04
    e.append(float( 2.3407700E-31))       #e05
    e.append(float(-1.59674E-38))         #e06
    e.append(float( 1.287E-45))           #e07
    e.append(float( 0.0))                 #e08

    f.append(float( 0.0))                 #f00
    f.append(float( 4.515344386039E1))    #f01
    f.append(float( 1.17912305209E-4))    #f02
    f.append(float( 5.8474201864E-10))    #f03
    f.append(float( 2.73386187E-15))      #f04
    f.append(float( 1.4308547E-20))       #f05
    f.append(float( 7.66562E-26))         #f06
    f.append(float( 4.2445E-31))          #f07
    f.append(float( 2.40E-36))            #f08

    rw.append(float( 0.0))
    iw.append(float( 0.0))

    p = 30.91849390613 * a
    q = c * f[8]

    for i in range(TOP-1,0,-1):
        q = c * (q + f[i])

    rw.append(float(q)) #rw[1]
    iw.append(float(p)) #iw[1]

    for i in range(2,TOP+1,+1):
        rw.append(float(q*rw[i-1]-p*iw[i-1]))
        iw.append(float(p*rw[i-1]+q*iw[i-1]))

    dx = d[TOP]*rw[TOP]
    dy = d[TOP]*iw[TOP]
    my = 0.0

    for i in range(TOP-1,0,-1):
        dx=dx+d[i]*rw[i]
        dy=dy+d[i]*iw[i]
        my=my+e[i]*iw[i]

    dx=dx+200000.0
    dy=dy+600000.0
    dy= dy/1000.0
    dx= dx/1000.0
    y =  dy
    x =  dx

    #print "y/x:",y,"/",x 
    return y,x

#------------------------------------------------------------------------------------------
def getevIDfromfile(fname):
    t1 = string.split(fname,"/",-1)
    t2 = string.split(t1[len(t1)-1],".",-1)
    return t2[0]
#------------------------------------------------------------------------------------------
def getevidx(evlist,sevID):
    #EventID is KP filename:
    if(sevID[0:2] == "KP"):
       for i in range(len(evlist)):
           if(evlist[i].evID1 == sevID):
              return i,evlist[i].evID1,evlist[i].evID2
    #Assume EventID is DD-ID:
    else:
       for i in range(len(evlist)):
           if(evlist[i].evID2 == long(sevID)):
              return i,evlist[i].evID1,evlist[i].evID2

    return -1,'None',-1  #event not included
#------------------------------------------------------------------------------------------
def getstationfrominventory(inv,netstat,channellist):

    pstat = "XXXXX"
    pcha  = "XXX"
    pnet  = "XX"
    ploc  = "NOTEXISTING"

    clist = string.split(channellist,",",-1)

    netstatsp = string.split(netstat,".",-1)

    if(len(netstatsp) == 2):
       case = 2
       netw = netstatsp[0]
       stat = netstatsp[1]
       loca = ""
    if(len(netstatsp) == 3):
       case = 3
       netw = netstatsp[0]
       stat = netstatsp[1]
       loca = netstatsp[2]
    if(len(netstatsp) != 2) and (len(netstatsp) != 3):
       case = 1
       netw = ""
       loca = ""
       stat = netstatsp[0]

    #Get information on network level:
    #print inv['CH']
    #Get information on station level:
    #print inv['CH.TORNY']
    #Get information on channel level:
    #print inv['CH.TORNY..HHZ']

    #get all keys:
    #inv.keys()
    #get specific key:
    #inv.keys()[0]

    #Print the code of a station:
    #print inv['CH.LIENZ'].code
    #print inv['CH.LIENZ'].lat

    #Get gain of channel:
    #print inv['CH.EMBD..HHN'][0].gain

    #Case 1: only station is known:
    if(case == 1):
       #Check if preferred channel is included:
       for i in range(len(clist)):
           for j in range(len(inv.keys())):
               keysp = string.split(inv.keys()[j],".",-1)
               if(len(keysp) != 4):
                  continue
               else:
                  if(keysp[1] == stat) and (keysp[3][0:2] == clist[i][0:2]):
                     pstat = keysp[1]
                     pcha  = keysp[3][0:2]
                     ploc  = keysp[2]
                     pnet  = keysp[0]
                     return pstat,pcha,pnet,ploc

    if(case == 2):
       #Check if preferred channel is included:
       for i in range(len(clist)):
           for j in range(len(inv.keys())):
               keysp = string.split(inv.keys()[j],".",-1)
               if(len(keysp) != 4):
                  continue
               else:
                  if(keysp[1] == stat) and (keysp[0] == netw) and (keysp[3][0:2] == clist[i][0:2]):
                     pstat = keysp[1]
                     pcha  = keysp[3][0:2]
                     ploc  = keysp[2]
                     pnet  = keysp[0]
                     return pstat,pcha,pnet,ploc

    if(case == 3):
       #Check if preferred channel is included:
       for i in range(len(clist)):
           for j in range(len(inv.keys())):
               keysp = string.split(inv.keys()[j],".",-1)
               if(len(keysp) != 4):
                  continue
               else:
                  if(keysp[1] == stat) and (keysp[0] == netw) and (keysp[2] == loca) and (keysp[3][0:2] == clist[i][0:2]):
                     pstat = keysp[1]
                     pcha  = keysp[3][0:2]
                     ploc  = keysp[2]
                     pnet  = keysp[0]
                     return pstat,pcha,pnet,ploc               

    #If exact stream was defined (case 3 + channel), use this definition even though it appears to be missing in inventory
    if(case == 3) and (len(clist)==1):
       pstat=stat
       pcha=clist[0][0:2]
       ploc=loca
       pnet=netw               

    return pstat,pcha,pnet,ploc

#------------------------------------------------------------------------------------------
def getSACfile(sacdir,netstat,comp,channellist):

    clist     = string.split(channellist,",",-1)
    netstatsp = string.split(netstat,".",-1)

    if(len(netstatsp) == 2):
       netw = netstatsp[0]
       stat = netstatsp[1]
       loca = ""
    if(len(netstatsp) == 3):
       netw = netstatsp[0]
       stat = netstatsp[1]
       loca = netstatsp[2]
    if(len(netstatsp) != 2) and (len(netstatsp) != 3):
       netw = ""
       loca = ""
       stat = netstatsp[0]

    #Check if preferred channel is included:
    for i in range(len(clist)):
        cch = clist[i]+comp
        sfile = "%s.%s.%s.%s.%s.SAC" % (sacdir,netw,stat,cch,loca)
        if(os.path.isfile(sfile)):
           return sfile

    return 'None'
#------------------------------------------------------------------------------------------
def gettraceposition(stream,netstat,comp,channellist):
    
    clist = string.split(channellist,",",-1)
    
    #Default:
    streamdef = "Full"    

    netstatsp = string.split(netstat,".",-1)
    if(len(netstatsp) == 2):
       netw = netstatsp[0]
       stat = netstatsp[1]
       loca = ""
       streamdef = "net-stat"
    if(len(netstatsp) == 3):
       netw = netstatsp[0]
       stat = netstatsp[1]
       loca = netstatsp[2]
       streamdef = "Full"
    if(len(netstatsp) != 2) and (len(netstatsp) != 3):
       netw = ""
       loca = ""
       stat = netstatsp[0]
       streamdef = "stat"

    #Check if preferred channel is included:
    for i in range(len(clist)):
        cch = clist[i]+comp
        for j in range(len(stream)):
            #print netw,"|",stat,"|",loca,"|",cch,"||",stream[j].stats.network,"|",stream[j].stats.station,"|",stream[j].stats.location,"|",stream[j].stats.channel 
            if(streamdef == "stat"):
               if(stream[j].stats.station == stat) and (stream[j].stats.channel == cch):
                  return j
            if(streamdef == "net-stat"):
               if(stream[j].stats.network == netw) and (stream[j].stats.station == stat) and (stream[j].stats.channel == cch):
                  return j
            if(streamdef == "Full"):
               if(stream[j].stats.network == netw) and (stream[j].stats.station == stat) and (stream[j].stats.location == loca) and (stream[j].stats.channel == cch):
                  return j

    return -1
#------------------------------------------------------------------------------------------
def getcoordinates4station(station,stationlist):

    lat = -999
    lon = -999

    for i in range(len(stationlist)):
        if(stationlist[i].staorg == station):
           return stationlist[i].lat,stationlist[i].lon

    return lat,lon

#------------------------------------------------------------------------------------------
def getpicks4station(station,channel,picklist):

    list = []

    #yy,mm,dd,hh,mi,ss,timestamp,phasetype,errlo,errup,quality,onsetqual,polar,used,station,network,location,channel,mode,author,pweight,arrres,arrazi,arrdist,arrtakeoff,statmag,statmagwgt,statamp,statamptimestamp
    for i in range(len(picklist)):
        #Check for match of station name and channel-type (if provided with pick):
        if(picklist[i].station == station) and (len(picklist[i].channel) >= 2) and (picklist[i].channel[0:2] == channel[0:2]):
           list.append(picklist[i])
        #Check only for match of station name (if no channel was provided with pick):
        if(picklist[i].station == station) and (len(picklist[i].channel) <  2):
           list.append(picklist[i])
    
    #Now sort with respect to arrivaltime:
    slist = sorted(list, key=lambda x: x.timestamp, reverse=False)

    return slist

#------------------------------------------------------------------------------------------
def loadstream2extract4Energy(ifile,evlist,dbgmode):

    list = []
    Elist = []

    #Open
    fp = open(ifile,"r")

    for line in fp:

        #Test if line is empty:
        if(len(line.split())<1):
           print "Line not considered:", line

        else:
           t = string.split(line,None,-1)

           #Check if line is commented:
           if(t[0].startswith("#")):
              continue

           format  = "AUTO"
           netstat = t[0]
           comlst  = "Z,N,E"
           #Use this for borehole stations like STIEG
           #comlst  = "Z,2,3"
           chalst  = t[1]

           #Check if filter is defined:
           f = string.split(line,"FILTER",-1)
           if(len(f)==2):
              s = string.split(f[1],None,-1)
              filter = s[0]
           else:
              filter = 'default'

           #Check if scale is defined:
           f = string.split(line,"SCALE",-1)
           if(len(f)==2):
              s = string.split(f[1],None,-1)
              scale = float(s[0])
           else:
              scale = 1.0

           #Check if integration is requested:
           f = string.split(line,"INTEG",-1)
           if(len(f)==2):
              s = string.split(f[1],None,-1)
              integ = int(s[0])
           else:
              integ = 0


           #Add station to Elist:
           Elist.append(netstat)

           #Now loop over events to be extracted:
           for k in range(len(evlist)):
               
               list.append(stream2extract(format,evlist[k].evID1,netstat,comlst,chalst,integ,filter,scale))
 
               if(dbgmode == 1):
                  print stream2extract(format,evlist[k].evID1,netstat,comlst,chalst,integ,filter,scale)


    #Close
    fp.close()

    return list,Elist

#------------------------------------------------------------------------------------------

def loadstream2extract(ifile,dbgmode):

    list = []

    #Open
    fp = open(ifile,"r")

    for line in fp:
        #Test if line is empty:
        if(len(line.split())<1):
           print "Line not considered:", line

        else:
           t = string.split(line,None,-1)

           #Check if line is commented:
           if(t[0].startswith("#")):
              continue

           format  = t[0]
           wafile  = t[1]
           netstat = t[2]
           comlst  = t[3]
           chalst  = t[4]

           #Check if filter is defined:
           f = string.split(line,"FILTER",-1)
           if(len(f)==2):
              s = string.split(f[1],None,-1)
              filter = s[0]
           else:
              filter = 'default'

           #Check if scale is defined:
           f = string.split(line,"SCALE",-1)
           if(len(f)==2):
              s = string.split(f[1],None,-1)
              scale = float(s[0])
           else:
              scale = 1.0

           #Check if integration is requested:
           f = string.split(line,"INTEG",-1)
           if(len(f)==2):
              s = string.split(f[1],None,-1)
              integ = int(s[0])
           else:
              integ = 0

           list.append(stream2extract(format,wafile,netstat,comlst,chalst,integ,filter,scale))

           if(dbgmode == 1):
              print stream2extract(format,wafile,netstat,comlst,chalst,integ,filter,scale)


    #Close
    fp.close()

    return list

#------------------------------------------------------------------------------------------
def loadMANUPICK(file,dbgmode):

    list = []

    if(dbgmode == 1):
       print "Load picks from MANUPICK file:",file

    #Open
    fp = open(file,"r")

    refflg = 0

    for line in fp:
        #Test if line contains reference time:
        ss1 = string.split(line,"/",-1)
        ss2 = string.split(line,":",-1)

        if(len(ss1) == 3) and (len(ss2) == 2) and (refflg == 0):
           refflg = 1
           yy = int(line[0:4])
           mm = int(line[5:7])
           dd = int(line[8:10])
           hh = int(line[11:13])
           mi = int(line[14:16])

           rdate = UTCDateTime("%04d-%02d-%02d %02d:%02d:%09.6f" % (yy,mm,dd,hh,mi,0.0))
        
           if(dbgmode == 1):   
              print line
              print "Reference time found: %04d %02d %02d %02d %02d %-s" % (yy,mm,dd,hh,mi,rdate.datetime)

        #Extract phases:
        if(refflg == 1 ) and (line[54:58] != "SKIP") and ((line[8:9] == "P") or (line[8:9] == "p") or (line[8:9] == "S") or (line[8:9] == "s")):
           sta = string.strip(line[0:7])
           pha = string.strip(line[8:15])
           wqa = string.strip(line[16:17])
           pol = string.strip(line[17:18])
           art = float(line[19:26])
           ius = int(line[27:28])

           #Get pickdate:
           pdate = UTCDateTime(rdate.timestamp + art)

           #Unknown parameters:
           net = ""
           loc = ""
           cha = ""
           mod = ""
           aut = ""
           pwg = 0.0
           res = 0.0
           azi = 0.0
           dis = 0.0
           toa = 0.0
           statmag = 0.0
           statmagwgt = 0.0
           statamp = 0.0
           statamptimestamp = 0.0
    
           #Add pick:
           list.append(arrivalpick(pdate.year,pdate.month,pdate.day,pdate.hour,pdate.minute,pdate.second+((pdate.microsecond)/1.0E+6),pdate.timestamp,pha,0.0,0.0,9,wqa,pol,ius,sta,net,loc,cha,mod,aut,pwg,res,azi,dis,toa,statmag,statmagwgt,statamp,statamptimestamp))
 
        #End of file reached:
        if(line[54:58] == "SKIP"):
           break

    #Close
    fp.close()

    return list
#------------------------------------------------------------------------------------------

def loadstations(ifile,dbgmode):

    list = []

    #Open
    fp = open(ifile,"r") 

    for line in fp:
        #Test if line is empty:
        if(len(line.split())<5):
           print "Line not considered:", line

        else:
           #Test if line is comment:
           t = string.split(line,None,-1)
           if(t[0].startswith("#") == False):
              staorg = string.strip(line[0:5])
              staali = string.strip(line[9:13])
              lat    = float(line[15:23])
              lon    = float(line[23:32])
              ele    = float(line[33:39])
              net    = ""
              loc    = ""

              list.append(stationobj(staorg,staali,net,loc,lat,lon,ele))
              if(dbgmode == 1):
                 print stationobj(staorg,staali,net,loc,lat,lon,ele)
              
    #Close
    fp.close()

    return list
#------------------------------------------------------------------------------------------

def loadrefevents(ifile,dbgmode):

    list = []

    #Open
    fp = open(ifile,"r")

    for line in fp:
        #Test if line is empty:
        if(len(line.split())<10):
           print "Line not considered:", line
        else:
           #Test if line is comment:
           #t=line.split()
           t = string.split(line,None,-1)
           s = string.split(line,"|",-1)

           if(t[0].startswith("#") == False):
              lon=float(line[0:8])
              lat=float(line[9:16])
              dep=float(line[17:22])
              mag=float(line[23:27])



              mat=string.strip(line[28:30])
              yy=int(line[31:35])
              mm=int(line[36:38])
              dd=int(line[39:41])
              hh=int(line[42:44])
              mi=int(line[45:47])
              ss=float(line[48:53])
              loc=string.strip(line[54:58])
              rms=float(line[60:65])
              gap=int(line[66:69])
              mdi=float(line[70:75])
              nob=int(line[76:79])
              ety=string.strip(line[80:81])
              lqa=string.strip(line[82:83])
              ID1=string.strip(line[84:104])    #File-ID (string)
              ID2=long(line[105:114])           #DD-evID (long)     

              #Get swiss coordinates:
              chx,chy = celleb(lon,lat)
              if(chy < 62.0) or (chy > 302.0) or (chx < 480.0) or (chx > 847.5):
                chy = 999
                chx = 999
              #Check if event was extracted from KP file structure:
              #if(ID1[0:2]=="KP"):
              if(loc != "SEDS"):
                 agy    = "SED_KP"
                 region = "None"
                 author = "KP"
                 ID3    = "None"          #Origin publicID (string)
                 ID4    = "None"          #Event  publicID (string)
                 fstat  = "REF"
                 laterr = -9
                 lonerr = -9
                 deperr = -9
                 mnobs  = -9
                 merr   = -9
                 modlID = '?'
                 magmeth = '?'

              if(loc == "SEDS"):
                 ID3    = string.strip(s[5])
                 ID4    = string.strip(s[6])
                 agy    = "SED"
                 region = string.strip(s[4])
                 author = string.strip(s[3])
                 fstat  = "REF"
                 laterr = -9
                 lonerr = -9
                 deperr = -9
                 mnobs  = -9
                 merr   = -9
                 modlID = '?'
                 magmeth = '?'

              #Check if event was extracted from SC3 DB structure:
              #STILL MISSING 
              #else:

              #Check if format is correct:
              if(ss < 60.0) and (mi < 60) and (hh < 24) and (dd < 32) and (mm < 13):
                 wdate = UTCDateTime("%04d-%02d-%02d %02d:%02d:%06.3f" % (yy,mm,dd,hh,mi,ss)) #Problems with seconds close to 60.0 e.g.: 59.99995 -> 60.0 Fix: increase precision
                 #wdate = UTCDateTime("%04d-%02d-%02d %02d:%02d:%09.6f" % (yy,mm,dd,hh,mi,ss))

                 #print ID1
                 list.append(hypo(yy,mm,dd,hh,mi,ss,wdate.timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mat,mnobs,merr,magmeth,rms,gap,mdi,nob,ety,lqa,chx,chy,agy,ID1,ID2,ID4,ID3,loc,modlID,author,region,fstat))
                 if(dbgmode == 1):
                    print hypo(yy,mm,dd,hh,mi,ss,wdate.timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mat,mnobs,merr,magmeth,rms,gap,mdi,nob,ety,lqa,chx,chy,agy,ID1,ID2,ID4,ID3,loc,modlID,author,region,fstat)
              else:
                 print "Warning: Corrupt origin time ",datetime2(yy,mm,dd,hh,mi,ss),ID1
           else:
              print "Line not considered:", line

    #Close
    fp.close()

    return list
#------------------------------------------------------------------------------------------
def readraw(ifile,minCC,CC_list,corrPha,corrSta,corrMaE):

    #Open
    fp = open(ifile,"r")

    for line in fp:
        #Test if line is empty:
        if(len(line.split()) != 27):
           print "Line not considered:", line
        else:
           #Test if line is comment:
           #t=line.split()
           t = string.split(line,None,-1)

           if(t[0].startswith("#") == False) and (float(t[7]) >= minCC) and ((long(t[0]) == corrMaE) or (long(t[1]) == corrMaE)) and (t[2] == corrSta) and (t[5][0:1] == corrPha):
              #OLD - Format of line  : id1,id2,stat,chan1,chan2,pha1,pha2,cc,lag,csf,edi,hdi,SN1,SN2,lat1,lon1,lat2,lon2,latS,lonS,epD1,epD2              
              #NEW - Format of line  : id1,id2,stat,chan1,chan2,pha1,pha2,cc,lag,csf,edi,hdi,SN1,SN2,lat1,lon1,lat2,lon2,latS,lonS,epD1,epD2,rlag,dep1,dep2,resE,resH
              #OLD - Format of object: id1,id2,stat,chan1,chan2,pha1,pha2,cc,lag,csf,edi,hdi,SN1,SN2,lat1,lon1,lat2,lon2,latS,lonS,epD1,epD2
              #NEW - Format of object: id1,id2,stat,chan1,chan2,pha1,pha2,cc,lag,csf,edi,hdi,SN1,SN2,lat1,lon1,lat2,lon2,latS,lonS,epD1,epD2,rlag,dep1,dep2,resE,resH
              CC_list.append(pairCC(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],t[14],t[15],t[16],t[17],t[18],t[19],t[20],t[21],t[22],t[23],t[24],t[25],t[26]))

    fp.close()

    return CC_list
#------------------------------------------------------------------------------------------

########################################################################
#Main code:

##################################################################################################################
#Get command line arguments:
oparser = optparse.OptionParser()
oparser.add_option('-f', '--file',               action='store', dest='WaveStatList',      help='List containing station-channel-events to plot; e.g. ./Sion_20130406.txt')
oparser.add_option('-e', '--eventlist',          action='store', dest='extractlist',       help='Reference list of hypocenters (SC3DB2KP.py format; default: /Users/tdiehl/lib/MANULOC_List_1984_2012_ECOS_short.log.latest)')
oparser.add_option('-s', '--stationlist',        action='store', dest='statlist',          help='Station list (MPX format; default: /Users/tdiehl/lib/sed_stations.GSE_SED.alias)')
oparser.add_option('-d', '--directoryKP',        action='store', dest='KPevDir',           help='Directory with KP MANUPICK/MANUPDE files; default: /Users/tdiehl/data/Swiss/events')
oparser.add_option('-y', '--y-axis-value',       action='store', dest='ytype',             help='What should be used for y-axes: trcnt: By trace-number, order of input file is kept (default); edist: use epicentral distance')
oparser.add_option('-r', '--y-axis-range',       action='store', dest='yrange',            help='y-range in km to plot (only used in ytype=edist): 0/200 (default) ')
oparser.add_option('-v', '--red-velocity',       action='store', dest='redvel',            help='reduction velocity in km/s (if <= 0, not used): 0 (default) ')
oparser.add_option('-a', '--align',              action='store', dest='alignflg',          help='Align with respect to (reference time): 0  = origin time; 10 = P-Pick (loaded from MANUPIC file) (default); 20 = Predicted P (from ak135)')
oparser.add_option('-z', '--scale',              action='store', dest='ascale',            help='Amplitude scale (normalized amplitudes * ascale), default: 1.0')
oparser.add_option('-o', '--output-eps',         action='store', dest='ofilefig',          help='Name of EPS file (without extensions; e.g. Sion_20130406); default: output')
oparser.add_option('-b', '--begin-window',       action='store', dest='plstwin',           help='Start of displayed window (with respect to reference time) in seconds; default: -0.5')
oparser.add_option('-n', '--end-window',         action='store', dest='plenwin',           help='End   of displayed window (with respect to reference time) in seconds; default: +1.0')
oparser.add_option('-E', '--Energy-Analysis',    action='store', dest='EnergyStationList', help='Perform Energy-Analysis for S-wave for stations specified in specified list (format is pyXCorr-stat-list); e.g. -E /Users/tdiehl/Documents/Swiss_Alps_Docu/2014_Diemtigen/Stations_Part_All.txt')
oparser.add_option('-p', '--arclink-profile',    action='store', dest='arclinkprofile',    help='define arclink profile: "SED-Internal"; "SED-Bhutan"; "SED-Internal-WS" (default)')
oparser.add_option('-w', '--waveform-archive',   action='store', dest='wavearc',           help='local archive with waveforms; default: /Users/tdiehl/data/waveforms/pl_Waveforms')
oparser.add_option('-S', '--Store-waveform',                action='store_true', dest='wavesto',      help='If set, waveforms are stored on disk in directory set with -w; default: False')
oparser.add_option('-L', '--Load-first-from-local-archive', action='store_true', dest='wavella',      help='If set, trying to load waveforms from local disk in directory set with -w, if waveform is not available try to load via arclink; default: False')
oparser.add_option('-l', '--Load-only-from-local-archive',  action='store_true', dest='wavello',      help='If set, trying to load waveforms from local disk in directory set with -w, no arclink request is done')

(options, arg) = oparser.parse_args()

#Check for arclinkprofile:
if(options.arclinkprofile != None):
   arclinkprofile = options.arclinkprofile
else:
   arclinkprofile = 'SED-Internal-WS'

#Check if waveforms should be stored on disk:
if(options.wavesto != None):
   wavesto = True
else:
   wavesto = False

#Check if waveforms should be first loaded from disk:
if(options.wavella != None):
   wavella = True
else:
   wavella = False

#Check if waveforms should be loaded only from disk:
if(options.wavello != None):
   wavello = True
else:
   wavello = False

#Check if waveform directory is specified:
if(options.wavearc != None):
   wavearc = options.wavearc
else:
   wavearc = 'pl_Waveforms'

#Check of Energy-analysis is requested:
if(options.EnergyStationList != None):
   EnergyStationList = options.EnergyStationList
   EnergyAnalysis = True
else:
   EnergyStationList = ""
   EnergyAnalysis = False

#Check if y-range is included as option:
if(options.yrange != None):
   yrangestr = options.yrange
   t = string.split(yrangestr,'/',-1)
   epimin=float(t[0])
   epimax=float(t[1])
else:
   epimin = 0
   epimax = 200

#Check if event list is included in list of input arguments:
if(options.extractlist != None):
   extractlist = options.extractlist
else:
   extractlist='lib/MANULOC_List_1984_2012_ECOS_short.log.latest'

#Waveform-Station list:
#Check if list is included in list of input arguments:
if(options.WaveStatList != None):
   WaveStatList = options.WaveStatList
else:
   #WaveStatList = "./Test_list.txt"
   #WaveStatList = "./Breithorn_20121128_list_SATI.txt"
   #WaveStatList = "./PizzoCastello_20130226_1157_FUSIO.txt"
   #WaveStatList = "./Rumisberg_2005.txt"
   WaveStatList = "./Sion_20130406.txt"

#List with station coordinates:
#Check if list is included in list of input arguments:
if(options.statlist != None):
   statlist = options.statlist
else:
   statlist = "lib/sed_stations.GSE_SED.alias"

#KP Event directory (directory where KP event files are stored (such as MANUPICK, MANUPDE), etc. 
#Check if directory is included in list of input arguments:
if(options.KPevDir != None):
   KPevDir = options.KPevDir
else:
   KPevDir = "/Users/tdiehl/data/Swiss/events"

#What should be used for y-axes:
# trcnt: By trace-number, order of input file is kept
# edist: use epicentral distance
#Check if option is included in list of input arguments:
if(options.ytype != None):
   ytype = options.ytype
else:
   ytype = "trcnt"

#Align switch:
#Align waveforms according to:
#0  = origin time
#10 = P-Pick (loaded from MANUPIC file)
#20 = Predicted P (from ak135)
#Other options (such as predicted by Nll) still missing
#Check if option is included in list of input arguments:
if(options.alignflg != None):
   alignflg = int(options.alignflg)
else:
   alignflg = 10

#Amplitude scale (normalized amplitudes * ascale):
#Check if option is included in list of input arguments:
if(options.ascale != None):
   ascale = float(options.ascale)
else:
   ascale = 1.00

#Save figure:
savefig = True
savefor = "eps"
#Check if option is included in list of input arguments:
if(options.ofilefig != None):
   ofilefig = options.ofilefig
else:
   ofilefig = 'output'
   #ofilefig = "PizzoCastello"
   #ofilefig = "Rumisberg_2005"
   #ofilefig = "Sierre_201304"

#Start-end of time window (with respect to reference time):
if(options.plstwin != None):
   plstwin = float(options.plstwin)
else:
   plstwin = -0.5
if(options.plenwin != None):
   plenwin = float(options.plenwin)
else:
   plenwin = 1.0

#Check if reduction velocity is set:
if(options.redvel != None):
   vred = float(options.redvel)
else:
   vred = -1.0

#Still Hardcoded:
#################
#Dump everything to ASCII File for GMT plot:
asciidump = True


#fontsize traces:
#fsizetrace = 8
fsizetrace = 4

#fontsize picks:
#fsizepick = 8
fsizepick = 4

#What window used for normalization:
#1: Entire trace
#2: Displayed window
normT = 2

#SC3 origins are preferred over MANULOC (reference) after:
SC3overKPdate = '2012-10-03 00:00:00'

#Default directory of waveforms (used for GSE2 files in automatic mode):
autodir_GSE2 = '/Volumes/share-sed-archive-3-$/events'

#Default directory of waveforms (used for SAC files in CH_SAC mode):
autodir_SAC = '/Users/tdiehl/data/HypoDD/Swiss/Waveforms'

#Arclink connection:
if(arclinkprofile == 'SED-Internal'):
   arcl_meth = 'arclink'
   #arcl_host = "eida.ethz.ch"
   arcl_host = "arclink.ethz.ch"
   #arcl_host = "rzseismo2.ethz.ch"
   #arcl_host = "arclinktest.ethz.ch"
   arcl_port = 18002
   #arcl_port = 18001
   #arcl_port = 18001
   arcl_user = "sc31007arc"
   ws_host   = 'http://arclink.ethz.ch'
   ws_user   = 'sed'
   ws_pasw   = 'sc31007arc'
   #Time window to etract (in seconds, still in respect to origin time, other definitions defined later):
   arcl_winst = -20.0
   arcl_winen = +240.0 

if(arclinkprofile == 'SED-Bhutan'):
   arcl_meth = 'arclink'
   #arcl_host = "arclink.ethz.ch"
   #arcl_host = "rzseismo2.ethz.ch"
   arcl_host = "arclinktest.ethz.ch"
   #arcl_port = 18002
   arcl_port = 18001
   arcl_user = "bhutan1156arc"
   ws_host   = 'http://arclink.ethz.ch'
   ws_user   = 'sed'
   ws_pasw   = 'sc31007arc'
   #Time window to etract (in seconds, still in respect to origin time, other definitions defined later):
   arcl_winst = -20.0
   arcl_winen = +240.0

if(arclinkprofile == 'SED-Internal-WS'):
   arcl_meth = 'FDSN-WS'
   arcl_host = 'arclink.ethz.ch'
   arcl_port = 18002
   arcl_user = "sc31007arc"
   ws_host   = 'http://arclink.ethz.ch'
   ws_user   = 'sed'
   ws_pasw   = 'sc31007arc'
   arcl_winst = -20.0
   arcl_winen = +240.0

#Load differentialtimes (differential time with respect to master event):
corrlag = False
#corrfil = '/Users/tdiehl/programs/python_codes/C-InterfaceTest/pyXCorr/combined_CC_run02_SomeOutliers.raw'
corrfil = '/Users/tdiehl/programs/python_codes/C-InterfaceTest/pyXCorr/combined_CC_run03_SomeOutliers.raw'
corrMaE = 200406030
corrSta = 'LKBD'
corrPha = 'P'
corrMin = 0.7

#Show figure in interactive tool:
showfig = True

#Debug-Mode (0 or 1):
dbgmode = 0

##################################################################################################################
#Load reference events:
print ("")
print ("--> Load MANULOC-List from: ",extractlist)
refevlist = loadrefevents(extractlist,dbgmode) 

#Load station coordinates:
print "--> Load Station-List from: ",statlist
stationlist = loadstations(statlist,dbgmode)

#Load Waveform-Station list:
if(EnergyAnalysis):
   streams2extract,EnergyList = loadstream2extract4Energy(EnergyStationList,refevlist,dbgmode)
else:
   streams2extract = loadstream2extract(WaveStatList,dbgmode)
   EnergyList=[]

#Load differential times for Phase and station:
cclag = []    #Empty list
if(corrlag):
   print ""
   print "--> Load CC-lags      from: ",corrfil," Master Event: ",corrMaE," for phase ",corrPha," station ",corrSta
   cclag = readraw(corrfil,corrMin,cclag,corrPha,corrSta,corrMaE)
   print "Number of CC-dts loaded:",len(cclag)
   for i in range(len(cclag)):
       print "%-s" % (cclag[i])

#Load inventory and initilize client via arclink:
print ""
print "--> Load inventory from",arcl_meth,arcl_host,arcl_port,arcl_user
arcl_client = Client(host=arcl_host, port=arcl_port, user=arcl_user)
inv = arcl_client.getInventory('*', '*', '*', '*',starttime=UTCDateTime("1984-01-01T00:00:00.0"),endtime=UTCDateTime("2020-01-01T00:00:00.0"))

#If waveforms should be loaded via WS, open WS-Client:
if(arcl_meth == 'FDSN-WS'):
   arcl_client = Client_WS(base_url=ws_host, user=ws_user, password=ws_pasw)

#Arc-link is preferred over GSE2 files after:
sc3date = UTCDateTime(SC3overKPdate)

#Initials:
cnttrc = 0
tridxlist = []
trlist = []
arcl_flg = 0

#Loop over streams to extract:
#format,wavefile,netstat,componentlist,channellist
pevID = ""
for i in range(len(streams2extract)):
    #Get current event ID from filename:
    cevID = getevIDfromfile(streams2extract[i].wavefile)

    #Check if column #2 is file or just a event ID
    if(os.path.isfile(streams2extract[i].wavefile)):
       IDisfile = True
       #print "ID is file"
    else:
       IDisfile = False
       #print "ID is ID"

    #Check if new event:
    if(cevID != pevID):
       #New event: Get index of event in refernce list:
       weidx,KPId,DDId = getevidx(refevlist,cevID)

       #Check if event was included:
       if(weidx < 0):
          print "WARNING: No event information found for file",streams2extract[i].wavefile,"--> skip file"
          continue

       #Check if event is in SC3 era:
       if(refevlist[weidx].timestamp < sc3date.timestamp):
          auto_era = 'KP'
       else:
          auto_era = 'SC'

    #Get components to extract:
    complist = string.split(streams2extract[i].componentlist,",",-1)

    #Get waveforms straight from SAC waveform file:
    if(streams2extract[i].format == "CH_SAC"):

       #Simply load one SAC file
       if(IDisfile):
          st = read(streams2extract[i].wavefile)

          #Number of traces to plot:
          cnttrc += 1

          #Append trace to list:
          trlist.append(st[0].copy())
          tridxlist.append(traceindex(cnttrc,weidx,streams2extract[i].integ,streams2extract[i].filter,streams2extract[i].scale))

          #Print information:
          print streams2extract[i],refevlist[weidx],0
       
       #Find matching SAC files (for multiple components):   
       else:
          #Loop over channels to extract:
          for j in range(len(complist)):
              #Setup searchpath:
              spath = "%s/%04d/%02d/%s" % (autodir_SAC,int(KPId[2:6]),int(KPId[6:8]),str(DDId))

              #Find matching file
              wfile = getSACfile(spath,streams2extract[i].netstat,complist[j],streams2extract[i].channellist)

              #File does not exist:
              if(wfile == 'None'):
                 print "Warning: No SAC file found in",spath,streams2extract[i].netstat,complist[j],streams2extract[i].channellist
                 continue
              #File found:
              else:
                 #Read SAC file
                 st = read(wfile)

                 #Number of traces to plot:
                 cnttrc += 1

                 #Append trace to list:
                 trlist.append(st[0].copy())
                 tridxlist.append(traceindex(cnttrc,weidx,streams2extract[i].integ,streams2extract[i].filter,streams2extract[i].scale))

                 #Print information:
                 print streams2extract[i],refevlist[weidx],0
            
    #Get waveforms straight from GSE2 waveform file:
    if(streams2extract[i].format == "KP_GSE2") or ((streams2extract[i].format == "AUTO") and (auto_era == 'KP')):
       #Setup filename:
       if(IDisfile):
          wfile = streams2extract[i].wavefile                                                     #Column 2 is Filename
       else:
          wfile = "%s/%04d/%02d/%s.GSE2" % (autodir_GSE2,int(KPId[2:6]),int(KPId[6:8]),KPId)      #Column 2 is event ID, use default GSE2 dierectory define in autodir_GSE2

       #Now load waveform file:
       st = read(wfile)

       #Now find position of channels to extract
       for j in range(len(complist)):
           tridx = gettraceposition(st,streams2extract[i].netstat,complist[j],streams2extract[i].channellist)

           #Check if event was included:
           if(tridx < 0):
              print "WARNING: No matching channels found for file",streams2extract[i],"--> skip file"
              continue

           #Number of traces to plot:
           cnttrc += 1

           #Append trace to list:
           trlist.append(st[tridx].copy())
           tridxlist.append(traceindex(cnttrc,weidx,streams2extract[i].integ,streams2extract[i].filter,streams2extract[i].scale))

           #Print information:
           print streams2extract[i],refevlist[weidx],tridx
 
    #Before extraction the data from arclink, check if they are available in the local archive:
    notinlocal = True

    #Split station name into network, code, location, and convert channel string to list
    lnet,lstat,lloc,lcha = getnetstatloc(streams2extract[i].netstat,streams2extract[i].channellist)

    if((streams2extract[i].format == "SC_ARCL") or ((streams2extract[i].format == "AUTO") and (auto_era == 'SC'))) and ((wavella) or (wavello)):

       #print "Try to read data"

       #Loop over components:
       files2read = []
       for j in range(len(complist)):

           #Setup filename of possible local file: warchive/eventID.net.sta.loc.cha+comp.Format
           files2read.append(wavearc + "/" + str(refevlist[weidx].evID2) + "." + lnet + "." + lstat + "." + lloc + "." + lcha[0] + complist[j] + ".SAC")

       #Now check if all components are available:
       for j in range(len(files2read)):
           #print "Test if file exists:",files2read[j]
           if(os.path.isfile(files2read[j])):
              notinlocal = False
           else:
              notinlocal = True
              break

       #If all data is available in local directory open and add them to stream:
       if(notinlocal == False):
          for j in range(len(files2read)):
              st = read(files2read[j])

              #Number of traces to plot:
              cnttrc += 1

              #Append trace to list:
              trlist.append(st[0].copy())
              tridxlist.append(traceindex(cnttrc,weidx,streams2extract[i].integ,streams2extract[i].filter,streams2extract[i].scale))

              #Print information:
              print streams2extract[i],refevlist[weidx],0

       
    #Get waveform via ARCLINK (not includded yet):
    if(notinlocal) and (wavello == False) and ((streams2extract[i].format == "SC_ARCL") or ((streams2extract[i].format == "AUTO") and (auto_era == 'SC'))): 
       #This is now done once at the beginning...
       #Define client (if not done earlier):
       #if(arcl_flg == 0):
          #arcl_client = Client(host=arcl_host, port=arcl_port, user=arcl_user)
          #Get entire inventory:
          #print ""
          #print "--> Load inventory from",arcl_host
          #inv = arcl_client.getInventory('*', '*', '*', '*')
          #arcl_flg = 1

       #Now get available channel-type + network-code + location-code:
       wstat,wcha,wnet,wloc = getstationfrominventory(inv,streams2extract[i].netstat,streams2extract[i].channellist)

       if(wloc  == "NOTEXISTING"):
          print "WARNING: Data not available in arclink:",streams2extract[i].netstat,streams2extract[i].channellist,cevID
          continue

       #Loop over channels to extract:
       for j in range(len(complist)):

           #Now extract:
           trc = []
           try:
              if(arcl_meth == 'arclink'):
                 trc = arcl_client.getWaveform(wnet, wstat, wloc, wcha+complist[j], UTCDateTime(refevlist[weidx].timestamp) + arcl_winst, UTCDateTime(refevlist[weidx].timestamp) + arcl_winen)
              if(arcl_meth == 'FDSN-WS'):
                 trc = arcl_client.get_waveforms(wnet, wstat, wloc, wcha+complist[j], UTCDateTime(refevlist[weidx].timestamp) + arcl_winst, UTCDateTime(refevlist[weidx].timestamp) + arcl_winen)
           except:
              #ARCLINK-request failed, return empty stream-object
              #print "obspy.arclink.client.ArcLinkException: No data available"
              trc = [] 

           #Check if trace could be extracted:
           if(len(trc) == 0):
              print "WARNING: Data not available in arclink:",wnet,wstat,wloc,wcha+complist[j],cevID
              continue

           #Number of traces to plot:
           cnttrc += 1

           #Store trace to local disk:
           if(wavesto):
              trc[0].write(wavearc + "/" + str(refevlist[weidx].evID2) + "." + lnet + "." + lstat + "." + lloc + "." + wcha+complist[j] + ".SAC", format="SAC")

           trlist.append(trc[0])
           tridxlist.append(traceindex(cnttrc,weidx,streams2extract[i].integ,streams2extract[i].filter,streams2extract[i].scale))

           #Print information:
           print streams2extract[i],refevlist[weidx]
       
    #Store current eventID
    pevID = cevID

########################################################################################################################################
#Now plot traces
if(EnergyAnalysis):
   #Loop over stations to analyze:
   for i in range(len(EnergyList)):
       #Get the channels for current station:
       clist =  getcomponents(EnergyList[i],trlist)

       #Loop over channels:
       for l in range(len(clist)):
           filter_trlist,filter_tridxlist,filter_cnttrc = filtertrlist(trlist,tridxlist,EnergyList[i],clist[l])
           ofilefig_new = ofilefig + "_" + EnergyList[i] + "_" + clist[l]
           print "Plot: Station-Channel",EnergyList[i],clist[l]
           plotstreams(refevlist,filter_trlist,filter_tridxlist,filter_cnttrc,alignflg,ytype,epimin,epimax,ofilefig_new,asciidump,dbgmode)

else:
   plotstreams(refevlist,trlist,tridxlist,cnttrc,alignflg,ytype,epimin,epimax,ofilefig,asciidump,dbgmode)


