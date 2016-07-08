#!/Users/timothy.lee/anaconda/envs/python2.7/bin/python

import re
import math
import time
import numpy as np
import cPickle
from obspy.core import UTCDateTime
from obspy.core.util import gps2DistAzimuth
#import pg
import os
import shutil
import glob
import string
import optparse

global VERY_SMALL_DOUBLE, SMALL_DOUBLE
VERY_SMALL_DOUBLE = 1.0e-30
SMALL_DOUBLE = 1.0e-8


###########################################################################################
# Classes:
#------------------------------------------------------------------------------------------
class hypo(object):

      def __init__(self,yy,mm,dd,hh,mi,ss,timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mtype,mnobs,merr,mmeth,rms,gap,mdist,nobs,etype,lqual,chx,chy,agency,evID1,evID2,evPID,orID1,methodID,earthmodelID,author,region,fstatus,DDx,DDy,DDz,DDxe,DDye,DDze,DDNCCP,DDNCCS,DDNCTP,DDNCTS,DDCCRMS,DDCTRMS,DDCID):
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
          self.DDx     = float(DDx)
          self.DDy     = float(DDy)
          self.DDz     = float(DDz)
          self.DDxe    = float(DDxe)
          self.DDye    = float(DDye)
          self.DDze    = float(DDze)
          self.DDNCCP  = int(DDNCCP)
          self.DDNCCS  = int(DDNCCS)
          self.DDNCTP  = int(DDNCTP)
          self.DDNCTS  = int(DDNCTS)
          self.DDCCRMS = float(DDCCRMS)
          self.DDCTRMS = float(DDCTRMS)
          self.DDCID   = int(DDCID)

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

###########################################################################################
#Subroutines/functions:

#------------------------------------------------------------------------------------------
def loadhypoDDevents(ifile,searchmethod,IDList,dbgmode):

    list = []
    list_blasts = []

    #Open
    fp = open(ifile,"r")

    for line in fp:


        #Test if line is empty:
        if(len(line.split())<24):
           print "Line not considered:", line
        else:
           #Test if line is comment:
           #t=line.split()
           t = string.split(line,None,-1)

           if(t[0].startswith("#") == False):
              lon=float(t[2])
              lat=float(t[1])
              dep=float(t[3])
              mag=float(t[16])
              yy=int(t[10])
              mm=int(t[11])
              dd=int(t[12])
              hh=int(t[13])
              mi=int(t[14])
              ss=float(t[15])
              ID2=long(t[0])

              #Get swiss coordinates:
              chx,chy = celleb(lon,lat)
              if(chy < 62.0) or (chy > 302.0) or (chx < 480.0) or (chx > 847.5):
                chy = 999
                chx = 999

              #hypoDD specific:
              DDx     = float(t[4])
              DDy     = float(t[5])
              DDz     = float(t[6])
              DDxe    = float(t[7])
              DDye    = float(t[8])
              DDze    = float(t[9])
              DDNCCP  = int(t[17])
              DDNCCS  = int(t[18])
              DDNCTP  = int(t[19])
              DDNCTS  = int(t[20])
              DDCCRMS = float(t[21])
              DDCTRMS = float(t[22])
              DDCID   = int(t[23])

              #Unknown:
              mat    = '?'
              loc    = '?'
              rms=-9
              gap=-9
              mdi=-9
              nob=-9
              ety='T'                  #Assume events are earthquakes...
              lqa='?'
              ID1='?'
              agy    = "SED_DD"
              region = "None"
              author = "DD"
              ID3    = "None"          #Origin publicID (string)
              ID4    = "None"          #Event  publicID (string)
              fstat  = "DD"
              laterr = -9
              lonerr = -9
              deperr = -9
              mnobs  = -9
              merr   = -9
              modlID = '?'
              magmeth = '?'

              #Check if ID list was provided:
              if(len(IDList) > 0):
                 excludeevent = True
                 for k in range(len(IDList)):
                     if(searchmethod == "ID-DD") and (IDList[k]==ID2):
                        excludeevent = False
                        break
                     if(searchmethod == "ID-KP") and (IDList[k]==ID1):
                        excludeevent = False
                        break
                 #Final check if event should be appended:
                 if(excludeevent):
                    continue

              #Check second:
              if(ss == 60.0):
                 print 'Second == 60.0 --> substract 0.00001'
                 ss = ss - 0.00001

              wdate = UTCDateTime("%04d-%02d-%02d %02d:%02d:%09.6f" % (yy,mm,dd,hh,mi,ss)) #Problems with seconds close to 60.0 e.g.: 59.999995 -> 60.0 Fix: increase precision

              list.append(hypo(yy,mm,dd,hh,mi,ss,wdate.timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mat,mnobs,merr,magmeth,rms,gap,mdi,nob,ety,lqa,chx,chy,agy,ID1,ID2,ID4,ID3,loc,modlID,author,region,fstat,DDx,DDy,DDz,DDxe,DDye,DDze,DDNCCP,DDNCCS,DDNCTP,DDNCTS,DDCCRMS,DDCTRMS,DDCID))
              if(dbgmode == 1):
                 print hypo(yy,mm,dd,hh,mi,ss,wdate.timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mat,mnobs,merr,magmeth,rms,gap,mdi,nob,ety,lqa,chx,chy,agy,ID1,ID2,ID4,ID3,loc,modlID,author,region,fstat,DDx,DDy,DDz,DDxe,DDye,DDze,DDNCCP,DDNCCS,DDNCTP,DDNCTS,DDCCRMS,DDCTRMS,DDCID)

    #Close
    fp.close()

    return list,list_blasts

#------------------------------------------------------------------------------------------
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False
#------------------------------------------------------------------------------------------
def loadpolygon(ifile):
    polylist = []

    #Open
    fp = open(ifile,"r")

    for line in fp:

        #Split string:
        t = string.split(line,None,-1)

        #Test if column 1 and 2 are numbers
        if(isfloat(t[0])) and (isfloat(t[1])):
           polylist.append((float(t[0]),float(t[1])))
           print "Polygon-Line     valid:",line
        else:
           print "Polygon-Line not valid:",line       
 
    #Close
    fp.close()
   
    print "Polygon: ",polylist 

    return polylist
#------------------------------------------------------------------------------------------
def loadIDlist(ifile,searchmethod):

    IDlist = []

    #Open
    fp = open(ifile,"r")

    for line in fp:
        #Test if line is empty:
        if(len(line.split())!=1):
           print "Line not considered:", line
        else:
           #Test if line is comment:
           #t=line.split()
           t = string.split(line,None,-1)

           if(t[0].startswith("#") == False):
              if(searchmethod == "ID-DD"):
                 IDlist.append(long(t[0]))
              if(searchmethod == "ID-KP"):
                 IDlist.append(str(t[0]))

    #Close
    fp.close()

    return IDlist
#------------------------------------------------------------------------------------------
def loadrefevents(ifile,searchmethod,IDList,dbgmode):

    list = []
    list_blasts = []

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

              #hypoDD specific (not available):
              DDx     = -9
              DDy     = -9
              DDz     = -9
              DDxe    = -9
              DDye    = -9
              DDze    = -9
              DDNCCP  = -9
              DDNCCS  = -9
              DDNCTP  = -9
              DDNCTS  = -9
              DDCCRMS = -9
              DDCTRMS = -9
              DDCID   = -9

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

              #Check if ID list was provided:
              if(len(IDList) > 0):
                 excludeevent = True
                 for k in range(len(IDList)):
                     if(searchmethod == "ID-DD") and (IDList[k]==ID2):
                        excludeevent = False
                        break
                     if(searchmethod == "ID-KP") and (IDList[k]==ID1):
                        excludeevent = False
                        break
                 #Final check if event should be appended:
                 if(excludeevent):
                    continue


              #Check if format is correct:
              if(ss < 60.0) and (mi < 60) and (hh < 24) and (dd < 32) and (mm < 13):
                 wdate = UTCDateTime("%04d-%02d-%02d %02d:%02d:%06.3f" % (yy,mm,dd,hh,mi,ss)) #Problems with seconds close to 60.0 e.g.: 59.99995 -> 60.0 Fix: increase precision
                 #wdate = UTCDateTime("%04d-%02d-%02d %02d:%02d:%09.6f" % (yy,mm,dd,hh,mi,ss))

                 #print ID1
                 list.append(hypo(yy,mm,dd,hh,mi,ss,wdate.timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mat,mnobs,merr,magmeth,rms,gap,mdi,nob,ety,lqa,chx,chy,agy,ID1,ID2,ID4,ID3,loc,modlID,author,region,fstat,DDx,DDy,DDz,DDxe,DDye,DDze,DDNCCP,DDNCCS,DDNCTP,DDNCTS,DDCCRMS,DDCTRMS,DDCID))

                 #Check if event is blast, if yes, add it to blast list:
                 if(ety == "E"):
                    list_blasts.append(hypo(yy,mm,dd,hh,mi,ss,wdate.timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mat,mnobs,merr,magmeth,rms,gap,mdi,nob,ety,lqa,chx,chy,agy,ID1,ID2,ID4,ID3,loc,modlID,author,region,fstat,DDx,DDy,DDz,DDxe,DDye,DDze,DDNCCP,DDNCCS,DDNCTP,DDNCTS,DDCCRMS,DDCTRMS,DDCID))

                 if(dbgmode == 1):
                    print hypo(yy,mm,dd,hh,mi,ss,wdate.timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mat,mnobs,merr,magmeth,rms,gap,mdi,nob,ety,lqa,chx,chy,agy,ID1,ID2,ID4,ID3,loc,modlID,author,region,fstat,DDx,DDy,DDz,DDxe,DDye,DDze,DDNCCP,DDNCCS,DDNCTP,DDNCTS,DDCCRMS,DDCTRMS,DDCID)
              else:
                 print "Warning: Corrupt origin time ",datetime2(yy,mm,dd,hh,mi,ss),ID1
           else:
              print "Line not considered:", line

    #Close
    fp.close()

    return list,list_blasts
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
# Determine if a point is inside a given polygon or not
# Polygon is a list of (x,y) pairs. This function
# returns True or False.  The algorithm is called
# the "Ray Casting Method".
# http://geospatialpython.com/2011/01/point-in-polygon.html

def point_in_poly(x,y,poly):

    n = len(poly)
    inside = False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside
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
def blast_check(kblist,slist,bcepimax,bcdayst,bcdayen,bcminct,fplog05,bcYYst,bcYYen):

    pblist = []


    print ""
    print "Blast Detector:"

    lS=len(slist)
    lB=len(kblist)
    print "Number of known blast-events: ",lB
    print ""
    print "Potential blast-events in extracted event list:"

    fplog05.write("Number of known blast-events: %8d\n" % (lB))
    fplog05.write("\n")
    fplog05.write("Potential blast-events in extracted event list:\n")

    #Loop over potential (known blast events):
    for i in range(lS):
        otim = UTCDateTime(slist[i].timestamp)      #Still missing, convert to local time...
        cntb = 0
        if(slist[i].etype != "E") and (slist[i].hh >= bcdayst) and (slist[i].hh <= bcdayen) and (otim.weekday < 6) and (otim.year >= bcYYst) and (otim.year < bcYYen):
           #Calculate distance to known blast and count number of blast events within radius:
           for j in range(lB):
               epid=gps2DistAzimuth(slist[i].lat,slist[i].lon,kblist[j].lat,kblist[j].lon)
               if(epid[0]/1000.0 <= bcepimax):
                  cntb += 1

           if(cntb >= bcminct):
              pblist.append(slist[i])
              print "#%03d | %s" % (cntb,slist[i])
              fplog05.write("#%03d | %s\n" % (cntb,slist[i]))

    return pblist
#------------------------------------------------------------------------------------------



###########################################################################################
# Main code:
###########################################################################################

#Get command line arguments:
oparser = optparse.OptionParser()
oparser.add_option('-f', '--file',               action='store', dest='extractlist', help='Reference list of hypocenters (SC3DB2KP.py format); default: /Users/tdiehl/lib/MANULOC_List_1984_2012_ECOS_short.log.latest')
oparser.add_option('-i', '--inputformat',        action='store', dest='informat',    help='Input format of hypocenter list (supported: SC3DB2KP, hypoDD); default: SC3DB2KP')
oparser.add_option('-o', '--output',             action='store', dest='outfi',       help='Output file with extracted events; default: episelfromlist.out')
oparser.add_option('-e', '--exportformat',       action='store', dest='exformat',    help='Export/Output format with extracted events (supported: SC3DB2KP, hypoDD); default: SC3DB2KP')
oparser.add_option('-t', '--eventtype',          action='store', dest='event2ex',    help='Eventtype to extract: T=Earthquakes; E=QuarryBlasts; T+E=Eqks+Blasts (default); A=All (including landslides etc.)')
oparser.add_option('-m', '--searchmethod',       action='store', dest='searchmethod',help='Searchmethod: "cir": Selects events within radius around point (see -c); "box" (default): Selects events within rectangular box (see -b); "ID-DD": Selects events which match hypoDD-IDs in a given list (see -l); "ID-KP": Selects events which match string-ID in a given list (see -l); "pol": Selects events within a polygon given in a file (see -p)')
oparser.add_option('-b', '--boxrange',           action='store', dest='boxrange',    help='if "-m box" this parameter defines the box by latmin/latmax/lonmin/lonmax; default: 45.0/49.0/5.0/11.0')
oparser.add_option('-c', '--circlerange',        action='store', dest='circlerange', help='if "-m cir" this parameter defines the circle by clat/clon/radius with radius in km; if clat > 90, clat is assumed to be DD-eventID, and lat/lon of corresponding event are used for clat/clon')
oparser.add_option('-l', '--IDlist',             action='store', dest='IDlist',      help='if "-m ID-DD" or "-m ID-KP" this parameter defines the file name of a event-ID list')
oparser.add_option('-p', '--polygonfile',        action='store', dest='Polyfi',      help='if "-m pol" this parameter defines the file name with the polygon definition lon lat pairs per line, not closed')
oparser.add_option('-M', '--MagnitudeRange',     action='store', dest='magrange',    help='Magnitude range: MagMin/MagMax; default -10/+10')
oparser.add_option('-D', '--DepthRange',         action='store', dest='deprange',    help='Depth range: DepMin/DepMax; default -999.9/+999.9')
oparser.add_option('-T', '--TimePeriod',         action='store', dest='TimePeriod',  help='Time Periode: YYYY-MM-DDTHH:MM:SS/YYYY-MM-DDTHH:MM:SS; default 1984-01-01T00:00:00/2099-12-31T23:59:00')
oparser.add_option('-B', '--Blast-Check',        action='store', dest='blastchk',    help='For every earthquake the entire catalog is searched for near-by blasts, argument is R/N/YYS/YYE where R: search radius in km (e.g. 4); N: Minimum number of known blast events within radius, YYS: Start-year for earthquake search, YYE: End-year for earthquake search')

(options, arg) = oparser.parse_args()

#Load a event list (created by SC3DB2KP.py) and select events based on different search criteria

#Input file:

#Check if blast check is required:
if(options.blastchk != None):
   blastchk = options.blastchk
else:
   blastchk = '-4.0/10/1984/2099'
sp = string.split(blastchk,'/',-1)              #Split string:
bcepimax = float(sp[0])                         #Radius to known blast events (km), if negative skip test
bcminct  = int(sp[1])                           #Minimum number of known blast events within radius
bcYYst   = int(sp[2])
bcYYen   = int(sp[3])
#Fixed:
bcdayst  = 6      #Start of working day (hour)
bcdayen  = 20     #End   of working day (hour)

#Check if event list is included in list of input arguments:
if(options.extractlist != None):
   extractlist = options.extractlist
else:
   extractlist='/Users/tdiehl/lib/MANULOC_List_1984_2012_ECOS_short.log.latest'
print "--> Input     File       :",extractlist

#Check if input format is specified:
if(options.informat != None):
   informat = options.informat
else:
   informat='SC3DB2KP'
print "--> Input     Format     :",informat

#Output file:
#Check if output file is included in list of input arguments:
if(options.outfi != None):
   outfi = options.outfi
else:
   outfi='episelfromlist.out'
print "--> Output    File       :",outfi

#Check if export/output format is specified:
if(options.exformat != None):
   exformat = options.exformat
else:
   exformat='SC3DB2KP'
print "--> Output/Ex Format     :",exformat

#Define event types to extract:
#T  : Tectonic (earthquake)
#E  : Explosion
#T+E: Tectonic + Explosion
#A  : All (also landslides, avalanches, etc.)
#Check if eventtype is included in list of input arguments:
if(options.event2ex != None):
   event2ex = options.event2ex
else:
   event2ex='T+E'
print "--> Event     Type(s)    :",event2ex

#Search mode:
# box : rectangular box defined by latmin,latmax,lonmin,lonmax
# cir : circular search defined by clat, clon, radius
#       radius is in km
#       if clat > 90 clat is assumed to be DD-eventID and clat,clon of corresponding event are used
# ID-DD : Select events from a list with hypoDD ID's
# ID-KP : Select events from a list with string ID (like SED's KP)
# pol   : Select events within a polygon
#Check if option is included in list of input arguments:
if(options.searchmethod != None):
   searchmethod = options.searchmethod
else:
   #searchmethod = 'cir'
   searchmethod = 'box'

#List Definition:
if(options.IDlist != None):
   IDlistfile = options.IDlist
else:
   IDlistfile = 'events_ID_DD'

#Polygon file:
if(options.Polyfi != None):
   Polyfile = options.Polyfi
else:
   Polyfile = 'Test.polygon'

#Box definition:
if(options.boxrange != None):
   boxrange = options.boxrange
else:
   boxrange = '45.0/49.0/5.0/11.0'
sp = string.split(boxrange,'/',-1)              #Split string:
latmin = float(sp[0])
latmax = float(sp[1])
lonmin = float(sp[2])
lonmax = float(sp[3])

#Circle defintion:
if(options.circlerange != None):
   circlerange = options.circlerange
else:
   circlerange = '201101011/0.0/2.5'            #Reference event is DD-evID = 201101011
sp = string.split(circlerange,'/',-1)           #Split string
ts1 = float(sp[0])
if(ts1 > 90.0):
   clat = long(sp[0])                           #clat is event ID -> long instead of float 
else:
   clat = ts1                                   #clat is latitude -> float
clon = float(sp[1])
radius = float(sp[2])

#Info to screen:
if(searchmethod == 'cir'):
   print "--> Search    Method     :",searchmethod,circlerange
else:
   print "--> Search    Method     :",searchmethod,boxrange

#Magnitude range:
if(options.magrange != None):
   magrange = options.magrange
else:
   magrange = '-10/+10'
print "--> Magnitude Range      :",magrange
sp = string.split(magrange,'/',-1)
magmin = float(sp[0])
magmax = float(sp[1])

#Depth range:
if(options.deprange != None):
   deprange = options.deprange
else:
   deprange = '-999.9/+999.9'
print "--> Depth     Range      :",deprange
sp = string.split(deprange,'/',-1)
depmin = float(sp[0])
depmax = float(sp[1])

#Time period:
if(options.TimePeriod != None):
   TimePeriod = options.TimePeriod
else:
   TimePeriod = '1984-01-01T00:00:00/2099-12-31T23:59:00'
print "--> Time      Period     :",TimePeriod
sp = string.split(TimePeriod,'/',-1)
timest = string.strip(sp[0])
timeen = string.strip(sp[1])

#####################################################################
#Other definition (fixed at the moment)

#List with station coordinates:
statlist = "/Users/tdiehl/lib/sed_stations.GSE_SED.alias"

#SC3 origins are preferred over MANULOC (reference) after:
SC3overKPdate = '2012-10-03 00:00:00'

#Debug-Mode (0 or 1):
dbgmode = 0

############
#start code:

stts = UTCDateTime(timest).timestamp
ents = UTCDateTime(timeen).timestamp

#Load polygon (if requested):
if(searchmethod == "pol"):
   geopolylist = loadpolygon(Polyfile)

#Load list with event IDs to extract:
if(searchmethod == "ID-DD") or (searchmethod == "ID-KP"):
   IDList = loadIDlist(IDlistfile,searchmethod)
else:
   IDList = []

#Load reference events:
if(informat == "SC3DB2KP"):
   print ""
   print "--> Load SC3DB2KP-List from: ",extractlist
   refevlist,refblasts = loadrefevents(extractlist,searchmethod,IDList,dbgmode)

if(informat == "hypoDD"):
   print ""
   print "--> Load hypoDD-List   from: ",extractlist
   refevlist,refblasts = loadhypoDDevents(extractlist,searchmethod,IDList,dbgmode)

#Load station coordinates:
print "--> Load Station-List  from: ",statlist
stationlist = loadstations(statlist,dbgmode)

#get reference event:
refev = 0
if(clat > 90.0):
   for i in range(len(refevlist)):
       if(refevlist[i].evID2 == clat):
          refev = i
          
          break
if(refev > 0):
   clat = refevlist[refev].lat
   clon = refevlist[refev].lon
   print ""
   print "Reference event:"
   print refevlist[refev] 

list = []

#yy,mm,dd,hh,mi,ss,timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mtype,mnobs,merr,mmeth,rms,gap,mdist,nobs,etype,lqual,chx,chy,agency,evID1,evID2,evPID,orID1,methodID,earthmodelID,author,region,fstatus
#Loop over events in list:
for i in range(len(refevlist)):

    #Check time:
    if(refevlist[i].timestamp < stts) or (refevlist[i].timestamp > ents):
       continue

    #Check event type:
    if(event2ex == "T") and (refevlist[i].etype != "T"):
       continue
    if(event2ex == "E") and (refevlist[i].etype != "E"):
       continue
    if(event2ex == "T+E") and ((refevlist[i].etype != "E") and (refevlist[i].etype != "T")):
       continue

    #Check depth:
    if(refevlist[i].dep < depmin) or (refevlist[i].dep > depmax):
       continue

    #Check magnitude:
    if(refevlist[i].mag < magmin) or (refevlist[i].mag > magmax):
       continue

    #Check epicenter:
    if(searchmethod == 'box') and ((refevlist[i].lat > latmax) or (refevlist[i].lat < latmin) or (refevlist[i].lon > lonmax) or (refevlist[i].lon < lonmin)):
       continue

    #Check epicenter:
    if(searchmethod == 'cir'):
       epid=gps2DistAzimuth(refevlist[i].lat,refevlist[i].lon,clat,clon)     
       if(epid[0]/1000.0 > radius): 
          continue

    #Check epicenter:
    if(searchmethod == 'pol'):
       if(point_in_poly(refevlist[i].lon,refevlist[i].lat,geopolylist) == False):
          continue

    #All tests passed: 
    list.append(refevlist[i])



#Loop over extracted events:
#open output file:
fp1 = open(outfi,"w")
print ""
print "Events matching search criteria and selected:"
print "============================================="
for i in range(len(list)):

    #Write SC3DB2KP format:
    if(exformat == "SC3DB2KP"):
       #Simple list to screen:
       #print list[i]
       #Detailed (full) list to screen:
       print "%8.3f %7.3f %5.1f %4.1f %-2s %04d %02d %02d %02d %02d %05.2f %-4s %6.3f %3d %5.1f %3d %1s %1s %-20s %9d | %4.0f | %4.0f | %-15s | %-25s | %s | %s" % (list[i].lon,list[i].lat,list[i].dep,list[i].mag,list[i].mtype[0:2],list[i].yy,list[i].mm,list[i].dd,list[i].hh,list[i].mi,list[i].ss,list[i].methodID,list[i].rms,list[i].gap,list[i].mdist,list[i].nobs,list[i].etype,list[i].lqual,list[i].evID1,list[i].evID2,list[i].chx,list[i].chy,list[i].author,list[i].region,list[i].orID1,list[i].evPID)
       #Detailed (full) list to file:
       fp1.write("%8.3f %7.3f %5.1f %4.1f %-2s %04d %02d %02d %02d %02d %05.2f %-4s %6.3f %3d %5.1f %3d %1s %1s %-20s %9d | %4.0f | %4.0f | %-15s | %-25s | %s | %s\n" % (list[i].lon,list[i].lat,list[i].dep,list[i].mag,list[i].mtype[0:2],list[i].yy,list[i].mm,list[i].dd,list[i].hh,list[i].mi,list[i].ss,list[i].methodID,list[i].rms,list[i].gap,list[i].mdist,list[i].nobs,list[i].etype,list[i].lqual,list[i].evID1,list[i].evID2,list[i].chx,list[i].chy,list[i].author,list[i].region,list[i].orID1,list[i].evPID))

    #Write hypoDD output:
    if(exformat == "hypoDD"):
       #Detailed (full) list to screen:
       print "%9d %10.6f %11.6f %9.3f %10.1f %10.1f %10.1f %8.1f %8.1f %8.1f %04d %02d %02d %02d %02d %06.3f %5.2f %5d %5d %5d %5d %6.3f %6.3f %3d" % (list[i].evID2,list[i].lat,list[i].lon,list[i].dep,list[i].DDx,list[i].DDy,list[i].DDz,list[i].DDxe,list[i].DDye,list[i].DDze,list[i].yy,list[i].mm,list[i].dd,list[i].hh,list[i].mi,list[i].ss,list[i].mag,list[i].DDNCCP,list[i].DDNCCS,list[i].DDNCTP,list[i].DDNCTS,list[i].DDCCRMS,list[i].DDCTRMS,list[i].DDCID)
       #Detailed (full) list to file:
       fp1.write("%9d %10.6f %11.6f %9.3f %10.1f %10.1f %10.1f %8.1f %8.1f %8.1f %04d %02d %02d %02d %02d %06.3f %5.2f %5d %5d %5d %5d %6.3f %6.3f %3d\n" % (list[i].evID2,list[i].lat,list[i].lon,list[i].dep,list[i].DDx,list[i].DDy,list[i].DDz,list[i].DDxe,list[i].DDye,list[i].DDze,list[i].yy,list[i].mm,list[i].dd,list[i].hh,list[i].mi,list[i].ss,list[i].mag,list[i].DDNCCP,list[i].DDNCCS,list[i].DDNCTP,list[i].DDNCTS,list[i].DDCCRMS,list[i].DDCTRMS,list[i].DDCID))

#Close file:
fp1.close()

###########################################
#Now check for mis-identified blasts (Blasts labeled as Earthquakes):
#Blast-Check (checks non-blast events for location close to known blasts and day-to-night time):
if(bcepimax > 0):
   blstlog = open(outfi+".PotentialBlasts","w")
   blist = blast_check(refblasts,list,bcepimax,bcdayst,bcdayen,bcminct,blstlog,bcYYst,bcYYen)
   blstlog.close()


print ""
print "Extracted events exported to                            :",outfi
if(bcepimax > 0):
   print "Log-file of exported earthquakes which migth be blasts  :",outfi+".PotentialBlasts"
