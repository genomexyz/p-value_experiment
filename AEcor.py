#!/usr/bin/python3

from datetime import datetime,timedelta
from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap

def parallax(height,lat,lon):
	satheight,satlat,satlon=35800,0,140.7

	dpi = 3.14159265

	# -- varius earth radius information
	radius_eq = 6378.077
	dheight = satheight
	radius_pole = 6356.577
	radius_ratio = radius_eq/radius_pole
	mean_radius = 0.5*(radius_eq+radius_pole)
	zdiff = 0.0

	# -- angle conversion to radians
	asatlat = satlat * dpi/180.0
	asatlon = satlon * dpi/180.0
	alat = lat*dpi/180.0
	alon = lon*dpi/180.0
	# -- cartesian coordinates for the satellite
	#	satlat_geod is the geodetic satellite latitude
	satlat_geod = ma.atan(tan(asatlat)*radius_ratio**2)
	xsat = dheight * ma.cos(satlat_geod) * ma.sin(asatlon)
	ysat = dheight * ma.sin(satlat_geod)
	zsat = dheight * ma.cos(satlat_geod) * ma.cos(asatlon)
	# -- cartesian coordinates of the surface point
	alat_geod = ma.atan(ma.tan(alat)*radius_ratio**2)
	radius_surf = radius_eq/sqrt(cos(alat_geod)**2 + radius_ratio**2 * ma.sin(alat_geod)**2)
	xsurf = radius_surf * ma.cos(alat_geod) * ma.sin(alon)
	ysurf = radius_surf * ma.sin(alat_geod)
	zsurf = radius_surf * ma.cos(alat_geod) * ma.cos(alon)
	# -- compute new radius ratio depending on height
	radius_ratio_local = ((radius_eq+height)/(radius_pole+height))**2
	# -- Satellite minus surface location
	xdiff = xsat - xsurf
	ydiff = ysat - ysurf
	zdiff = zsat - zsurf
	# -- compute local zenith angle
	xfact = ma.sqrt(xdiff**2 + ydiff**2 + zdiff**2)
	zen = (xdiff*xsurf+ydiff*ysurf+zdiff*zsurf)/(mean_radius*xfact)
	zen = ma.acos(zen)
	zen = zen*180.0/dpi
	# -- equation to solve for the line of sight at height Z
	e1 = xdiff**2 + radius_ratio_local*ydiff**2 + zdiff**2
	e2 = 2.0 * (xsurf*xdiff + radius_ratio_local*ysurf*ydiff + zsurf*zdiff)
	e3 = xsurf**2 + zsurf**2 + radius_ratio_local*ysurf**2 - (radius_eq+height)**2
	corr = (ma.sqrt(e2**2 - 4.0*e1*e3) - e2)/2.0/e1
	#	 corrected surface coordinates
	xcorr = xsurf + corr*xdiff
	ycorr = ysurf + corr*ydiff
	zcorr = zsurf + corr*zdiff
	#	 convert back to latitude and longitude
	latcorr = ma.atan(ycorr/sqrt(xcorr**2 + zcorr**2))
	latcorr = ma.atan(ma.tan(latcorr)/radius_ratio**2) * 180.0/dpi
	loncorr = ma.atan2(xcorr,zcorr) * 180.0/dpi

	return latcorr,loncorr

ef getheight(IR):
	tsurf=25+273.15
	lapserate=5.95
	height=(tsurf-IR)/5.95
	return height

def corrparallaxrr(lats,lons,rr,IR):
	rrcorrected=zeros(shape(rr))	
	for i in range (len(rr)):
		for j in range(len(rr[0])):
			if rr[i,j]>0.2: #kondisi minimum dianggap hujan 
				height = getheight(IR[i,j])
				latcorr,loncorr = parallax(height,lats[i],lons[j])
				y = int((latcorr-lats[0])/(lats[1]-lats[0]))#Posisi lintang baru
				x = int((loncorr-lons[0])/(lons[1]-lons[0]))#posisi bujur baru
				try: #jika ada yang error langsung dilewati saja 
					rrcorrected[y,x]=rr[i,j]
				except:
					pass
	return rrcorrected
