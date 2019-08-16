#!/usr/bin/python3

from datetime import datetime,timedelta
from netCDF4 import Dataset
import numpy as np
import gdal
import os
import sys
import re

#setting
bulan = [['Jan', '01'], ['Feb', '02'], ['Mar', '03'], ['Apr', '04'],
['May', '05'], ['Jun', '06'], ['Jul', '07'], ['Aug', '08'],
['Sep', '09'], ['Oct', '10'], ['Nov', '11'], ['Dec', '12']]
searchradius = 80
resol = 0.05
startlon = 70.
endlon = 160.
startlat = -20.
endlat = 70.
#-6.127171, 106.652756
targetlat = LATVAL
targetlon = LONVAL

#spasial, start from bottom left!!!
latitude = np.arange(startlat+resol,endlat+resol/2,resol)
longitude = np.arange(startlon+resol,endlon+resol/2,resol)

def autoestimator(ir):
	return 1.1183 * 10**(11.) * np.exp(-3.6382 * 10**(-2.) * ir**(1.2))

def findidx(lat, lon, matlat, matlon):
	for k in range(len(matlat)):
		if lat < matlat[k]:
			idxlat1 = abs(lat - matlat[k])
			idxlat2 = abs(lat - matlat[k-1])
			latidx = k if idxlat1 < idxlat2 else k-1
			break
	for k in range(len(matlon)):
		if lon < matlon[k]:
			idxlon1 = abs(lon - matlon[k])
			idxlon2 = abs(lon - matlon[k-1])
			lonidx = k if idxlon1 < idxlon2 else k-1
			break
	return latidx, lonidx

def findidxunstructured(lat, lon, matlat2d, matlon2d):
	newmatunstruct = np.zeros(np.shape(matlat2d))
	for k in range(len(newmatunstruct)):
		for l in range(len(newmatunstruct[0])):
			distlat = lat - matlat2d[k,l]
			distlon = lon - matlon2d[k,l]
			newmatunstruct[k,l] = np.sqrt(distlat**2. + distlat**2.)
	latidx, lonidx = np.unravel_index(newmatunstruct.argmin(), newmatunstruct.shape)
	return latidx, lonidx

def minimizemat(latsearch, lonsearch, matlat, matlon, data):
	idxlat, idxlon = findidx(latsearch, lonsearch, matlat, matlon)
	newlatmat = matlat[idxlat-searchradius:idxlat+searchradius]
	newlonmat = matlon[idxlon-searchradius:idxlon+searchradius]
	newdata = data[idxlat-searchradius:idxlat+searchradius, idxlon-searchradius:idxlon+searchradius]
	return newlatmat, newlonmat, newdata

def parallax(height,lat,lon):
	satheight,satlat,satlon=35800,0,140.7

	# -- varius earth radius information
	radius_eq = 6378.077
	dheight = satheight
	radius_pole = 6356.577
	radius_ratio = radius_eq/radius_pole
	mean_radius = 0.5*(radius_eq+radius_pole)
	zdiff = 0.0

	# -- angle conversion to radians
	asatlat = satlat * np.pi/180.0
	asatlon = satlon * np.pi/180.0
	alat = lat*np.pi/180.0
	alon = lon*np.pi/180.0
	# -- cartesian coordinates for the satellite
	#	satlat_geod is the geodetic satellite latitude
	satlat_geod = np.arctan(np.tan(asatlat)*radius_ratio**2)
	xsat = dheight * np.cos(satlat_geod) * np.sin(asatlon)
	ysat = dheight * np.sin(satlat_geod)
	zsat = dheight * np.cos(satlat_geod) * np.cos(asatlon)
	# -- cartesian coordinates of the surface point
	alat_geod = np.arctan(np.tan(alat)*radius_ratio**2)
	radius_surf = radius_eq/np.sqrt(np.cos(alat_geod)**2 + radius_ratio**2 * np.sin(alat_geod)**2)
	xsurf = radius_surf * np.cos(alat_geod) * np.sin(alon)
	ysurf = radius_surf * np.sin(alat_geod)
	zsurf = radius_surf * np.cos(alat_geod) * np.cos(alon)
	# -- compute new radius ratio depending on height
	radius_ratio_local = ((radius_eq+height)/(radius_pole+height))**2
	# -- Satellite minus surface location
	xdiff = xsat - xsurf
	ydiff = ysat - ysurf
	zdiff = zsat - zsurf
	# -- compute local zenith angle
	xfact = np.sqrt(xdiff**2 + ydiff**2 + zdiff**2)
	zen = (xdiff*xsurf+ydiff*ysurf+zdiff*zsurf)/(mean_radius*xfact)
	zen = np.arccos(zen)
	zen = zen*180.0/np.pi
	# -- equation to solve for the line of sight at height Z
	e1 = xdiff**2 + radius_ratio_local*ydiff**2 + zdiff**2
	e2 = 2.0 * (xsurf*xdiff + radius_ratio_local*ysurf*ydiff + zsurf*zdiff)
	e3 = xsurf**2 + zsurf**2 + radius_ratio_local*ysurf**2 - (radius_eq+height)**2
	corr = (np.sqrt(e2**2 - 4.0*e1*e3) - e2)/2.0/e1
	#	 corrected surface coordinates
	xcorr = xsurf + corr*xdiff
	ycorr = ysurf + corr*ydiff
	zcorr = zsurf + corr*zdiff
	#	 convert back to latitude and longitude
	latcorr = np.arctan(ycorr/np.sqrt(xcorr**2 + zcorr**2))
	latcorr = np.arctan(np.tan(latcorr)/radius_ratio**2) * 180.0/np.pi
	loncorr = np.arctan2(xcorr,zcorr) * 180.0/np.pi

	return latcorr,loncorr

def getheight(IR):
	tsurf=25+273.15
	lapserate=5.95
	height=(tsurf-IR)/5.95
	return height

def corrparallaxrr(lats,lons,rr,IR):
	rrcorrected=zeros(shape(rr))	
	for i in range (len(rr)):
		for j in range(len(rr[0])):
			if rr[i,j]>0.2: #kondisi minimum dianggap rain 
				height = getheight(IR[i,j])
				latcorr,loncorr = parallax(height,lats[i],lons[j])
				y = int((latcorr-lats[0])/(lats[1]-lats[0]))#Posisi lintang baru
				x = int((loncorr-lons[0])/(lons[1]-lons[0]))#posisi bujur baru
				try: #jika ada yang error langsung dilewati saja 
					rrcorrected[y,x]=rr[i,j]
				except:
					pass
	return rrcorrected

def recreatelatlon(lat, lon, matpgm):
	newcor = np.zeros((len(lat), len(lon), 2))
	for k in range(len(lat)):
		for l in range(len(lon)):
			height = getheight(matpgm[k,l])
			latcor1, loncor1 = parallax(height,lat[k],lon[l])
			newcor[k,l,0] = latcor1
			newcor[k,l,1] = loncor1
	return newcor

if len(sys.argv) < 2:
	print('please provide the month name')
	exit()
else:
	none = True
	for i in range(len(bulan)):
		if sys.argv[1] == bulan[i][0]:
			bulanexecangka = bulan[i][1]
			bulanexec = bulan[i][0]
			none = False
			break
	if none == True:
		print('month name not match')
		print('the option is Jan, Feb, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec')
		exit()


rainall = []
raincorall = []
dataall = []
for hari in range(1,32):
	if hari < 10:
		haristr = '0'+str(hari)
	else:
		haristr = str(hari)
	for jam in range(24):
		if jam < 10:
			jamstr = '0'+str(jam)
		else:
			jamstr = str(jam)
		pgmfile = 'GAME/2018/%s/IR1/HMW818%s%s%sIR1.pgm'%(bulanexec, bulanexecangka, haristr, jamstr)
		datfile = 'CAL/2018/%s/HMW818%s%s%sCAL.dat'%(bulanexec, bulanexecangka, haristr, jamstr)
		print('processing %s...'%(pgmfile))
		try:
			ds = gdal.Open(pgmfile)
			irdat = np.array(ds.GetRasterBand(1).ReadAsArray())
			opendat = open(datfile)
		except AttributeError:
			print('broken or not found pgm file')
			continue
		except FileNotFoundError:
			print('broken or not found dat file')
			continue
		datall = opendat.read()
		datall = datall.split('\n')
		irdatsuhu = np.zeros(np.shape(irdat))
		#create cor matrix IR1
		corIR1 = np.zeros((256))
		for i in range(256):
			if i < 10:
				istr = '  '+str(i)
			elif i < 100:
				istr = ' '+str(i)
			else:
				istr = str(i)
			for j in range(len(datall)):
				if 'IR1 Temperature of %s pixval:'%(istr) in datall[j]:
					corIR1[i] = float(datall[j][len('IR1 Temperature of %s pixval:'%(istr)):])
		for i in range(len(irdat)):
			for j in range(len(irdat[0])):
				idx = irdat[i,j]
				irdatsuhu[i,j] = corIR1[idx]
		minilatmat, minilonmat, miniirdatsuhu = minimizemat(targetlat, targetlon, latitude, longitude, irdatsuhu)
		#without parallax correction
		rain = autoestimator(miniirdatsuhu[searchradius, searchradius])
		rainall.append(rain)
#		print(rain)
		
		#with parallax correction
		newcor = recreatelatlon(minilatmat, minilonmat, miniirdatsuhu)
		newminilat = newcor[:,:,0]
		newminilon = newcor[:,:,1]
#		print(np.min(newminilat), np.max(newminilat))
#		print(np.min(newminilon), np.max(newminilon))
#		print(np.min(minilatmat), np.max(minilatmat))
#		print(np.min(minilonmat), np.max(minilonmat))
#		print(minilonmat)
		rllatidx, rllonidx = findidxunstructured(targetlat, targetlon, newminilat, newminilon)
		raincor = autoestimator(miniirdatsuhu[rllatidx, rllonidx])
		raincorall.append(raincor)
		dataall.append(pgmfile)
#		print(raincor)

outdata = open(bulanexec+'.csv', 'w')
outdata.write('data, rain,rain correction\n')
for i in range(len(rainall)):
	if i == len(rainall) - 1:
		closing = ''
	else:
		closing = '\n'
	outdata.write('%s,%s,%s%s'%(dataall[i], rainall[i], raincorall[i], closing))
