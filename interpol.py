"""
class for preforming interpolation between grids defined by lats/lons either in netCDF4 or grib files
"""
#for array manipulation
import numpy as np
import numpy.ma as ma

#for interpolation
from scipy.spatial import cKDTree

#for IO
import pygrib
from gribapi import *
from datetime import date

#plot
import matplotlib.pyplot as plt

class KdInterpol:

	def __init__(self,source,target):
		
		def lon_lat_to_cartesian(lon,lat,R=1):
			"""
			calculates lon, lat coordinates of a point on a sphere with
			radius R
			"""
			
			lon_r=np.radians(lon)
			lat_r=np.radians(lat)

			x=R*np.cos(lat_r)*np.cos(lon_r)
			y=R*np.cos(lat_r)*np.sin(lon_r)
			z=R*np.sin(lat_r)
			return x,y,z
	
		"""
		read data from files
		"""
		self.source=source
		self.target=target

		self.D1=pygrib.open(self.source)
		self.latSOURCE,self.lonSOURCE=self.D1[8].latlons()
				
		self.D2=pygrib.open(self.target)
		self.latTARGET,self.lonTARGET=self.D2[1].latlons()	

		"""
		flatten grid coordinates into cartesian x,y,z
		"""
		self.xs,self.ys,self.zs=lon_lat_to_cartesian(self.lonSOURCE.flatten(),self.latSOURCE.flatten())
		self.xt,self.yt,self.zt=lon_lat_to_cartesian(self.lonTARGET.flatten(),self.latTARGET.flatten())

	def __call__(self):
		
		def write_to_file(variable):
			"""
			writing to grib file using low level C bindings
			"""
			
			OldTemplate=self.D1[8]		

			NewTemplate=self.D2[1]		
			
			grbout=open('test.grb','wb')

			OldTemplate.dataDate=20140728	
			OldTemplate.values=variable
			
			OldTemplate.latitudeOfFirstGridPoint=NewTemplate.latitudeOfFirstGridPoint
			OldTemplate.latitudeOfFirstGridPointInDegrees=NewTemplate.latitudeOfFirstGridPointInDegrees
			OldTemplate.longitudeOfFirstGridPoint=NewTemplate.longitudeOfFirstGridPoint
			OldTemplate.longitudeOfFirstGridPointInDegrees=NewTemplate.longitudeOfFirstGridPointInDegrees

			OldTemplate.latitudeOfLastGridPoint=NewTemplate.latitudeOfLastGridPoint
			OldTemplate.latitudeOfLastGridPointInDegrees=NewTemplate.latitudeOfLastGridPointInDegrees
			OldTemplate.longitudeOfLastGridPoint=NewTemplate.longitudeOfLastGridPoint
			OldTemplate.longitudeOfLastGridPointInDegrees=NewTemplate.longitudeOfLastGridPointInDegrees

			msg=OldTemplate.tostring()	
			grbout.write(msg)
			grbout.close()
			
			test=pygrib.open('test.grb')
			test=test[1]
					
		def extrapolate_simpel(TEMPvarTMP):
			"""
			covers the coast with water ~10km inland then reapplies the coastline with a mask from
			the new grid. This is done to remove things such as the "london pier"
			"""
			mask_var=self.D2[3].values

			for shift in (-1,1):
				for axis in (0,1):        
					a_shifted=np.roll(TEMPvarTMP,shift=shift,axis=axis)
					idx=~a_shifted.mask * TEMPvarTMP.mask
					TEMPvarTMP[idx]=a_shifted[idx]
			TEMPvarTMP=ma.masked_array(TEMPvarTMP,mask=np.logical_not(mask_var))
			return TEMPvarTMP	
		
		"""
		build KDtree based in source grid
		"""
		tree = cKDTree(zip(self.xs, self.ys, self.zs))

		"""
		interpolates data from source to target with the nearest neighbor
		or the IDW of square method using 10 neighbors

		this needs to be done iteralivly over all variables in the source!
		"""
		
		self.temp_var=self.D1[8].values #set variable, this needs to be looped!
		self.temp_var=extrapolate_simpel(self.temp_var)
		
		#nearest neighbor interpolation algorithm
		d, inds = tree.query(zip(self.xt, self.yt, self.zt), k = 1)
		
		self.temp_var_nearest = self.temp_var.flatten()[inds].reshape(self.lonTARGET.shape)
		
		#Inverse distance weighting with 10 neighbors (k=10) algorithm
		d, inds = tree.query(zip(self.xt, self.yt, self.zt), k = 10)
		w = 1.0 / d**2
		self.temp_var_idw = np.sum(w * self.temp_var.flatten()[inds], axis=1) / np.sum(w, axis=1)
		self.temp_var_idw.shape = self.lonTARGET.shape
			
		write_to_file(self.temp_var_nearest)

	
	


	def plot_algorithms(self):	
		test=pygrib.open('test.grb')
		self.temp_var_test=test[1].values

		print self.temp_var_test.shape
		print self.temp_var_nearest.shape
		
		plt.figure(figsize=(10,5))
		plt.subplot(121)
		plt.pcolormesh(self.temp_var_test)
		plt.xlim([0, self.temp_var_test.shape[0]])
		plt.ylim([0, self.temp_var_test.shape[1]])
		plt.colorbar()
		plt.title("from grib")

		plt.subplot(122)
		plt.pcolormesh(self.temp_var_nearest)
		plt.colorbar()
		plt.xlim([0, self.temp_var_nearest.shape[0]])
		plt.ylim([0, self.temp_var_nearest.shape[1]])
		plt.title("from numpy array");
	
		plt.show()
	
#.........................................................................
def main():#test
	file1='NS02_201109240000+024H00M'
	file2='NS02_SURF_201407241200+012H00M'
	interpol=KdInterpol(file1,file2)
	interpol()
	interpol.plot_algorithms()

if __name__== "__main__":
	main()
