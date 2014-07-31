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
import time

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
		
		def write_to_file(variable,iter,filename):
			"""
			writing to grib file
			"""	
			filename=filename[49:] #filename without path
			variable=variable.filled(9999) #grib_apis interpretation of 9999 are '--' 	

			OldTemplate=self.D1[iter]		
			
			grbout=open(filename,'ab')

			OldTemplate.dataDate=time.strftime("%Y%m%d")
			OldTemplate.values=variable	
			
			msg=OldTemplate.tostring()	
			grbout.write(msg)
			grbout.close()
					
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
			TEMPvarTMP=ma.masked_array(TEMPvarTMP,mask=np.logical_not(mask_var),fill_value=np.nan)	
			
			return TEMPvarTMP	
		
		"""
		build KDtree based in source grid
		"""
		tree = cKDTree(zip(self.xs, self.ys, self.zs))

		"""
		interpolates data from source to target with the nearest neighbor

		this needs to be done iteralivly over all variables in the source!
		"""
	
		self.i=1
		for MSG in self.D1:

			print "Interpolating and writing msg nr: ",self.i
			
			self.temp_var=self.D1[self.i].values #read current grib messege
			self.temp_var=extrapolate_simpel(self.temp_var) #makes the data fit the new grid
			
			"""
			nearest neighbor interpolation algorithm
			"""
			d, inds = tree.query(zip(self.xt, self.yt, self.zt), k = 1)	
			self.temp_var_nearest = self.temp_var.flatten()[inds].reshape(self.lonTARGET.shape)
			write_to_file(self.temp_var_nearest,self.i,self.source)
			
			self.i+=1 #iterator, goes from 1 to last grib messege
	
	@staticmethod
	def plot_algorithms(file1,file2):	
		file1=pygrib.open(file1)
		temp_var_file1=file1[8].values

		file2=pygrib.open(file2)	
		temp_var_file2=file2[3].values
		
		file1.close()
		file2.close()
		
		plt.figure(figsize=(10,5))
		plt.subplot(121)
		plt.pcolormesh(temp_var_file1)
		plt.xlim([0, temp_var_file1.shape[0]])
		plt.ylim([0, temp_var_file1.shape[1]])
		plt.colorbar()
		plt.title("from file1")

		plt.subplot(122)
		plt.pcolormesh(temp_var_file2)
		plt.colorbar()
		plt.xlim([0, temp_var_file2.shape[0]])
		plt.ylim([0, temp_var_file2.shape[1]])
		plt.title("from file2");
	
		plt.show()
	
#.........................................................................
import glob
import os
def main():
	MODE='plot'# or 'run'

	if MODE=='plot':
		file1='NS02_201101010000+096H00M'
		file2='NS02_SURF_201407241200+012H00M'
		KdInterpol.plot_algorithms(file1,file2)
	if MODE=='run':
		INPUTpath=''
		OUTPUTpath=''
		file2=''
		os.chdir(OUTPUTpath)
		for fname in sorted(glob.iglob(INPUTpath)):
				print fname[49:] #only filename without path
				file1=fname
				interpol=KdInterpol(file1,file2)
				interpol()

if __name__== "__main__":
	main()
