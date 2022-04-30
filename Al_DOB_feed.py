##########
'''
Aim= get the date of observation in the formate of yyyy:Month name:dd HH:MM:SS.S
logic: used the the common date of observation files and all files related and tracde back to its location.

'''







import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from astropy.io import fits
from astropy import visualization as av
import scipy.misc 
import math as mt
from jdcal import gcal2jd,jd2gcal
import datetime
import os
import timeit
import scipy as sp
import cv2
import pathlib
import imageio
import numpy.ma as ms
from scipy.ndimage import label, generate_binary_structure,find_objects,measurements,map_coordinates
import scipy.stats as si

startTime = timeit.default_timer()
totelIm=0


a=np.loadtxt('Datesheet.dat')
b=np.loadtxt('all_ti_date.dat')
c=np.loadtxt('all_al_date.dat')
d=np.loadtxt('all_al_date_n_fname.dat',dtype='str')
e=np.loadtxt('all_ti_date_n_fname.dat',dtype='str')
f=np.loadtxt('Common.dat')
#print(a[0],b[0],c[0],d[0],e[0],f[0])

al_fn_list=[]
ti_fn_list=[]

CHarray=[]
BParray=[]
ARarray=[]
BGarray=[]
FDarray=[]

CHarea=[]
BParea=[]
ARarea=[]
BGarea=[]

CHa=[]
BPa=[]
ARa=[]
BGa=[]
FDa=[]

CHi=[]
BPi=[]
ARi=[]
BGi=[]
Fdi=[]

n_AR=[]
n_BP=[]
n_CH=[]
l_DOB=[]

feed_dates=[]


for i in range (len(f)):  #almesh images
  com=f[i]
  al_f=np.where(c==com)
  al_fn_list.append(d[al_f[0][0]])
  
Length1=len(al_fn_list)


for i in range (len(f)): #tipoly images
  com=f[i]
  ti_f=np.where(b==com)
  ti_fn_list.append(e[ti_f[0][0]])
Length=len(ti_fn_list)
#print(ti_fn_list, al_fn_list)
al_fn_list.sort()
ti_fn_list.sort()

print(ti_fn_list[55],al_fn_list[55])

for l in range(Length):
 try:
   img1=fits.open(al_fn_list[l])
   scidata1=img1[0].data
   DOB1=img1[0].header['DATE_OBS']
   #print('al mesh',DOB1)
   cen=img1[0].header['XCEN']
   ycen=img1[0].header['YCEN']
   b=np.matrix(scidata1)
   
   size=scidata1.shape
   m=scidata1.mean()
 
   dob_str=DOB1
   dob_obj = datetime.datetime.strptime(dob_str, '%Y-%m-%dT%H:%M:%S.%f')
   Mon_nam=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
   M=dob_obj.month
   #contam_time='2007-Mar-01 09:00:00'2008-01-04T11:03:43.883

   Y=dob_obj.year
   D=dob_obj.day
   Hr=dob_obj.hour
   Min=dob_obj.minute
   day_minute=round((Min/5))*5
   img_DnM=[]
   H_Y_M=('%i%i' % (Y,M))
 
   d=((Hr/24)+(Min/(24*60)))
   day=D+d
   dY=((M-1)*30.5+day)/365.25  
   obsDate=(Y+dY)
   date_array=dob_str.replace('T', ' ')
   mon=Mon_nam[int(date_array[5:7])-1] #extracted month 
   #date_array=date_array.replace(date_array[5:7],mon)
   date_array=date_array[:5]+mon+date_array[7:]
   print(date_array)
   feed_dates.append(date_array)
   #print()
 
 except:
   pass
 
np.savetxt('feed_dates.dat',feed_dates,fmt='%s') 
 
 
 
