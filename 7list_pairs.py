
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
#pathlib.Path("Temp_Output").mkdir(parents=True, exist_ok=True) 

newlist=[]
GOESdata=[]

DM=[]
tiDate=[]
tiFname=[]
alDate=[]
alFname=[]
imgsDate=[]
path = "/media/adithyahn/Adi_Backup_drive/XRT_data/solar.physics.montana.edu/HINODE/XRT/SCIA/synop_images/comp_part1"

fname = []
Fname=[]

for root,d_names,f_names in os.walk(path):
  for f in f_names:
    if f.endswith(".fits"):
      Fname.append(f) #with only file names
      fname.append(os.path.join(root, f)) #file name along with path
fs=[] #name string

fname.sort()
Fname.sort()
#print(Fname)
Length=len(fname)
print(Length)
for l in range(Length):
 try:
   img=fits.open(fname[l])
   scidata=img[0].data
   DOB=img[0].header['DATE_OBS']
   cen=img[0].header['XCEN']
   ycen=img[0].header['YCEN']
   b=np.matrix(scidata)
 except:
   pass

 start=0
 Filter1=img[0].header['EC_FW1_']
 Filter2=img[0].header['EC_FW2_']
 
 if Filter1== 'Ti_poly':
   start=1
 elif Filter1== 'Open':
   if Filter2== 'Ti_poly':
     start=1 
 
 if Filter1== 'Al_mesh':
   start=2
 elif Filter1== 'Open':
   if Filter2== 'Al_mesh':
     start=2
     
 else:
   start=0
 
 if start==1:
   #print("A",Fname[l] )  
   size=scidata.shape
   m=scidata.mean()
   dateA=[]
 
   dob_str=DOB
   dob_obj = datetime.datetime.strptime(dob_str, '%Y-%m-%dT%H:%M:%S.%f')

   M=dob_obj.month
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
   dateA.append(Y)
   dateA.append(M)
   dateA.append(D)
   dateA.append(Hr)
   dateA.append(Min)
   dateTi=int(int(str(Y)+str(M)+str(D)+str(Hr)))
   tiDate.append(dateTi)
   tiFname.append(fname[l])
   
 if start==2:
   #print("T",Fname[l] ) 
   size=scidata.shape
   m=scidata.mean()
   dateA=[]
 
   dob_str=DOB
   dob_obj = datetime.datetime.strptime(dob_str, '%Y-%m-%dT%H:%M:%S.%f')

   M=dob_obj.month
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
   
   obsDate=(Y+dY)
   dateA.append(Y)
   dateA.append(M)
   dateA.append(D)
   dateA.append(Hr)
   dateA.append(Min)
   dateAl=int(int(str(Y)+str(M)+str(D)+str(Hr)))
   DM.append(dateAl) 
   alDate.append(dateAl)
   imgsDate.append(Fname[l])
   alFname.append(fname[l])
 print('[',l+1,'/',Length,']')
alDate=(np.array(alDate)).round(5)
tiDate=(np.array(tiDate)).round(5)
common=list(set(alDate)&set(tiDate)) 
#print(tiDate,alDate,common) 
#print(alFname)
np.savetxt('Datesheet.dat',DM,fmt='%4.0f')
np.savetxt('Common.dat',common,fmt='%4.0f')
np.savetxt('all_al_date_n_fname.dat',alFname,fmt='%s')
np.savetxt('all_imgs_date_scmp.dat',imgsDate,fmt='%s')
np.savetxt('all_al_date.dat',alDate,fmt='%s')
np.savetxt('all_ti_date.dat',tiDate,fmt='%s')
np.savetxt('all_ti_date_n_fname.dat',tiFname,fmt='%s')

stopTime = timeit.default_timer()
runtime=stopTime-startTime
print('')
print('......COMPLEETED.....')






