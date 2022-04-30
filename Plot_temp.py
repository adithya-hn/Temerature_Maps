import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pylab import *
from astropy.io import fits
import scipy.misc
import math as mt
from jdcal import gcal2jd,jd2gcal
import os
import timeit
from scipy import stats as S
import scipy as sp
import pathlib
import pandas as pd


raw=np.loadtxt('Temp_data.dat')
df = pd.DataFrame(raw)
df = df.apply (pd.to_numeric, errors='coerce')
df = df.dropna()

np.savetxt('Temp_data_F.dat',df)
#print(raw.shape,df.shape)


data=(np.loadtxt('Temp_data_F.dat')).transpose()
data1=(np.loadtxt('ssn_mod_1.dat')).transpose()
Data2=data#(np.loadtxt('Al_X-ray_data_2008-12.dat')).transpose()
#print(data.shape)
pos=np.where(data<0)
data = np.delete(data, pos, axis=1)
#print(data.shape)
pathlib.Path("plots").mkdir(parents=True, exist_ok=True) 

ax1=Data2[11]
ay1=Data2[16]

x1=data[11]
y1=data[0] #CH temp
chM=np.mean(y1)
chE=np.std(y1)#S.sem(y1)


x2=data[11]
y2=data[1]  #xBP temp
bpM=np.mean(y2)
bpE=np.std(y2)

x3=data[11]
y3=data[2]  #AR temp
arM=np.mean(y3)
arE=np.std(y3)

x4=data[11]
y4=data[3]  #BG temp
bgM=np.mean(y4)
bgE=np.std(y4)

x5=data[11]
y5=data[4] #%CH

x6=data[11]
y6=data[5] #%xbps

x7=data[11]
y7=data[6] #%AR

x8=data[11]
y8=data[7] #%BG

x9=data[11]
y9=data[8] #nXBP

x10=data[11]
y10=data[9] #nAR

x11=data[11]
y11=data[10] #nCH

x12=data[11]
y12=data[12] #CHa

x13=data[11]
y13=data[13] #BPa


x14=data[11]
y14=data[14] #ARa

x15=data[11]
y15=data[15] #BGa


x16=data1[3]
y16=data1[4] #ssn

FDarea=y15+y14+y13+y12
x17=data[11]#FDT
y17=data[16]

x18=data[11]#FDintR
y18=data[17]

x19=data[11]#CHintR
y19=data[18]

x20=data[11]#BPintR
y20=data[19]

x21=data[11]#ARintR
y21=data[20]

x22=data[11]#BGintR
y22=data[21]

y23=data[23]# shape


cl=len(y5)
xl=len(y6)
al=len(y7)
bl=len(y8)


avgFDT=y17.sum()/(len(y17))
AgCp=y5.sum()/cl
AgLp=y6.sum()/xl
AgAp=y7.sum()/al
AgBp=y8.sum()/bl
Efdt=np.std(y17)
print(avgFDT/1000000)

print(chM/1000000,'CH',bpM/1000000,'XBP',arM/1000000,'AR',bgM/1000000,'BG',avgFDT/1000000,'FDT')
print(chE/1000000,'CHe',bpE/1000000,'XBPe',arE/1000000,'ARe',bgE/1000000,'BGe',Efdt/1000000)
print('------Percentage-------')
print(AgCp,AgLp,AgAp,AgBp)
rc('axes', linewidth=1)
#plt.style.use('dark_background')
plt.rcParams["xtick.major.size"] = 10
fig= plt.figure(dpi=300) #figsize=(9,5),
ax = fig.subplots()
ax.xaxis.set_tick_params(size=0.5)
ax.yaxis.set_tick_params(size=0.5)
ax.tick_params(axis='both',direction='in', length=6, width=1)
ax.tick_params(which='minor',direction='in', length=3, width=1)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.minorticks_on()
#plot=int(input())2
for plot in range(25):
##1
 if plot==1:
  #plt.text(0.99995*(max(x1)),0.94*(min(y1)+max(y1)),'(b)')
  plt.title('Variation of Temperature of Coronal Holes', fontsize=12)
  plt.ylabel('Average Temperature of CHs',fontsize=11)
  plt.xlabel('Date of Observation',fontsize=11)
  plt.plot(x1,y1, 'k',linewidth=1)
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  #plt.savefig('plots/Al_mesh-totInCH.eps')
  plt.savefig('plots/AvgTp_CH.png')
  plt.show()#block=False)
  plt.show(block=False)
  plt.pause(2)
  plt.close()
##2
 if plot==2:
  fig1= plt.figure(dpi=300)

  ax = fig1.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
    #plt.rcParams["xtick.major.size"] = 15
  #plt.text(0.9999*(max(x2)),0.94*(min(y2)+max(y2)),'(d)' )
  plt.title('Variation of Average Temperature of XBPs',fontsize=12)
  plt.ylabel('Average Temperature of XBPs',fontsize=11)
  plt.xlabel('Date of Observation',fontsize=11)
  plt.plot(x2,y2, 'k',linewidth=1)
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  #plt.savefig('plots/Al_mesh-totInXBPs.eps')
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.savefig('plots/AvgTp_XBP.png')
  plt.show()#block=False)
  #plt.pause(2)
  #plt.close()

##3
 if plot==3:
  #fig= plt.figure(figsize=(9,5),dpi=300)
  fig= plt.figure(dpi=300)
  ax = fig.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
  #plt.text(0.9999*(max(x3)),0.94*(min(y3)+max(y3)),'(a)' )
  plt.rcParams["xtick.major.size"] = 15
  plt.title('Variation of Average Temperature of Active Regions',fontsize=12)
  plt.ylabel('Average Temperature of ARs',fontsize=11)
  plt.xlabel('Date of Observation',fontsize=11)
  plt.plot(x3,y3, 'k',linewidth=1)
  #plt.savefig('plots/Al_mesh-totInAR.eps')
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.savefig('plots/AvgTp_AR.png')
  plt.show(block=False)
  plt.pause(2)
  plt.close()
  
##4
 if plot==4:
  #fig= plt.figure(figsize=(9,5),dpi=300)
  plt.rcParams["xtick.major.size"] = 15
  fig= plt.figure(dpi=300)
  ax = fig.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
  #plt.style.use('dark_background')
  #plt.text(0.9999*(max(x4)),0.92*(min(y4)+max(y4)),'(c)' )
  plt.title('Variation of Average Temperature of Background',fontsize=12)
  plt.ylabel('Average Temperature of BG ',fontsize=11)
  plt.xlabel('Date of observation',fontsize=11)
  plt.plot(x4,y4, 'k',linewidth=1)
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.savefig('plots/AvgTp_BG.png')
  plt.show(block=False)
  plt.pause(2)
  plt.close()
 
 '''
##5
 if plot==5:
  #fig= plt.figure(figsize=(9,5),dpi=300)
  plt.rcParams["xtick.major.size"] = 15
  fig= plt.figure(dpi=300)
  ax = fig.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
  plt.text(0.9999*(max(x5)),0.94*(min(y5)+max(y5)),'(b)')
  plt.title('Temperature Contribution of CH',fontsize=13)
  plt.ylabel('Percentage',fontsize=11)
  plt.xlabel('Date of observation',fontsize=11)
  plt.plot(x5,y5, 'k',linewidth=1)
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.savefig('plots/Tp_C_CH.png')
  plt.show(block=False)
  plt.pause(2)
  plt.close()

##6
 if plot==6:
  #fig= plt.figure(figsize=(9,5),dpi=300)
  plt.rcParams["xtick.major.size"] = 15
  fig= plt.figure(dpi=300)
  ax = fig.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
  plt.text(0.9999*(max(x6)),0.94*(min(y6)+max(y6)),'(d)')
  plt.title('Temperature Contribution of XBPs',fontsize=13)
  plt.ylabel('Percentage',fontsize=11)
  plt.xlabel('Date of observation',fontsize=11)
  plt.plot(x6,y6, 'k',linewidth=1)
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.savefig('plots/Tp_C_XBPs.png')
  plt.show(block=False)
  plt.pause(2)
  plt.close()

##7
 if plot==7:
  #fig= plt.figure(figsize=(9,5),dpi=300)
  plt.rcParams["xtick.major.size"] = 15
  fig= plt.figure(dpi=300)
  ax = fig.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
  plt.text(0.9999*(max(x7)),0.94*(min(y7)+max(y7)),'(a)')
  plt.title('Temperature Contribution of ARs',fontsize=13)
  plt.ylabel('Percentage',fontsize=11)
  plt.xlabel('Date of observation',fontsize=11)
  plt.plot(x7,y7, 'k',linewidth=1)
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.savefig('plots/Tp_C_AR.png')
  plt.show(block=False)
  plt.pause(2)
  plt.close()
  
##8
 if plot==8:
  #fig= plt.figure(figsize=(9,5),dpi=300)
  plt.rcParams["xtick.major.size"] = 15
  fig= plt.figure(dpi=300)   
  ax = fig.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
  plt.text(0.9999*(max(x8)),0.8*(min(y8)+max(y8)),'(c)' )
  plt.title('Temperature Contribution of BG',fontsize=13)
  plt.ylabel('Percentage',fontsize=11)
  plt.xlabel('Date of observation',fontsize=11)
  plt.plot(x8,y8, 'k',linewidth=1)
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.savefig('plots/Tp_C_BG.png')
  plt.show(block=False)
  plt.pause(2)
  plt.close()
 '''
##15
 if plot==15:
  #y18=np.multiply(y18,np.where(y23==2048, 0.25,1))
  #fig= plt.figure(figsize=(9,5),dpi=300)
  plt.rcParams["xtick.major.size"] = 15
  fig= plt.figure(dpi=300)   
  ax = fig.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
  plt.text(0.9999*(max(x8)),0.8*(min(y8)+max(y8)),'(c)' )
  plt.title('Variation of Mean of FDI ratio',fontsize=12)
  plt.ylabel('Mean of FDI ratio',fontsize=11)
  plt.xlabel('Date of observation',fontsize=11)
  plt.plot(x18,(y18/FDarea), 'k',linewidth=1)
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.savefig('plots/FDintRat.eps')
  plt.show(block=False)
  plt.pause(2)
  plt.close()
  
 if plot==16:
  #fig= plt.figure(figsize=(9,5),dpi=300)
  plt.rcParams["xtick.major.size"] = 15
  fig= plt.figure(dpi=300)
  ax = fig.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
  #plt.text(0.9999*(max(x16)),0.94*(min(y16)+max(y16)),'(g)')
  plt.title('Variation of Sunspot Numbers',fontsize=12)
  plt.ylabel('SSN',fontsize=11)
  plt.xlabel('Date of observation',fontsize=11)
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.plot(x16,y16, 'k',linewidth=1)
  plt.savefig('plots/ssn.png')
  plt.show(block=False)
  plt.pause(2)
  plt.close()
  
##16
 if plot==17:
  #fig= plt.figure(figsize=(9,5),dpi=300)
  plt.rcParams["xtick.major.size"] = 15
  fig= plt.figure(dpi=300)
  ax = fig.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
  #plt.text(0.9999*(max(x17)),0.94*(min(y17)+max(y17)),'(e)')
  #plt.text(0.9999*(max(x17)),0.94*(min(y17)+max(y17)),'(d)')
  plt.title('Variation of Full Disk Temperature',fontsize=12)
  plt.ylabel('FDT',fontsize=11)
  plt.xlabel('Date of observation',fontsize=11)
  plt.plot(x17,y17, 'k',linewidth=1)
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.savefig('plots/FDT.png')
  plt.show(block=False)
  plt.pause(2)
  plt.close()
 if plot==18:
   #fig= plt.figure(figsize=(9,5),dpi=300)
  plt.rcParams["xtick.major.size"] = 15
  fig= plt.figure(dpi=300)
  ax = fig.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
  #plt.text(0.9999*(max(x16)),0.94*(min(y16)+max(y16)),'(g)')
  plt.title('Variation of Al-mesh Intensity ',fontsize=12)
  plt.ylabel('Al-mesh Total Intensity',fontsize=11)
  plt.xlabel('Date of observation',fontsize=11)
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.plot(ax1,ay1, 'k',linewidth=1)
  plt.savefig('plots/al-mesh_tot_int.png')
  plt.show(block=False)
  plt.pause(2)
  plt.close()

####################

 if plot==19:
  ARdt=np.array(["%.3f" % i for i in x3],dtype=float64)
  ARDl=list(map(float,ARdt))
  ssnDl=list(map(float,x16)) 
  Dl=list((set(ARDl)&set(ssnDl)))
  Dli=sorted(Dl)
  Dlist=np.array(Dli)
  AR_T=[]
  XBP_T=[]
  CH_T=[]
  BG_T=[]
  FDT_T=[]
  SSN=[]
  for k in range (len(Dlist)):
    xk=np.where(ARdt== Dlist[k])
    sk = np.where(x16==Dlist[k])
    Xk=xk[0][0]
    Sk=sk[0][0]
    AR_T.append(y3[Xk])
    XBP_T.append(y2[Xk])
    CH_T.append(y1[Xk])
    BG_T.append(y4[Xk])
    FDT_T.append(y17[Xk])
    SSN.append(y16[Sk])
  #print(len(AR_T),len(SSN))  
  
  plt.rcParams["xtick.major.size"] = 15
  fig= plt.figure(dpi=300)
  ax = fig.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
  
  corr=S.spearmanr(AR_T,SSN)
  CR=round(corr[0],2)
  plt.text(0.28*(min(AR_T)+max(AR_T)),0.88*(min(SSN)+max(SSN)),(str('R = ')+str(CR)) )#SpC,bbox=dict(facecolor='none', edgecolor='black', boxstyle='square,pad=1'))
  #plt.text(0.92*(max(AR_T)),0.94*(min(SSN)+max(SSN)),'(h)')
  plt.title('AR Temperature vs Sunspot Numbers',fontsize=12)
  plt.ylabel('SSN',fontsize=11)
  plt.xlabel('AR Temperature',fontsize=11)
  #plt.xlim((0,800))
  plt.plot(AR_T,SSN, 'k*',markersize=4)
  plt.savefig('plots/ARtp vs ssn.png')
  plt.show(block=False)
  plt.pause(2)
  plt.close()

#if plot==20:
  plt.rcParams["xtick.major.size"] = 15
  fig= plt.figure(dpi=300)
  ax = fig.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
  
  corr=S.spearmanr(CH_T,SSN)
  CR=round(corr[0],2)
  plt.text(0.3*(min(CH_T)+max(CH_T)),0.88*(min(SSN)+max(SSN)),(str('R = ')+str(CR)) )#SpC,bbox=dict(facecolor='none', edgecolor='black', boxstyle='square,pad=1'))
  #plt.text(0.92*(max(CH_T)),0.94*(min(SSN)+max(SSN)),'(h)')
  plt.title('CH Temperature vs Sunspot Numbers',fontsize=12)
  plt.ylabel('SSN',fontsize=11)
  plt.xlabel('CH Temperature',fontsize=11)
  #plt.xlim((0,800))
  plt.plot(CH_T,SSN, 'k*',markersize=4)
  plt.savefig('plots/CHtp vs ssn.png')
  plt.show(block=False)
  plt.pause(2)
  plt.close()
  
#if plot==21:
  plt.rcParams["xtick.major.size"] = 15
  fig= plt.figure(dpi=300)
  ax = fig.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
  
  corr=S.spearmanr(XBP_T,SSN)
  CR=round(corr[0],2)
  plt.text(0.3*(min(XBP_T)+max(XBP_T)),0.88*(min(SSN)+max(SSN)),(str('R = ')+str(CR)) )#SpC,bbox=dict(facecolor='none', edgecolor='black', boxstyle='square,pad=1'))
  #plt.text(0.92*(max(XBP_T)),0.94*(min(SSN)+max(SSN)),'(h)')
  plt.title('XBP Temperature vs Sunspot Numbers',fontsize=12)
  plt.ylabel('SSN',fontsize=11)
  plt.xlabel('XBP Temperature',fontsize=11)
  #plt.xlim((0,800))
  plt.plot(XBP_T,SSN, 'k*',markersize=4)
  plt.savefig('plots/XBPtp vs ssn.png')
  plt.show(block=False)
  plt.pause(2)
  plt.close()

#if plot==22:
  plt.rcParams["xtick.major.size"] = 15
  fig= plt.figure(dpi=300)
  ax = fig.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
  
  corr=S.spearmanr(BG_T,SSN)
  CR1=round(corr[0],2)
  plt.text(0.35*(min(BG_T)+max(BG_T)),0.88*(min(SSN)+max(SSN)),(str('R = ')+str(CR1)) )#SpC,bbox=dict(facecolor='none', edgecolor='black', boxstyle='square,pad=1'))
  #plt.text(0.92*(max(BG_T)),0.94*(min(SSN)+max(SSN)),'(h)')
  plt.title('BG Temperature vs Sunspot Numbers',fontsize=12)
  plt.ylabel('SSN',fontsize=11)
  plt.xlabel('BG Temperature',fontsize=11)
  plt.plot(BG_T,SSN, 'k*',markersize=4)
  plt.savefig('plots/BGtp vs ssn.png')
  plt.show(block=False)
  plt.pause(2)
  plt.close()
  
#if plot==23:
  plt.rcParams["xtick.major.size"] = 15
  fig= plt.figure(dpi=300)
  ax = fig.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
  
  corr=S.spearmanr(FDT_T,SSN)
  CR=round(corr[0],2)
  plt.text(0.35*(min(FDT_T)+max(FDT_T)),0.88*(min(SSN)+max(SSN)),(str('R = ')+str(CR)) )#SpC,bbox=dict(facecolor='none', edgecolor='black', boxstyle='square,pad=1'))
  #plt.text(0.92*(max(FDT_T)),0.94*(min(SSN)+max(SSN)),'(h)')
  plt.title('Full Disk Temperature vs Sunspot Numbers',fontsize=12)
  plt.ylabel('SSN',fontsize=11)
  plt.xlabel('FD Temperature',fontsize=11)
  #plt.xlim((0,800))
  plt.plot(FDT_T,SSN, 'k*',markersize=4)
  plt.savefig('plots/FDTtp vs ssn.png')
  plt.show(block=False)
  plt.pause(2)
  plt.close()

 if plot==20: #CH
   #fig= plt.figure(figsize=(9,5),dpi=300)
  #y19 = np.multiply(y19, np.where(y23 == 2048, 0.25, 1))
  plt.rcParams["xtick.major.size"] = 15
  fig= plt.figure(dpi=300)
  ax = fig.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
  #plt.text(0.9999*(max(x16)),0.94*(min(y16)+max(y16)),'(g)')
  plt.title('Variation of CH Mean Intensity Ratio ',fontsize=12)
  plt.ylabel('CHs Mean Intensity Ratio',fontsize=11)
  plt.xlabel('Date of observation',fontsize=11)
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.plot(x19,(y19/y12), 'k',linewidth=1)
  plt.savefig('plots/CHintRat.eps')
  plt.show(block=False)
  plt.pause(2)
  plt.close()

 if plot==21: #BP
   #fig= plt.figure(figsize=(9,5),dpi=300)
  #y20 = np.multiply(y20, np.where(y23 == 2048, 0.25, 1))
  plt.rcParams["xtick.major.size"] = 15
  fig= plt.figure(dpi=300)
  ax = fig.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
  #plt.text(0.9999*(max(x16)),0.94*(min(y16)+max(y16)),'(g)')
  plt.title('Variation of XBP Mean Intensity Ratio ',fontsize=12)
  plt.ylabel('XBPs Mean Intensity Ratio',fontsize=11)
  plt.xlabel('Date of observation',fontsize=11)
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.plot(x20,(y20/y13), 'k',linewidth=1)
  plt.savefig('plots/XBPintRat.eps')
  plt.show(block=False)
  plt.pause(2)
  plt.close()

 if plot==21: #AR
   #fig= plt.figure(figsize=(9,5),dpi=300)
  #y21 = np.multiply(y21, np.where(y23 == 2048, 0.25, 1))
  plt.rcParams["xtick.major.size"] = 15
  fig= plt.figure(dpi=300)
  ax = fig.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
  #plt.text(0.9999*(max(x16)),0.94*(min(y16)+max(y16)),'(g)')
  plt.title('Variation of AR Mean Intensity Ratio ',fontsize=12)
  plt.ylabel('AR Mean Intensity Ratio',fontsize=11)
  plt.xlabel('Date of observation',fontsize=11)
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.plot(x21,(y21/y14), 'k',linewidth=1)
  plt.savefig('plots/ARintRat.eps')
  plt.show()#block=False)
  #plt.pause(2)
  #plt.close()

 if plot==22: #BG
   #fig= plt.figure(figsize=(9,5),dpi=300)
  #y22 = np.multiply(y22, np.where(y23 == 2048, 0.25, 1))
  plt.rcParams["xtick.major.size"] = 15
  fig= plt.figure(dpi=300)
  ax = fig.subplots()
  ax.xaxis.set_tick_params(size=0.5)
  ax.yaxis.set_tick_params(size=0.5)
  ax.tick_params(axis='both',direction='in', length=6, width=1)
  ax.tick_params(which='minor',direction='in', length=3, width=1)
  ax.yaxis.set_ticks_position('both')
  ax.xaxis.set_ticks_position('both')
  ax.minorticks_on()
  #plt.text(0.9999*(max(x16)),0.94*(min(y16)+max(y16)),'(g)')
  plt.title('Variation of BG Mean Intensity Ratio ',fontsize=12)
  plt.ylabel('BG Mean Intensity Ratio',fontsize=11)
  plt.xlabel('Date of observation',fontsize=11)
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.plot(x22,(y22/y15), 'k',linewidth=1)
  plt.savefig('plots/BGintRat.eps')
  plt.show(block=False)
  plt.pause(2)
  plt.close()



