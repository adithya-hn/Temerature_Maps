import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from astropy.io import fits
from astropy import visualization as av
import scipy.misc
import math as mt
from jdcal import gcal2jd, jd2gcal
import datetime
import os
import timeit
import scipy as sp
import cv2
import pathlib
import imageio
import numpy.ma as ms
from scipy.ndimage import label, generate_binary_structure, find_objects, measurements, map_coordinates,shift
import scipy.stats as si
import zipfile


from skimage import data, color
from skimage.transform import hough_circle, hough_circle_peaks
from skimage.feature import canny
from skimage.draw import circle_perimeter
from skimage.util import img_as_ubyte

startTime = timeit.default_timer()
totelIm = 0
filetime=0
prvTime=startTime

#a = np.loadtxt('Datesheet.dat')
b = np.loadtxt('all_ti_date.dat')
c = np.loadtxt('all_al_date.dat')
d = np.loadtxt('all_al_date_n_fname.dat', dtype='str')
e = np.loadtxt('all_ti_date_n_fname.dat', dtype='str')
f = np.loadtxt('Common.dat')

r_data = (np.loadtxt('Resp_data_offcial_with_date.txt')).transpose()
temp = r_data[0]
Tiresp = r_data[1]
Alresp = r_data[2]
pathlib.Path("TempMaps").mkdir(parents=True, exist_ok=True)
pathlib.Path("IRmaps").mkdir(parents=True, exist_ok=True)
#pathlib.Path("imgs2").mkdir(parents=True, exist_ok=True)
pathlib.Path("TempFits").mkdir(parents=True, exist_ok=True)
#pathlib.Path("IR_fits").mkdir(parents=True, exist_ok=True)

# print(r_data[1])


al_fn_list = []
ti_fn_list = []

CHarray = []
BParray = []
ARarray = []
BGarray = []
FDarray = []

CHarea = []
BParea = []
ARarea = []
BGarea = []

CHa = []
BPa = []
ARa = []
BGa = []
FDa = []

CHi = []
BPi = []
ARi = []
BGi = []
Fdi = []

n_AR = []
n_BP = []
n_CH = []
l_DOB = []
conv_rat=[]

IR_data = []
bpIRsum=[]
chIRsum=[]
bgIRsum=[]
ArIRsum=[]

shapeArry=[]


for i in range(len(f)):  # almesh images
    com = f[i]
    al_f = np.where(c == com)
    al_fn_list.append(d[al_f[0][0]])

Length1 = len(al_fn_list)

for i in range(len(f)):  # tipoly images
    com = f[i]
    ti_f = np.where(b == com)
    ti_fn_list.append(e[ti_f[0][0]])

Length = len(ti_fn_list)
al_fn_list.sort()
ti_fn_list.sort()
pair_imgs=(np.vstack((al_fn_list,ti_fn_list))).transpose()
#print(pair_imgs)
for l in range(Length):
    try:
        img1 = fits.open(pair_imgs[l][0])#(al_fn_list[l])
        scidata1 = img1[0].data
        DOB1 = img1[0].header['DATE_OBS']
        dob_str = DOB1
        dob_obj = datetime.datetime.strptime(dob_str, '%Y-%m-%dT%H:%M:%S.%f')
        M = dob_obj.month
        Y = dob_obj.year
        D = dob_obj.day
        Hr = dob_obj.hour
        Min = dob_obj.minute
        day_minute = round((Min / 5)) * 5
        H_Y_M = ('%i%i' % (Y, M))
        dateAl = int(int(str(Y) + str(M) + str(D) + str(Hr)))
        x1cen = img1[0].header['XCEN']
        y1cen = img1[0].header['YCEN']
        x1scale = img1[0].header['XSCALE']
        y1scale = img1[0].header['YSCALE']
        scidata1=shift(scidata1,((x1cen/x1scale),(y1cen/y1scale)),cval=1)
        scidata = scidata1
        size1 = scidata1.shape
        Ms = str(M).zfill(2)
        Ds = str(D).zfill(2)
        HrS = str(Hr).zfill(2)
        str1='comp_XRT'
        ti_fl = np.where(b == dateAl)

        if size1[0] == 1024:
            kcanny = canny(scidata1, sigma=2, low_threshold=8, high_threshold=100)
            hough_radii = np.arange(466, 475, 2)
            hough_res = hough_circle(kcanny, hough_radii)
            accums, cx, cy, rad = hough_circle_peaks(hough_res, hough_radii,total_num_peaks=1)
            Rad = int(rad[0])
            center = (int(cx[0]),int(cy[0]))
            scidata1 = shift(scidata1, (-(512-center[0]), -(512-center[1])), cval=1)
            center=(512,512)

        else:
            kcanny = canny(scidata1, sigma=2, low_threshold=5, high_threshold=30)
            hough_radii = np.arange(924, 955, 2)
            hough_res = hough_circle(kcanny, hough_radii)
            accums, cx, cy, rad = hough_circle_peaks(hough_res, hough_radii, total_num_peaks=1)
            Rad = int(rad[0])
            center = (int(cx[0]), int(cy[0]))
            #print(center)
            #print(1024-center[0],1024-center[1])
            scidata1 = shift(scidata1, (-(1024-center[0]), -(1024-center[1])), cval=1)
            center=(1024,1024)
            #ti_fl = np.where(b == dateAl)
            #print('got tipoly loc')
            #print(dateAl, b[ti_fl[0][0]], e[ti_fl[0][0]])


    except:
        pass

    try:
        #print(dateAl,b[ti_fl[0][0]])
        #print(b.count(dateAl))
        #print(ti_fl = np.where(b == dateAl))  ######Ti poly #to avaid mismatching
        #print(dateAl,b[ti_fl[0][0]],e[ti_fl[0][0]] )
        f_name=e[ti_fl[0][0]]
        img2 = fits.open(e[ti_fl[0][0]])
        #img2 = fits.open(ti_fn_list[l])
        scidata2 = img2[0].data
        DOB2 = img2[0].header['DATE_OBS']
        #print('ti poly',DOB2)
        x2cen = img2[0].header['XCEN']
        y2cen = img2[0].header['YCEN']
        x2scale = img2[0].header['XSCALE']
        y2scale = img2[0].header['YSCALE']
        scidata2 = shift(scidata2, ((x2cen / x2scale), (y2cen / y2scale)),cval=1)
        #b = np.matrix(scidata2)
        size2 = scidata2.shape
        dob_str2 = DOB2
        dob_obj2 = datetime.datetime.strptime(dob_str2, '%Y-%m-%dT%H:%M:%S.%f')
        M2 = dob_obj2.month
        Y2 = dob_obj2.year
        D2 = dob_obj2.day
        Hr2 = dob_obj2.hour
        Min2 = dob_obj2.minute
        day_minute2 = round((Min2 / 5)) * 5
        H_Y_M2 = ('%i%i' % (Y2, M2))
        dateTi = int(int(str(Y2) + str(M2) + str(D2) + str(Hr2)))
        #print('ok')
        if size2[0] == 1024:
            kcanny2 = canny(scidata2, sigma=2, low_threshold=20, high_threshold=80)
            hough_rad2 = np.arange(466, 475, 2)
            hough_res2 = hough_circle(kcanny2, hough_rad2)
            accums2, cx2, cy2, rad2 = hough_circle_peaks(hough_res2, hough_rad2,total_num_peaks=1)
            Rad2 = int(rad2[0])
            center2 = (int(cx2[0]), int(cy2[0]))
            scidata2= shift(scidata2, (-(512-center2[0]), -(512-center2[1])), cval=1)
            center2=(512,512)

        else:
            kcanny2 = canny(scidata2, sigma=2, low_threshold=2, high_threshold=20)
            hough_rad2 = np.arange(924, 955, 2)
            hough_res2 = hough_circle(kcanny2, hough_rad2)
            accums2, cx2, cy2, rad2 = hough_circle_peaks(hough_res2, hough_rad2,total_num_peaks=1)
            Rad2 = int(rad2[0])
            center2 = (int(cx2[0]), int(cy2[0]))
            scidata2 = shift(scidata2, (-(1024-center2[0]), -(1024-center2[1])), cval=1)
            center2=(1024,1024)
	
    except:
        print('skipped')
        pass

    #print(size1,size2)
    if size1 == size2:
        start = 1
        if dateAl == dateTi:
            start = 1
        else:
            print('Time is not matched')
            start = 0
    else:
        start = 0
        print('size not matched',DOB1, size1,size2)


    if start == 1:
        size1 = scidata1.shape
        shapeArry.append(size1[0])
        m = scidata1.mean()
        dob_str = DOB1
        dob_obj = datetime.datetime.strptime(dob_str, '%Y-%m-%dT%H:%M:%S.%f')
        M = dob_obj.month
        Y = dob_obj.year
        D = dob_obj.day
        Hr = dob_obj.hour
        Min = dob_obj.minute
        day_minute = round((Min / 5)) * 5
        img_DnM = []
        H_Y_M = ('%i%i' % (Y, M))

        d = ((Hr / 24) + (Min / (24 * 60)))
        day = D + d
        dY = ((M - 1) * 30.5 + day) / 365.25
        obsDate = (Y + dY)

        # print(DOB1)
        l_DOB.append(obsDate)
        aa = np.zeros((size1))
        CHmask = np.zeros(size1, np.uint8)
        BPmask = np.zeros(size1, np.uint8)
        ARmask = np.zeros(size1, np.uint8)
        ARmask1 = np.zeros(size1, np.uint8)
        BGmask = np.zeros(size1, np.uint8)
        #print(center)
        h = center[0]
        k = center[1]
        R = Rad - 20

        circ = cv2.circle(aa, (h, k), R, (255, 0, 0), -1)  # disk
        circle = circ.astype(np.bool_)
        Circle = np.invert(circle)  # hole
        mask = ms.array(scidata1, mask=Circle)
        disk = ms.array(scidata1, mask=circle)  # hided disk
        DD = disk * 0
        dd = ms.array(DD, mask=Circle)

        sun = DD.data  # only solar disk
        index = np.nonzero(sun)
        mtxData = sun[index[0], index[1]]
        s = np.std(mtxData)
        mn = mtxData.mean()
        md = np.median(mtxData)

        brightTh = mn * 2.5  # active regin thresh 1.7 #mn*.3
        threshInt = md * 0.3
        # print(md,mn,s)

        CHt = np.where(sun < threshInt, 255, 0)
        BPt = np.where(sun > brightTh, 255, 0)

        CHm = ms.array(CHt, mask=circle)
        BPm = ms.array(BPt, mask=circle)
        BPd = BPm * 0
        CHd = CHm * 0
        # chSum=(CHt.sum())/255

        CHs = CHd.astype(np.uint8)  # Converting images to 8bit
        BPs = BPd.astype(np.uint8)

        # morphological closing
        kernel = np.ones((15, 15), np.uint8)
        kernel1 = np.ones((5, 5), np.uint8)
        CH = cv2.morphologyEx(CHs, cv2.MORPH_CLOSE, kernel)
        BP = cv2.morphologyEx(BPs, cv2.MORPH_CLOSE, kernel1)

        CHcont, hierarchy = cv2.findContours(CH, cv2.RETR_EXTERNAL,
                                             cv2.CHAIN_APPROX_NONE)  # ret _ext to detect outer most cont , approx_none -no need to apporoximate give full cont pts
        BPcont, hierarchy = cv2.findContours(BP, cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)

        zimage = av.ZScaleInterval(nsamples=600, contrast=0.25, max_reject=0.5, min_npixels=5, krej=2.5,
                                   max_iterations=5)
        z = zimage.get_limits(scidata1)
        ZI = np.clip(scidata1, z[0], z[1])  # zimage(scidata)
        ZS = (ZI / ZI.max()) * 255  #
        SUN = ((ZS).astype(np.uint8))

        imgCont = cv2.cvtColor(SUN, cv2.COLOR_GRAY2BGR)
        Org_img = cv2.cvtColor(SUN, cv2.COLOR_GRAY2BGR)
        CH_img = cv2.cvtColor(SUN, cv2.COLOR_GRAY2BGR)
        BP_img = cv2.cvtColor(SUN, cv2.COLOR_GRAY2BGR)
        AR_img = cv2.cvtColor(SUN, cv2.COLOR_GRAY2BGR)
        AR_img1 = cv2.cvtColor(SUN, cv2.COLOR_GRAY2BGR)
        BG_img = cv2.cvtColor(SUN, cv2.COLOR_GRAY2BGR)
        comb_img = imgCont
        ar = []
        BP_cont = []  # AR excluded
        AR_cont = []
        AR_area = 0
        nBP = len(BPcont)
        nCH = len(CHcont)
        region = np.zeros(size1, np.uint8)
        # ar second operation
        ARTh = mn * 2.5  # mn*2.5
        ARt = np.where(sun > ARTh, 255, 0)
        ARs = ARt.astype(np.uint8)
        AR = cv2.morphologyEx(ARs, cv2.MORPH_CLOSE, kernel1)
        ARcont, hierarchy = cv2.findContours(AR, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        nAR = len(ARcont)
        for j in range(nAR):
            b1 = cv2.contourArea(ARcont[j])
            cv2.drawContours(ARmask1, [ARcont[j]], -1, (255, 0, 0), cv2.FILLED)
            if b1 > 1000:
                AR_cont.append(ARcont[j])
                ARarea.append(b1)
                cv2.drawContours(imgCont, ARcont[j], 1, (0, 0, 255), cv2.FILLED)
                cv2.drawContours(region, ARcont[j], -1, (255, 255, 255), cv2.FILLED)
                cv2.drawContours(AR_img1, [ARcont[j]], -1, (255, 255, 0), 2)  # ar image
                cv2.drawContours(comb_img, [ARcont[j]], -1, (255, 255, 0), 2)
                # cv2.drawContours(ARmask1, [ARcont[j]],-1 ,(255,0,0), cv2.FILLED)
            else:
                BP_cont.append(ARcont[j])
        Bo_ARmask2 = ARmask1.astype(np.bool_)
        In_ARmask2 = np.invert(Bo_ARmask2)
        AR_masked2_sun = ms.array(sun, mask=ARmask1)
        index2 = np.nonzero(AR_masked2_sun)
        mtxData2 = AR_masked2_sun[index2[0], index2[1]]
        s1 = np.std(mtxData2)
        MeanII = mtxData2.mean()
        _BPTh = MeanII * 1.7  # mn*2.5
        _BPt = np.where(AR_masked2_sun > _BPTh, 255, 0)
        _BPs = _BPt.astype(np.uint8)
        _BP = cv2.morphologyEx(_BPs, cv2.MORPH_CLOSE, kernel)
        _BPcont, hierarchy = cv2.findContours(_BP, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        n_BP_ = len(_BPcont)
        counter = 0
        for p in range(n_BP_):
            _c = cv2.contourArea(_BPcont[p])
            # Mom=cv2.moments(_BPcont[p])
            # cx = int(Mom['m10']/Mom['m00'])
            # moment=Mom['m10']
            # print("---")
            if _c > 1000:
                # BP_cont.append(ARcont[j])
                # cv2.drawContours(comb_img, [_BPcont[p]],-1,(255,255,255), 2)
                counter = 1
            elif _c > 10 and _c < 1000:
                BP_cont.append(_BPcont[p])



        color = (255, 255, 10)
        color1 = (25, 10, 100)
        sun_img = cv2.cvtColor(SUN, cv2.COLOR_GRAY2BGR)
        #cv2.drawContours(sun_img, BR, -1, (115, 100, 50), 2)

        #sun_img = cv2.circle(sun_img,center ,Rad, color1, 2)
        Sun_img = cv2.circle(sun_img, center, Rad, color, 2)
        #cv2.drawContours(Sun_img, AR_cont, -1, (255, 255, 0), 2)
        #imageio.imwrite('C-fit/cirT3_fit{}.jpg'.format(dob_str), Sun_img[::-1])
        #########~~~~~~~~~##########

        ###################--Temp map--#############################################################################-

        IR = abs(scidata2) / abs(scidata1)  # Ti-Poly/Al-Mesh
        IR = IR.round(3)  # smoothning

        # reponse data interpolation
        sld2 = (l + 1) * 26  # selected data array
        sld1 = l * 26
        ratio = np.array(Tiresp[sld1:sld2] / Alresp[sld1:sld2])
        # print(ratio,sld1,sld2)
        # tp=np.linspace(MinT,MaxT,200)

        Fr = []
        x1 = []
        y1 = []

        for k in range(25):
            yi = np.linspace(ratio[k], ratio[k + 1], 100)
            xi = (temp[k + 1] - temp[k]) / (ratio[k + 1] - ratio[k]) * (yi - ratio[k]) + temp[k]
            x1.append(xi[0:99])
            y1.append(yi[0:99])
            #print(ratio[k])

        x1 = (np.array(x1)).flatten()
        y1 = (np.array(y1)).flatten()
        y1 = y1.round(3)
        TempMap = np.zeros(IR.shape)

        for m in range(len(y1)):
            TempMap[IR == y1[m]] = x1[m]
        # print(TempMap[555,555], scidata1[555,555])

        disk = ms.array(TempMap, mask=circle)  # hided disk
        DD = disk * 0
        TempMap = DD.data
        #hdu = fits.PrimaryHDU(TempMap)
        #hdu.writeto("TempFits/Temp_Map{}.fits".format(DOB1),overwrite=True)
        #zip_file=zipfile.ZipFile("TempFits/Temp_Map{}.zip".format(DOB1),'w')
        #zip_file.write("TempFits/Temp_Map{}.fits".format(DOB1),compress_type= zipfile.ZIP_DEFLATED)
        #zip_file.close()
        #os.remove("TempFits/Temp_Map{}.fits".format(DOB1))
        #/media/adithyahn/Adi_Backup_drive/FitsTemp_map

        IRdisk = ms.array(IR, mask=circle)  # hided disk
        IRDD = IRdisk * 0
        IRMap = IRDD.data
        IRsum = IRMap.sum()
        #hdu = fits.PrimaryHDU(IRMap)
        #hdu.writeto("IR_fits/IR_Map{}.fits".format(DOB1),overwrite=True)

        index = np.nonzero(TempMap)
        mtxData1 = TempMap[index[0], index[1]]
        s1 = np.std(mtxData1)
        mn1 = mtxData1.mean()
        md1 = np.median(mtxData1)
        mn1_mk = mn1 / 1000000
        print('[', l+1, '/', Length, '] ', DOB1, 'Avg. FDT', mn1_mk.round(2),'PRV_TIME',filetime)
        zimage = av.ZScaleInterval(nsamples=600, contrast=0.25, max_reject=0.5, min_npixels=5, krej=2.5,
                                   max_iterations=5)
        z1 = zimage.get_limits(TempMap)
        ZI1 = np.clip(TempMap, z1[0], z1[1])  # zimage(scidata)
        ZS1 = (ZI1 / ZI1.max()) * 255
        tsunz = ((ZS1).astype(np.uint8))


        #########################################################################################///////////////////

        ######--IR contribution---#######
        # scaling IR map
        z2 = zimage.get_limits(IRMap)
        ZI2 = np.clip(IRMap, z2[0], z2[1])  # zimage(scidata)
        ZS2 = (ZI2 / ZI2.max()) * 255
        IRz = ((ZS2).astype(np.uint8))

        CH_IR = cv2.cvtColor(IRz, cv2.COLOR_GRAY2BGR)
        BP_IR = cv2.cvtColor(IRz, cv2.COLOR_GRAY2BGR)
        AR_IR = cv2.cvtColor(IRz, cv2.COLOR_GRAY2BGR)
        BG_IR = cv2.cvtColor(IRz, cv2.COLOR_GRAY2BGR)
        BG_Ir = cv2.cvtColor(IRz, cv2.COLOR_GRAY2BGR) #why?
        comb_IR = cv2.cvtColor(IRz, cv2.COLOR_GRAY2BGR)

        CHirMask = np.zeros(size1, np.uint8)
        BPirMask = np.zeros(size1, np.uint8)
        ARirMask = np.zeros(size1, np.uint8)
        ARirMask1 = np.zeros(size1, np.uint8)
        BGirMask = np.zeros(size1, np.uint8)

        cv2.drawContours(BG_IR, CHcont, -1, (255, 0, 0), cv2.FILLED)
        cv2.drawContours(BG_IR, BP_cont, -1, (255, 0, 0), cv2.FILLED)
        cv2.drawContours(BG_IR, AR_cont, -1, (255, 0, 0), cv2.FILLED)
        cv2.drawContours(ARirMask1, AR_cont, -1, (255, 0, 0), cv2.FILLED)
        cv2.drawContours(CHirMask, CHcont, -1, (255, 0, 0), cv2.FILLED)
        cv2.drawContours(BPirMask, BP_cont, -1, (255, 0, 0), cv2.FILLED)
        #cv2.drawContours(CH_IR, CHcont, -1, (255, 0, 255), 2)
        #cv2.drawContours(BP_IR, BP_cont, -1, (255, 0, 0), 2)
        cv2.drawContours(comb_IR, CHcont, -1, (255, 0, 255), 2)
        cv2.drawContours(comb_IR, BP_cont, -1, (255, 0, 0), 2)
        cv2.drawContours(comb_IR, AR_cont, -1, (255, 255, 0), 2)

        Bo_CHirMask = CHirMask.astype(np.bool_)
        Bo_BPirMask = BPirMask.astype(np.bool_)
        Bo_BGirMask = BGirMask.astype(np.bool_)
        Bo_ARirMask1 = ARirMask1.astype(np.bool_)
        In_CHirMask = np.invert(Bo_CHirMask)
        In_BGirMask = np.invert(Bo_BGirMask)
        In_BPirMask = np.invert(Bo_BPirMask)
        In_ARirMask1 = np.invert(Bo_ARirMask1)

        # masking
        CH_IRmasked_sun = ms.array(IRMap, mask=In_CHirMask)
        BP_IRmasked_sun = ms.array(IRMap, mask=In_BPirMask)
        BG_IRmasked_sun = ms.array(IRMap, mask=Bo_BGirMask)
        AR_IRmasked1_sun = ms.array(IRMap, mask=In_ARirMask1)

        bgIRimage = ms.array(IRMap, mask=In_BGirMask)
        BGIr = ((bgIRimage * 0).data).astype(np.uint8)
        BPIr = ((BP_IRmasked_sun * 0).data).astype(np.uint8)
        BP_IRmasked1_sun = ms.array(IRMap, mask=In_BPirMask)
        bpIRsum1 = BP_IRmasked1_sun.sum()# both are same
        #print(bpIRsum1,BP_IRmasked_sun.sum())

        bpIRsum.append(BP_IRmasked_sun.sum())
        chIRsum.append(CH_IRmasked_sun.sum())
        bgIRsum.append(BG_IRmasked_sun.sum())
        ArIRsum.append(AR_IRmasked1_sun.sum())

        ######-------|||-------###########

        CH_Temp = cv2.cvtColor(tsunz, cv2.COLOR_GRAY2BGR)
        BP_Temp = cv2.cvtColor(tsunz, cv2.COLOR_GRAY2BGR)
        AR_Temp = cv2.cvtColor(tsunz, cv2.COLOR_GRAY2BGR)
        BG_Temp = cv2.cvtColor(tsunz, cv2.COLOR_GRAY2BGR)
        BG_T = cv2.cvtColor(tsunz, cv2.COLOR_GRAY2BGR)
        comb_Temp = cv2.cvtColor(tsunz, cv2.COLOR_GRAY2BGR)

        CHTmask = np.zeros(size1, np.uint8)
        BPTmask = np.zeros(size1, np.uint8)
        ARTmask = np.zeros(size1, np.uint8)
        ARTmask1 = np.zeros(size1, np.uint8)
        BGTmask = np.zeros(size1, np.uint8)

        cv2.drawContours(BGTmask, CHcont, -1, (255, 0, 0), cv2.FILLED)
        cv2.drawContours(BGTmask, BP_cont, -1, (255, 0, 0), cv2.FILLED)
        cv2.drawContours(BGTmask, AR_cont, -1, (255, 0, 0), cv2.FILLED)
        cv2.drawContours(ARTmask1, AR_cont, -1, (255, 0, 0), cv2.FILLED)
        #cv2.drawContours(CH_Temp, CHcont, -1, (255, 0, 255), 2)
        #cv2.drawContours(BP_Temp, BP_cont, -1, (255, 0, 0), 2)
        cv2.drawContours(comb_Temp, CHcont, -1, (255, 0, 255), 2)
        cv2.drawContours(comb_Temp, BP_cont, -1, (255, 0, 0), 2)
        cv2.drawContours(comb_Temp, AR_cont, -1, (255, 255, 0), 2)
        cv2.drawContours(CHTmask, CHcont, -1, (255, 0, 0), cv2.FILLED)
        cv2.drawContours(BPTmask, BP_cont, -1, (255, 0, 0), cv2.FILLED)

        Bo_CHmask = CHTmask.astype(np.bool_)
        Bo_BPmask = BPTmask.astype(np.bool_)
        Bo_BGmask = BGTmask.astype(np.bool_)
        Bo_ARmask1 = ARTmask1.astype(np.bool_)
        In_CHmask = np.invert(Bo_CHmask)
        In_BGmask = np.invert(Bo_BGmask)
        In_BPmask = np.invert(Bo_BPmask)
        In_ARmask1 = np.invert(Bo_ARmask1)

        # masking
        CH_masked_sun = ms.array(TempMap, mask=In_CHmask)
        BP_masked_sun = ms.array(TempMap, mask=In_BPmask)
        BG_masked_sun = ms.array(TempMap, mask=Bo_BGmask)
        AR_masked1_sun = ms.array(TempMap, mask=In_ARmask1)

        bgimage = ms.array(TempMap, mask=In_BGmask)
        BGI = ((bgimage * 0).data).astype(np.uint8)
        BPI = ((BP_masked_sun * 0).data).astype(np.uint8)
        BP_masked1_sun = ms.array(TempMap, mask=In_BPmask)
        bpsum1 = BP_masked1_sun.sum()

        # Numebrs
        no_of_AR = (len(AR_cont))
        no_of_BP = (len(BP_cont))
        no_of_CH = len(CHcont)
        # print (no_of_BP, no_of_AR)
        n_AR.append(no_of_AR)
        n_BP.append(no_of_BP)
        n_CH.append(no_of_CH)

        # Totel area
        cha = np.count_nonzero(CHTmask)
        bpa = np.count_nonzero(BPTmask)
        ara = np.count_nonzero(ARTmask1)
        A5 = np.count_nonzero(TempMap)  # fuldisk size
        crc_area= np.pi*R*R
        c_rat=(A5/crc_area)*100 #convertion ratio
        A4 = np.count_nonzero(In_BGmask)  # Bg + oudside disk
        A6 = size1[0] * size1[1]
        A7 = A6 - A5
        A8 = A4 - A7#bg count
        #no need to correct area since it does itself since we are dealing with average values
        '''
        if size1[0] == 2048: #area corrected for big size
            cha = cha / 2
            bpa = bpa / 2
            ara = ara / 2
            A8=A8/2
        '''

        CHa.append(cha)
        BPa.append(bpa)
        ARa.append(ara)
        BGa.append(A8)
        FDa.append(A5)
        conv_rat.append(c_rat)

        # Totel Temperature
        bpsum = BP_masked_sun.sum()
        csum = CH_masked_sun.sum()
        bgsum = BG_masked_sun.sum()
        Arsum = AR_masked1_sun.sum()
        tsum = TempMap.sum()
        '''
        if size1[0]==2048:
            bpsum = BP_masked_sun.sum()/2
            csum = CH_masked_sun.sum()/2
            bgsum = BG_masked_sun.sum()/2
            Arsum = AR_masked1_sun.sum()/2
            tsum = TempMap.sum()/2
        '''

        # Totel=csum+bgsum+bpsum+Arsum
        XBP = bpsum  ##

        # Average temperature
        ar_temp = Arsum / ara
        ch_temp = csum / cha
        bp_temp = bpsum / bpa
        bg_temp = bgsum / A8
        FD_temp = tsum / A5
        Totel = ar_temp + bg_temp + bp_temp + ch_temp

        IR_data.append(IRsum)#area not corrected, will be done in plotter
        CHarray.append(ch_temp)
        BParray.append(bp_temp)
        ARarray.append(ar_temp)
        BGarray.append(bg_temp)
        FDarray.append(FD_temp)
        CHi.append((((ch_temp) / Totel) * 100))
        BPi.append((((bp_temp) / Totel) * 100))
        ARi.append((((ar_temp) / Totel) * 100))
        BGi.append((((bg_temp) / Totel) * 100))


        dname=dob_str+' | '+str(DOB2)
        #print(dname)
        imageio.imwrite('TempMaps/T_map_{}.jpg'.format(dname), comb_Temp[::-1])
        #imageio.imwrite('imgs2/IR_{}.jpg'.format(dob_str), Org_img[::-1])
        imageio.imwrite('IRmaps/COMB_img{}.jpg'.format(dob_str), IRz[::-1])
        f = open('Temp_data.dat', 'a')
        np.savetxt('Temp_data.dat', np.c_[
            CHarray, BParray, ARarray, BGarray, CHi, BPi, ARi, BGi, n_BP, n_AR, n_CH, l_DOB, CHa, BPa, ARa, BGa, FDarray, IR_data,chIRsum,bpIRsum,ArIRsum,bgIRsum,conv_rat,shapeArry],
                   fmt='%11.5f',
                   header=' CH Int,   XBP Int,     AR Int,  Background | %intensity-CH   XBP          AR         BG		nXBP		nAR		nCH		l_DOB	    CHa       BPa       ARa        BGa     FD-int    Int-rat     chIR    xbpIR    arIR    bgIR   conversion%     shape')
        f.close()
        tempstopTime = timeit.default_timer()
        filetime = tempstopTime - prvTime
        prvTime=tempstopTime

stopTime = timeit.default_timer()
runtime = (stopTime - startTime)
TotTime=runtime/3600 #in Hours

print('')
print('......COMPLEETED.....')
print('Time taken', TotTime,'Avg. Time per image', runtime/Length)
print('--------------------------------------------------')