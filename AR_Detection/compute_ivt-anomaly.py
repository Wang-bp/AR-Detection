# -*- coding: utf-8 -*-
'''Compute 3d THR on IVT. Process a single data file.

Input:
    IVT (integrated Vapor Transport) in netCDF format.

    Data are assumed to be in the format: (time, level, latitude, longitude)
    dimensions (level dimension is optional, if present, should be a
    singleton axis of length 1).

    NOTE: data should have proper time, latitude and longitude axes!

Usage:

    Change global parameters in the Globals section to point to the storage
    location of IVT data, and specificy an output folder to save results.

    Specify the latitudinal domain in LAT1, LAT2.

    The KERNEL parameter specifies the $t$ and $s$ parameters of the
    structuring element size. See paper for more details.
    $t$ is in number time steps, $s$ is number of grid cells.

    SHIFT_LON shifts the longitude by a given degree of longitudes, so
    that the Pacific and Atlantic basins can be centered.

    Run the script as:
        ```
        python compute_thr_singlefile.py
        ```

Author: guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
Update time: 2019-12-06 22:55:33.
'''

from __future__ import print_function


#--------------Globals------------------------------------------
#-----------IVT data----------------------
VARIN='ivt'          # data id in nc file

LAT1=-90; LAT2=90      # degree, latitude domain
LON1=0; LON2=359.75;

#-------Structuring element for erosion (E)-------
KERNEL=[4]   # half length of time (time steps), and half length of spatial (number of grids)

SHIFT_LON=0  # shift longitudinally to center Pacific(0) and Altantic(80)

#------------------Output folder------------------
OUTPUTDIR='/home/jingzhao/wangbp/era5/ivt250_summer_mon/'




#--------Import modules-------------------------
import os
import cdms2 as cdms
import MV2 as MV
import numpy as np
from netCDF4 import Dataset
from skimage import morphology
from cdms2.selectors import Selector
from utils import funcs
import time
import glob
import multiprocessing


def readVar(abpath_in, varid):
    '''Read in netcdf variable

    Args:
        abpath_in (str): absolute file path to nc file.
        varid (str): id of variable to read.
    '''

    print('\n# <readVar>: Read in file:\n',abpath_in)
    fin=cdms.open(abpath_in,'r')
    var=fin(varid)

    #-----------------Try get latitude-----------------
    try:
        latax=var.getLatitude()
        if latax is None:
            raise Exception("latax is None")
    except:
        try:
            lat=fin('lat')[:,0]
        except:
            raise Exception("Can't find latitude axis in data.")
        else:
            latax=cdms.createAxis(lat)
            latax.designateLatitude()
            latax.id='y'
            latax.name='latitude'
            latax.unit='degree_north'

    #-----------------Try get longitude-----------------
    try:
        lonax=var.getLongitude()
        if lonax is None:
            raise Exception("lonax is None")
    except:
        try:
            lon=fin('lon')[0,:]
        except:
            raise Exception("Can't find longitude axis in data.")
        else:
            lonax=cdms.createAxis(lon)
            lonax.designateLongitude()
            lonax.id='x'
            lonax.name='longitude'
            lonax.unit='degree_east'

    axislist=var.getAxisList()
    axislist[-2]=latax
    axislist[-1]=lonax

    var.setAxisList(axislist)
    fin.close()

    return var


def filterData_250(ivt,ivt_mean):

    ndim=np.ndim(ivt)
    ivt=ivt(squeeze=1)

    ivtrec = ivt - ivt_mean
    ivtano = ivt - ivt_mean

    if ndim==4:
        levax=cdms.createAxis([0,])
        levax.designateLevel()
        levax.id='z'
        levax.name='level'
        levax.units=''

        ivt=funcs.addExtraAxis(ivt,levax,1)
        ivtano=funcs.addExtraAxis(ivtano,levax,1)
        ivtrec=funcs.addExtraAxis(ivtrec,levax,1)

    axislist=ivt.getAxisList()
    ivtano.setAxisList(axislist)
    ivtrec.setAxisList(axislist)

    ivtrec.id='ivt_rec'
    ivtrec.long_name='Integreated moisture transport, anomaly wrt minimal reconstruction (same as ivt_ano)'
    ivtrec.standard_name=ivtrec.long_name
    ivtrec.title=ivtrec.long_name
    ivtrec.units=ivt.units

    ivtano.id='ivt_ano'
    ivtano.long_name='Integreated moisture transport, anomaly wrt minimal reconstruction'
    ivtano.standard_name=ivtano.long_name
    ivtano.title=ivtano.long_name
    ivtano.units=ivt.units

    return ivt, ivtrec, ivtano


def filterData(ivt,kernel,verbose=True):

    ndim=np.ndim(ivt)
    ivt=ivt(squeeze=1)

    #-------------------3d ellipsoid-------------------
    ele=funcs.get3DEllipse(*kernel)
    #dt=kernel[0] # half length in time dimesion

    #################### use a cube to speed up ##############
    if kernel[0]>=10 or kernel[1]>=6:
        ele=np.ones(ele.shape)
    ##########################################################

    if verbose:
        print('\n# <filterData>: Computing erosion ...')

    lm=morphology.erosion(ivt.data,selem=ele)

    if verbose:
        print('\n# <filterData>: Computing reconstruction ...')

    ivtrec=morphology.reconstruction(lm,ivt,method='dilation')
    ivtrec=MV.array(ivtrec)
    ivtano=ivt-ivtrec
    ivtano=MV.array(ivtano)
    ivtrec=MV.array(ivtrec)

    if ndim==4:
        levax=cdms.createAxis([0,])
        levax.designateLevel()
        levax.id='z'
        levax.name='level'
        levax.units=''

        ivt=funcs.addExtraAxis(ivt,levax,1)
        ivtano=funcs.addExtraAxis(ivtano,levax,1)
        ivtrec=funcs.addExtraAxis(ivtrec,levax,1)

    axislist=ivt.getAxisList()
    ivtano.setAxisList(axislist)
    ivtrec.setAxisList(axislist)

    ivtrec.id='ivt_rec'
    ivtrec.long_name='Integrated moisture transport, minimal reconstruction'
    ivtrec.standard_name=ivtrec.long_name
    ivtrec.title=ivtrec.long_name
    ivtrec.units=ivt.units

    ivtano.id='ivt_ano'
    ivtano.long_name='Integreated moisture transport, anomaly wrt minimal reconstruction'
    ivtano.standard_name=ivtano.long_name
    ivtano.title=ivtano.long_name
    ivtano.units=ivt.units

    return ivt, ivtrec, ivtano


def shift_lon_varid(varid,lon,SHIFT_LON):
    #转某个度数

    lon = lon.tolist()
    shift_lon_index = lon.index(SHIFT_LON)
    varid_right = varid[:,shift_lon_index:]
    varid_left = varid[:,:shift_lon_index]
    varid_new = np.concatenate((varid_right,varid_left),axis=1)

    return varid_new


def main(IVT_FILE):

    #-----------Read in data----------------------
    var=readVar(IVT_FILE, 'ivt')

    #-----------------Shift longitude-----------------
    var=var(latitude=(LAT1, LAT2))
    var=var(longitude=(LON1, LON2))
    var=var(longitude=(SHIFT_LON+LON1,SHIFT_LON+LON2))

    #----------------------Do THR----------------------
    #ivt, ivtrec, ivtano=filterData(var, KERNEL)
    ivt, ivtrec, ivtano=filterData_250(var, IVT_mean)

    #--------Save------------------------------------
    fname=os.path.split(IVT_FILE)[1]
    #file_out_name='%s-minimal-rec-ano-kernel-t%d-s%d.nc'\
    #        %(os.path.splitext(fname)[0], KERNEL[0], KERNEL[1])
    file_out_name='%s-ano-250.nc' %(os.path.splitext(fname)[0])

    abpath_out=os.path.join(OUTPUTDIR,file_out_name)
    print('\n### <testrotatingfilter>: Saving output to:\n',abpath_out)
    fout=cdms.open(abpath_out,'w')
    fout.write(ivt,typecode='f')
    fout.write(ivtrec,typecode='f')
    fout.write(ivtano,typecode='f')
    fout.close()



#-------------Main---------------------------------
if __name__=='__main__':

    # 程序开始
    time_start=time.time()

    if not os.path.exists(OUTPUTDIR):
        os.makedirs(OUTPUTDIR)

    #-----------Read in lon & lat---------------------
    #abpath_in='/home/jingzhao/wangbp/era5/data_summer_mon/era5_ivt_197904.nc'
    #abpath_in=Dataset(abpath_in)
    #lat=(abpath_in.variables['latitude'][:])
    #lon=(abpath_in.variables['longitude'][:])
    #abpath_in.close()

    #------------climatological mean-------------
    IVT_mean = np.loadtxt("/home/jingzhao/wangbp/era5/ivt_summer_mon/ivt_seasonal_mean_1979-2005.txt")
    #IVT_mean = shift_lon_varid(IVT_mean,lon,SHIFT_LON)

    #-----------------ivt reconstruction and anomalies-----------------
    #SOURCEDIR='/home/jingzhao/wangbp/era5/ivt_summer_mon/'
    #RECORD_FILE='era5_ivt_20*.nc'
    #abpath_in=os.path.join(SOURCEDIR,RECORD_FILE)
    #IVT_FILE = glob.glob(abpath_in)
    #IVT_FILE.sort()
    #IVT_FILE = IVT_FILE[36:]
    
    IVT_FILE = ['/home/jingzhao/wangbp/era5/ivt_summer_mon/era5_ivt_201108.nc','/home/jingzhao/wangbp/era5/ivt_summer_mon/era5_ivt_201508.nc']

    filename = []
    for file in IVT_FILE:
        name = file[-18:-3]+'-ano-250.nc'  
        abpath = os.path.join(OUTPUTDIR, name)
        if os.path.exists(abpath):
            continue
        else:
            filename.append(file)

    pool=multiprocessing.Pool(processes=2)  # 创建一个线程池
    pool.map(main, filename)  # 往线程池中填线程
    pool.close()  # 关闭线程池，不再接受线程
    pool.join()  # 等待线程池中线程全部执行完

    # 程序结束
    time_end=time.time()
    print('totally cost',time_end-time_start,'S')




