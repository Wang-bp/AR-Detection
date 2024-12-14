'''Compute Integrated Vapor Transport (IVT) from u- and v- vapor flux components.

u-vapor flux component's standard_name:
    "eastward_atmosphere_water_transport_across_unit_distance".
v-vapor flux component's standard_name:
    "northward_atmosphere_water_transport_across_unit_distance".

IVT = sqrt(uflux^2 + vflux^2), in kg/m/s.

Author: guangzhi XU (xugzhi1987@gmail.com)
Update time: 2020-07-22 10:08:31.
'''

#--------------Globals------------------------------------------
OUTPUTFILE='/home/jingzhao/wangbp/era5/ivt_summer_mon/'

#--------Import modules-------------------------
import numpy as np
from ipart.utils import funcs
from ipart.utils import plot

import os
import time
import glob
import multiprocessing



def main(UFLUX_FILE):

    VFLUX_FILE = UFLUX_FILE
    #-----------Read in data----------------------
    ufluxNV=funcs.readNC(UFLUX_FILE, UFLUX_VARID)
    vfluxNV=funcs.readNC(VFLUX_FILE, VFLUX_VARID)

    ivtdata=np.ma.sqrt(ufluxNV.data**2+vfluxNV.data**2)
    ivtNV=funcs.NCVAR(ivtdata, 'ivt', ufluxNV.axislist, {'name': 'ivt',
        'long_name': 'integrated vapor transport (IVT)',
        'standard_name': 'integrated vapor transport (IVT)',
        'title': 'integrated vapor transport (IVT)',
        'units': getattr(ufluxNV, 'units', '')})

    #--------Save------------------------------------
    IVT_FILE='era5_ivt_%s.nc'   %UFLUX_FILE[-9:-3]
    abpath_out=os.path.join(OUTPUTFILE,IVT_FILE)
    print('\n### <compute_ivt>: Saving output to:\n',abpath_out)
    funcs.saveNC(abpath_out, ivtNV, 'w')


#-------------Main---------------------------------
if __name__=='__main__':

    # 程序开始
    time_start=time.time()

    #-----------uflux----------------------
    UFLUX_VARID='p71.162'
    VFLUX_VARID='p72.162'
    '''
    SOURCEDIR='/home/jingzhao/wangbp/era5/additional/'
    RECORD_FILE='era5_ivt_*.nc'
    abpath_in=os.path.join(SOURCEDIR,RECORD_FILE)
    UFLUX_FILE = glob.glob(abpath_in)
    UFLUX_FILE.sort()
    
    pool=multiprocessing.Pool(processes=24)  # 创建一个线程池
    pool.map(main, UFLUX_FILE)  # 往线程池中填线程
    pool.close()  # 关闭线程池，不再接受线程
    pool.join()  # 等待线程池中线程全部执行完
    '''
    fn = "/home/jingzhao/xiaohui/era5/era5_vq_200609.nc"
    for ii in range(2006,2010):
        fn = "/home/jingzhao/xiaohui/era5/era5_vq_" + str(ii) + "09.nc"
        main(fn)
    
    # 程序结束
    time_end=time.time()
    print('totally cost',time_end-time_start,'S')













