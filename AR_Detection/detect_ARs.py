# -*- coding: utf-8 -*-
"""Detect atmospheric rivers (ARs) from integrated water vapor transports
(IVTs) using Top-hat by reconstruction (THR) algorithm.

In this script, read in some input data, detect ARs, and the detection results
are yielded once for each time step, and the results are saved to disk as long
as they are computed.

See also detect_ARs.py, where the detection results are collected and saved
to disk in one go.

# Input data

1. uflux, vflux:

    Instantaneous vertically integrated moisture flux, in kg/(m*s).
    u-moisture flux component's standard_name:
        "eastward_atmosphere_water_transport_across_unit_distance".
    v-moisture flux component's standard_name:
        "northward_atmosphere_water_transport_across_unit_distance".

    Data should be formatted into 4D (time, singleton_z, latitude, longitude),
    or 3D (time, latitude, longitude).

2. IVT:

    This is the magnitude of vertically integrated moisture flux, i.e.
    IVT^2 = uflux^2 + vflux^2.

3. THR results:

    Instantaneous THR reconstruction and anomalies of IVT, in kg/(m*s).
    This is the outcome of the THR process. See compute_thr_singlefile.py,
    compute_thr_multifile.py for scripts to perform this process.
    The THR algorithm is implemented in ipart.thr.THR.

uflux, vflux, IVT, and THR result data should be formatted into 4D
(time, singleton_z, latitude, longitude) or 3D (time, latitude, longitude),
and have compatible shapes.

The user also needs to provide time, latitude and longitude axes metadata.
In this script, these data are read in from the netCDF files using the
netCDF4 package. If you are using some other package, e.g.
CDAT, xarray, iris or something else, please adjust the relevant code
accordingly.

# Domain

Take only northern hemisphere, shift longitude to 80 E.

# Output data

1. labels:

    Labels of detected ARs. This is a 3D ndarray with dimension
    (time, lat, lon). At each time point, a unique integer label is assign
    to each detected AR, and the AR region is filled with the label value in
    the (lat, lon) map.

2. angles:

    Orientation differences between AR axes and horizontal moisture fluxes,
    measured in degrees.

3. crossfluxes:

    Cross-sectional moisture fluxes in all ARs, in kg/(m*s), computed as
    the projection of the total moisture flux onto the local AR axis.

labels, angles and crossfluxes have the same dimension and are saved into a
netCDF file.

4. result_df:

    Table of detected AR records. The columns of the table includes:

        id, time, centroid_x, centroid_y, axis_x, axis_y, ... etc.

    This table is saved to a .csv file.

5. AR detection result plots (optional):

    If set `PLOT=True`, will also plot out the IVT, THR reconstruction and
    THR anomaly distributions at each time point when any AR is detected.
    The boundary of all detected ARs are also marked out.

Author: guangzhi XU (xugzhi1987@gmail.com)
Update time: 2020-07-22 10:13:31.
"""


from __future__ import print_function

#######################################################################
#                               Globals                               #
#######################################################################
#------------------Output folder------------------
OUTPUTDIR='/home/jingzhao/wangbp/era5/arinfo_asia_250/'


PLOT=True          # create maps of found ARs or not

LAT1=-90; LAT2=90      # degree, latitude domain
SHIFT_LON=0          # degree, shift left bound to longitude.

PARAM_DICT={
    # kg/m/s, define AR candidates as regions >= than this anomalous ivt.
    # If None is given, compute a threshold based on anomalous ivt data. See
    # the docstring of ipart.AR_detector.determineThresLow() for details.
    'thres_low' : 250,
    # km^2, drop AR candidates smaller than this area.
    'min_area': 50*1e4,
    # km^2, drop AR candidates larger than this area.
    'max_area': 1800*1e4,
    # float, min length/width ratio.
    'min_LW': 0,
    # degree, exclude systems whose centroids are lower than this latitude.
    # NOTE this is the absolute latitude for both NH and SH. For SH, systems
    # with centroid latitude north of -20 will be excluded.
    'min_lat': 0,
    # degree, exclude systems whose centroids are higher than this latitude.
    # NOTE this is the absolute latitude for both NH and SH. For SH, systems
    # with centroid latitude south of -80 will be excluded.
    'max_lat': 80,
    # km, ARs shorter than this length is treated as relaxed.
    'min_length': 1000,
    # km, ARs shorter than this length is discarded.
    'min_length_hard': 800,     #1500
    # degree lat/lon, error when simplifying axis using rdp algorithm.
    'rdp_thres': 2,
    # grids. Remove small holes in AR contour.
    'fill_radius': None,
    # do peak partition or not, used to separate systems that are merged
    # together with an outer contour.
    'single_dome': True,
    # max prominence/height ratio of a local peak. Only used when single_dome=True
    'max_ph_ratio': 0.7,
    # minimal proportion of flux component in a direction to total flux to
    # allow edge building in that direction
    'edge_eps': 0,  # IMPORTANT: set to 0 when no uq and vq
    # bool, if True, treat the data as zonally cyclic (e.g. entire hemisphere
    # or global). ARs covering regions across the longitude bounds will be
    # correctly treated as one. If your data is not zonally cyclic, or a zonal
    # shift of the data can put the domain of interest to the center, consider
    # doing the shift and setting this to False, as it will save computations.
    'zonal_cyclic': True,
    }



#--------Import modules-------------------------
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import cartopy.crs as ccrs
from netCDF4 import Dataset
import time
import glob
import multiprocessing

from ipart.utils import funcs
from ipart.utils import plot
from ipart.AR_detector import plotAR, findARsGen



def main(IVT_FILE_NAME):

    YEAR = IVT_FILE_NAME[-17:-13]
    MONTH = IVT_FILE_NAME[-13:-11]
    #TIME_START='%d-01-01 00:00:00' %YEAR
    #TIME_END='%d-02-03 18:00:00' %YEAR

    LABEL_FILE_OUT_NAME='ar_label-angle-flux_%s%s.nc' %(YEAR,MONTH)
    RECORD_FILE_OUT_NAME='ar_records_%s%s.csv' %(YEAR,MONTH)

    #-----------Read in flux data----------------------
    #quNV=funcs.readNC(UQ_FILE_NAME, UQ_VAR)
    #qvNV=funcs.readNC(VQ_FILE_NAME, VQ_VAR)

    #-----------------Shift longitude-----------------
    #quNV=quNV.shiftLon(SHIFT_LON)
    #qvNV=qvNV.shiftLon(SHIFT_LON)

    #-------------------Read in ivt and THR results-------------------
    #��dataset��
    ivtNV=funcs.readNC(IVT_FILE_NAME, 'ivt')
    ivtrecNV=funcs.readNC(IVT_FILE_NAME, 'ivt_rec')
    ivtanoNV=funcs.readNC(IVT_FILE_NAME, 'ivt_ano')

    #--------------------Slice data--------------------
    #quNV=quNV.sliceData(TIME_START, TIME_END).sliceData(LAT1, LAT2, axis=2).squeeze()
    #qvNV=qvNV.sliceData(TIME_START, TIME_END).sliceData(LAT1, LAT2, axis=2).squeeze()
    #ivtNV=ivtNV.sliceData(TIME_START, TIME_END).sliceData(LAT1, LAT2, axis=1).squeeze()
    #ivtrecNV=ivtrecNV.sliceData(TIME_START, TIME_END).sliceData(LAT1, LAT2, axis=1).squeeze()
    #ivtanoNV=ivtanoNV.sliceData(TIME_START, TIME_END).sliceData(LAT1, LAT2, axis=1).squeeze()
    ivtNV=ivtNV.sliceData(LAT1, LAT2, axis=1).squeeze()
    ivtrecNV=ivtrecNV.sliceData(LAT1, LAT2, axis=1).squeeze()
    ivtanoNV=ivtanoNV.sliceData(LAT1, LAT2, axis=1).squeeze()

    #---------------Create dummy uq, vq---------------
    quNV=funcs.NCVAR(np.zeros(ivtNV.shape), 'uflux', ivtNV.axislist, ivtNV.attributes)
    qvNV=funcs.NCVAR(np.zeros(ivtNV.shape), 'vflux', ivtNV.axislist, ivtNV.attributes)

    #--------------------Data shape check--------------------
    if np.ndim(quNV)!=3 or np.ndim(qvNV)!=3:
        raise Exception("<qu> and <qv> should be 3D data.")
    if quNV.shape!=qvNV.shape or ivtNV.shape!=quNV.shape:
        raise Exception("Data shape dismatch: qu.shape=%s; qv.shape=%s; ivt.shape=%s"\
                %(quNV.shape, qvNV.shape, ivtNV.shape))

    #-----------------Get coordinates-----------------
    latax=quNV.getLatitude()
    lonax=quNV.getLongitude()
    timeax=ivtNV.getTime()

    #-----------------Prepare outputs-----------------
    if not os.path.exists(OUTPUTDIR):
        os.makedirs(OUTPUTDIR)

    if PLOT:
        plot_dir=os.path.join(OUTPUTDIR, 'plots_%s%s' %(YEAR,MONTH))
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

    # nc file to save AR location labels
    abpath_out=os.path.join(OUTPUTDIR, LABEL_FILE_OUT_NAME)
    print('\n### <detect_ARs2>: Saving output to:\n',abpath_out)
    ncfout=Dataset(abpath_out, 'w')

    # csv file to save AR record table
    abpath_out=os.path.join(OUTPUTDIR, RECORD_FILE_OUT_NAME)
    print('\n### <detect_ARs2>: Saving output to:\n',abpath_out)
    # Necessary: to remove ... in csv file
    if sys.version_info.major==2:
        np.set_printoptions(threshold=np.inf)
    elif sys.version_info.major==3:
        np.set_printoptions(threshold=sys.maxsize)

    with open(abpath_out, 'a') as dfout:

        #############################################################
        #                     Start processing                      #
        #############################################################
        finder_gen = findARsGen(ivtNV.data, ivtrecNV.data, ivtanoNV.data,
                quNV.data, qvNV.data, latax, lonax, times=timeax, **PARAM_DICT)
      
        next(finder_gen)  # prime the generator to prepare metadata

        for (tidx, timett, label, result_df) in finder_gen:

            #------------------Save record to csv file------------------
            result_df.to_csv(dfout, header=dfout.tell()==0, index=False)

            #-------------------Save labels to nc file-------------------
            funcs.saveNCDims(ncfout, label.axislist)
            funcs._saveNCVAR(ncfout, label, 'int')
            funcs._saveNCVAR(ncfout, angle)
            funcs._saveNCVAR(ncfout, cross)

            #-------------------Plot------------------------
            if PLOT:

                timett_str=str(timett)

                slab=ivtNV.data[tidx]
                #slabrec=ivtrecNV.data[tidx]
                slabano=ivtanoNV.data[tidx]

                plot_vars=[slab,slabano]
                titles=['IVT', 'IVT-250ano']
                iso=plot.Isofill(plot_vars,12,1,1,min_level=0,max_level=800)

                figure=plt.figure(figsize=(12,10),dpi=100)

                for jj in range(len(plot_vars)):
                    ax=figure.add_subplot(2,1,jj+1)
                    pobj=plot.plot2(plot_vars[jj],iso,ax,projection='cyl',
                            xarray=lonax, yarray=latax,
                            title='%s %s' %(timett_str, titles[jj]),
                            fix_aspect=False)

                    bmap=pobj.bmap
                    plotAR(result_df,ax,bmap,lonax)

                #----------------- Save plot------------
                plot_save_name='ar_%s' %(timett_str)
                plot_save_name=os.path.join(plot_dir,plot_save_name)
                print('\n# <detect_ARs2>: Save figure to', plot_save_name)
                figure.savefig(plot_save_name+'.png',dpi=100,bbox_inches='tight')

                plt.close('all')

    #----------------Close the nc file----------------
    ncfout.close()




#-------------Main---------------------------------
if __name__=='__main__':

    # ����ʼ
    time_start=time.time()
    num = xxx
    
    
    for mon in range(5,10):
        #-----------------ivt summer data file-----------------
        SOURCEDIR='/home/jingzhao/wangbp/era5/ivt250_summer_mon/'
        RECORD_FILE='era5_ivt_*' + str(mon).zfill(2) + '-ano-250.nc'  
        abpath_in=os.path.join(SOURCEDIR,RECORD_FILE)
        if mon==5:
          IVT_FILE_NAME = glob.glob(abpath_in)
        else:
          IVT_FILE_NAME_mon = glob.glob(abpath_in)
          IVT_FILE_NAME = IVT_FILE_NAME + IVT_FILE_NAME_mon
      
    IVT_FILE_NAME.sort(key = lambda x: x[-17:-11])
        
    #IVT_FILE_NAME = IVT_FILE_NAME[-24:]
    if num<6:
        IVT_FILE_NAME = IVT_FILE_NAME[num*24:(num+1)*24]
    else:
        IVT_FILE_NAME = IVT_FILE_NAME[num*24:]
    print('-----IVT_FILE_NAME:', IVT_FILE_NAME)
    pool=multiprocessing.Pool(processes=24)  # ����һ���̳߳�
    pool.map(main, IVT_FILE_NAME)  # ���̳߳������߳�
    pool.close()  # �ر��̳߳أ����ٽ����߳�
    pool.join()  # �ȴ��̳߳����߳�ȫ��ִ����
    
    # �������
    time_end=time.time()
    print('totally cost',time_end-time_start,'S')
            





