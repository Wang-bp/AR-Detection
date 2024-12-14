'''Link ARs found in consecutive time slices to form tracks

# Input data

AR records at individual time steps. This should be the output from
detect_ARs.py or detect_ARs_generator_version.py in csv format.

# Output data

1. trackdf:

    Table of AR tracks. Columns of this table are the same as the input
    AR record table, plus one additional column containing the track id.

    This table is saved to a .csv file.

2. track movement and Hausdorff linkage schematic plots (optional):

    If set `PLOT=True`, will also plot out the track movements.
    If set `SCHEMATIC=True`, will also plot out the schematic plots showing
    the tracking procedure at consecutive time points (t=t and t=t+1).


Author: guangzhi XU (xugzhi1987@gmail.com)
Update time: 2020-06-05 23:22:25.
'''

from __future__ import print_function

#--------------Globals------------------------------------------
RECORD_FILE_IN_NAME='/home/jingzhao/wangbp/era5/AR_TEST/THR/arinfo_250/ar_records_*00.csv'
OUTPUTDIR='/home/jingzhao/wangbp/era5/AR_TEST/THR/tracks/'


PLOT=True         # plot track movements or not
SCHEMATIC=False   # plot schematic or not
LAT1=-90; LAT2=90; LON1=0; LON2=360;         # domain to plot

# Int, hours, gap allowed to link 2 records. Should be the time resolution of
# the data.
TIME_GAP_ALLOW=6

# tracking scheme. 'simple': all tracks are simple paths.
# 'full': use the network scheme, tracks are connected by their joint points.
TRACK_SCHEME='simple'  # 'simple' | 'full'

# int, max Hausdorff distance in km to define a neighborhood relationship
MAX_DIST_ALLOW=1200  # km

# int, min duration in hrs to keep a track.
MIN_DURATION=24

# int, min number of non-relaxed records in a track to keep a track.
MIN_NONRELAX=2




#--------Import modules-------------------------
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import cartopy.crs as ccrs
import time
import glob
import multiprocessing

from ipart.AR_tracer import readCSVRecord, trackARs, filterTracks, plotAR



def main(filename):

    name = filename[-18:-4]
    RECORD_FILE_OUT_NAME='ar_tracks_'+name+'.csv'

    #-----------Read in records---------------------
    print('\n### <trace_ARs>: Read in file:\n', filename)
    ardf=readCSVRecord(filename)

    if not os.path.exists(OUTPUTDIR):
        os.makedirs(OUTPUTDIR)

    if PLOT or SCHEMATIC:
        plot_dir=os.path.join(OUTPUTDIR, 'plots_%s' %name)
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

    #-------------------Track ars-------------------
    track_list=trackARs(ardf, TIME_GAP_ALLOW, MAX_DIST_ALLOW,
            track_scheme=TRACK_SCHEME, isplot=SCHEMATIC, plot_dir=plot_dir)

    #------------------Filter tracks------------------
    track_list=filterTracks(track_list, MIN_DURATION, MIN_NONRELAX)

    #-------------------Save output-------------------
    latax=np.arange(LAT1, LAT2)
    lonax=np.arange(LON1, LON2)

    for ii in range(len(track_list)):
        tii=track_list[ii]
        trackidii='%d%d' %(tii.data.loc[0,'time'].year, ii+1)
        tii.data.loc[:,'trackid']=trackidii
        for jj in range(len(tii.data.loc[0:,'time'])):
          tii.data.loc[jj,'year']=tii.data.loc[jj,'time'].year
          tii.data.loc[jj,'mon']=tii.data.loc[jj,'time'].month
          tii.data.loc[jj,'day']=tii.data.loc[jj,'time'].day
          tii.data.loc[jj,'hour']=tii.data.loc[jj,'time'].hour

        tii.trackid=trackidii

        if ii==0:
            trackdf=tii.data
        else:
            trackdf=pd.concat([trackdf,tii.data],ignore_index=True)

        if PLOT:
            figure=plt.figure(figsize=(12,6),dpi=100)
            ax=figure.add_subplot(111)
            plotAR(tii,latax,lonax,ax,True)

            #----------------- Save plot------------
            plot_save_name='ar_track_%s' %trackidii
            plot_save_name=os.path.join(plot_dir,plot_save_name)
            print('\n# <trace_ARs>: Save figure to', plot_save_name)
            figure.savefig(plot_save_name+'.png',dpi=100,bbox_inches='tight')

            plt.close(figure)

    #--------Save------------------------------------
    abpath_out=os.path.join(OUTPUTDIR, RECORD_FILE_OUT_NAME)
    print('\n### <trace_ARs>: Saving output to:\n',abpath_out)
    if sys.version_info.major==2:
        np.set_printoptions(threshold=np.inf)
    elif sys.version_info.major==3:
        np.set_printoptions(threshold=sys.maxsize)
    trackdf.to_csv(abpath_out,index=False)    


#-------------Main---------------------------------
if __name__=='__main__':

    # 程序开始
    time_start=time.time()
    
    filelist=glob.glob(RECORD_FILE_IN_NAME)
    filelist.sort()



    filename = []
    for file in filelist:
        name = 'ar_tracks_'+file[-18:]
        abpath = os.path.join(OUTPUTDIR, name)
        if os.path.exists(abpath):
            continue
        else:
            filename.append(file)
    

    pool=multiprocessing.Pool(processes=27)  # 创建一个线程池
    pool.map(main, filename)  # 往线程池中填线程
    pool.close()  # 关闭线程池，不再接受线程
    pool.join()  # 等待线程池中线程全部执行完
    
    # 程序结束
    time_end=time.time()
    print('totally cost',time_end-time_start,'S')


















