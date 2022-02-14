#!/usr/bin/env python
# coding: utf-8

# #  Study of signal correlation in raft 
# Goal : This Notebook will generate plots allowing to perform a diagnostic on the signal Covariance between pixels inside a ccd channels 
#        This is a complement to CTE and full well studies  in BOT data .
# 
# Author : P.Antilogus
# 
# Version : 9th February 2022 01:50 
# 
# New : To be run at CCIN2P3 
################################## 

from eochar.bot_frame_op import *
from eochar.frame_cte_noise import *

# system imports
import os
import time
from sys import exit 
import glob

#  Specific package (display , pyfits ..) 
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
# Firefly client imports   
#from IPython.display import IFrame
import matplotlib
matplotlib.rcParams['figure.dpi'] = 120

# LSST stack imports
from lsst.daf.persistence import Butler
import lsst.afw.display as afwDisplay
from lsst.ip.isr import IsrTask
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath

# do we use the eotest resultes for the gain
display_in_electron=True
#
if display_in_electron :
    try : 
        from get_EO_analysis_results import get_EO_analysis_results
        from get_EO_analysis_files import get_EO_analysis_files
        from exploreFocalPlane import exploreFocalPlane
        from exploreRaft import exploreRaft
        from eTraveler.clientAPI.connection import Connection
        g = get_EO_analysis_results(db=DB_key[data_to_select['Data_location']])
        eotest_db=True
    except:
        print('No access to eotest DB , so the gain of all devices will be set to 1 ( 1 ADU counted as 1 e- )')
        eotest_db=False
else :
    eotest_db=False
    
# CONFIGURATION FOR THE CURRENT EXECUTION  ========================================
# ---- raft and associated run ============ To be updated if needed 
# 'raft' :  list of raft 
# 'run' :  list of run , ex : ['9876','9874']  (remark : run is a string )   for all run of this raft : '*' 
# 'sensor' :  list of sensor , ex ['S00','S01','S02'] 
# amplifier : list of amplifier , ex [1,2]  , for 1 to 16 you can use [-1] instead 
run_all=['13144']
#,'12849','12851','12857','12860']
#,'12606']
#
raft_itl=['R01', 'R02', 'R03', 'R10', 'R20', 'R41', 'R42', 'R43']
raft_e2v=['R11', 'R12', 'R13', 'R14', 'R21', 'R22', 'R23', 'R24', 'R30', 'R31', 'R32', 'R33', 'R34']
sensors_raft=['S00','S01','S02','S10','S11','S12','S20','S21','S22']
raft_corner=['R00','R04','R40','R44']
sensors_corner=['SG0','SG1','SW0','SW1']
sensors_8ch=['SW0','SW1']
segments=['10','11','12','13','14','15','16','17','07','06','05','04','03','02','01','00']

#
all_sensors={}
all_good_sensors={}
#
for raft in raft_itl+raft_e2v:
    all_sensors[raft]=sensors_raft
    all_good_sensors[raft]={}
    for sensor in sensors_raft:
        all_good_sensors[raft][sensor]=True
for raft in raft_corner:
    all_sensors[raft]=sensors_corner
    all_good_sensors[raft]={}
    for sensor in sensors_corner:
        all_good_sensors[raft][sensor]=True
#
amplifier=[-1]
#
#raft=raft_itl+raft_corner+raft_e2v
#raft=raft_e2v
#raft=raft_corner
raft=['R14']
#raft=['R10','R11','R12','R20','R21','R22','R30','R31','R32']
#raft=['R00']
##

# output file directory
output_data='/sps/lsst/users/tguillem/web/debug/'

print('Configuration arguments: ', str(sys.argv))
str_run_all = str(sys.argv[1])
str_raft = str(sys.argv[2])
str_all_sensors = str(sys.argv[3])
repo_path=str(sys.argv[4])
output_data=str(sys.argv[5])

run_all=[str_run_all]
raft=[str_raft]
all_sensors[str_raft]=[str_all_sensors]
#butler = Butler(repo_path)

#
# colors for each channel setup
cmap=plt.get_cmap('gist_ncar')
colors=[cmap(j)[:3] for j in np.linspace(0,1,17)]

number_of_pair_max=3000
run_visits=np.zeros((len(run_all),number_of_pair_max,2),dtype=np.int64)
time_visits=np.zeros((len(run_all),number_of_pair_max))
number_of_pair_per_run=np.zeros((len(run_all)),dtype=np.int16)
run=[]
irun=0
for run_cur in run_all :
    # activate the butler
    BOT_REPO_DIR = '/sps/lsst/groups/FocalPlane/SLAC/run5/butler/'+run_cur
    butler = Butler(BOT_REPO_DIR)
    # select only FLAT inages from PTC 
    visits=butler.queryMetadata('raw', ['visit','ExpTime','filter'], dataId={'run':run_cur,'testType': 'FLAT','imageType': 'FLAT'})
    print('We found ',len(visits),' flats in run ',run_cur)
    if len(visits) < 2 : 
        print('No enough good flat in run ',run_cur,' to proceed' )
        continue
    # select exposure without density filters (too faint for cov calculation )  
    visits_nnd=[(visits[i][0],visits[i][1]) for i in range(len(visits)) if not('ND_' in visits[i][2]) ]
    print('We found ',len(visits_nnd),' flats not using Neutral Density filters (too low flux) in run ',run_cur)
    if len(visits_nnd) < 2 : 
        print('No enough good flat in run ',run_cur,' to proceed' )
        continue
    #restrict to 20 visits to debug
    #visits_nnd = visits_nnd[:20]
    # order them by exposure time & visit  order
    visits_nd=np.array(visits_nnd, dtype=[('visit', np.int64), ('time', 'f')])
    visits_nd.sort(order=['time','visit'])
    #visits=visits_nd['visit']
    good_visits = []
    ptim=[]
    visit1=-1
    for ivisit in range(len(visits_nd)) :
        if visit1==-1 : 
            visit1=ivisit
            continue
        visit2=ivisit
        # to be keept the 2 exposure should one after the other in visit number and have the same exposure time
        if (visits_nd['visit'][visit2]-visits_nd['visit'][visit1]!=1) or  (np.abs(visits_nd['time'][visit2]-visits_nd['time'][visit1])>0.001 ) :
            print("Mismatch in sequence of exposure",visits_nd['visit'][visit1],visits_nd['time'][visit1],visits_nd['visit'][visit2],visits_nd['time'][visit2])
            visit1=ivisit
            continue
        else :
            good_visits.append((visits_nd['visit'][visit1],visits_nd['visit'][visit2]))
            ptim.append(visits_nd['time'][visit1])
            visit1=-1
    ##restrict to 10 good visits to debug        
    ##good_visits = good_visits[:10]
    ##ptim=ptim[:10]
    if len(ptim)==0 :
        print('No enough good flat pair in this run ', run_cur )
        continue
    run.append(run_cur)
    number_of_pair_per_run[irun]=len(ptim)
    run_visits[irun,:number_of_pair_per_run[irun],:]=np.array(good_visits)
    time_visits[irun,:number_of_pair_per_run[irun]]=np.array(ptim)
    print('Run=',run_cur,', time distribution of the ',number_of_pair_per_run[irun],' selected flat pair ',time_visits[irun,:number_of_pair_per_run[irun]])
    irun+=1
    #number_of_pair_per_run=2
#
number_of_raft=len(raft)
number_of_run=len(run)
# set the max value for the tables creation
number_of_sensor=len(sensors_raft)
#how many amplifiers 
if -1 in amplifier :
    amplifier_list=np.array(range(1,17))
else :
    amplifier_list=np.array(amplifier)
amplifier_list.sort()
number_of_amplifier=len(amplifier_list)
if number_of_run < 1 : 
    print('No data found')
    raise 


def find_pair(raft_cur,all_good_sensors):
    # compute the noise per image
    # Do all the fft plots if asked 
    # intialisation 
    number_of_raft_sensor=len(all_sensors[raft_cur])
    #
    all_file=np.zeros((number_of_run,number_of_raft_sensor,number_of_pair_max,2),dtype=np.object_)
    #
    first_ccd_run0=-1
    #
    for irun in range(len(run)) : 
        run_cur=run[irun]
        #
        BOT_REPO_DIR = '/sps/lsst/groups/FocalPlane/SLAC/run5/butler/'+run[irun]
        butler = Butler(BOT_REPO_DIR)
        #
        t0=0
        #
        for iccd in range(number_of_raft_sensor) : 
            #
            for it in range(number_of_pair_per_run[irun]) : 
                for ifile in range(2) : 
                    dataId2 = {'visit': int(run_visits[irun,it,ifile]), 'detectorName': all_sensors[raft_cur][iccd], 'raftName': raft_cur}
                    try : 
                        file=butler.get('raw_filename', dataId2)[0][:-3]
                    except : 
                        print('Error ??? file ',ifile,' of the pair ',it,' not found for : ',dataId2)
                        all_good_sensors[raft_cur][all_sensors[raft_cur][iccd]]=False
                        break
                    if irun==0 and first_ccd_run0<0: first_ccd_run0=iccd
                    all_file[irun,iccd,it,ifile]=file
    return  all_file , first_ccd_run0        


def cov_in_raft(raft_cur,all_file, first_ccd_run0 ,show=True):
    #
    if raft_cur in raft_corner :
        nb_amp_in_raft=48
    else :
        nb_amp_in_raft=144
    #
    number_of_raft_sensor=len(all_sensors[raft_cur])
    #  
    # we grant the image size to split in 5 in // and serial direction 
    line_split=5
    serial_split=5 
    # get the amplifier shape for this raft
    fits=pyfits.open(all_file[0,first_ccd_run0 ,0,0])
    first_line,first_p_over,first_col,first_s_over=image_area(fits)
    fits.close()
    #ncol=len(fits[1].data[100,:])
    #nline==len(fits[1].data[:,100])
    if first_p_over-first_line == 2000 :
        # itl 
        xf=first_col+19
        xl=first_s_over-20
        yf=first_line+50
        yl=first_p_over-20
    else :
        #e2v 
        xf=first_col+21
        xl=first_s_over-21
        yf=first_line+50
        yl=first_p_over-22

#
    yl1=yl+1
    xl1=xl+1
    ncol=xl-xf
    nline=yl-yf
    #
    flux=np.zeros((number_of_run,9,number_of_pair_max,16))
    #
    nflux=np.zeros((number_of_run,9,16),dtype=np.int16)
    # covar 10 and 01 measured in 5 image slices 
    covars=np.zeros((number_of_run,9,number_of_pair_max,16,serial_split))
    covarp=np.zeros((number_of_run,9,number_of_pair_max,16,line_split))
    # covar 10 and 01 measured on the full  image
    covarst=np.zeros((number_of_run,9,number_of_pair_max,16))
    covarpt=np.zeros((number_of_run,9,number_of_pair_max,16))
    # ratio of covar 10 and 01 
    covarps=np.zeros((number_of_run,9,number_of_pair_max,16))
    #
    for irun in range(number_of_run) :
        #
        for iccd in range(number_of_raft_sensor) :
#        for iccd in [6] :
            if not(all_good_sensors[raft_cur][all_sensors[raft_cur][iccd]]) :
                continue
            if ( all_sensors[raft_cur][iccd] in sensors_8ch ):
                number_of_amp_sensor=8
            else : 
                    number_of_amp_sensor=16
            for ipair in range(number_of_pair_per_run[irun]) :
                #
                fname=[all_file[irun,iccd,ipair,0],all_file[irun,iccd,ipair,1]]
                all=InFile(dirall=fname,Slow=False,verbose=False)
                #
                for k in range(number_of_amp_sensor) :
                    # did we pass the ultra saturation state : stop
                    if ipair >  nflux[irun,iccd,k] : 
                        # set it anyway , to have a simple way to estimate min-max flux per sensor
                        flux[irun,iccd,ipair,k]=flux[irun,iccd,ipair-1,k] 
                        continue
                    # select-truncate  the image for the channel ch_cur
                    Im0=np.median(all.all_file[0].Image[k][yf:yl1,xf:xl1])
                    Im1=np.median(all.all_file[1].Image[k][yf:yl1,xf:xl1])
                    ImDiff=((all.all_file[0].Image[k][yf:yl1,xf:xl1]/Im0-all.all_file[1].Image[k][yf:yl1,xf:xl1]/Im1)*(Im0+Im1)/2.)
                    flux[irun,iccd,ipair,k]=(Im0+Im1)/2.
                    # we are in a ultra saturation : stop 
                    if ipair>0 and flux[irun,iccd,ipair,k]<flux[irun,iccd,ipair-1,k] : 
                        flux[irun,iccd,ipair,k]=flux[irun,iccd,ipair-1,k] 
                        continue
                    nflux[irun,iccd,k]+=1
                    av=np.mean((ImDiff[:,:]))
                    std=np.std((ImDiff[:,:]))
                    gooda= (np.abs(ImDiff-av)<4.*std)
                    nb_gooda=np.sum(gooda)
                    diff_tmean=np.sum(ImDiff*gooda)/nb_gooda
                    # serial
                    good=gooda[:,:-1]*gooda[:,1:]
                # 
                    count=np.sum(np.sum(np.split(good,serial_split,axis=1),axis=2),axis=1)
                    countt=np.sum(good)
                #
                #   as it's computed from an image diff , we have to  /2 , 
                #   and we will plot in .1% unit , so we *1000
                    covars[irun,iccd,ipair,k,:]=np.sum(np.sum(np.split((ImDiff[:,:-1]-diff_tmean)
                                                *(ImDiff[:,1:]-diff_tmean)*good,serial_split,axis=1)
                                                ,axis=2),axis=1)/count/flux[irun,iccd,ipair,k]/2.*1000.
                    covarst[irun,iccd,ipair,k]=np.sum((ImDiff[:,:-1]-diff_tmean)
                                                *(ImDiff[:,1:]-diff_tmean)*good
                                                )/countt/flux[irun,iccd,ipair,k]/2.*1000.
                    
                #   //
                    good=gooda[:-1,:]*gooda[1:,:]
                # 
                    count=np.sum(np.sum(np.split(good,line_split,axis=0),axis=2),axis=1)
                #
                #   as it's computed from an image diff , we have to  /2 , 
                #   and we will plot in .1% unit , so we *1000
                    covarp[irun,iccd,ipair,k,:]=np.sum(np.sum(np.split((ImDiff[:-1,:]-diff_tmean)
                                                *(ImDiff[1:,:]-diff_tmean)*good,line_split,axis=0)
                                                ,axis=2),axis=1)/count/flux[irun,iccd,ipair,k]/2.*1000.
                    covarpt[irun,iccd,ipair,k]=np.sum((ImDiff[:-1,:]-diff_tmean)
                                                *(ImDiff[1:,:]-diff_tmean)*good
                                                )/countt/flux[irun,iccd,ipair,k]/2.*1000.
                # ratio of covar 
                    covarps[irun,iccd,ipair,k]= covarpt[irun,iccd,ipair,k]/covarst[irun,iccd,ipair,k]  
            for k in range(number_of_amp_sensor) :
                # per amplifier  plot , serial , // & ratio 
                # flux max min per sensor 
                flux_max=np.max(flux[irun,iccd,:nflux[irun,iccd,k],:number_of_amp_sensor])*1.1 
                #flux_min=min(1100,np.min(flux[irun,iccd,:nflux[irun,iccd,k],:number_of_amp_sensor]))*.9
                #
                fig=plt.figure(figsize=[15,15])
                title='covariance  vs flux for : raft=%s  ccd=%s Seg=%s (Hdu=%d) for run %s' %(raft_cur,all_sensors[raft_cur][iccd],segments[k],k+1,run[irun])
                fig.suptitle(title)
                title='cov(1,0) '
                fig.add_subplot(2,3,1,title=title)
                for isplit in range(serial_split) :
                    label='%d/%d serial fraction' %(isplit+1,serial_split)
                    plt.plot(flux[irun,iccd,:nflux[irun,iccd,k],k],covars[irun,iccd,:nflux[irun,iccd,k],k,isplit],label=label)
                label='full image' 
                plt.plot(flux[irun,iccd,:nflux[irun,iccd,k],k],covarst[irun,iccd,:nflux[irun,iccd,k],k],label=label)
                plt.ylabel('cov(1,0)/flux * 1000. ')
                plt.xlabel('flux in ADU')
                plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
                plt.yscale('symlog')
                plt.xlim(0.,flux_max)
                plt.ylim(0.)
                plt.ticklabel_format(axis='x',style='sci', scilimits=(0,0))
                #plt.xscale('symlog')
                plt.grid(color='grey',linestyle='dotted',which='major',axis='both')
                title='cov(0,1)'
                fig.add_subplot(2,3,3,title=title)
                for isplit in range(serial_split) :
                    label='%d/%d // fraction' %(isplit+1,serial_split)
                    plt.plot(flux[irun,iccd,:nflux[irun,iccd,k],k],covarp[irun,iccd,:nflux[irun,iccd,k],k,isplit],label=label)
                label='full image' 
                plt.plot(flux[irun,iccd,:nflux[irun,iccd,k],k],covarpt[irun,iccd,:nflux[irun,iccd,k],k],label=label)
                plt.ylabel('cov(0,1)/flux * 1000. ')
                plt.xlabel('flux in ADU')
                plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
                plt.yscale('symlog')
                plt.xlim(0.,flux_max)
                plt.ylim(0.)
                plt.ticklabel_format(axis='x',style='sci', scilimits=(0,0))
                #plt.xscale('symlog')
                plt.grid(color='grey',linestyle='dotted',which='major',axis='both')
                title='cov(0,1)/cov(1,0)'
                fig.add_subplot(2,3,5,title=title)
                plt.plot(flux[irun,iccd,:nflux[irun,iccd,k],k],covarps[irun,iccd,:nflux[irun,iccd,k],k])
                plt.ylabel('cov(0,1)/cov(1,0) ')
                plt.xlabel('flux in ADU')
                plt.xlim(0.,flux_max)
                plt.ylim(0.)
                plt.ticklabel_format(axis='x',style='sci', scilimits=(0,0))
                #plt.xscale('symlog')
                plt.grid(color='grey',linestyle='dotted',which='major',axis='both')
                #plt.show()
                if show : plt.show()
                hdu='%d' % (k+1)
                root_plt=os.path.join(output_data,run[irun],raft_cur,all_sensors[raft_cur][iccd],hdu)
                os.makedirs(root_plt,exist_ok=True)
                plotfile=os.path.join(root_plt,'Cov_vs_Flux.png')
                print (plotfile)
                fig.savefig(plotfile,bbox_inches='tight')
                plt.close(fig) 
        # flux max min per sensor 
            flux_max=np.max(flux[irun,iccd,:number_of_pair_per_run[irun],:number_of_amp_sensor])*1.1 
            #flux_min=min(1100,np.min(flux[irun,iccd,:number_of_pair_per_run[irun],:number_of_amp_sensor]))*.9
        # per CCD plot , serial 
            fig=plt.figure(figsize=[15,15])
            fig.suptitle('cov(1,0)/Flux vs flux per slice of serial register : '+raft_cur+' ccd '+all_sensors[raft_cur][iccd]+' for run '+run[irun])
            #
            for k in range(number_of_amp_sensor) :
                title=' Seg=%s (Hdu=%d)' %(segments[k],k+1)
                fig.add_subplot(4,4,k+1,title=title)
                for isplit in range(serial_split) :
                    label='%d/%d serial fraction' %(isplit+1,serial_split)
                    plt.plot(flux[irun,iccd,:nflux[irun,iccd,k],k],covars[irun,iccd,:nflux[irun,iccd,k],k,isplit],label=label)
                label='full image' 
                plt.plot(flux[irun,iccd,:nflux[irun,iccd,k],k],covarst[irun,iccd,:nflux[irun,iccd,k],k],label=label)
                if k%4==0 :
                    plt.ylabel('cov(1,0)/flux * 1000. ')
                if k-number_of_amp_sensor>-5 :
                    plt.xlabel('flux in ADU')
                if (k+1)%4 ==0 :
                    plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
                plt.yscale('symlog')
                plt.xlim(0,flux_max)
                plt.ylim(0.)
                plt.ticklabel_format(axis='x',style='sci', scilimits=(0,0))
                #plt.xscale('symlog')
                plt.grid(color='grey',linestyle='dotted',which='major',axis='both')
                #plt.xlim(-1.,1000.)
            #plt.show()
            if show : plt.show()
            root_plt=os.path.join(output_data,run[irun],raft_cur,all_sensors[raft_cur][iccd])
            os.makedirs(root_plt,exist_ok=True)
            plotfile=os.path.join(root_plt,'Cov10_vs_Flux.png')
            print (plotfile)
            fig.savefig(plotfile,bbox_inches='tight')
            plt.close(fig) 
        # per CCD plot , // 
            fig=plt.figure(figsize=[15,15])
            fig.suptitle('cov(0,1)/Flux vs flux per slice of column : '+raft_cur+' ccd '+all_sensors[raft_cur][iccd]+' for run '+run[irun])
            #
            for k in range(number_of_amp_sensor) :
                title=' Seg=%s (Hdu=%d)' %(segments[k],k+1)
                fig.add_subplot(4,4,k+1,title=title)
                for isplit in range(serial_split) :
                    label='%d/%d // fraction' %(isplit+1,serial_split)
                    plt.plot(flux[irun,iccd,:nflux[irun,iccd,k],k],covarp[irun,iccd,:nflux[irun,iccd,k],k,isplit],label=label)
                label='full image' 
                plt.plot(flux[irun,iccd,:nflux[irun,iccd,k],k],covarpt[irun,iccd,:nflux[irun,iccd,k],k],label=label)
                if k%4==0 :
                    plt.ylabel('cov(0,1)/flux * 1000. ')
                if k-number_of_amp_sensor>-5 :
                    plt.xlabel('flux in ADU')
                if (k+1)%4 ==0 :
                    plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
                plt.yscale('symlog')
                plt.xlim(0.,flux_max)
                plt.ylim(0.)
                plt.ticklabel_format(axis='x',style='sci', scilimits=(0,0))
                #plt.xscale('symlog')
                plt.grid(color='grey',linestyle='dotted',which='major',axis='both')
                 #plt.xlim(-1.,1000.)
            #plt.show()
            if show : plt.show()
            root_plt=os.path.join(output_data,run[irun],raft_cur,all_sensors[raft_cur][iccd])
            os.makedirs(root_plt,exist_ok=True)
            plotfile=os.path.join(root_plt,'Cov01_vs_Flux.png')
            print (plotfile)
            fig.savefig(plotfile,bbox_inches='tight')
            plt.close(fig) 
        # per CCD plot , // / serial 
            fig=plt.figure(figsize=[15,15])
            fig.suptitle('cov(0,1)/cov(1,0) vs flux  : '+raft_cur+' ccd '+all_sensors[raft_cur][iccd]+' for run '+run[irun])
            #
            for k in range(number_of_amp_sensor) :
                title=' Seg=%s (Hdu=%d)' %(segments[k],k+1)
                fig.add_subplot(4,4,k+1,title=title)
                plt.plot(flux[irun,iccd,:nflux[irun,iccd,k],k],covarps[irun,iccd,:nflux[irun,iccd,k],k])
                if k%4==0 :
                    plt.ylabel('cov(0,1)/cov(1,0) ')
                if k-number_of_amp_sensor>-5 :
                    plt.xlabel('flux in ADU')
                plt.xlim(0.,flux_max)
                ymax=min(20.,np.max(covarps[irun,iccd,:nflux[irun,iccd,k],k]))*1.1
                plt.ylim(0.,ymax)
                plt.ticklabel_format(axis='x',style='sci', scilimits=(0,0))
                #plt.xscale('symlog')
                plt.grid(color='grey',linestyle='dotted',which='major',axis='both')
                #plt.xlim(-1.,1000.)
            #plt.show()
            if show : plt.show()
            root_plt=os.path.join(output_data,run[irun],raft_cur,all_sensors[raft_cur][iccd])
            os.makedirs(root_plt,exist_ok=True)
            plotfile=os.path.join(root_plt,'Cov_01div10_vs_Flux.png')
            print (plotfile)
            fig.savefig(plotfile,bbox_inches='tight')
            plt.close(fig) 

        # per raft plot 
        if raft_cur in raft_corner :
            corner=True
        else :
            corner=False
        # per CCD plot , // / serial 
        fig=plt.figure(figsize=[15,15])
        fig.suptitle('cov(0,1)/cov(1,0) vs flux  : '+raft_cur+' for run '+run[irun])
        #
        ifig=0
        for iccd in range(number_of_raft_sensor) :
#
            ifig+=1
            if corner and iccd==2 :
                ifig+=1
            if not(all_good_sensors[raft_cur][all_sensors[raft_cur][iccd]]) :
                continue
            ax=fig.add_subplot(3,3,ifig,title=all_sensors[raft_cur][iccd])        
#        for iccd in [6] :
            if ( all_sensors[raft_cur][iccd] in sensors_8ch ):
                number_of_amp_sensor=8
            else : 
                    number_of_amp_sensor=16  
            for k in range(number_of_amp_sensor) :
                label='Seg. %s (Hdu=%d)' % (segments[k],k+1)    
                plt.plot(flux[irun,iccd,:nflux[irun,iccd,k],k],covarps[irun,iccd,:nflux[irun,iccd,k],k],color=colors[k],label=label)
        # flux max min per sensor 
            flux_max=np.max(flux[irun,iccd,:number_of_pair_per_run[irun],:number_of_amp_sensor])*1.1 
            #flux_min=min(1100,np.min(flux[irun,iccd,:number_of_pair_per_run[irun],:number_of_amp_sensor]))*.9
        #
            if iccd%3==0 or corner :
                plt.ylabel('cov(0,1)/cov(1,0) ')
            if iccd-number_of_raft_sensor>-4 or corner :
                plt.xlabel('flux in ADU')
            plt.xlim(0.,flux_max)
            ymax=min(20.,np.max(covarps[irun,iccd,:number_of_pair_per_run[irun],:]))*1.1
            plt.ylim(0.,ymax)
            plt.ticklabel_format(axis='x',style='sci', scilimits=(0,0))
            #plt.xscale('symlog')
            plt.grid(color='grey',linestyle='dotted',which='major',axis='both')
            if corner :
                if iccd%2 ==1 : plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
            elif iccd==2 :
                plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
        if show : plt.show()
        root_plt=os.path.join(output_data,run[irun],raft_cur)
        os.makedirs(root_plt,exist_ok=True)
        plotfile=os.path.join(root_plt,'Cov_01div10_vs_Flux.png')
        print (plotfile)
        fig.savefig(plotfile,bbox_inches='tight')
        plt.close(fig) 

    return 
#    return     flux,covars,covarp

#noise=np.zeros((number_of_raft,number_of_run,number_of_sensor*16)) 
#noise_std=np.zeros((number_of_raft,number_of_run,number_of_sensor*16)) 
#cor=np.zeros((number_of_raft,number_of_run,144,144))
#
for iraft in range(len(raft)) :
    raft_cur=raft[iraft]
    # find all the files to compute the correlation between channels
    all_file, first_ccd_run0=find_pair(raft_cur,all_good_sensors)
    # do all the cov plots for this raft and all run 
#    flux,covars,covarp=cov_in_raft(raft_cur,all_file,show=False)
    if first_ccd_run0>-1 : cov_in_raft(raft_cur,all_file, first_ccd_run0,show=False)
    #
