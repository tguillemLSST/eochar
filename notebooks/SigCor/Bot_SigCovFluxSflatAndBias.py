#!/usr/bin/env python
# coding: utf-8

# Study of signal cavariance between channels in raft 
# Goal : This Notebook will generate plots allowing to perform a diagnostic on the signal covariance between channels inside a raft in BOT data .
# Author : P.Antilogus
# Version : 21st January 2021 
# New : Trying to run at CC-IN2P3 (T. Guillemin)

# Method: 
# 
# Signal covariance,  is estimated using flat pair differences at 2 fluxes (SFLAT data) , and in BIAS 
# the correlation  is estimated for each pair at 2 # scale : 
# - at the pixel level 
# - at the line level ( covariance between mean signal per line) 
# 
# In practice the covariance in BIAS is an offset present in all covriances , at high flux some new covariance pop up that are in most of the case an indication of gain common change between amplifier. 
# 
# The code plots the correlation per sensor and per raft . 

# system imports
import os
import time
from sys import exit 
import glob
#import pdb 
from astropy.time import Time

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
#import lsst.afw.display as afwDisplay
#from lsst.ip.isr import IsrTask
#import lsst.afw.geom as afwGeom
#import lsst.afw.math as afwMath

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
    
# load the frame analysis code 
#get_ipython().run_line_magic('run', '-i  /home/antilog/repos/eochar/python/lsst/eochar/bot_frame_op.py')
# this seems correct
#os.system('python /sps/lsst/users/tguillem/Rubin/Focal_Plane/lsst_distrib/w_2022_01/eochar/python/lsst/eochar/bot_frame_op.py')
#print('bot_frame_op loaded')
#from ...python/lsst/eochar import bot_frame_op
from eochar.bot_frame_op import *

# activate the butler
#repo_path = '/sps/lsst/users/tguillem/Rubin/Focal_Plane/lsst_distrib/w_2022_01/data_sflat_test'
#butler = Butler(repo_path)

#check the butler content
#registry = butler.registry
#for col in registry.queryCollections():
#        print(col)
#print('butler configured')

# CONFIGURATION FOR THE CURRENT EXECUTION  ========================================
# ---- raft and associated run ============ To be updated if needed 
# 'raft' :  list of raft 
# 'run' :  list of run , ex : ['9876','9874']  (remark : run is a string )   for all run of this raft : '*' 
# 'sensor' :  list of sensor , ex ['S00','S01','S02'] 
# amplifier : list of amplifier , ex [1,2]  , for 1 to 16 you can use [-1] instead 

run_all=['13151']

raft_itl=['R01', 'R02', 'R03', 'R10', 'R20', 'R41', 'R42', 'R43']
raft_e2v=['R11', 'R12', 'R13', 'R14', 'R21', 'R22', 'R23', 'R24', 'R30', 'R31', 'R32', 'R33', 'R34']
sensors_raft=['S00','S01','S02','S10','S11','S12','S20','S21','S22']
raft_corner=['R00','R04','R40','R44']
sensors_corner=['SG0','SG1','SW0','SW1']
sensors_8ch=['SW0','SW1']
#
all_sensors={}
#
for raft in raft_itl+raft_e2v:
    all_sensors[raft]=sensors_raft
for raft in raft_corner:
    all_sensors[raft]=sensors_corner
#
amplifier=[-1]
#
#raft=raft_itl+raft_e2v+raft_corner
#raft=raft_corner
#raft=['R23']
#raft=['R10','R11','R12','R20','R21','R22','R30','R31','R32']
#raft=['R00']
##
raft=['R14']  

# output file directory
#output_data='/home/antilog/DATA/eochar6'
output_data='/sps/lsst/users/tguillem/web'

print('Configuration arguments: ', str(sys.argv))
str_run_all = str(sys.argv[1])
str_raft = str(sys.argv[2])
str_all_sensors = str(sys.argv[3])
repo_path=str(sys.argv[4])
output_data==str(sys.argv[5])

run_all=[str_run_all]
raft=[str_raft]
all_sensors[str_raft]=[str_all_sensors]
butler = Butler(repo_path)

number_of_pair_per_run_max=40
# number of # flux (=3 1 BIAS , 2 super flat  here , for high and low sflat exposure fluxes )
nb_flux=3
run_visits=np.zeros((len(run_all),number_of_pair_per_run_max,2))
time_visits=np.zeros((len(run_all),nb_flux))
number_of_pair_per_run_ct=np.zeros((len(run_all),nb_flux),dtype=np.int16)
run=[]
for run_cur in run_all :
    # select 10 bias inages 
    visits_bias=butler.queryMetadata('raw', ['visit'], dataId={'run':run_cur,'testType': 'BIAS','imageType': 'BIAS'})
    if len(visits_bias) < 4 : 
        print('No good-enough BIAS in run ',run_cur)
        continue
    visits_bias.sort()
    nb_bias=min(16,len(visits_bias[2:]))
    #
    visits_flat=butler.queryMetadata('raw', ['visit'], dataId={'run':run_cur,'testType': 'SFLAT','imageType': 'FLAT'})
    if len(visits_flat) < 1 : 
        print('No good flat in run ',run_cur)
        continue
    #
    print('--------visits_flat')
    print(visits_flat)
    print('--------')
    visits_flat.sort()
    visits=visits_bias[2:2+nb_bias]+visits_flat
    good_visits = []
    ptim=[]
    for visit1, visit2 in zip(visits[:-1:2], visits[1::2]):
        #dId = {'visit': visit1, 'detector': 2}
        #correct value of 'detector'?
        dId = {'visit': visit1, 'detector': 66}
        raw1 = butler.get('raw', **dId)
        time1 = raw1.getInfo().getVisitInfo().getExposureTime()

        # Get ISR data for second image
        dId = {'visit': visit2, 'detector': 66}
        raw2 = butler.get('raw', **dId)
        time2 = raw2.getInfo().getVisitInfo().getExposureTime()
        if abs(time1 - time2) > 0.01:
            print("Mismatched exptimes")
        else :
            good_visits.append((visit1,visit2))
            ptim.append(time1)
    print('--------End loop visits')        
    run.append(run_cur)
    irun=len(run)-1
    iflux=nb_flux
    print(ptim)
    previous_time=ptim[-1]+1
    for ipair in range(len(ptim)-1,-1,-1) :
        #if ptim[ipair] < previous_time :
        if ptim[ipair] != previous_time :   
            iflux-=1
            if iflux<0 :
                print('Error , the exposure sequence is not the expected one (increasing exposure time with ',nb_flux,' different exposures time). List of time :',ptim)
                exit('Execution stopped on Error')            
            number_of_pair_per_run_ct[irun,iflux]= ipair+1
            time_visits[irun,iflux]= ptim[ipair]
        run_visits[irun,ipair,0]=good_visits[ipair][0]
        run_visits[irun,ipair,1]=good_visits[ipair][1]
        previous_time=ptim[ipair]
    if iflux!=0 :
        print('Error , the exposure sequence is not the expected one (increasing exposure time with ',nb_flux,' different exposures time). List of time :',ptim)
        exit('Execution stopped on Error')            
#
#number_of_pair_per_run=2
print('Time distribution of the #pair ',ptim)
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

print('Config OK')

def find_pair(raft_cur):
    # compute the noise per image
    # Do all the fft plots if asked 
    # intialisation 
    number_of_raft_sensor=len(all_sensors[raft_cur])
    #
    all_file=np.zeros((number_of_run,number_of_raft_sensor,number_of_pair_per_run_max,2),dtype=np.object_)
    #         
    for irun in range(len(run)) : 
        run_cur=run[irun]
        #
        #
        t0=0
        #
        for iccd in range(number_of_raft_sensor) : 
            #
            for it in range(number_of_pair_per_run_ct[irun,-1]) : 
                for ifile in range(2) : 
                    dataId2 = {'visit': run_visits[irun,it,ifile], 'detectorName': all_sensors[raft_cur][iccd], 'raftName': raft_cur}
                    try : 
                        file=butler.get('raw_filename', dataId2)[0][:-3]
                    except : 
                        print('Error , data not found : ',dataId2)
                        break
                    all_file[irun,iccd,it,ifile]=file
    return  all_file          

def sig_cor_in_raft(raft_cur,all_file,show_raft=True,show_ccd=False,per_sensor=False):
    #
    if raft_cur in raft_corner :
        nb_amp_in_raft=48
    else :
        nb_amp_in_raft=144
    #
    number_of_raft_sensor=len(all_sensors[raft_cur])
    # labels for the plots
    label_txt=np.zeros((number_of_raft_sensor),dtype=np.object_)
    label_pos=np.zeros((number_of_raft_sensor))
    if number_of_raft_sensor==9 : 
        label_chan=np.zeros((18),dtype=np.object_)
        label_chan_pos=np.zeros((18))
    else : 
        label_chan=np.zeros((6),dtype=np.object_)
        label_chan_pos=np.zeros((6))    
    ch_cur=0
    ichan_lab=0
    for  iccd in range(number_of_raft_sensor) : 
            label_chan[ichan_lab]='%d' % (1)
            label_chan_pos[ichan_lab]=ch_cur
            ichan_lab+=1
            if ( all_sensors[raft_cur][iccd] in sensors_8ch ):
                number_of_amp_sensor=8
            else : 
                number_of_amp_sensor=16
                label_chan[ichan_lab]='%d' % (9)
                label_chan_pos[ichan_lab]=ch_cur+8
                ichan_lab+=1
            #
            label_txt[iccd]=all_sensors[raft_cur][iccd]
            label_pos[iccd]=ch_cur+number_of_amp_sensor/2-.5
            ch_cur+=number_of_amp_sensor
    # number of plots : number_of_pair_per_run * 2 flux per pair , so for each flux 1 line with 2 correlation plots , the first linee is reeserved for the flux per channel  
    plot_nb_line=1+2
    plot_nb_col=nb_flux
    #  
    # we will compute 2 covar between channel, 1 between pixel and one for a fraction of line 
    # we can split a line  in 1/1 (mean per line) and for e2v :  , 1/2 , 1/4 1/8 1/59 
    # line_split=1
    # Remark : this is a simplified version of the code, for eeochar use  
    #  we only compute covar between averaged signal per line (equive to line_split=1)
    split_text='per line'
    # get the amplifier shape for this raft
    fits=pyfits.open(all_file[0,0,0,0])
    first_line,first_p_over,first_col,first_s_over=image_area(fits)
    fits.close()
    #ncol=len(fits[1].data[100,:])
    #nline==len(fits[1].data[:,100])
    xf=first_col+20
    xl=first_s_over-20
    yf=first_line+50
    yl=first_p_over-20
    #
    ncol=xl-xf
    nline=yl-yf
    #
 
    #
    for irun in range(number_of_run) :
        # 
        #
        ImDiff=np.zeros((nb_amp_in_raft,nline,ncol))
        Flux=np.zeros((nb_flux,nb_amp_in_raft))
        good=np.zeros((nb_amp_in_raft,nline,ncol),dtype=bool)
        nb_good=np.zeros((nb_amp_in_raft))
        diff_tmean=np.zeros((nb_amp_in_raft))
        # 
        nb_b=(nline)   
        b_count=np.zeros((nb_amp_in_raft,nline))
        #
        #
        covarb=np.zeros((nb_amp_in_raft,nb_amp_in_raft))
        corb=np.zeros((nb_flux,nb_amp_in_raft,nb_amp_in_raft))
        covar=np.zeros((nb_amp_in_raft,nb_amp_in_raft))
        cor=np.zeros((nb_flux,nb_amp_in_raft,nb_amp_in_raft))
        cor_max=np.zeros((nb_flux,nb_amp_in_raft,nb_amp_in_raft))
        corb_max=np.zeros((nb_flux,nb_amp_in_raft,nb_amp_in_raft))
        for ipair in range(number_of_pair_per_run_ct[irun,-1]) :
            offset=0
            for jflux in range(nb_flux):
                if ipair < number_of_pair_per_run_ct[irun,jflux] :
                    iflux=jflux
                    break
                else :
                    offset=number_of_pair_per_run_ct[irun,jflux]
            nb_iflux=number_of_pair_per_run_ct[irun,iflux]-offset
            #
            ch_cur=0
            for iccd in range(number_of_raft_sensor) :
                if ( all_sensors[raft_cur][iccd] in sensors_8ch ):
                    number_of_amp_sensor=8
                else : 
                    number_of_amp_sensor=16
                #
                fname=[all_file[irun,iccd,ipair,0],all_file[irun,iccd,ipair,1]]
                if (ipair+1)%10==0 and (iccd+1)%5==0 :
                    print(Time.now(),' Run:',irun,' Raft:',raft_cur,'; pair of file number:',ipair+1,' among:',number_of_pair_per_run_ct[irun,-1],' current raft channels ',ch_cur+1,' processed among ' , nb_amp_in_raft)
                    all=InFile(dirall=fname,Slow=False,verbose=True)
                else :
                    all=InFile(dirall=fname,Slow=False,verbose=False)                    
                #
                for k in range(number_of_amp_sensor) :
                    # select-truncate  the image for the channel ch_cur
                    Im0=np.median(all.all_file[0].Image[k][yf:yl,xf:xl])
                    Im1=np.median(all.all_file[1].Image[k][yf:yl,xf:xl])
                    if iflux==0 : 
                        # don't re-scale in flux for bias 
                        ImDiff[ch_cur]=(all.all_file[0].Image[k][yf:yl,xf:xl]-all.all_file[1].Image[k][yf:yl,xf:xl])
                    else :     
                        ImDiff[ch_cur]=((all.all_file[0].Image[k][yf:yl,xf:xl]/Im0-all.all_file[1].Image[k][yf:yl,xf:xl]/Im1)*(Im0+Im1)/2.)
                    Flux[iflux,ch_cur]+=(Im0+Im1)/2./nb_iflux
                    av=np.mean((ImDiff[ch_cur][:,:]))
                    std=np.std((ImDiff[ch_cur][:,:]))
                    # All measurements are done with signal within 4 sigma to the <> . 
                    good[ch_cur,:,:] = (np.abs(ImDiff[ch_cur]-av)<4.*std) 
                    nb_good[ch_cur]=np.sum(good[ch_cur,:,:])
                    diff_tmean[ch_cur]=np.sum(ImDiff[ch_cur]*good[ch_cur,:,:])/nb_good[ch_cur]
                    # 
                    b_count[ch_cur,:]=np.sum(good[ch_cur,:,:],axis=1)
                    #
                    ch_cur1=0
                    for iccd1 in range(number_of_raft_sensor) :
                        if ch_cur1 > ch_cur : break
                        if ( all_sensors[raft_cur][iccd1] in sensors_8ch ):
                            number_of_amp_sensor_1=8
                        else : 
                            number_of_amp_sensor_1=16
                        #
                        for l in range(number_of_amp_sensor_1) :
                            if ch_cur1 > ch_cur : break
                            # compute covar betwen pixels of # channels (ch_cur and ch_cur1 ) 
                            count=np.sum(good[ch_cur,:,:]*good[ch_cur1,:,:])
                            covar[ch_cur,ch_cur1]=np.sum((ImDiff[ch_cur,:,:]-diff_tmean[ch_cur])*good[ch_cur,:,:]
                                             *(ImDiff[ch_cur1,:,:]-diff_tmean[ch_cur1])*good[ch_cur1,:,:]/count)
                            covar[ch_cur1,ch_cur]=covar[ch_cur,ch_cur1]
                            # compute covar between <pixel> per line between # channels (ch_cur and ch_cur1 ) 
                            covarb[ch_cur,ch_cur1]=np.sum(
        (np.sum(ImDiff[ch_cur,:,:]*good[ch_cur,:,:],axis=1)/b_count[ch_cur,:]-diff_tmean[ch_cur])
       *(np.sum(ImDiff[ch_cur1,:,:]*good[ch_cur1,:,:],axis=1)/b_count[ch_cur1,:]-diff_tmean[ch_cur1])
            )/nb_b
                            covarb[ch_cur1,ch_cur]=covarb[ch_cur,ch_cur1]
                            ch_cur1+=1
                    ch_cur+=1
            ch_cur=0
            for iccd in range(number_of_raft_sensor) :
                if ( all_sensors[raft_cur][iccd] in sensors_8ch ):
                    number_of_amp_sensor=8
                else : 
                    number_of_amp_sensor=16
                #
                for k in range(number_of_amp_sensor) :
                    ch_cur1=0
                    for iccd1 in range(number_of_raft_sensor) :
                        if ch_cur1 > ch_cur : break
                        if ( all_sensors[raft_cur][iccd1] in sensors_8ch ):
                            number_of_amp_sensor_1=8
                        else : 
                            number_of_amp_sensor_1=16
                        #
                        for l in range(number_of_amp_sensor_1) :
                            if ch_cur1 > ch_cur : break
                            # we /2 to take into account that we computed this with image dif , but are interested by single image stat 
                            cor[iflux,ch_cur1,ch_cur]+=covar[ch_cur1,ch_cur]/nb_iflux/2.
                            corb[iflux,ch_cur1,ch_cur]+=covarb[ch_cur1,ch_cur]/nb_iflux/2.
#                            cor[iflux,ch_cur1,ch_cur]+=covar[ch_cur1,ch_cur]/(Flux[iflux,ch_cur]+Flux[iflux,ch_cur1])*2./nb_iflux
#                            corb[iflux,ch_cur1,ch_cur]+=covarb[ch_cur1,ch_cur]/(Flux[iflux,ch_cur]+Flux[iflux,ch_cur1])*2./nb_iflux
                            ch_cur1+=1
                    ch_cur+=1
            #
        for iflux in range(nb_flux):
            ch_cur=0
            for iccd in range(number_of_raft_sensor) :
                if ( all_sensors[raft_cur][iccd] in sensors_8ch ):
                    number_of_amp_sensor=8
                else : 
                    number_of_amp_sensor=16
                #
                for k in range(number_of_amp_sensor) :
                    ch_cur1=0
                    for iccd1 in range(number_of_raft_sensor) :
                        if ch_cur1 > ch_cur : break
                        if ( all_sensors[raft_cur][iccd1] in sensors_8ch ):
                            number_of_amp_sensor_1=8
                        else : 
                            number_of_amp_sensor_1=16
                        #
                        for l in range(number_of_amp_sensor_1) :
                            if ch_cur1 > ch_cur : break
                            cor[iflux,ch_cur,ch_cur1]=cor[iflux,ch_cur1,ch_cur]
                            corb[iflux,ch_cur,ch_cur1]=corb[iflux,ch_cur1,ch_cur]
                            if ch_cur != ch_cur1 :
                                cor_max[iflux,ch_cur,ch_cur1]=cor[iflux,ch_cur1,ch_cur]
                                corb_max[iflux,ch_cur,ch_cur1]=corb[iflux,ch_cur1,ch_cur]
                                cor_max[iflux,ch_cur1,ch_cur]=cor[iflux,ch_cur1,ch_cur]
                                corb_max[iflux,ch_cur1,ch_cur]=corb[iflux,ch_cur1,ch_cur]
                            ch_cur1+=1
                    ch_cur+=1
        
        # per CCD plot 
        if per_sensor :
            chs=0
            for iccd in range(number_of_raft_sensor) :
                if ( all_sensors[raft_cur][iccd] in sensors_8ch ):
                        number_of_amp_sensor=8
                else : 
                        number_of_amp_sensor=16
                fig=plt.figure(figsize=[15,15])
                fig.suptitle('Signal channels cov in raft '+raft_cur+' ccd '+all_sensors[raft_cur][iccd]+' for run '+run[irun])
                #
                chl=chs+number_of_amp_sensor
                for iflux in range(nb_flux):
                    loc=plot_nb_col+1+iflux
                    if iflux==0 : 
                        title='BIAS channels cov(pixel)' 
                    else:
                        title='Channels cov(pixel), exp. %fs' %(time_visits[irun,iflux])
                    fig.add_subplot(plot_nb_line,plot_nb_col,loc,title=title)
                    cmin=np.min(cor[iflux,chs:chl,chs:chl])
                    #if cmin<-.1 : cmin=-1.
                    cmax=np.max(cor_max[iflux,chs:chl,chs:chl])
                    #cmin=min(-.1,np.min(cor[ipair,chs:chl,chs:chl]))
                    #if cmin<-.1 : cmin=-1.
                    #cmax=int(cor_max[iccd,ipair]*10)/10.+.1
                    plt.imshow(cor[iflux,chs:chl,chs:chl],cmap=matplotlib.cm.nipy_spectral,origin='lower',vmin=cmin,vmax=cmax)
                    #plt.imshow(cor[irun,:,:],cmap=matplotlib.cm.nipy_spectral,shape=(144,144),origin='lower',vmin=-.15,vmax=1.)
                    #plt.xticks(label_pos,label_txt)
                    #plt.yticks(label_pos,label_txt)
                    plt.colorbar()
                    #
                    loc=(plot_nb_col)*2+1+iflux
                    if iflux==0 : 
                        title='BIAS channels cov(line)' 
                    else:
                        title='Channels cov(line), exp. %fs' %(time_visits[irun,iflux])
                    fig.add_subplot(plot_nb_line,plot_nb_col,loc,title=title)
                    #cmin=min(-.1,np.min(corb[ipair,chs:chl,chs:chl]))
                    #if cmin<-.1 : cmin=-1.
                    #cmax=int(corb_max[iccd,ipair]*10)/10.+.1
                    cmin=np.min(corb[iflux,chs:chl,chs:chl])
                    #if cmin<-.1 : cmin=-1.
                    cmax=np.max(corb_max[iflux,chs:chl,chs:chl])
                    plt.imshow(corb[iflux,chs:chl,chs:chl],cmap=matplotlib.cm.nipy_spectral,origin='lower',vmin=cmin,vmax=cmax)
                    #plt.imshow(cor[irun,:,:],cmap=matplotlib.cm.nipy_spectral,shape=(144,144),origin='lower',vmin=-.15,vmax=1.)
                    #plt.xticks(label_pos,label_txt)
                    #plt.yticks(label_pos,label_txt)
                    plt.colorbar()
                #
                ax1=fig.add_subplot(plot_nb_line,1,1)
                #ax2=ax1.twiny()
                #plt.tick_params(axis='x',top='on',labeltop='on')
                #,bottom='off',labelbottom='off')
                #
                for iflux in range(nb_flux) :
                    if iflux==0 :
                        label='<BIAS>' 
                        ax1.plot(range(number_of_amp_sensor),Flux[iflux,chs:chl],label=label)
                        label='Var(Signal_pixel) BIAS' 
                        ax1.plot(range(number_of_amp_sensor),np.diag(cor[iflux,chs:chl,chs:chl]),label=label)
                        label='Var(Signal_line)  BIAS' 
                        ax1.plot(range(number_of_amp_sensor),np.diag(corb[iflux,chs:chl,chs:chl]),label=label)
                    else :    
                        label='<flux> exp. %fs' %(time_visits[irun,iflux])
                        ax1.plot(range(number_of_amp_sensor),Flux[iflux,chs:chl],label=label)
                        label='Var(flux_pixel) exp. %fs' %(time_visits[irun,iflux])
                        ax1.plot(range(number_of_amp_sensor),np.diag(cor[iflux,chs:chl,chs:chl]),label=label)
                        label='Var(flux_line) exp. %fs' %(time_visits[irun,iflux])
                        ax1.plot(range(number_of_amp_sensor),np.diag(corb[iflux,chs:chl,chs:chl]),label=label)
                ax1.set_yscale('log')
                ax1.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
                ax1.set_ylabel('Signal mean (var) in image in ADU (ADU**2)')
                ax1.set_xlabel('Channel (hdu-1)')
                #plt.show()
                if show_ccd : plt.show()
                root_plt=os.path.join(output_data,run[irun],raft_cur,all_sensors[raft_cur][iccd])
                os.makedirs(root_plt,exist_ok=True)
                plotfile=os.path.join(root_plt,'cor_in_signal.png')
                print (plotfile)
                fig.savefig(plotfile,bbox_inches='tight')
                plt.close(fig) 
                #
                chs=chl
        # per raft plot 
        fig=plt.figure(figsize=[15,15])
        fig.suptitle('Signal channels cov  in raft :'+raft_cur+' for run '+run[irun])
        #
        for iflux in range(nb_flux):
            loc=plot_nb_col+1+iflux
            if iflux==0 : 
                title='BIAS channels cov(pixel)' 
            else:
                title='Channels cov(pixel), exp. %fs' %(time_visits[irun,iflux])
            fig.add_subplot(plot_nb_line,plot_nb_col,loc,title=title)
            cmin=np.min(cor[iflux,:,:])
            #if cmin<-.1 : cmin=-1.
            cmax=np.max(cor_max[iflux,:,:])
            #cmin=min(-.1,np.min(cor[ipair,:,:]))
            #if cmin<-.1 : cmin=-1.
            #cmax=int(np.max(cor_max[:,ipair])*10)/10.+.1
            plt.imshow(cor[iflux,:,:],cmap=matplotlib.cm.nipy_spectral,origin='lower',vmin=cmin,vmax=cmax)
            #plt.imshow(cor[irun,:,:],cmap=matplotlib.cm.nipy_spectral,shape=(144,144),origin='lower',vmin=-.15,vmax=1.)
            plt.xticks(label_pos,label_txt)
            plt.yticks(label_pos,label_txt)
            plt.colorbar()
            #
            loc=(plot_nb_col)*2+1+iflux
            if iflux==0 : 
                title='BIAS channels cov(line)' 
            else:
                title='Channels cov(line), exp. %fs' %(time_visits[irun,iflux])
            fig.add_subplot(plot_nb_line,plot_nb_col,loc,title=title)
            #cmin=min(-.1,np.min(corb[ipair,:,:]))
            #if cmin<-.1 : cmin=-1.
            #cmax=int(np.max(corb_max[:,ipair])*10)/10.+.1
            cmin=np.min(corb[iflux,:,:])
            #if cmin<-.1 : cmin=-1.
            cmax=np.max(corb_max[iflux,:,:])
            #
            plt.imshow(corb[iflux,:,:],cmap=matplotlib.cm.nipy_spectral,origin='lower',vmin=cmin,vmax=cmax)
            #plt.imshow(cor[irun,:,:],cmap=matplotlib.cm.nipy_spectral,shape=(144,144),origin='lower',vmin=-.15,vmax=1.)
            plt.xticks(label_pos,label_txt)
            plt.yticks(label_pos,label_txt)
            plt.colorbar()
        #
        ax1=fig.add_subplot(plot_nb_line,1,1)
        ax2=ax1.twiny()
        #plt.tick_params(axis='x',top='on',labeltop='on')
        #,bottom='off',labelbottom='off')
        #
        for iflux in range(nb_flux) :
            if iflux==0 :
                label='<BIAS>' 
                ax1.plot(range(nb_amp_in_raft),Flux[iflux,:],label=label)
                label='Var(Signal_pixel) BIAS' 
                ax1.plot(range(nb_amp_in_raft),np.diag(cor[iflux,:,:]),label=label)
                label='Var(Signal_line)  BIAS' 
                ax1.plot(range(nb_amp_in_raft),np.diag(corb[iflux,:,:]),label=label)
            else :    
                label='<flux> exp. %fs' %(time_visits[irun,iflux])
                ax1.plot(range(nb_amp_in_raft),Flux[iflux,:],label=label)
                label='Var(flux_pixel) exp. %fs' %(time_visits[irun,iflux])
                ax1.plot(range(nb_amp_in_raft),np.diag(cor[iflux,:,:]),label=label)
                label='Var(flux_line) exp. %fs' %(time_visits[irun,iflux])
                ax1.plot(range(nb_amp_in_raft),np.diag(corb[iflux,:,:]),label=label)
        ax1.set_yscale('log')
        ax1.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
        ax1.set_ylabel('Signal mean (var) in image in ADU (ADU**2)')
        ax1.set_xlabel('Sensor')
        ax1.set_xticks(label_pos)
        ax1.set_xticklabels(label_txt)
        ax2.set_xlim(ax1.get_xlim())
        ax2.tick_params(labeltop='on')
        ax2.set_xticks(label_chan_pos)
        ax2.set_xticklabels(label_chan)
        plt.xlabel('Channel')
        title=''
        #plt.show()
        if show_raft : plt.show()
        root_plt=os.path.join(output_data,run[irun],raft_cur)
        os.makedirs(root_plt,exist_ok=True)
        plotfile=os.path.join(root_plt,'cor_in_signal.png')
        print (plotfile)
        fig.savefig(plotfile,bbox_inches='tight')
        plt.close(fig) 
#    return noise,noise_std,cor
    return 

#noise=np.zeros((number_of_raft,number_of_run,number_of_sensor*16)) 
#noise_std=np.zeros((number_of_raft,number_of_run,number_of_sensor*16)) 
#cor=np.zeros((number_of_raft,number_of_run,144,144))
#
for iraft in range(len(raft)) :
    raft_cur=raft[iraft]
    # find all the files to compute the correlation between channels
    all_file=find_pair(raft_cur)
    # do all the correlation plots for this raft and all run 
    sig_cor_in_raft(raft_cur,all_file,show_raft=True,show_ccd=False,per_sensor=True)
    #
