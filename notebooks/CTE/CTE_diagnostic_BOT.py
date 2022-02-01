#!/usr/bin/env python
# coding: utf-8

# # Notebook : CTE_diagnostic for full BOT 
# 
# Goal : This Notebook will generate plots allowing to perform  a diagnostic for channels with a bad CTE .  
# 
# Author : P.Antilogus
# 
# Version : 23rd Nov 2020 
# 
#     The data used are the PTC flats ( large data set at many fluxes), the study of the evolution of the first overscan pixels with the flux give an unique indication on : 
#     - is it a trap ? :  number of defered charges ~ a few e- ( 1 to 12) ~ constant for flux between 100 e- 10^4 
#     - is there issue in the overscan shape ( various plot related to the overscan shape including : 
#         - evolution of the overscan in function of the flux
#         - comparison of the noise estimated with the overscan associated to the image area or to the non-image corner
#         - direct plot of the overscan for all the fluxes 
# 
#     The data for different runs can be overploted to see the evolution / differences 
#     The raft , dataset ( BNL , SLAC ..) has to be provided , and a subselection can be done 
#     on the run , devices and amplifiers . 
#     
#     You can descide to only do the summary plots (for all amplifiers ) , or display the diagnostic plots
#     for all or a fraction of the amplifiers . See box 3 for configuration parameters ( serial_polt ...) 
#     
#     
# 
# Remark  :
# 
#      To provide a flux scale in e- a gain is needed , the one computed from the Fe55 data from the same run is used
#      by doing a query of the eotest DB . If you haven't installed the needed software the gain will be set to 1 
#      for all channels ... but you can change this default value if you want  (0.7 is better for e2v and 0.9 for itl )  
# 


from eochar.bot_frame_op import * 
from eochar.frame_cte_noise import *

# system imports
import os
import time
import sys
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
print(os.environ['PYTHONPATH'])
from eTraveler.clientAPI.connection import Connection
print('OK1')
#from importlib import import_module
#datacat_utilities = import_module('datacat-utilities')
#from datacat_utilities import datacat_utilities
#from datacat-utilites import datacat-utilites
#.get_EO_analysis_results import get_EO_analysis_results
#datacat_utilities = __import__('datacat-utilities')
sys.path.append("/sps/lsst/users/tguillem/Rubin/Focal_Plane/lsst_distrib/w_2022_01/database_eotest/datacat-utilities")
#datacat_utilities = __import__('datacat-utilities')

print('OK2')
if display_in_electron :
    try : 
        from get_EO_analysis_results import get_EO_analysis_results
        from get_EO_analysis_files import get_EO_analysis_files
        from exploreFocalPlane import exploreFocalPlane
        from exploreRaft import exploreRaft
        from eTraveler.clientAPI.connection import Connection
        g = get_EO_analysis_results(db=DB_key[data_to_select['Data_location']])
        eotest_db=True
        unit='e-'
    except:
        print('No access to eotest DB , so the gain of all devices will be set to 1 ( 1 ADU counted as 1 e- )')
        eotest_db=False
        unit='ADU'
else :
    eotest_db=False
    unit='ADU'
    
# activate the butler
repo_path = '/sps/lsst/users/tguillem/Rubin/Focal_Plane/lsst_distrib/w_2022_01/data_PTC'
butler = Butler(repo_path)

# CONFIGURATION FOR THE CURRENT EXECUTION  ========================================
# ---- raft and associated run ============ To be updated if needed 
# 'raft' :  list of raft 
# 'run' :  list of run , ex : ['9876','9874']  (remark : run is a string )   for all run of this raft : '*' 
# 'sensor' :  list of sensor , ex ['S00','S01','S02'] 
# amplifier : list of amplifier , ex [1,2]  , for 1 to 16 you can use [-1] instead 
run_all=['*']
#run_all=['12680']
#
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
raft=['R14']
#raft=['*']
all_sensors['R14']=['S22']
#all_sensors['R33']=['S00']
#
#directory to output data
#root_dir='/home/antilog/DATA/6'
root_dir='/sps/lsst/users/tguillem/web/PTC/'

visits_all=[]
run=[]
for run_cur in run_all : 
    # select only bias inages 
    #visits=butler.queryMetadata('raw', ['visit'],dataId={'run':run_cur,'imageType': 'FLAT','testtype' : 'FLAT'})
    #visits+=butler.queryMetadata('raw', ['visit'],dataId={'run':run_cur,'imageType': 'FLAT','testtype' : 'SFLAT'})
    visits=butler.queryMetadata('raw', ['visit'],dataId={'imageType': 'FLAT','testtype' : 'FLAT'})
    visits+=butler.queryMetadata('raw', ['visit'],dataId={'imageType': 'FLAT','testtype' : 'SFLAT'})
    if len(visits) < 1 : 
        print('No image in run ',run_cur)
    else :
        visits.sort()
        visits_all.append(visits)
        run.append(run_cur)
        print('For run ',run_cur,' we identified ',len(visits_all[-1]),' FLAT images from FLAT and SFLAT testtype')
#
nb_visits_all=np.zeros((len(run_all)),dtype=np.int16)
for irun in range(len(run_all)) :
    nb_visits_all[irun]=len(visits_all[irun])
#
nb_visits_max=np.int(np.max(nb_visits_all))
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
    print('No data found for ',data_to_select)
    raise 


# CONFIGURATION FOR THE CURRENT EXECUTION  ========================================
# Do we look for serial CTE ?
serial=True
# do we plot per amplifiers serial CTE diagnostic plots  
serial_plot=True
# Do we look for // CTE ?
parallel=True
# do we plot per amplifiers // CTE diagnostic plots  
parallel_plot=True
# Do we put the plots on screen ? ( or just in file ) 
plot_on_screen=False
# default gain if eotest data not available ( 0.7 is better of e2v , 0.9 for itl ) 
default_gain=1.

# create the data table with the right dimension 
if parallel :
    cte_p_data=np.zeros((number_of_run,number_of_raft,number_of_sensor),dtype=np.object_)
if serial : 
    cte_s_data=np.zeros((number_of_run,number_of_raft,number_of_sensor),dtype=np.object_)
file_name=np.zeros((number_of_run,number_of_raft,number_of_sensor),dtype=np.object_) 
#
for irun in range(len(run)) : 
    run_cur=run[irun]
    print('CTE Analysis of run ',run_cur)
    #
    output_dir=os.path.join(root_dir,run_cur)
    #
    for iraft in range(len(raft)) :
        raft_cur=raft[iraft]
        print('- Analysis for raft ',raft_cur)  
        #
        if eotest_db : 
            raft_list, data = g.get_tests(site_type="I&T-Raft", test_type="gain", run =run_cur)  # get the data for I&T-Raft
            res = g.get_results(test_type='gain', data=data, device=raft_cur)  # get the data for a raft
        #
        sensors=all_sensors[raft_cur]
        for iccd in range(len(sensors)) : 
            # 
            gain=True
            # get the gain for this sensor :
            if eotest_db :
                sensor_id=[]
                for d in data['steps']['fe55_raft_analysis']['fe55_raft_analysis']:
                    if 'slot' in d.keys():
                        if d['slot']==sensors[iccd]:
                            sensor_id=d['sensor_id']
                if len(sensor_id)==0 :
                            all_gain=np.ones((16))*default_gain 
                            gain=False
                else :
                    all_gain=np.array(res[sensor_id])
                    # protection against badly measured gain : we accept only a gain within 15% of the median
                    medgain=np.median(all_gain)
                    if medgain < 0.5 : 
                        # sensor not operational or what ? 
                        print('Gain found in eotest for this sensor looks bad (median gain=',medgain,'), sensor dead  ? We set all gain to ',default_gain)
                        all_gain=np.ones((16))*default_gain 
                    else : 
                        for igain in range(len(all_gain)):
                            if np.abs(all_gain[igain]-medgain)>.15 :
                                print('Gain channel ',igain+1,' found out of range ',all_gain[igain],' Channel could be dead or its a bad measurement . We set it to the median gain =',medgain)
                                all_gain[igain]=medgain    
                    print('sensor =',sensor_id) 
            else : 
                # if no DB info set the gain to 1 
                all_gain=np.ones((16))*default_gain
                gain=False
            # 
            new=True
            for ifile in range(nb_visits_all[irun]) :
                dataId2 = {'visit': visits_all[irun][ifile], 'detectorName': sensors[iccd], 'raftName': raft_cur}
                try : 
                    fileN=butler.get('raw_filename', dataId2)[0][:-3]
                except : 
                    print('Error , data not found : ',dataId2)
                    continue
                if new : 
                    file_name[irun,iraft,iccd]=[]
                    new=False 
                file_name[irun,iraft,iccd].append(fileN)
            print(file_name[irun,iraft,iccd])    
            print('- Analysis for CCD ',sensors[iccd],' using ',len(file_name[irun,iraft,iccd]),' files')  
            if len(file_name[irun,iraft,iccd])<1 : 
                break 
            file=Ifile(dirall=file_name[irun,iraft,iccd],fkey={},Slow=False)
            if serial :
                cte_s_data[irun,iraft,iccd]=cte(all_file=file.all_file,gain=all_gain,serie=True)
            if parallel: 
                cte_p_data[irun,iraft,iccd]=cte(all_file=file.all_file,gain=all_gain,serie=False)
            nf=0
            # fill a text with the info to print on each plot TBD
            if gain :
                ccd='gain from Fe55; run='+run_cur
            else : 
                ccd='1e- set to 1 ADU; run='+run_cur
            for ichan in amplifier_list  :
                # when we give the list to plot it's in general for raft , not corner raft , so we limit to corner raft channels if needed 
                if ichan>file.all_file[-1].HduMax : break
                ch=ichan-1
                if serial and  serial_plot :
                    cte_s_data[irun,iraft,iccd].plot_cte(ch=ch,ccd_name=ccd,nf=nf,on_screen=plot_on_screen,root_dir=output_dir,unit=unit)
                if parallel and parallel_plot :
                    cte_p_data[irun,iraft,iccd].plot_cte(ch=ch,ccd_name=ccd,nf=nf,on_screen=plot_on_screen,root_dir=output_dir,unit=unit)
    


# In[7]:


# Summary plots (cte vs flux ) for the # selected sensors , overlaping the # selected run 
# 
sel_color=['r','g','b','c','m','y','k']
nb_color=len(sel_color)
sel_mark=['<','>','*','^','v','+']
nb_mark=len(sel_mark)
#
yv=5.0e-6
yy=[yv,yv]
xx=[10.,300000.]
#
for iraft in range(len(raft)) :
    raft_cur=raft[iraft]
    sensors=all_sensors[raft_cur] 
    # we set the output to the first run ...but at the end we will the link the produced file in all directories 
    if serial :
        for iccd in range(len(sensors)) :
            fig=plt.figure(figsize=[20,20])
            fig.suptitle('Serial CTE :'+raft_cur+' '+sensors[iccd],fontsize=14)
            for ichan in amplifier_list  :
                # when we give the list to plot it's in general for raft , not corner raft , so we limit to corner raft channels if needed 
                if ichan>cte_s_data[0,iraft,iccd].first_file.HduMax : break
                ch=ichan-1
                #
                ax=fig.add_subplot(4,4,ch+1)
                for irun in range(number_of_run) :
                      plt.errorbar(cte_s_data[irun,iraft,iccd].cte_flux_s[ch,0:cte_s_data[irun,iraft,iccd].lmax[ch]],cte_s_data[irun,iraft,iccd].ylev[ch,0:cte_s_data[irun,iraft,iccd].lmax[ch]], yerr=cte_s_data[irun,iraft,iccd].ylev_std[ch,0:cte_s_data[irun,iraft,iccd].lmax[ch]],fmt=sel_mark[irun%nb_mark],color=sel_color[irun%nb_color],label=run[irun])
                # 
                plt.plot(xx,yy,'g')
                label='channel Hdu %d  : <flux> in %s' % (ichan,unit)
                plt.xlabel(label)
                if ch%4 == 0 :
                    plt.ylabel('1-CTE (serial CTE)')
                plt.xscale('log')
                plt.yscale('log')
                if (ch+1)%4==0 :
                    plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
                plt.xlim(xx[0],xx[1])
                plt.ylim(1e-8,5e-5)
            if plot_on_screen : plt.show()
            root_plt=os.path.join(root_dir,run[0],raft_cur,sensors[iccd])
            # create the directorty
            os.makedirs(root_plt,exist_ok=True)
            #
            plotfile=os.path.join(root_plt,'Summary_serialCTE.png')            
            fig.savefig(plotfile,bbox_inches='tight')
            if not(plot_on_screen) : plt.close(fig)
            # link in all other run directory 
            for irun in range(1,number_of_run) :
                root_plt=os.path.join(root_dir,run[irun],raft_cur,sensors[iccd])
                # create the directorty
                os.makedirs(root_plt,exist_ok=True)
                # make the link
                os.symlink(plotfile,root_plt+'/Summary_serialCTE.png')
    if parallel :
        for iccd in range(len(sensors)) :
            fig=plt.figure(figsize=[20,20])
            fig.suptitle('// CTE :'+raft_cur+' '+sensors[iccd],fontsize=14)
            for ichan in amplifier_list  :
                # when we give the list to plot it's in general for raft , not corner raft , so we limit to corner raft channels if needed 
                if ichan>cte_s_data[0,iraft,iccd].first_file.HduMax : break
                ch=ichan-1
                #
                ax=fig.add_subplot(4,4,ch+1)
                for irun in range(number_of_run) :
                    plt.errorbar(cte_p_data[irun,iraft,iccd].cte_flux_s[ch,0:cte_p_data[irun,iraft,iccd].lmax[ch]],cte_p_data[irun,iraft,iccd].ylev[ch,0:cte_p_data[irun,iraft,iccd].lmax[ch]], yerr=cte_p_data[irun,iraft,iccd].ylev_std[ch,0:cte_p_data[irun,iraft,iccd].lmax[ch]],fmt=sel_mark[irun%nb_mark],color=sel_color[irun%nb_color],label=run[irun])
                # 
                plt.plot(xx,yy,'g')
                label='channel Hdu %d  : <flux> in %s' % (ichan,unit)
                plt.xlabel(label)
                if ch%4 == 0 :
                    plt.ylabel('1-CTE (serial CTE)')
                plt.xscale('log')
                plt.yscale('log')
                if (ch+1)%4==0 :
                    plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
                plt.xlim(xx[0],xx[1])
                plt.ylim(1e-8,5e-5)
            if plot_on_screen : plt.show()
            root_plt=os.path.join(root_dir,run[0],raft_cur,sensors[iccd])
            # create the directorty
            os.makedirs(root_plt,exist_ok=True)
            #
            plotfile=os.path.join(root_plt,'Summary_parallelCTE.png')            
            fig.savefig(plotfile,bbox_inches='tight')
            if not(plot_on_screen) : plt.close(fig)
            # link in all other run directory 
            for irun in range(1,number_of_run) :
                root_plt=os.path.join(root_dir,run[0],raft_cur,sensors[iccd])
                # create the directorty
                os.makedirs(root_plt,exist_ok=True)
                # make the link
                os.symlink(plotfile,root_plt+'/Summary_parallelCTE.png')


# In[8]:


# Summary plots (cte vs flux ) for the # selected sensors , overlaping the # selected run 
# colors for each channel setup
cmap=plt.get_cmap('gist_ncar')
colors=[cmap(j)[:3] for j in np.linspace(0,1,17)]
# 
sel_color=['r','g','b','c','m','y','k']
nb_color=len(sel_color)
sel_mark=['<','>','*','^','v','+']
nb_mark=len(sel_mark)
#
yv=5.0e-6
yy=[yv,yv]
xx=[10.,300000.]
#
for iraft in range(len(raft)) :
    raft_cur=raft[iraft]
    if raft_cur in raft_corner :
        corner=True
    else :
        corner=False
    sensors=all_sensors[raft_cur] 
    # we set the output to the first run ...but at the end we will the link the produced file in all directories 
    output_dir=os.path.join(root_dir,run[0],raft_cur)
    if serial :
        fig=plt.figure(figsize=[20,20])
        fig.suptitle('Serial CTE :'+raft_cur,fontsize=14)
        ifig=0
        for iccd in range(len(sensors)) :  
            ifig+=1
            if corner and iccd==2 :
                ifig+=1
            ax=fig.add_subplot(3,3,ifig)
            for ichan in amplifier_list  :
                # when we give the list to plot it's in general for raft , not corner raft , so we limit to corner raft channels if needed 
                if ichan>cte_s_data[0,iraft,iccd].first_file.HduMax : break
                ch=ichan-1
                #
                for irun in range(number_of_run) :
                    if ch==0 or irun==0 :
                        labelx='%d (%s)' % (ch+1,run[irun])
                        plt.errorbar(cte_s_data[irun,iraft,iccd].cte_flux_s[ch,0:cte_s_data[irun,iraft,iccd].lmax[ch]],cte_s_data[irun,iraft,iccd].ylev[ch,0:cte_s_data[irun,iraft,iccd].lmax[ch]], yerr=cte_s_data[irun,iraft,iccd].ylev_std[ch,0:cte_s_data[irun,iraft,iccd].lmax[ch]],fmt=sel_mark[irun%nb_mark],color=colors[ch],label=labelx)
                    else :    
                        plt.errorbar(cte_s_data[irun,iraft,iccd].cte_flux_s[ch,0:cte_s_data[irun,iraft,iccd].lmax[ch]],cte_s_data[irun,iraft,iccd].ylev[ch,0:cte_s_data[irun,iraft,iccd].lmax[ch]], yerr=cte_s_data[irun,iraft,iccd].ylev_std[ch,0:cte_s_data[irun,iraft,iccd].lmax[ch]],fmt=sel_mark[irun%nb_mark],color=colors[ch])
                # 
            plt.plot(xx,yy,'g')
            label='CCD %s' % sensors[iccd]
            plt.xlabel(label)
            if ifig%3 == 1 :
                plt.ylabel('1-CTE (serial CTE)')
            plt.xscale('log')
            plt.yscale('log')
            if corner :
                if iccd%2 ==1 : plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
            elif iccd==2 :
                plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
            plt.xlim(xx[0],xx[1])
            plt.ylim(1e-8,5e-5)
        if plot_on_screen : plt.show()
        root_plt=os.path.join(root_dir,run[0],raft_cur)
        # create the directorty
        os.makedirs(root_plt,exist_ok=True)
        #
        plotfile=os.path.join(root_plt,'SerialCTE.png')            
        fig.savefig(plotfile,bbox_inches='tight')
        if not(plot_on_screen) : plt.close(fig)
        # link in all other run directory 
        for irun in range(1,number_of_run) :
            root_plt=os.path.join(root_dir,run[irun],raft_cur)
            # create the directorty
            os.makedirs(root_plt,exist_ok=True)
            # make the link
            os.symlink(plotfile,root_plt+'/SerialCTE.png')
    if parallel :
        fig=plt.figure(figsize=[20,20])
        fig.suptitle('// CTE :'+raft_cur,fontsize=14)
        ifig=0
        for iccd in range(len(sensors)) :
            ifig+=1
            if corner and iccd==2 :
                ifig+=1
            ax=fig.add_subplot(3,3,ifig)
            for ichan in amplifier_list  :
                # when we give the list to plot it's in general for raft , not corner raft , so we limit to corner raft channels if needed 
                if ichan>cte_s_data[0,iraft,iccd].first_file.HduMax : break
                ch=ichan-1
                #
                for irun in range(number_of_run) :
                    if ch==0 or irun==0 :
                        labelx='%d (%s)' % (ch+1,run[irun])
                        plt.errorbar(cte_p_data[irun,iraft,iccd].cte_flux_s[ch,0:cte_p_data[irun,iraft,iccd].lmax[ch]],cte_p_data[irun,iraft,iccd].ylev[ch,0:cte_p_data[irun,iraft,iccd].lmax[ch]], yerr=cte_p_data[irun,iraft,iccd].ylev_std[ch,0:cte_p_data[irun,iraft,iccd].lmax[ch]],fmt=sel_mark[irun%nb_mark],color=colors[ch],label=labelx)
                    else :    
                        plt.errorbar(cte_p_data[irun,iraft,iccd].cte_flux_s[ch,0:cte_p_data[irun,iraft,iccd].lmax[ch]],cte_p_data[irun,iraft,iccd].ylev[ch,0:cte_p_data[irun,iraft,iccd].lmax[ch]], yerr=cte_p_data[irun,iraft,iccd].ylev_std[ch,0:cte_p_data[irun,iraft,iccd].lmax[ch]],fmt=sel_mark[irun%nb_mark],color=colors[ch])
                # 
            plt.plot(xx,yy,'g')
            label='CCD %s' % sensors[iccd]
            plt.xlabel(label)
            if iccd%3 == 0 :
                plt.ylabel('1-CTE (// CTE)')
            plt.xscale('log')
            plt.yscale('log')
            if corner :
                if iccd%2 ==1 : plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
            elif iccd==2 :
                plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
            plt.xlim(xx[0],xx[1])
            plt.ylim(1e-8,5e-5)
        if plot_on_screen : plt.show()
        root_plt=os.path.join(root_dir,run[0],raft_cur)
        # create the directorty
        os.makedirs(root_plt,exist_ok=True)
        #
        plotfile=os.path.join(root_plt,'ParallelCTE.png')            
        fig.savefig(plotfile,bbox_inches='tight')
        if not(plot_on_screen) : plt.close(fig)
        # link in all other run directory 
        for irun in range(1,number_of_run) :
            root_plt=os.path.join(root_dir,run[irun],raft_cur)
            # create the directorty
            os.makedirs(root_plt,exist_ok=True)
            # make the link
            os.symlink(plotfile,root_plt+'/ParallelCTE.png')


# In[ ]:





# In[ ]:




