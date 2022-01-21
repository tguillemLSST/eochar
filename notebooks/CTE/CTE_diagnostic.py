#!/usr/bin/env python
# coding: utf-8

# # Notebook : CTE_diagnostic
# 
# Goal : This Notebook will generate plots allowing to perform  a diagnostic for channels with a bad CTE .  
# 
# Author : P.Antilogus
# 
# Version : 4th July 2019 13:30 
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

import matplotlib
from eochar.bot_frame_op import *

# CONFIGURATION FOR THE CURRENT EXECUTION  ========================================
# ---- raft and associated run ============ To be updated if needed 
# 'RTM' :  full RTM name 
# 'Data_location'  :  0 = BNL_RAFT_ROOT ,  1 = SLAC_RAFT_ROOT 
# 'run' :  list of run , ex : ['9876','9874']  (remark : run is a string )   for all run of this raft : '*' 
# 'sensor' :  list of sensor , ex ['S00','S01','S02']  ,  for all sensor ['S*']
# amplifier : list of amplifier , ex [1,2]  , for 1 to 16 you can use [-1] instead 
#data_to_select={'RTM':'RTM-004','Data_location':1,'run':['*'],'sensor':['S*'],'amplifier':[-1]}
#data_to_select={'RTM':'RTM-004','Data_location':1,'run':['9012'],'sensor':['S01'],'amplifier':[4]}
#data_to_select={'RTM':'RTM-004','Data_location':1,'run':['9012'],'sensor':['S*'],'amplifier':[-1]}
data_to_select={'RTM':'RTM-017','Data_location':1,'run':['11166'],'sensor':['S*'],'amplifier':[-1]}
# Do we look for serial CTE ?
serial=True
# do we plot per amplifiers serial CTE diagnostic plots  
serial_plot=True
# Do we look for // CTE ?
parallel=True
# do we plot per amplifiers // CTE diagnostic plots  
parallel_plot=True
# default gain if eotest data not available ( 0.7 is better of e2v , 0.9 for itl ) 
default_gain=1.

import numpy as np
import glob  
try : 
    from get_EO_analysis_results import get_EO_analysis_results
    from get_EO_analysis_files import get_EO_analysis_files
    from exploreFocalPlane import exploreFocalPlane
    from exploreRaft import exploreRaft
    from eTraveler.clientAPI.connection import Connection
    g = get_EO_analysis_results(db='Prod')
    eotest_db=True
except:
    print('No access to eotest DB , so the gain of all devices will be set to 1 ')
    eotest_db=False

# LSST stack imports
from lsst.daf.persistence import Butler

# the list of super flat for a raft - run : BNL_RAFT_ROOT+raft+'/'+run+'/cte_raft/*/*/'+'Sxx'+'/*_superflat_high.fits' 
##
if False:
    count = get_ipython().getoutput('set | grep  .ncsa.illinois.edu | wc')
    if int(count[0].split()[0])> 1 : 
        NCSA=True
    else : 
        NCSA=False
        if NCSA : 
            SLACmirror='/project/rgruendl/SLACmirror'
            # the list of super flat for a raft - run : BNL_RAFT_ROOT+raft+'/'+run+'/cte_raft/*/*/'+'Sxx'+'/*_superflat_high.fits' 
            BNL_RAFT_ROOT=SLACmirror+'/BNLmirror/mirror/BNL-prod/prod/LCA-11021_RTM/LCA-11021_'
            #
            SLAC_RAFT_ROOT=SLACmirror+'/fs3/jh_archive/LCA-11021_RTM/LCA-11021_'
        else :
            # the list of super flat for a raft - run : BNL_RAFT_ROOT+raft+'/'+run+'/cte_raft/*/*/'+'Sxx'+'/*_superflat_high.fits' 
            BNL_RAFT_ROOT='/nfs/farm/g/lsst/u1/mirror/BNL-prod/prod/LCA-11021_RTM/LCA-11021_'
            #
            SLAC_RAFT_ROOT='/gpfs/slac/lsst/fs1/g/data/jobHarness/jh_archive/LCA-11021_RTM/LCA-11021_'
            #
            #
            # order in which they are added to ROOT should follow what will be given to 'data_location' 
            ROOT=[BNL_RAFT_ROOT,SLAC_RAFT_ROOT]

print('read butler')
# activate the butler
repo_path = '/sps/lsst/users/tguillem/Rubin/Focal_Plane/lsst_distrib/w_2022_01/data_PTC'
butler = Butler(repo_path)

# how many run & sensor 
sensor_dir=[]
run_loc=[]
# get all directory  for all sensor and all run 
for run_cur in data_to_select['run'] :
    for sensor_cur  in data_to_select['sensor'] :
        root_dir=ROOT[data_to_select['Data_location']]+data_to_select['RTM']+'/'+run_cur+'/flat_pair_raft_acq/v0/*/'+sensor_cur
        run_loc+=glob.glob(root_dir)
# order them all 
run_loc.sort()
# filter : keep only 1 proceesing for a given run and sensor (the first seen )
run_name_list=[]
sensor_list=[]
sensor_all=[]
sensor_dir=[]
for run_cur in run_loc :
    name_list=run_cur.split('/')
    run_number=name_list[-5]
    if not(run_number in run_name_list ) :
        # new run , let see which sensor are in it 
        run_name_list.append(run_number)
        if len(sensor_all) > len(sensor_list) :
            sensor_list=sensor_all
        sensor_all=[]
    sensor_name=name_list[-1]
    if not(sensor_name in sensor_all) :
        sensor_all.append(sensor_name)
        # we keep the first processed data for this run and sensor
        sensor_dir.append({'dir':run_cur,'run':run_number,'sensor':sensor_name})
if len(sensor_all) > len(sensor_list) :
    sensor_list=sensor_all
sensor_list.sort()    
#
number_of_run=len(run_name_list)
number_of_sensor=len(sensor_list)
#
if number_of_run < 1 : 
    print('No data found for ',data_to_select)
    raise 
#
txt='We will analyse the '
if serial :
    txt+='serial '
    if parallel :
        txt+='and parallel '
elif parallel :
    txt+='parallel '
else :
    print('None of the possible option (serial and/or parallel) have been selected ....')
    raise
print(txt,'CTE for ',number_of_run,' run(s) (',run_name_list,') from the raft ',data_to_select['RTM'])
#how many amplifiers 
if -1 in data_to_select['amplifier'] :
    amplifier_list=np.array(range(1,17))
else :
    amplifier_list=np.array(data_to_select['amplifier'])
amplifier_list.sort()
number_of_amplifier=len(amplifier_list)
if number_of_amplifier ==16 : 
    print('For the 16 amplifiers of ',number_of_sensor,' sensors (',sensor_list,')')
else :
    print('including ',number_of_sensor,' sensor(s) (',sensor_list,') and ',number_of_amplifier,'amplifier(s) (',amplifier_list,')')
# create the data table with the right dimension 
if parallel :
    cte_p_data=np.zeros((number_of_run,number_of_sensor),dtype=np.object_)
if serial : 
    cte_s_data=np.zeros((number_of_run,number_of_sensor),dtype=np.object_)


# In[ ]:


# fill the cte data  and do the key diagnotic plots for each selected channels 
for dir_cur in sensor_dir : 
    irun=run_name_list.index(dir_cur['run'])
    iccd=sensor_list.index(dir_cur['sensor'])
    #
    all_file=glob.glob(dir_cur['dir']+'/*flat*flat*.fits')
    if len(all_file)<1 : 
        continue
    print('Processing ',len(all_file),' flats for run',dir_cur['run'],' and sensor ',dir_cur['sensor'])
    #
    if eotest_db : 
        raft_list, data = g.get_tests(site_type="I&T-Raft", test_type="gain", run =run_name_list[irun])  # get the data for I&T-Raft
        res = g.get_results(test_type='gain', data=data, device=raft_list[0])  # get the data for a raft
     # get the gain for this sensor :
    if eotest_db :
        sensor_id=[]
        for d in data['steps']['fe55_raft_analysis']['fe55_raft_analysis']:
            if 'slot' in d.keys():
                if d['slot']==sensor_list[iccd]:
                    sensor_id=d['sensor_id']
        if len(sensor_id)==0 :
                    all_gain=np.ones((16))*default_gain 
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
        sensor_id=''
        all_gain=np.ones((16))*default_gain 
    #
    all_file.sort()
    file=Ifile(dirall=all_file,fkey={},Slow=False)
    if serial :
        cte_s_data[irun,iccd]=cte(all_file=file.all_file,gain=all_gain,serie=True)
    if parallel: 
        cte_p_data[irun,iccd]=cte(all_file=file.all_file,gain=all_gain,serie=False)
    nf=0
    ccd=data_to_select['RTM']+','+run_name_list[irun]+':'+sensor_list[iccd]
    for ichan in range(number_of_amplifier)  :
        ch=amplifier_list[ichan]-1
        if serial and  serial_plot :
            cte_s_data[irun,iccd].plot_cte(ch=ch,ccd_name=ccd,nf=nf,on_screen=True)
        #
        if parallel and parallel_plot :
            cte_p_data[irun,iccd].plot_cte(ch=ch,ccd_name=ccd,nf=nf,on_screen=True)


# In[ ]:


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
if serial :
    for iccd in range(number_of_sensor) :
        fig=plt.figure(figsize=[20,20])
        fig.suptitle('Serial CTE :'+data_to_select['RTM']+':'+sensor_list[iccd],fontsize=14)
        for ch in range(16):
            ax=fig.add_subplot(4,4,ch+1)
            for irun in range(number_of_run) :
                  plt.errorbar(cte_s_data[irun,iccd].cte_flux_s[ch,0:cte_s_data[irun,iccd].lmax[ch]],cte_s_data[irun,iccd].ylev[ch,0:cte_s_data[irun,iccd].lmax[ch]], yerr=cte_s_data[irun,iccd].ylev_std[ch,0:cte_s_data[irun,iccd].lmax[ch]],fmt=sel_mark[irun%nb_mark],color=sel_color[irun%nb_color],label=run_name_list[irun])
            # 
            plt.plot(xx,yy,'g')
            label='channel %d  : <flux> in e-' % (ch+1)
            plt.xlabel(label)
            if ch%4 == 0 :
                plt.ylabel('1-CTE (serial CTE)')
            plt.xscale('log')
            plt.yscale('log')
            if (ch+1)%4==0 :
                plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
            plt.xlim(xx[0],xx[1])
            plt.ylim(1e-8,5e-5)
        plt.show()
        plotfile=data_to_select['RTM']+':'+sensor_list[iccd]+'_serialCTE.png'
        fig.savefig(plotfile,bbox_inches='tight')
if parallel :
    for iccd in range(number_of_sensor) :
        fig=plt.figure(figsize=[20,20])
        fig.suptitle('// CTE :'+data_to_select['RTM']+':'+sensor_list[iccd],fontsize=14)
        for ch in range(16):
            ax=fig.add_subplot(4,4,ch+1)
            for irun in range(number_of_run) :
                  plt.errorbar(cte_p_data[irun,iccd].cte_flux_s[ch,0:cte_p_data[irun,iccd].lmax[ch]],cte_p_data[irun,iccd].ylev[ch,0:cte_p_data[irun,iccd].lmax[ch]], yerr=cte_p_data[irun,iccd].ylev_std[ch,0:cte_p_data[irun,iccd].lmax[ch]],fmt=sel_mark[irun%nb_mark],color=sel_color[irun%nb_color],label=run_name_list[irun])
            # 
            plt.plot(xx,yy,'g')
            label='channel %d  : <flux> in e-' % (ch+1)
            plt.xlabel(label)
            if ch%4 == 0 :
                plt.ylabel('1-CTE (// CTE)')
            plt.xscale('log')
            plt.yscale('log')
            if (ch+1)%4==0 :
                plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
            plt.xlim(xx[0],xx[1])
            plt.ylim(1e-8,5e-5)
        plt.show()
        plotfile=data_to_select['RTM']+':'+sensor_list[iccd]+'_parallelCTE.png'
        fig.savefig(plotfile,bbox_inches='tight')


# In[ ]:


# Summary plots (cte vs flux ) for the # selected sensors , overlaping the # selected run 
#
# colors for each channel setup
cmap=plt.get_cmap('gist_ncar')
colors=[cmap(j)[:3] for j in np.linspace(0,1,17)]
#
sel_color=['r','g','b','c','m','y','k']
nb_color=len(sel_color)
sel_mark=['<-','>-','*-','^-','v-','+-']
nb_mark=len(sel_mark)
#
yv=5.0e-6
yy=[yv,yv]
xx=[10.,300000.]
#
if serial :
    fig=plt.figure(figsize=[20,20])
    fig.suptitle('Serial CTE :'+data_to_select['RTM'],fontsize=14)
    for iccd in range(number_of_sensor) :
        ax=fig.add_subplot(3,3,iccd+1)
        for ch in range(16):
            for irun in range(number_of_run) :
                if ch==0 or irun==0 :
                    labelx='%d (%s)' % (ch+1,run_name_list[irun])
                    plt.errorbar(cte_s_data[irun,iccd].cte_flux_s[ch,0:cte_s_data[irun,iccd].lmax[ch]],cte_s_data[irun,iccd].ylev[ch,0:cte_s_data[irun,iccd].lmax[ch]], yerr=cte_s_data[irun,iccd].ylev_std[ch,0:cte_s_data[irun,iccd].lmax[ch]],fmt=sel_mark[irun%nb_mark],color=colors[ch],label=labelx)
                else :    
                    plt.errorbar(cte_s_data[irun,iccd].cte_flux_s[ch,0:cte_s_data[irun,iccd].lmax[ch]],cte_s_data[irun,iccd].ylev[ch,0:cte_s_data[irun,iccd].lmax[ch]], yerr=cte_s_data[irun,iccd].ylev_std[ch,0:cte_s_data[irun,iccd].lmax[ch]],fmt=sel_mark[irun%nb_mark],color=colors[ch])
            # 
        plt.plot(xx,yy,'g')
        label='CCD %s' % sensor_list[iccd]
        plt.xlabel(label)
        if iccd%3 == 0 :
            plt.ylabel('1-CTE (serial CTE)')
        plt.xscale('log')
        plt.yscale('log')
        if (iccd+1)%3==0 :
                plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
        plt.xlim(xx[0],xx[1])
        plt.ylim(1e-8,5e-5)
    plt.show()
    plotfile=data_to_select['RTM']+':all_CCD_serialCTE.png'
    fig.savefig(plotfile,bbox_inches='tight')
if parallel :
    fig=plt.figure(figsize=[20,20])
    fig.suptitle('// CTE :'+data_to_select['RTM'],fontsize=14)
    for iccd in range(number_of_sensor) :
        ax=fig.add_subplot(3,3,iccd+1)
        for ch in range(16):
            for irun in range(number_of_run) :
                if ch==0 or irun==0 :
                    labelx='%d (%s)' % (ch+1,run_name_list[irun])
                    plt.errorbar(cte_p_data[irun,iccd].cte_flux_s[ch,0:cte_p_data[irun,iccd].lmax[ch]],cte_p_data[irun,iccd].ylev[ch,0:cte_p_data[irun,iccd].lmax[ch]], yerr=cte_p_data[irun,iccd].ylev_std[ch,0:cte_p_data[irun,iccd].lmax[ch]],fmt=sel_mark[irun%nb_mark],color=colors[ch],label=labelx)
                else :
                    plt.errorbar(cte_p_data[irun,iccd].cte_flux_s[ch,0:cte_p_data[irun,iccd].lmax[ch]],cte_p_data[irun,iccd].ylev[ch,0:cte_p_data[irun,iccd].lmax[ch]], yerr=cte_p_data[irun,iccd].ylev_std[ch,0:cte_p_data[irun,iccd].lmax[ch]],fmt=sel_mark[irun%nb_mark],color=colors[ch])
        # 
        plt.plot(xx,yy,'g')
        label='CCD %s' % sensor_list[iccd]
        plt.xlabel(label)
        if iccd%3 == 0 :
            plt.ylabel('1-CTE (// CTE)')
        plt.xscale('log')
        plt.yscale('log')
        if (iccd+1)%4==0 :
            plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
        plt.xlim(xx[0],xx[1])
        plt.ylim(1e-8,5e-5)
    plt.show()
    plotfile=data_to_select['RTM']+':all_CCD_parallelCTE.png'
    fig.savefig(plotfile,bbox_inches='tight')


# In[ ]:


#           
yv=5.0e-6
yy=[yv,yv]
xx=[10.,300000.]
#sel_color=['r','g','b','c','m','y','k']
nb_color=len(sel_color)
sel_mark=['<','>','*','^','v','+']
nb_mark=len(sel_mark)
flux_sel=[1000,50000]
if serial :
    x=range(cte_s_data[0,0].first,cte_s_data[0,0].first+28)
    for flux in flux_sel:
        for iccd in range(number_of_sensor) :
            fig=plt.figure(figsize=[20,20])
            title='Serial overscan at %d e- %s:%s' % (flux,data_to_select['RTM'],sensor_list[iccd])
            fig.suptitle(title,fontsize=14)
            for ch in range(16):
                ax=fig.add_subplot(4,4,ch+1)
                for irun in range(number_of_run) :
                    ii=np.argmin(np.abs(cte_s_data[irun,iccd].cte_flux_s[ch,:cte_s_data[irun,iccd].lmax[ch]]-flux))
                    yplt=cte_s_data[irun,iccd].cte_y_s[ch,ii,:]
                    #y_min=min(max(min(np.min(yplt)*1.2,0.),-10.),y_min)
                    #y_max=max(min(np.max(yplt)*1.2,100.),y_max)
                    plt.plot(x,yplt,color=sel_color[irun%nb_color],label=run_name_list[irun])
                label='channel %d  : serial overscan ' % (ch+1)
                plt.xlabel(label)
                if ch%4 == 0 :
                    plt.ylabel('flux in e-')
                #plt.xscale('log')
                plt.yscale('symlog')
                if (ch+1)%4==0 :
                    plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
                #plt.xlim(xx[0],xx[1])
                #plt.ylim(,5e-5)
            plt.show()
            plotfile='%s:%sat%d_serialoverscan.png' % (data_to_select['RTM'],sensor_list[iccd],flux)
            fig.savefig(plotfile,bbox_inches='tight')
if parallel :
    x=range(cte_p_data[0,0].first,cte_p_data[0,0].first+28)
    for flux in flux_sel:
        for iccd in range(number_of_sensor) :
            fig=plt.figure(figsize=[20,20])
            title='// overscan at %d e- %s:%s' % (flux,data_to_select['RTM'],sensor_list[iccd])
            fig.suptitle(title,fontsize=14)
            for ch in range(16):
                ax=fig.add_subplot(4,4,ch+1)
                for irun in range(number_of_run) :
                    ii=np.argmin(np.abs(cte_p_data[irun,iccd].cte_flux_s[ch,:cte_s_data[irun,iccd].lmax[ch]]-flux))
                    yplt=cte_p_data[irun,iccd].cte_y_s[ch,ii,:]
                    #y_min=min(max(min(np.min(yplt)*1.2,0.),-10.),y_min)
                    #y_max=max(min(np.max(yplt)*1.2,100.),y_max)
                    plt.plot(x,yplt,color=sel_color[irun%nb_color],label=run_name_list[irun])
                label='channel %d  : // overscan ' % (ch+1)
                plt.xlabel(label)
                if ch%4 == 0 :
                    plt.ylabel('flux in e-')
                #plt.xscale('log')
                plt.yscale('symlog')
                if (ch+1)%4==0 :
                    plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
                #plt.xlim(xx[0],xx[1])
                #plt.ylim(,5e-5)
            plt.show()
            plotfile='%s:%sat%d_paralleloverscan.png' % (data_to_select['RTM'],sensor_list[iccd],flux)
            fig.savefig(plotfile,bbox_inches='tight')


# In[ ]:


# summary for the full raft at low and high flux overlaing the # selected run 
# this box is executed only if all channels of the raft are selected 
#
if number_of_sensor < 9 : 
    print('You should not run this box with less than 9 sensors selected for a given raft')
    raise
# create the associated table with the right dimension 
if parallel :
    cte_p_1000=np.zeros((number_of_run,number_of_sensor,16))
    cte_p_1000_std=np.zeros((number_of_run,number_of_sensor,16))
    cte_p_50000=np.zeros((number_of_run,number_of_sensor,16))
    cte_p_50000_std=np.zeros((number_of_run,number_of_sensor,16))
if serial : 
    cte_s_1000=np.zeros((number_of_run,number_of_sensor,16))
    cte_s_1000_std=np.zeros((number_of_run,number_of_sensor,16))
    cte_s_50000=np.zeros((number_of_run,number_of_sensor,16))
    cte_s_50000_std=np.zeros((number_of_run,number_of_sensor,16))
#
# fill the cte data 
for irun in range(number_of_run) :
    for iccd in range(number_of_sensor) :
         for ch in range(16) :
            # this parameter is to avoid to look for low flux in saturated image
            alf_icur=int(cte_s_data[irun,iccd].lmax[ch]/2)
            #
            if serial :
                ii=np.argmin(np.abs(np.array(cte_s_data[irun,iccd].cte_flux_s[ch,:alf_icur])-1000)) 
                 # as error are large at 1000 I extrapolate from points above 1000 ???? hum
                cte_s_1000[irun,iccd,ch]=cte_s_data[irun,iccd].ylev[ch,ii]
                cte_s_1000_std[irun,iccd,ch]=cte_s_data[irun,iccd].ylev_std[ch,ii]
                ii=np.argmin(np.abs(np.array(cte_s_data[irun,iccd].cte_flux_s[ch,:])-50000)) 
                cte_s_50000[irun,iccd,ch]=cte_s_data[irun,iccd].ylev[ch,ii]
                cte_s_50000_std[irun,iccd,ch]=cte_s_data[irun,iccd].ylev_std[ch,ii]
            #
            if parallel :
                ii=np.argmin(np.abs(np.array(cte_p_data[irun,iccd].cte_flux_s[ch,:alf_icur])-1000)) 
                # as error are large at 1000 I extrapolate from points above 1000 ???? hum
                cte_p_1000[irun,iccd,ch]=cte_p_data[irun,iccd].ylev[ch,ii]
                cte_p_1000_std[irun,iccd,ch]=cte_p_data[irun,iccd].ylev_std[ch,ii]
                ii=np.argmin(np.abs(np.array(cte_p_data[irun,iccd].cte_flux_s[ch,:])-50000)) 
                cte_p_50000[irun,iccd,ch]=cte_p_data[irun,iccd].ylev[ch,ii]
                cte_p_50000_std[irun,iccd,ch]=cte_p_data[irun,iccd].ylev_std[ch,ii]
            
# labels for the plots
label_txt=np.zeros((9),dtype=np.object_)
label_pos=np.zeros((9))
label_chan=np.zeros((18),dtype=np.object_)
label_chan_pos=np.zeros((18))
for jcd in range(3) : 
    for icd in range(3) :
        s='S%d%d' % (jcd,icd)
        iccd=icd+3*jcd
        label_txt[iccd]=s
        label_pos[iccd]=iccd*16+7.5
        label_chan[iccd*2]='%d' % (1)
        label_chan[iccd*2+1]='%d' % (9)
        label_chan_pos[iccd*2]=iccd*16
        label_chan_pos[iccd*2+1]=iccd*16+8
xx1=[0.,144.]
yv1=5.0e-6
yy1=[yv1,yv1]
sel_mark=['<','^','>','v','*','+']
nb_mark=len(sel_mark)
sel_color=['r','g','b','c','m','y','k']
nb_color=len(sel_color)
#
if serial :
    fig=plt.figure(figsize=[15,10])
    fig.suptitle(data_to_select['RTM']+' serial CTE ',fontsize=14)
    ax1=fig.add_subplot(1,1,1)
    ax2=ax1.twiny()
    #
    for irun in range(number_of_run) :
        ax1.plot(np.ravel(cte_s_1000[irun,:,:]),sel_mark[(2*irun)%nb_mark],color=sel_color[(2*irun)%nb_mark],label='low flux'+':'+run_name_list[irun])
        ax1.plot(np.ravel(cte_s_50000[irun,:,:]),sel_mark[(2*irun+1)%nb_mark],color=sel_color[(2*irun+1)%nb_mark],label='high flux'+':'+run_name_list[irun])
    ax1.plot(xx1,yy1,'g')
    ax1.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
    ax1.set_yscale('log')
    ax1.set_ylabel('1-CTE')
    ax1.set_xlabel('Sensor')
    ax1.set_xlim(0,144)
    ax1.set_xticks(label_pos)
    ax1.set_xticklabels(label_txt)
    ax2.set_xlim(ax1.get_xlim())
    ax2.tick_params(labeltop='on')
    ax2.set_xticks(label_chan_pos)
    ax2.set_xticklabels(label_chan)
    ax2.set_xlabel(data_to_select['RTM']+' Channels')
    plt.show()
    plotfile=data_to_select['RTM']+'_serialCTE.png'
    fig.savefig(plotfile,bbox_inches='tight')
if parallel : 
    fig=plt.figure(figsize=[15,10])
    fig.suptitle(data_to_select['RTM']+' // CTE ',fontsize=14)
    ax1=fig.add_subplot(1,1,1)
    ax2=ax1.twiny()
    for irun in range(number_of_run) :
        ax1.plot(np.ravel(cte_p_1000[irun,:,:]),sel_mark[(irun*2)%nb_mark],color=sel_color[(2*irun)%nb_mark],label='low flux'+':'+run_name_list[irun])
        ax1.plot(np.ravel(cte_p_50000[irun,:,:]),sel_mark[(irun*2+1)%nb_mark],color=sel_color[(2*irun+1)%nb_mark],label='high flux'+':'+run_name_list[irun])
    ax1.plot(xx1,yy1,'g')
    ax1.set_xlim(0,144)
    ax1.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
    ax1.set_yscale('log')
    ax1.set_ylabel('1-CTE')
    ax1.set_xlabel('Sensor')
    ax1.set_xticks(label_pos)
    ax1.set_xticklabels(label_txt)
    ax2.set_xlim(ax1.get_xlim())
    ax2.tick_params(labeltop='on')
    ax2.set_xticks(label_chan_pos)
    ax2.set_xticklabels(label_chan)
    ax2.set_xlabel(data_to_select['RTM']+' Channels')
    plt.show()
    plotfile=data_to_select['RTM']+'_paralCTE.png'
    fig.savefig(plotfile,bbox_inches='tight')


# In[ ]:


# serial trap overview (manely usefull for e2v device)
flux_cte_s=[100.,1000.,5000.]
#flux_cte_p=[1000.,50000.]
flux_flag=['<','>','*','v']
sel_color=['r','g','b','c','m','y','k']
nb_color=len(sel_color)
flux_label=['100e-','1ke-','5ke-','1k-100e-']
trap_s_flux=np.zeros((number_of_run,len(flux_cte_s)+1,9,16))
#
xx1=[0.,17.]
yv1=15
yy1=[yv1,yv1]
yv2=1
yy2=[yv2,yv2]

#
if serial : 
    for iccd in range(number_of_sensor) :
        fig=plt.figure(figsize=[15,10])
        fig.suptitle('Serial trap identification :'+data_to_select['RTM']+':'+sensor_list[iccd],fontsize=14)
        ax1=fig.add_subplot(1,1,1)
        for irun in range(number_of_run) :
            for i_flux in range(len(flux_cte_s)) :
                for ch in range(16) :
                    ii=np.argmin(np.abs(cte_s_data[irun,iccd].cte_flux_s[ch,:cte_s_data[irun,iccd].lmax[ch]]-flux_cte_s[i_flux]))
                    trap_s_flux[irun,i_flux,iccd,ch]=cte_s_data[irun,iccd].cte_y_s[ch,ii,0]+cte_s_data[irun,iccd].cte_y_s[ch,ii,1]
                    if i_flux==1 :
                        trap_s_flux[irun,-1,iccd,ch]= trap_s_flux[irun,i_flux,iccd,ch]-trap_s_flux[irun,i_flux-1,iccd,ch]
                ax1.plot(range(1,17),np.ravel(trap_s_flux[irun,i_flux,iccd,:]),flux_flag[i_flux],color=sel_color[irun%nb_color],label=flux_label[i_flux]+':'+run_name_list[irun])
                if i_flux==1 :
                    ax1.plot(range(1,17),np.ravel(trap_s_flux[irun,-1,iccd,:]),flux_flag[-1],color=sel_color[irun%nb_color],label=flux_label[-1]+':'+run_name_list[irun])

        ax1.plot(xx1,yy1,'g')
        ax1.plot(xx1,yy2,'g--')
        ax1.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
        #ax1.set_yscale('log')
        ax1.set_ylabel('Overscan pix1 + pix2 content in e-',fontsize=14)
        ax1.set_xlabel('Channel')
        ax1.set_xlim(0,16)
        plt.show()
        plotfile=data_to_select['RTM']+':'+sensor_list[iccd]+'_trapserial.png'
        fig.savefig(plotfile,bbox_inches='tight')


# In[ ]:





# In[ ]:




