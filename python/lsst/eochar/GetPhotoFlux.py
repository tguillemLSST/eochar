import numpy as np

def GetPhotoFlux(filename,ExpTime):
    """
    Estimate from the monitoring photodiode current at plateau and the EXPTIME ( from the fits header of the exposure )
    an integrated flux and error . This  integrated flux is expected to be linearly related to the flux received by the CCDs .
    If the exposure duration estimated from the photo_diode data is ~2.5 'photodiode time slice' away from the requested EXPTIME
    , then this 'photodiode time' is used instead of EXPTIME to compute the flux , and the estimated error on the 
    measured flux is enlarged.
     Parameters
    ----------
     filename: str  , name of photodiode readings file. 
      ExpTime: float, Exposure lenght as read from the EXPTIME fits header keyword  
    Returns
    -------
    float: The integrated flux  
    float: and an estimate of its sigma .   
    """


    try :
        # read the data 
        photodiode=np.loadtxt(filename)
        # number of measurement 
        nb_measure=len(photodiode[:,0])
        # spacing between each photo diode measurement 
        delta=photodiode[1:,0]-photodiode[:-1,0]
        # Estimate the location  of photodiode data  =================
        # step 1 : estimate the photodiode noize in area without signal : easy just take photodiode flux <0  
        neg_photo=np.array([photodiode[i,1] for i in range(nb_measure) if photodiode[i,1]<0.])
        if len(neg_photo)<5 :
            # conservative guess if not enough data without flux 
            zero=5E-12 
            zero_sig=1E-12
        else : 
            # remark that we will not cut at 4 sigma , but at 4 time the largest excurtion : 
            # the sample is small it's safer ...we are just looking for a scale here  
            zero=max(abs(neg_photo))*4
            zero_sig=neg_photo.std()
            if zero_sig > 1e-12 : 
                print('WARNING : in file %s photo diode noise in area without signal is large %f ' % (filename,zero_sig))   
        # Step 2: estimate the start and  end of the zone with flux : 
        # Step 2.1 : start : step between no significat signal to significant signal 
        #  4 sigma cut , see above in the def of "zero"  , but in case  no step found look for less significant step 
        for j in range(1,4) :
            firstl=[i for i in range(nb_measure-1) if photodiode[i,1]<=zero/j and photodiode[i+1,1]>zero/j ]
            if len(firstl)>0 : break
        # Step 2.2 : end : step between significant signal to no significat signal 
        for j in range(1,4) :                
            lastl=[i+1 for i in range(nb_measure-1) if photodiode[i,1]>zero/j and photodiode[i+1,1]<=zero/j ]
            if len(lastl)>0 : break
        # step 2.3 : No step found : photodiode recording doesn't includes "no-signal area" ( at least 1 example found for the start  in run 12877 ),
        # force the start or end
        if len(firstl)<1 :
            firstl=[0]
        if len(lastl)<1 :
            lastl=[nb_measure-1]
        # step 2.4 : if more than one "start" or "end" step , slect the  one with the largest physical signal 
        first=firstl[0]+1
        last=lastl[0]-1
        # fix issues in case of photodiode glitch : select edges with largest signal 
        if len(firstl)>1 :
            #print('WARNING : photo diode fluctuation ends to a possible bad edge detection (file %s)'%(filename),firstl,lastl)
            for first_cur in firstl[1:] : 
                if photodiode[first_cur+1,1] > photodiode[first,1] :
                    first=first_cur+1
            #print("selected start signal ",first)
        if len(lastl)>1 :
            #print('WARNING : photo diode fluctuation ends to a possible bad edge detection (file %s)'%(filename),firstl,lastl)
            for last_cur in lastl[1:] : 
                if photodiode[last_cur-1,1] > photodiode[last,1] :
                    last=last_cur-1
            #print("selected end signal ",last)
        # step 2.5 : relax the requirement on significant signal for the  start-1 and end+1 pixel : they could contain signal  
        # compute the noise in signal level data removing first and last pixel ( rise and fall area) 
        photo_noise_signal=photodiode[first+1:last,1].std()
        #
        if first>0 and abs(photodiode[first-1,1])>5.*photo_noise_signal :
            first-=1
        if last<nb_measure-1 and abs(photodiode[last+1,1])>5.*photo_noise_signal :
            last+=1
        # step 3 : compute the flux as measured by photo diode , and the expected precision in this measurement 
        # step 3.1 : compute  the photo_diode mean flux , first and last pixel  with signal excluded , and the associated signal dispeersion
        photo_mean=photodiode[first+1:last,1].mean()
        photo_noise=photodiode[first+1:last,1].std()
        # step 3.2 : do a crude estimate of the  exposure time , using an renomalized length  for the first and last pixel 
        # ( not correct due to the start&ed expectd photodiode signal shape ... but the best  that can be done with the available information ) 
        tot_time_observed=delta[first]*photodiode[first,1]/photo_mean+delta[last]*photodiode[last,1]/photo_mean+photodiode[last,0]-photodiode[first+1,0]
        # step 3.3 : We do expecte that the  requested ExpTime is the best estimate of the exposure time  (= it has an offset but it as the lowest dispersion on the real ExpTime)
        #            Still we observed a few ( =1 in run 12877 : flat_ND_OD1.0_SDSSi_426.0_flat0_021 ) where the real exposure time  was significatly different from the requested one
        #            in this case we consider  that the real exposure time is a function of tot_time_observed (the function  is such that tot_time_observed = ExpTime ...as the ExpTime  is offseted , and the time from the photo-diode,tot_time_observed, is not in "second" unit ... yes yes )
        if delta[last]<.1 :
            Corected_Observed=tot_time_observed*0.99883708+0.03216984
        else :
            Corected_Observed=tot_time_observed*0.99931203+0.11776485           
        OutOfTime=abs((ExpTime-Corected_Observed)/delta[last])            
        if OutOfTime>1.3 :
            # The real exposure Time is significantly # for the exposure time measured by the photodiode 
            # we will fix it , but the precision in the estimated photodiode flux will be degraded ( and biased as the delta[last] above is not perfect as it changes with the sampling )
            # by the worse precision  in the exposure length : ~ 1 photodiode time slice ?   
            PhotoFlux=photo_mean*Corected_Observed
            PhotoFluxStd=np.sqrt((photo_noise/np.sqrt(last-first-1.)*ExpTime)**2+(photo_mean*delta[last])**2)
            print('WARNING : from file %s the flux is estimated with the exposure time estimated from the photodiode data (%f s) instead of the given EXPTIME (%f s)' % (filename,tot_time_observed-delta[last],ExpTime))  
        else :
            PhotoFlux=photo_mean*ExpTime
            PhotoFluxStd=photo_noise/np.sqrt(last-first-1.)*ExpTime
        return PhotoFlux,PhotoFluxStd
    except IOError :
        print('Could not find %s, just hoping it is usesless'%filename)
        return -1
