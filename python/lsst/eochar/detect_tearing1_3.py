import numpy as np

#
def no_tearing(fits,channels=range(1,17),verbose=False):
    # =========================
    # name     :  no_tearing   function
    # author   :  P.Antilogus
    # vesrion  :  1.3  2018/05/20  add a detection using the ratio between hedges pixel at 3 location ( line ~ 0, 1000 , 2000) instead of 2 ( 0 , 2000)
    # version  :  1.2  2018/05/20  extend the overscan considered      
    # version  :  1.1  2018/05/07  fix a typo in the parameter names ( channel ==> channels ) 
    # version  :  1.0  2018/05/05  initial realease 
    # Description  :  Routine to indentify tearing in e2v images
    # Input :
    # fits     = pointer to a pyfits object - file , with the image to test
    # channels = list of ccd channels to test for tearing , default = range(1,17)
    # verbose  = level of message to print , default = 0 (=none ) , option : 0,1,2
    #
    # output :
    # tearing_flag : 0 tearing , 1 no tearing , 2 error as processing data/cannot test tearing on this image( too low flux ? ) 
    #
    #
    # interesting variable inside the code ( and that could be returned by a simple edition of it , see last line ) :
    # k_ratio      : np array ((len(channels),2))   with the tearing "ratios" results for each channels edges
    # k_error      : np array ((len(channels),2))   with the error on k_ratio values
    # above_nb     : number of amplifers edges above tearing cuts 
    #
    # method :
    #  In case of tearing the content of the first and last serial column of each amplifier drop between the line affected by
    # the tearing to the first line of the image .
    # Comparing the ratio of this pixel to the nearest physical pixel , for the  ~ first 100 and ~ last 100 lines
    # is sensitive to any jump due to the tearing .
    #==========================
    #
    # tearing selection cuts 
    #  number of sigma cut to reject non-tearing ( 2 sigma cut is enough as we ask to pass it 10 time )
    nsigma=1.
    # cut on ratio difference : > 5%  , requesting a minimal ratio diff on the other amp side > 2%
    rdif_cut_up=0.05 
    rdif_cut_min=-0.01
    # accpetable number of hedge above ration difference, for a non tearing image 
    above_cut=10
    #
    nb_channels=len(channels)
    above_nb=0
    k_ratio=np.zeros((nb_channels,2))
    k_error=np.zeros((nb_channels,2))
    i_ch=-1
    try:
        for ch in channels :
            i_ch+=1
            ov=fits[ch].data[0:2000,530:].mean(axis=1)
            ratio=[(fits[ch].data[0:2000,11]-ov)/(fits[ch].data[0:2000,10]-ov),(fits[ch].data[0:2000,520]-ov)/(fits[ch].data[0:2000,521]-ov)]
            bad=False
            for i in range(2) :
                # compute the ratio of the first hedge at 3 location
                r_val=np.array([ratio[i][100:200].mean(),ratio[i][950:1050].mean(),ratio[i][1900:2000].mean()])
                # avoid saturation : to be confirmed
                if r_val[0]<1.0 :
                    bad=True
                    break
                e_val=np.array([ratio[i][100:200].std(),ratio[i][950:1050].std(),ratio[i][1900:2000].std()])
                r_max=np.argmax(r_val)
                r_min=np.argmin(r_val)
                # get the difference from the mean and the max
                k_ratio[i_ch,i]=r_val[r_max]-r_val[r_min]
                k_error[i_ch,i]=np.sqrt(e_val[r_max]**2/100.+e_val[r_min]**2/100)
            if bad : continue
           # cut at nsigma the non-tearing
            if (k_ratio[i_ch,0]-rdif_cut_up> nsigma*k_error[i_ch,0] and k_ratio[i_ch,1]-rdif_cut_min> nsigma*k_error[i_ch,1] ) : above_nb+=1
            if (k_ratio[i_ch,1]-rdif_cut_up> nsigma*k_error[i_ch,1] and k_ratio[i_ch,0]-rdif_cut_min> nsigma*k_error[i_ch,0] ) : above_nb+=1
        if above_nb > above_cut :
            tearing_flag=0
        else:
            tearing_flag=1
    except :
        if verbose :
            print('Warning : tearing id failed . Remark : tearing identification cannot run on too low flux images or biases')  
        tearing_flag=2
    if not(tearing_flag) and verbose :
        print('Info :  this image has tearing (detected from %d amplifier edges ) ' % (above_cut))
    # debug output
    return tearing_flag,above_nb,k_ratio,k_error
    # minimal output 
    #return tearing_flag

