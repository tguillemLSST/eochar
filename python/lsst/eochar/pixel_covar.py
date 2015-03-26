
# coding: utf-8

#get_ipython().magic(u'pylab inline')

from lsst.eotest.sensor.MaskedCCD import MaskedCCD
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import numpy as np
from math import sqrt
import os

# Ideally, we make use of the standard eotest interface.
# This requires at this stage that the 16 amps be present 
# in the files
file1="test_flat1.fits"
file2="test_flat2.fits"
if not os.path.exists(file1): simulateFlat(file1, 15000, 5, hdus=16)
if not os.path.exists(file2): simulateFlat(file2, 15000, 5, hdus=16)

ccd1=MaskedCCD(file1)
ccd2=MaskedCCD(file2)

if ccd1.md.get('EXPTIME') != ccd2.md.get('EXPTIME'):
    raise RuntimeError("Exposure times for files %s, %s do not match"%(file1, file2))


#proof of concept: look only at 1 amp
amp=5

image1 = ccd1.unbiased_and_trimmed_image(amp)
image2 = ccd2.unbiased_and_trimmed_image(amp)

#subtract the two images
fdiff = afwImage.MaskedImageF(image1, True)
fdiff -= image2
ncols,nrows=fdiff.getDimensions()
fdiff_vals=fdiff.getArrays()[0]


# Create a simple mask based on thresholding the variance
Nsigma=3
fdiff_array = fdiff.getArrays()[0]
diffmed  = np.median(fdiff_array)
varclip = lambda x : afwMath.makeStatistics(x, afwMath.VARIANCECLIP).getValue()
diffvar = varclip(fdiff)
thresh = Nsigma*sqrt(diffvar)
fdiff_mask = np.ma.masked_where(abs( fdiff_array-diffmed ) > thresh,fdiff_array,copy=True)

def pixel_covar():
    for k in range(0,5):
        for l in range(0,5):
            if k==0 and l==0:
                print np.var(fdiff_mask)
            else:
                temp1=fdiff_mask[l:nrows  ,   k:ncols]
                data1=temp1.data
                mask1=temp1.mask
                temp2=fdiff_mask[0:nrows-l  , 0:ncols-k]
                data2=temp2.data
                mask2=temp2.mask
                or_mask=np.logical_or(mask1,mask2)
                npixused1=or_mask.size-or_mask.sum()
                sub1=np.ma.MaskedArray(data1,or_mask)
                sub2=np.ma.MaskedArray(data2,or_mask)
                sum11=sub1.sum()
                sum21=sub2.sum()
                sum121=(sub1*sub2).sum()
                #print (sum121 - sum11*sum21/npixused1)/npixused1


from numba import double
from numba.decorators import jit, autojit
numba_version = autojit(pixel_covar)

#numba_version()
pixel_covar()

# Finally compute the pixel covariance
# Clearly the code is not good here, as the nested loop take time.
# Proof of concept : set k and l to fixed values:
k=1
l=1
i_range = range(k,ncols)
j_range = range(l,nrows)

if False:
    sum1 = 0.; sum2 = 0.; sum12 = 0.
    npixused = 0
    ntotal = 0
    for j in j_range:
        for i in i_range:
            ntotal += 1
            val1 = fdiff_mask[j,i]; 
            val2 = fdiff_mask[j-l,i-k]
            if (val1 is not np.ma.masked) and (val2 is not np.ma.masked) :
                #print i,j,val1,sub1[j-l,i-k],val2,sub2[j-l,i-k], sub1[j-l,i-k] is not np.ma.masked, sub2[j-l,i-k] is not np.ma.masked
                sum1 += val1; 
                sum2 += val2
                sum12 += val1*val2
                npixused += 1


    print  (sum12 - sum1*sum2/npixused)/npixused
    print np.var(fdiff_mask)
    print np.var(fdiff_array)


