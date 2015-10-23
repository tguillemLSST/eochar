"""
Pixel covariance calculator

Compute the pixel covariance matrix on a pair of
flat fields, taking masks into account.
If (k,l) is the (row,col) position of a pixel,
its pixel covariance is
Mean(Diff(j,i)*Diff(j-l,i-k)) - Mean(Diff(j,i))*Mean(Diff(j-l,i-k))
where Diff stands for the (flat1-flat2) array, and "*" in the
first term is an element-wise product.
"""

from lsst.eotest.sensor.MaskedCCD import MaskedCCD
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import numpy as np
import os

def create_diff(ccd_1, ccd_2, amp=5):
    """
    returns the difference of two flat images from the sensor
    at position amp, after standard unbiasing and trimming
    """
    image1 = ccd_1.unbiased_and_trimmed_image(amp)
    image2 = ccd_2.unbiased_and_trimmed_image(amp)

    #subtract the two images
    diff_image = afwImage.MaskedImageF(image1, True)
    diff_image -= image2
    return diff_image

def variance_threshold_masking(fdiff, sigma):
    """
    simple masking based on the condition variance>sigma
    """
    arr = fdiff.getImage().getArray()
    stats = afwMath.makeStatistics(fdiff, afwMath.MEDIAN | afwMath.VARIANCECLIP)
    diffmed = stats.getValue(afwMath.MEDIAN)
    diffvar = stats.getValue(afwMath.VARIANCECLIP)
    thresh = sigma*np.sqrt(diffvar)
    masked_arr = np.ma.masked_where\
        (abs(arr-diffmed) > thresh, arr, copy=True)
    return masked_arr

def pixel_covar(arr, k_idx, l_idx, verbose):
    """
    compute the (k,l)-pixel covariance value
    """
    if k_idx == 0 and l_idx == 0:
        return np.var(arr)

    nrow, ncol = arr.shape
    temp1 = arr[l_idx:nrow, k_idx:ncol]
    data1 = temp1.data
    mask1 = temp1.mask
    temp2 = arr[0:nrow-l_idx, 0:ncol-k_idx]
    data2 = temp2.data
    mask2 = temp2.mask
    or_mask = np.logical_or(mask1, mask2)
    sub1 = np.ma.MaskedArray(data1, or_mask)
    sub2 = np.ma.MaskedArray(data2, or_mask)
    res = (sub1*sub2).mean() - sub1.mean()*sub2.mean()
    if verbose:
        print res
    return res

def run_pixel_covar(ccd_1, ccd_2, amp, npix, sigma):
    """
    computes the npix x npix array of pixel covariance.
    """
    if ccd_1.md.get('EXPTIME') != ccd_2.md.get('EXPTIME'):
        raise RuntimeError("Exposure times for files %s, %s do not match"\
                               %(file1, file2))

    #diff the 2 flats, for the given amp
    fdiff = create_diff(ccd_1, ccd_2, amp)

    #create a mask based on a simple thresholding variance
    fdiff_mask = variance_threshold_masking(fdiff, sigma)

    res_array = np.zeros((npix, npix))
    for k in range(0, npix):
        for l in range(0, npix):
            res_array[k, l] = pixel_covar(fdiff_mask, k, l, verbose=False)
    return res_array

if __name__ == "__main__":

    from lsst.eotest.sensor.sim_tools import simulateFlat

    # create simulated flats, or use local ones if alerady created
    file1 = "test_flat1.fits"
    file2 = "test_flat2.fits"
    if not os.path.exists(file1):
        simulateFlat(file1, 15000, 5, hdus=16)
    if not os.path.exists(file2):
        simulateFlat(file2, 15000, 5, hdus=16)

    ccd1 = MaskedCCD(file1)
    ccd2 = MaskedCCD(file2)


    #run the pixel covariance for 25 pixels : 5x5 (k,l) pairs
    results = run_pixel_covar(ccd1, ccd2, amp=5, npix=5, sigma=3)
    print results

