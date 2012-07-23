# overrides for pipe_tasks OutlierRejectedCoaddTask.ConfigClass
from lsst.obs.sdss.selectSdssImages import SelectSdssImagesTask

root.select.retarget(SelectSdssImagesTask)

# Not psf-matching for Summer2012 as no deblender and the noise is better behaved
root.desiredFwhm=None

# For Summer2012, this is needed to fix the background striping due to insufficient background modelling
root.doSigmaClip=True
root.sigmaClip=3.0

# Use smaller warping of mask plane
root.warp.maskWarpingKernelName='bilinear'

# Optimization for memory
root.subregionSize=[1000,1000]

# Use quality 2 = acceptable (3 = good, 1 = bad) or better for coadd
root.select.quality=2
