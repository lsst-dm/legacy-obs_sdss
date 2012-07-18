# overrides for pipe_tasks OutlierRejectedCoaddTask.ConfigClass
from lsst.obs.sdss.selectSdssImages import SelectSdssImagesTask

root.select.retarget(SelectSdssImagesTask)

root.desiredFwhm=None
root.doSigmaClip=True
root.sigmaClip=3.0
root.warp.maskWarpingKernelName='bilinear'
root.subregionSize=[1000,1000]
