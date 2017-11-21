# overrides for pipe_tasks CoaddTask.ConfigClass
from lsst.obs.sdss.selectSdssImages import SelectSdssImagesTask
from lsst.obs.sdss.scaleSdssZeroPoint import ScaleSdssZeroPointTask

config.select.retarget(SelectSdssImagesTask)
config.scaleZeroPoint.retarget(ScaleSdssZeroPointTask)

config.subregionSize = (2500, 2500)
config.sigmaClip = 5

# Configs for deep coadd
config.coaddName = 'deep'
config.select.maxFwhm = 2.0
config.select.quality = 2
