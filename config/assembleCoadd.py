# overrides for pipe_tasks CoaddTask.ConfigClass
from lsst.obs.sdss.selectSdssImages import SelectSdssImagesTask
from lsst.obs.sdss.scaleSdssZeroPoint import ScaleSdssZeroPointTask

config.select.retarget(SelectSdssImagesTask)
config.scaleZeroPoint.retarget(ScaleSdssZeroPointTask)

config.matchBackgrounds.usePolynomial = True
config.matchBackgrounds.binSize = 128
config.matchBackgrounds.order = 4
config.subregionSize = (2500, 2500)
config.sigmaClip = 5
config.maxMatchResidualRatio = 1.7
config.maxMatchResidualRMS = 1.0
config.autoReference = False

# Configs for deep coadd
config.coaddName = 'deep'
config.select.maxFwhm = 2.0
config.select.quality = 2
