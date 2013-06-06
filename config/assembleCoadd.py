# overrides for pipe_tasks CoaddTask.ConfigClass
from lsst.obs.sdss.selectSdssImages import SelectSdssImagesTask
from lsst.obs.sdss.scaleSdssZeroPoint import ScaleSdssZeroPointTask

root.select.retarget(SelectSdssImagesTask)
root.scaleZeroPoint.retarget(ScaleSdssZeroPointTask)

root.matchBackgrounds.usePolynomial=True
root.matchBackgrounds.binSize=128
root.matchBackgrounds.order=4
root.subregionSize=(2500, 2500)
root.sigmaClip=5
root.maxMatchResidualRatio=1.7
root.maxMatchResidualRMS=1.0
root.autoReference=False

#Configs for deep coadd
root.coaddName='deep'
root.select.maxFwhm=2.0
root.select.quality=2

