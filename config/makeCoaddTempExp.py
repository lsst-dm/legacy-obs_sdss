# overrides for pipe_tasks CoaddTask.ConfigClass
from lsst.obs.sdss.selectSdssImages import SelectSdssImagesTask

root.select.retarget(SelectSdssImagesTask)

root.doOverwrite=True

#Configs for deep coadd
root.coaddName='deep'
root.select.maxFwhm=2.0
root.select.quality=2
