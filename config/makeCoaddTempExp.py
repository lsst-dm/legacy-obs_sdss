# overrides for pipe_tasks CoaddTask.ConfigClass
from lsst.obs.sdss.selectSdssImages import SelectSdssImagesTask

root.select.retarget(SelectSdssImagesTask)

#configs for deep coadd
root.doOverwrite=True
root.coaddName='deep'
root.select.maxFwhm=2.0
root.select.quality=2
