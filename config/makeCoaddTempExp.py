# overrides for pipe_tasks CoaddTask.ConfigClass
from lsst.obs.sdss.selectSdssImages import SelectSdssImagesTask

config.select.retarget(SelectSdssImagesTask)

config.doOverwrite = True

# Configs for deep coadd
config.coaddName = 'deep'
config.select.maxFwhm = 2.0
config.select.quality = 2
