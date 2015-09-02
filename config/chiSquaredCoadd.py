# overrides for pipe_tasks ChiSquaredCoaddTask.ConfigClass
from lsst.obs.sdss.selectSdssImages import SelectSdssImagesTask

config.select.retarget(SelectSdssImagesTask)
