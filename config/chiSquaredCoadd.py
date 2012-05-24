# overrides for pipe_tasks ChiSquaredCoaddTask.ConfigClass
from lsst.obs.sdss.selectSdssImages import SelectSdssImagesTask

root.select.retarget(SelectSdssImagesTask)
