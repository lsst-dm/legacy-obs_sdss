# overrides for pipe_tasks CoaddTask.ConfigClass
from lsst.obs.sdss.selectSdssImages import SelectSdssImagesTask

root.select.retarget(SelectSdssImagesTask)
