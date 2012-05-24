# overrides for pipe_tasks OutlierRejectedCoaddTask.ConfigClass
from lsst.obs.sdss.selectSdssImages import SelectSdssImagesTask

root.select.retarget(SelectSdssImagesTask)
