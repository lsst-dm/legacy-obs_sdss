# overrides for pipe_tasks CoaddTask.ConfigClass
from lsst.obs.sdss.selectSdssImages import SelectSdssImagesTask
from lsst.obs.sdss.scaleSdssZeroPoint import ScaleSdssZeroPointTask

root.select.retarget(SelectSdssImagesTask)
root.scaleZeroPoint.retarget(ScaleSdssZeroPointTask)
