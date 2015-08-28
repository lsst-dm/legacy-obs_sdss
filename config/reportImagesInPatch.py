# overrides for pipe_tasks ReportImagesToCoaddTask.ConfigClass
from lsst.obs.sdss.selectSdssImages import SelectSdssImagesTask

config.select.retarget(SelectSdssImagesTask)
