"""
SDSS-specific overrides for RunIsrTask
"""
from lsst.obs.sdss.sdssNullIsr import SdssNullIsrTask

config.isr.retarget(SdssNullIsrTask)
