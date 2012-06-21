from lsst.obs.sdss.calibrate import SdssCalibrateTask
root.calibrate.retarget(SdssCalibrateTask)

import lsst.meas.astrom.catalogStarSelector
self.u.starSelector.name = "catalog"
self.z.starSelector.name = "catalog"
