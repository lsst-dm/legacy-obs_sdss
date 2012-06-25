from lsst.obs.sdss.calibrate import SdssCalibrateTask
root.calibrate.retarget(SdssCalibrateTask)

import lsst.meas.astrom.catalogStarSelector
root.calibrate.repair.doInterpolate = False
root.calibrate.repair.doCosmicRay = False
root.calibrate.u.starSelector.name = "catalog"
root.calibrate.z.starSelector.name = "catalog"
