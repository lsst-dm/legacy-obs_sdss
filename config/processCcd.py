from lsst.obs.sdss.calibrate import SdssCalibrateTask
root.calibrate.retarget(SdssCalibrateTask)

import lsst.meas.astrom.catalogStarSelector
root.calibrate.u.starSelector.name = "catalog"
root.calibrate.g.starSelector.name = "catalog"
root.calibrate.r.starSelector.name = "catalog"
root.calibrate.i.starSelector.name = "catalog"
root.calibrate.z.starSelector.name = "catalog"

# Use the PSF determined by SDSS for calibration, without trying to fit anything better,
# because the exposures cover too small a region to do a good job of PSF fitting
root.calibrate.useExposurePsf = True
root.calibrate.doPsf = False

root.calibrate.background.binSize = 512
root.calibrate.detection.background.binSize = 512
root.detection.background.binSize = 512

root.loadSdssWcs = True
root.calibrate.astrometry.forceKnownWcs = True
root.calibrate.astrometry.solver.calculateSip = False

