from lsst.obs.sdss.calibrate import SdssCalibrateTask
root.calibrate.retarget(SdssCalibrateTask)

import lsst.meas.astrom.catalogStarSelector
root.calibrate.repair.doInterpolate = False
root.calibrate.repair.doCosmicRay = False
root.calibrate.u.starSelector.name = "catalog"
root.calibrate.z.starSelector.name = "catalog"
root.detection.background.binSize = 512

import lsst.meas.extensions.multiShapelet
root.measurement.algorithms.names += ("multishapelet.psf", "multishapelet.exp", "multishapelet.dev", 
                                      "multishapelet.combo")
root.measurement.apCorrFluxes.names += ("multishapelet.exp", "multishapelet.dev", "multishapelet.combo")
root.measurement.slots.modelFlux = "multishapelet.combo.flux"
