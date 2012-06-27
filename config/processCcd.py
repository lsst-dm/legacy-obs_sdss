from lsst.obs.sdss.calibrate import SdssCalibrateTask
root.calibrate.retarget(SdssCalibrateTask)

import lsst.meas.astrom.catalogStarSelector
root.calibrate.u.starSelector.name = "catalog"
root.calibrate.z.starSelector.name = "catalog"
root.calibrate.u.useInputPsf = True
root.calibrate.z.useInputPsf = True
root.detection.background.binSize = 512

try:
    import lsst.meas.extensions.multiShapelet
    root.measurement.algorithms.names += ("multishapelet.psf", "multishapelet.exp", "multishapelet.dev", 
                                          "multishapelet.combo")
    root.measurement.apCorrFluxes += ("multishapelet.exp.flux", "multishapelet.dev.flux",
                                      "multishapelet.combo.flux")
    root.measurement.slots.modelFlux = "multishapelet.combo.flux"
except ImportError:
    # TODO: find a better way to log this
    print "WARNING: Could not import lsst.meas.extensions.multiShapelet; model fluxes not enabled!"
    
