from lsst.obs.sdss.calibrate import SdssCalibrateTask
root.calibrate.retarget(SdssCalibrateTask)

# Inputs to coadd have already been background-subtracted; do we need to do it again?
# (answer obviously changes if we implement background-matching).
# Defaults here just do one, final step to tweak it up, assuming that there's no big
# DC value to begin with.
root.calibrate.doBackground = False
root.calibrate.detection.reestimateBackground = False
root.detection.reestimateBackground = True

root.calibrate.u.useInitialPsf = True
root.calibrate.g.useInitialPsf = True
root.calibrate.r.useInitialPsf = True
root.calibrate.i.useInitialPsf = True
root.calibrate.z.useInitialPsf = True

root.calibrate.astrometry.forceKnownWcs = True
root.calibrate.astrometry.solver.calculateSip = False
