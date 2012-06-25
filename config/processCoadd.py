from lsst.obs.sdss.calibrate import SdssCalibrateTask
root.calibrate.retarget(SdssCalibrateTask)

# Inputs to coadd have already been background-subtracted; do we need to do it again?
# (answer obviously changes if we implement background-matching).
# Defaults here just do one, final step to tweak it up, assuming that there's no big
# DC value to begin with.
root.calibrate.doBackground = False
root.calibrate.detection.reestimateBackground = False
root.detection.reestimateBackground = True

useMatchedToPsf = True

for filterName in ("u", "g", "r", "i", "z"):
    subConfig = getattr(root.calibrate, filterName)
    subConfig.psfDeterminer["pca"].spatialOrder    = 1  # Should be spatially invariant
    subConfig.psfDeterminer["pca"].kernelSizeMin   = 31 # Larger Psfs
    subConfig.starSelector["secondMoment"].fluxLim = 3000.0
    # If useMatchedToPsf=True, the previous three lines only get used for determining the sources that
    # go into aperture correction; we determine the PSF then throw it away
    subConfig.useInputPsf = useMatchedToPsf

root.calibrate.repair.doInterpolate = False
root.calibrate.repair.doCosmicRay = False

root.calibrate.astrometry.forceKnownWcs = True
root.calibrate.astrometry.solver.calculateSip = False

 # Remove flags.pixel.interpolated.any
root.calibrate.computeApCorr.badFlags = ("flags.pixel.edge", "flags.pixel.saturated.any")
root.calibrate.photocal.badFlags = ('flags.pixel.edge','flags.pixel.saturated.any')

# JFB: this wasn't being set before #2188, but it probably should have been changed when the other detection
# threshold was.
root.calibrate.detection.thresholdType = "pixel_stdev"

root.detection.thresholdType = "pixel_stdev"
