from lsst.obs.sdss.calibrate import SdssCalibrateTask
root.calibrate.retarget(SdssCalibrateTask)

# Inputs to coadd have already been background-subtracted; do we need to do it again?
# (answer obviously changes if we implement background-matching).
# Defaults here just do one, final step to tweak it up, assuming that there's no big
# DC value to begin with.
root.calibrate.doBackground = True
root.calibrate.detection.reEstimateBackground = True
root.detection.reEstimateBackground = True
# Official Summer2012 background binSize
root.calibrate.background.binSize = 512
root.calibrate.detection.background.binSize = 512
root.detection.background.binSize = 512

# Setting this to False will result in a simple Gaussian model PSF 
root.calibrate.useExposurePsf = False
root.calibrate.doPsf = False
root.calibrate.initialPsf.model='DoubleGaussian'

# The settings below only matter for determining stars to pass to the aperture corrections unless
# useInputPsf is not True.
for filterName in ("u", "g", "r", "i", "z"):
    subConfig = getattr(root.calibrate, filterName)
    subConfig.psfDeterminer["pca"].spatialOrder    = 2  
    subConfig.psfDeterminer["pca"].kernelSizeMin   = 31 # Larger Psfs
    subConfig.starSelector["secondMoment"].fluxLim = 3000.0
    subConfig.starSelector.name = "catalog"

root.calibrate.astrometry.forceKnownWcs = True
root.calibrate.astrometry.solver.calculateSip = False

# Remove flags.pixel.interpolated.any
root.calibrate.computeApCorr.badFlags = ("flags.pixel.edge", "flags.pixel.saturated.any")
root.calibrate.computeApCorr.order = 0
root.calibrate.photocal.badFlags = ('flags.pixel.edge','flags.pixel.saturated.any')

# Official config for Summer 2012
root.calibrate.detection.thresholdType = "pixel_stdev"
root.detection.thresholdType = "pixel_stdev"

# For detection
root.calibrate.initialPsf.fwhm=1.7

try:
    # Enable multiShapelet for model mags.
    import lsst.meas.extensions.multiShapelet
    root.measurement.algorithms.names += ("multishapelet.psf", "multishapelet.exp", "multishapelet.dev", 
                                          "multishapelet.combo")
    root.measurement.slots.modelFlux = "multishapelet.combo.flux"
    # PSF should be exactly double-Gaussian (zeroth-order shapelet)
    root.measurement.algorithms["multishapelet.psf"].innerOrder = 0
    root.measurement.algorithms["multishapelet.psf"].outerOrder = 0
    # too many INTERP pixels on coadds, so we relax the masking in modeling
    for name in ("exp", "dev", "combo"):
        root.measurement.algorithms["multishapelet." + name].badMaskPlanes = ["EDGE", "SAT"]
except ImportError:
    # TODO: find a better way to log this
    print "WARNING: Could not import lsst.meas.extensions.multiShapelet; model fluxes not enabled!"
