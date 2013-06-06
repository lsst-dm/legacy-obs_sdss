root.calibrate.doBackground = True
root.calibrate.detection.reEstimateBackground = True
root.detection.reEstimateBackground = True

# Official Summer2012
root.calibrate.background.binSize = 512
root.calibrate.detection.background.binSize = 512
root.detection.background.binSize = 512
root.calibrate.detection.thresholdType = "pixel_stdev"
root.detection.thresholdType = "pixel_stdev"

#set in #2800 default already
root.calibrate.astrometry.forceKnownWcs = True
root.calibrate.astrometry.solver.calculateSip = False

# Remove flags.pixel.interpolated.any
root.calibrate.photocal.badFlags = ('flags.pixel.edge','flags.pixel.saturated.any')

"""Ap correction defaults to true happens in the measurement algorithm: 'correctFluxes'
SDSS the standard aperture correction is quoted as out to a radius of 7.43.
According to http://www.sdss.org/dr7/algorithms/photometry.html#photo_profile this corresponds to 18.58.
Note that 7.43/0.3961270 = 18.7566 <> 18.58.Why?"""

#psf flux = ap flux at this radius. Will also be applied to galaxies Same everywhere'
root.measurement.algorithms['correctfluxes'].apCorrRadius = 18.58 #pixels

root.detection.thresholdValue = 3.0
root.doDeblend=True
root.deblend.maxNumberOfPeaks=40


try:
    # Enable multiShapelet for model mags.
    import lsst.meas.extensions.multiShapelet
    root.measurement.algorithms.names |= lsst.meas.extensions.multiShapelet.algorithms
    root.measurement.slots.modelFlux = "multishapelet.combo.flux"
    root.measurement.algorithms["multishapelet.psf"].innerOrder = 0
    root.measurement.algorithms["multishapelet.psf"].outerOrder = 0
    # too many INTERP pixels on coadds, so we relax the masking in modeling
    for name in ("exp", "dev", "combo"):
        root.measurement.algorithms["multishapelet." + name].badMaskPlanes = ["EDGE", "SAT"]
except ImportError:
    # TODO: find a better way to log this
    print "WARNING: Could not import lsst.meas.extensions.multiShapelet; model fluxes not enabled!"
