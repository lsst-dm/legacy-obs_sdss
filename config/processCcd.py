from lsst.obs.sdss.sdssNullIsr import SdssNullIsrTask

# Read post-ISR data from fpC, fpM, asTrans, etc. files and remove overlap
config.isr.retarget(SdssNullIsrTask)

# Use the PSF determined by SDSS, without trying to fit anything better,
# because the exposures cover too small a region to do a good job of PSF fitting
config.charImage.doMeasurePsf = False

config.charImage.repair.doInterpolate = False
config.charImage.repair.doCosmicRay = False

config.charImage.background.binSize = 512
config.charImage.detectAndMeasure.detection.includeThresholdMultiplier = 10.0
config.charImage.detectAndMeasure.detection.background.binSize = 512
config.charImage.detectAndMeasure.detection.background.binSize = 512
config.charImage.detectAndMeasure.measurement.algorithms.names = [
    "base_SdssCentroid",
    "base_PixelFlags",
    "base_SdssShape",
    "base_PsfFlux",
    "base_CircularApertureFlux",
]
config.charImage.detectAndMeasure.measurement.algorithms["base_CircularApertureFlux"].radii = [7.0]
config.charImage.detectAndMeasure.measurement.slots.centroid = "base_SdssCentroid"
config.charImage.detectAndMeasure.measurement.slots.apFlux = "base_CircularApertureFlux_7_0"
config.charImage.detectAndMeasure.measurement.slots.modelFlux = None
config.charImage.detectAndMeasure.measurement.slots.instFlux = None
config.charImage.detectAndMeasure.measurement.slots.calibFlux = "base_CircularApertureFlux_7_0"
config.charImage.detectAndMeasure.measureApCorr.refFluxName = "slot_CalibFlux"

config.calibrate.detectAndMeasure.detection.background.binSize = 512
config.calibrate.detectAndMeasure.detection.background.binSize = 512
# we rarely run PSF determination on SDSS data, so use the output of the star selector instead
config.calibrate.detectAndMeasure.measureApCorr.inputFilterFlag = "calib_psfCandidate"

# use the WCS determined by SDSS (why?)
config.calibrate.astrometry.forceKnownWcs = True

# Ap correction defaults to true happens in the measurement algorithm: 'correctFluxes'
# SDSS the standard aperture correction is quoted as out to a radius of 7.43.
# According to http://www.sdss.org/dr7/algorithms/photometry.html#photo_profile this corresponds to 18.58.
# Note that 7.43/0.3961270 = 18.7566 <> 18.58. Why?

#psf flux = ap flux at this radius. Will also be applied to galaxies Same everywhere'
#config.calibrate.detectAndMeasure.measurement.algorithms['correctfluxes'].apCorrRadius = 18.58 #pixels
