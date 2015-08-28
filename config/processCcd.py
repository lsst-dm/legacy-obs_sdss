from lsst.obs.sdss.calibrate import SdssCalibrateTask
config.calibrate.retarget(SdssCalibrateTask)

import lsst.meas.astrom.catalogStarSelector
config.calibrate.u.starSelector.name = "catalog"
config.calibrate.g.starSelector.name = "catalog"
config.calibrate.r.starSelector.name = "catalog"
config.calibrate.i.starSelector.name = "catalog"
config.calibrate.z.starSelector.name = "catalog"

# Use the PSF determined by SDSS for calibration, without trying to fit anything better,
# because the exposures cover too small a region to do a good job of PSF fitting
config.calibrate.useExposurePsf = True
config.calibrate.doPsf = False

config.calibrate.background.binSize = 512
config.calibrate.detection.background.binSize = 512
config.detection.background.binSize = 512

config.loadSdssWcs = True
config.calibrate.astrometry.forceKnownWcs = True

#Important for producing coaddTempExps without nans
config.removeOverlap=False

"""Ap correction defaults to true happens in the measurement algorithm: 'correctFluxes'
SDSS the standard aperture correction is quoted as out to a radius of 7.43.
According to http://www.sdss.org/dr7/algorithms/photometry.html#photo_profile this corresponds to 18.58.
Note that 7.43/0.3961270 = 18.7566 <> 18.58.Why?"""

#psf flux = ap flux at this radius. Will also be applied to galaxies Same everywhere'
#config.measurement.algorithms['correctfluxes'].apCorrRadius = 18.58 #pixels
