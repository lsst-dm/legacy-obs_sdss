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

#Important for producing coaddTempExps without nans
root.removeOverlap=False

"""Ap correction defaults to true happens in the measurement algorithm: 'correctFluxes'
SDSS the standard aperture correction is quoted as out to a radius of 7.43.
According to http://www.sdss.org/dr7/algorithms/photometry.html#photo_profile this corresponds to 18.58.
Note that 7.43/0.3961270 = 18.7566 <> 18.58.Why?"""

#psf flux = ap flux at this radius. Will also be applied to galaxies Same everywhere'
root.measurement.algorithms['correctfluxes'].apCorrRadius = 18.58 #pixels
