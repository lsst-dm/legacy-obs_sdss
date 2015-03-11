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
