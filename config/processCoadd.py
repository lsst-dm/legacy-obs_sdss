config.calibrate.doBackground = True
config.calibrate.detection.reEstimateBackground = True
config.detection.reEstimateBackground = True

# Official Summer2012
config.calibrate.background.binSize = 512
config.calibrate.detection.background.binSize = 512
config.detection.background.binSize = 512
config.calibrate.detection.thresholdType = "pixel_stdev"
config.detection.thresholdType = "pixel_stdev"

#set in #2800 default already
config.calibrate.astrometry.forceKnownWcs = True

# Remove flags.pixel.interpolated.any
config.calibrate.photocal.badFlags = ('flags.pixel.edge','flags.pixel.saturated.any')

"""Ap correction defaults to true happens in the measurement algorithm: 'correctFluxes'
SDSS the standard aperture correction is quoted as out to a radius of 7.43.
According to http://www.sdss.org/dr7/algorithms/photometry.html#photo_profile this corresponds to 18.58.
Note that 7.43/0.3961270 = 18.7566 <> 18.58.Why?"""

#psf flux = ap flux at this radius. Will also be applied to galaxies Same everywhere'
config.measurement.algorithms['correctfluxes'].apCorrRadius = 18.58 #pixels

config.detection.thresholdValue = 3.0
config.doDeblend=True
config.deblend.maxNumberOfPeaks=40
