# 
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import lsst.daf.base as dafBase
import lsst.pipe.base as pipeBase
import lsst.pex.config as pexConfig
import lsst.meas.algorithms as measAlg
import lsst.meas.astrom as measAstrom
import lsst.afw.table as afwTable
from lsst.pipe.tasks.calibrate import InitialPsfConfig, CalibrateConfig, CalibrateTask

from lsst.meas.photocal import PhotoCalTask
from lsst.pipe.tasks.astrometry import AstrometryTask
from lsst.pipe.tasks.repair import RepairTask

import math

class SdssCalibratePerFilterConfig(pexConfig.Config):
    starSelector = measAlg.starSelectorRegistry.makeField(
        doc="algorithm to do star/galaxy classification for PSF and aperture corrections",
        default="secondMoment"
    )
    psfDeterminer = measAlg.psfDeterminerRegistry.makeField(
        doc=("algorithm to determine PSF (always run even if we discard the PSF, "
             "to get inputs for aperture correction)"),
        default="pca"
    )
    useInputPsf = pexConfig.Field(
        dtype=bool,
        default=False,
        doc=("If True, discard the PSF produced by the psfDeterminer and use the one "
             "attached to the exposure (this would be the SDSS PSF for SFM, and the analytic "
             "matched-to PSF on coadds).  This also means we skip initial measurement with a "
             "simple guess-PSF.")
    )

    def setDefaults(self):
        self.starSelector["secondMoment"].clumpNSigma = 2.0
        self.psfDeterminer["pca"].nEigenComponents = 4
        self.psfDeterminer["pca"].kernelSize = 7
        self.psfDeterminer["pca"].spatialOrder = 2
        self.psfDeterminer["pca"].kernelSizeMin = 25

# Note that this class does not inherit from CalibrateConfig; that has lots of things we don't want,
# and duck-typing means we don't need the inheritance for its own sake.
class SdssCalibrateConfig(pexConfig.Config):
    initialPsf = pexConfig.ConfigField(dtype=InitialPsfConfig, doc=InitialPsfConfig.__doc__)
    doBackground = pexConfig.Field(
        dtype = bool,
        doc = ("Estimate and subtract background? (note that this only affects the initial background,"
               " not post-detection re-estimation and subtraction)"),
        default = True,
    )
    doComputeApCorr = pexConfig.Field(
        dtype = bool,
        doc = "Compute aperture correction?",
        default = True,
    )    
    doPhotoCal = pexConfig.Field(
        dtype = bool,
        doc = "Compute photometric zeropoint?",
        default = True,
    )
    repair = pexConfig.ConfigurableField(target = RepairTask, doc = "")
    background = pexConfig.ConfigField(
        dtype = measAlg.estimateBackground.ConfigClass,
        doc = "Background estimation configuration"
        )
    detection = pexConfig.ConfigurableField(
        target = measAlg.SourceDetectionTask,
        doc = "Initial (high-threshold) detection phase for calibration",
    )
    initialMeasurement = pexConfig.ConfigurableField(
        target = measAlg.SourceMeasurementTask,
        doc = "Initial measurements used to feed PSF determination and astrometry",
    )
    measurement = pexConfig.ConfigurableField(
        target = measAlg.SourceMeasurementTask,
        doc = "Post-PSF-determination measurements used to feed other calibrations",
    )
    computeApCorr = pexConfig.ConfigField(dtype = measAlg.ApertureCorrectionConfig,
                                          doc = measAlg.ApertureCorrectionConfig.__doc__)
    astrometry    = pexConfig.ConfigurableField(target = AstrometryTask, doc = "")
    photocal      = pexConfig.ConfigurableField(target = PhotoCalTask, doc="")

    u = pexConfig.ConfigField(dtype=SdssCalibratePerFilterConfig, doc="u-band specific config fields")
    g = pexConfig.ConfigField(dtype=SdssCalibratePerFilterConfig, doc="g-band specific config fields")
    r = pexConfig.ConfigField(dtype=SdssCalibratePerFilterConfig, doc="r-band specific config fields")
    i = pexConfig.ConfigField(dtype=SdssCalibratePerFilterConfig, doc="i-band specific config fields")
    z = pexConfig.ConfigField(dtype=SdssCalibratePerFilterConfig, doc="z-band specific config fields")
    
    def validate(self):
        pexConfig.Config.validate(self)
        if self.initialMeasurement.prefix == self.measurement.prefix:
            raise ValueError("SdssCalibrateConfig.initialMeasurement and SdssCalibrateConfig.measurement "
                             "have the same prefix; field names may clash.")
        if self.initialMeasurement.doApplyApCorr:
            raise ValueError("Cannot apply aperture corrections to pre-PSF initial measurements.")
        if self.measurement.doApplyApCorr and not self.doComputeApCorr:
            raise ValueError("Cannot apply aperture correction without computing it")

    def setDefaults(self):
        self.detection.includeThresholdMultiplier = 10.0
        self.initialMeasurement.prefix = "initial."
        self.initialMeasurement.doApplyApCorr = False
        self.initialMeasurement.algorithms.names = ["flags.pixel", "shape.sdss", "flux.psf", "flux.sinc"]
        self.initialMeasurement.slots.apFlux = "flux.sinc"
        self.initialMeasurement.slots.modelFlux = None
        self.initialMeasurement.slots.instFlux = None
        self.background.binSize = 512
        self.detection.background.binSize = 512
        self.computeApCorr.alg1.name = "flux.psf"
        self.computeApCorr.alg2.name = "flux.sinc"
        

class SdssCalibrateTask(CalibrateTask):
    """SDSS-specific version of lsst.pipe.tasks.calibrate.CalibrateTask
    """
    ConfigClass = SdssCalibrateConfig

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.algMetadata = dafBase.PropertyList()
        self.makeSubtask("repair")
        self.makeSubtask("detection", schema=self.schema)
        self.makeSubtask("initialMeasurement", schema=self.schema, algMetadata=self.algMetadata)
        self.makeSubtask("astrometry", schema=self.schema)
        self.starSelectors = {}
        self.psfDeterminers = {}
        for filterName in ("u", "g", "r", "i", "z"):
            # We don't pass a schema to the star selectors and PSF determiners (it's optional) because we
            # don't currently have a way to make them all share the same flag field.
            subConfig = getattr(self.config, filterName)
            self.starSelectors[filterName] = subConfig.starSelector.apply()
            if subConfig.psfDeterminer.active is None:
                self.psfDeterminers[filterName] = None
            else:
                self.psfDeterminers[filterName] = subConfig.psfDeterminer.apply()
        self.makeSubtask("measurement", schema=self.schema, algMetadata=self.algMetadata)
        self.makeSubtask("photocal", schema=self.schema)

    @pipeBase.timeMethod
    def run(self, exposure, defects=None, idFactory=None):
        """Calibrate an exposure: measure PSF, subtract background, measure astrometry and photometry

        @param[in,out]  exposure   Exposure to calibrate; measured PSF will be installed there as well
        @param[in]      defects    List of defects on exposure
        @param[in]      idFactory  afw.table.IdFactory to use for source catalog.
        @return a pipeBase.Struct with fields:
        - psf: Point spread function
        - apCorr: Aperture correction
        - sources: Sources used in calibration
        - matches: Astrometric matches
        - matchMeta: Metadata for astrometric matches
        """
        assert exposure is not None, "No exposure provided"

        inputPsf = exposure.getPsf()

        if idFactory is None:
            idFactory = afwTable.IdFactory.makeSimple()

        filterName = exposure.getFilter().getName()

        filterConfig = getattr(self.config, filterName) 

        if filterConfig.useInputPsf and inputPsf is None:
            raise ValueError("Cannot run with useInputPsf=True when there is no input PSF.")

        if not filterConfig.useInputPsf:
            self.installInitialPsf(exposure)
            self.log.log(self.log.INFO, "Ignoring input (psField or matched-to) PSF.")
            keepCRs = True   # we'll remove them when we have a better PSF
        else:
            self.log.log(self.log.INFO, "Running with input (psField or matched-to) PSF.")
            keepCRs = None  # this is the last repair we need to run; defer to config values

        self.repair.run(exposure, defects=defects, keepCRs=keepCRs)
        self.display('repair', exposure=exposure)

        if self.config.doBackground:
            with self.timer("background"):
                bg, exposure = measAlg.estimateBackground(exposure, self.config.background, subtract=True)
                del bg
            self.display('background', exposure=exposure)

        table = afwTable.SourceTable.make(self.schema, idFactory)
        table.setMetadata(self.algMetadata)
        detRet = self.detection.makeSourceCatalog(table, exposure)
        sources = detRet.sources

        # If we're using the input PSF, we only need to do one measurement step, and we do that now.
        # If not, we do the initial measurement with the fake PSF in a prefixed part of the schema.

        if not filterConfig.useInputPsf:
            self.initialMeasurement.measure(exposure, sources)
        else:
            self.measurement.measure(exposure, sources)

        # We always run astrometry; if you want to effectively turn it off, set
        # "forceKnownWcs=True" and "calibrateSip=False" in config.astrometry.
        # In that case, we'll still match to the reference catalog, but we won't update the Wcs.
        astromRet = self.astrometry.run(exposure, sources)
        matches = astromRet.matches
        matchMeta = astromRet.matchMeta

        # We always run star selection, but you can set 'starSelector.name = "catalog"'
        # to use stars from the astrometry.net catalog.
        psfCandidateList = self.starSelectors[filterName].selectStars(exposure, sources)
        self.log.log(self.log.INFO, "'%s' PSF star selector found %d candidates" 
                     % (filterConfig.starSelector.name, len(psfCandidateList)))

        # We always run PSF determination
        psf, cellSet = self.psfDeterminers[filterName].determinePsf(exposure, psfCandidateList, self.metadata)
        self.log.log(self.log.INFO, "PSF determination using %d/%d stars." % 
                     (self.metadata.get("numGoodStars"), self.metadata.get("numAvailStars")))
        if  filterConfig.useInputPsf:
            psf = inputPsf
        exposure.setPsf(psf)
        
        # If we aren't using the input PSF, we need to re-repair and re-measure before doing
        # aperture corrections and photometric calibration.
        if not filterConfig.useInputPsf:
            self.repair.run(exposure, defects=defects, keepCRs=None)
            self.display('repair', exposure=exposure)
            self.measurement.measure(exposure, sources)   # don't use run, because we don't have apCorr yet

        if self.config.doComputeApCorr:
            apCorr = self.computeApCorr(exposure, cellSet)
        else:
            apCorr = None

        if self.measurement.config.doApplyApCorr:
            assert(apCorr is not None)
            self.measurement.applyApCorr(sources, apCorr)

        if self.config.doPhotoCal:
            photocalRet = self.photocal.run(matches, exposure.getFilter().getName())
            zp = photocalRet.photocal
            self.log.log(self.log.INFO, "Photometric zero-point: %f" % zp.getMag(1.0))
            exposure.getCalib().setFluxMag0(zp.getFlux(0))

        self.display('calibrate', exposure=exposure, sources=sources, matches=matches)

        return pipeBase.Struct(
            exposure = exposure,
            psf = psf,
            apCorr = apCorr,
            sources = sources,
            matches = matches,
            matchMeta = matchMeta,
        )
