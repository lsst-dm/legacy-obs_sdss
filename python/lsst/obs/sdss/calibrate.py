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
from lsst.meas.astrom.catalogStarSelector import CatalogStarSelector
import lsst.afw.table as afwTable
import lsst.afw.math as afwMath
from lsst.meas.base.sfm import SingleFrameMeasurementTask
from lsst.pipe.tasks.calibrate import InitialPsfConfig, CalibrateTask

from lsst.meas.photocal import PhotoCalTask
from lsst.pipe.tasks.astrometry import AstrometryTask
from lsst.pipe.tasks.repair import RepairTask

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
    doPsf = pexConfig.Field(
        dtype = bool,
        doc = "Perform PSF fitting?",
        default = True,
    )
    useExposurePsf = pexConfig.Field(
        dtype=bool,
        default=False,
        doc=("Specify the source of initial PSF."
            "If True, use the PSF in the input exposure (and ignore initialPsf)."
            "If False, use the PSF defined by initialPsf (and ignore the PSF in the exposure)."
            "See also doPsf, which controls whether a better PSF is determined.")
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
        target = SingleFrameMeasurementTask,
        doc = "Initial measurements used to feed PSF determination and astrometry",
    )
    measurement = pexConfig.ConfigurableField(
        target = SingleFrameMeasurementTask,
        doc = "Post-PSF-determination measurements used to feed other calibrations",
    )
    astrometry    = pexConfig.ConfigurableField(target = AstrometryTask, doc = "")
    photocal      = pexConfig.ConfigurableField(target = PhotoCalTask, doc="")

    u = pexConfig.ConfigField(dtype=SdssCalibratePerFilterConfig, doc="u-band specific config fields")
    g = pexConfig.ConfigField(dtype=SdssCalibratePerFilterConfig, doc="g-band specific config fields")
    r = pexConfig.ConfigField(dtype=SdssCalibratePerFilterConfig, doc="r-band specific config fields")
    i = pexConfig.ConfigField(dtype=SdssCalibratePerFilterConfig, doc="i-band specific config fields")
    z = pexConfig.ConfigField(dtype=SdssCalibratePerFilterConfig, doc="z-band specific config fields")

    useExposurePsf = pexConfig.Field(
        dtype=bool,
        optional = True,
        doc=("Specify the source of initial PSF."
            "If True, use the PSF in the input exposure (and ignore initialPsf)."
            "If False, use the PSF defined by initialPsf (and ignore the PSF in the exposure)."
            "See also doPsf, which controls whether a better PSF is determined.")
    )
    doPsf = pexConfig.Field(
        dtype=bool,
        optional = True,
        doc=("Do PSF fitting?")
    )
    minPsfCandidates = pexConfig.Field(
        dtype=int, default=1,
        doc=("If the number of candidates returned by the star selector is less than this amount, "
             "retry with the catalog star selector")
    )
    
    def validate(self):
        pexConfig.Config.validate(self)
        if self.initialMeasurement.target.tableVersion != self.measurement.target.tableVersion:
            raise ValueError("tableVersions of measurement subtasks do not match")

    def setDefaults(self):
        self.detection.includeThresholdMultiplier = 10.0
        #self.initialMeasurement.algorithms.names = ["flags.pixel", "shape.sdss", "flux.psf", "flux.sinc"]
        #self.initialMeasurement.slots.apFlux = "flux.sinc"
        self.initialMeasurement.algorithms.names = ["base_SdssCentroid", "base_PixelFlags", "base_SdssShape", "base_PsfFlux", "base_SincFlux"]
        self.initialMeasurement.slots.centroid = "base_SdssCentroid"
        self.initialMeasurement.slots.apFlux = "base_SincFlux"
        self.initialMeasurement.slots.modelFlux = None
        self.initialMeasurement.slots.instFlux = None
        self.repair.doInterpolate = False
        self.repair.doCosmicRay = False        

class SdssCalibrateTask(CalibrateTask):
    """SDSS-specific version of lsst.pipe.tasks.calibrate.CalibrateTask
    """
    ConfigClass = SdssCalibrateConfig

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.schema1 = afwTable.SourceTable.makeMinimalSchema()
        minimalCount = self.schema1.getFieldCount()
        self.tableVersion = self.config.measurement.target.tableVersion
        self.schema1.setVersion(self.tableVersion)
        self.algMetadata = dafBase.PropertyList()
        self.makeSubtask("repair")
        self.makeSubtask("detection", schema=self.schema1)

        beginInitial = self.schema1.getFieldCount()
        self.makeSubtask("initialMeasurement", schema=self.schema1, algMetadata=self.algMetadata)
        endInitial = self.schema1.getFieldCount()

        self.makeSubtask("astrometry", schema=self.schema1)
        self.starSelectors = {}
        self.psfDeterminers = {}
        if self.tableVersion == 0:
            self.psfCandidateKey = self.schema1.addField(
                "calib.psf.candidate", type="Flag", 
                doc="Set if the source was selected by the star selector algorithm"
            )
            self.psfUsedKey = self.schema1.addField(
                "calib.psf.used", type="Flag",
                doc="Set if the source was used in PSF determination"
            ) 
        else:                                         
            self.psfCandidateKey = self.schema1.addField(
                "calib_psf_candidate", type="Flag", 
                doc="Set if the source was selected by the star selector algorithm"
            )
            self.psfUsedKey = self.schema1.addField(
                "calib_psf_used", type="Flag",
                doc="Set if the source was used in PSF determination"
            ) 

        for filterName in ("u", "g", "r", "i", "z"):
            # We don't pass a schema to the star selectors and PSF determiners (it's optional) because we
            # don't currently have a way to make them all share the same flag field.
            subConfig = getattr(self.config, filterName)
            self.starSelectors[filterName] = subConfig.starSelector.apply()
            self.psfDeterminers[filterName] = subConfig.psfDeterminer.apply()

        self.makeSubtask("photocal", schema=self.schema1)
        # create a schemaMapper to map schema1 into schema2
        self.schemaMapper = afwTable.SchemaMapper(self.schema1)
        if self.tableVersion == 0:
            separator = "."
        else:
            separator =  "_"
        count = 0
        for item in self.schema1:
            count = count + 1
            field = item.getField()
            name = field.getName()
            if count > beginInitial and count <= endInitial:
                name = "initial" + separator + name
            self.schemaMapper.addMapping(item.key, name)

        # measurements fo the second measurement step done with a second schema
        schema = self.schemaMapper.editOutputSchema()
        self.makeSubtask("measurement", schema=schema, algMetadata=self.algMetadata)
        # the final schema is the same as the schemaMapper output
        self.schema = self.schemaMapper.getOutputSchema()

    def getCalibKeys(self):
        return (self.psfCandidateKey, self.psfUsedKey)

    @pipeBase.timeMethod
    def run(self, exposure, defects=None, idFactory=None):
        """Calibrate an exposure: measure PSF, subtract background, measure astrometry and photometry

        @param[in,out]  exposure   Exposure to calibrate
        @param[in]      defects    List of defects on exposure
        @param[in]      idFactory  afw.table.IdFactory to use for source catalog.
        @return a pipeBase.Struct with fields:
        - backgrounds: Array of background objects that were subtracted from the exposure
        - psf: Point spread function
        - sources: Sources used in calibration
        - matches: Astrometric matches
        - matchMeta: Metadata for astrometric matches
        """

        psf = None
        matches = None
        matchMeta = None
        cellSet = None
        backgrounds = []

        if idFactory is None:
            idFactory = afwTable.IdFactory.makeSimple()

        filterName = exposure.getFilter().getName()
        filterConfig = getattr(self.config, filterName)

        if self.config.useExposurePsf is None:
            useExposurePsf = filterConfig.useExposurePsf
        else:
            useExposurePsf = self.config.useExposurePsf

        if useExposurePsf:
            if exposure.getPsf() is None:
                raise ValueError("Cannot run with useExposurePsf=True when there is no input PSF.")        
            self.log.log(self.log.INFO, "Initial PSF is from input exposure; ignoring config.initialPsf.")
        else:
            self.installInitialPsf(exposure)
            self.log.log(self.log.INFO, "Initial PSF is from config.initialPsf; ignoring PSF in exposure.")

        initialPsf = exposure.getPsf()
        if initialPsf is None:
            raise RuntimeError("Internal error: exposure has no PSF when it should")
        
        if self.config.doPsf:
            keepCRs = True   # we'll remove them when we have a better PSF
        else:
            keepCRs = None  # this is the last repair we need to run; defer to config values

        self.repair.run(exposure, defects=defects, keepCRs=keepCRs)
        self.display('repair', exposure=exposure)

        if self.config.doBackground:
            with self.timer("background"):
                bg, exposure = measAlg.estimateBackground(exposure, self.config.background, subtract=True)
                backgrounds.append(bg)
            self.display('background', exposure=exposure)

        table1 = afwTable.SourceTable.make(self.schema1, idFactory)
        table1.setMetadata(self.algMetadata)
        detRet = self.detection.makeSourceCatalog(table1, exposure)
        sources1 = detRet.sources
        if detRet.fpSets.background:
            backgrounds.append(detRet.fpSets.background)

        # If we're using the input PSF, we only need to do one measurement step, and we do that now.
        # If not, we do the initial measurement with the fake PSF in the intial schema1

        if self.config.doPsf:
            self.initialMeasurement.measure(exposure, sources1)
     
        # make a second table with which to do the second measurement
        # the schemaMapper will copy the footprints and ids, which is all we need.
        # Note that the old measurements fields are copied to "initial_" names
        table2 = afwTable.SourceTable.make(self.schema, idFactory)
        table2.setMetadata(self.algMetadata)
        sources = afwTable.SourceCatalog(table2)
        # transfer to a second table
        sources.extend(sources1, self.schemaMapper)
        if self.measurement.tableVersion == 0:
            separator = "."
        else:
            separator = "_"
        if sources1.hasCentroidSlot():
            sources.defineCentroid("initial" + separator + sources1.getCentroidDefinition())
        if sources1.hasShapeSlot():
            sources.defineShape("initial" + separator + sources1.getShapeDefinition())
        if sources1.hasPsfFluxSlot():
            sources.definePsfFlux("initial" + separator + sources1.getPsfFluxDefinition())
        if sources1.hasInstFluxSlot():
            sources.defineInstFlux("initial" + separator + sources1.getInstFluxDefinition())
        if sources1.hasModelFluxSlot():
            sources.defineModelFlux("initial" + separator + sources1.getModelFluxDefinition())
        if sources1.hasApFluxSlot():
            sources.defineApFlux("initial" + separator + sources1.getApFluxDefinition())
     
        if not self.config.doPsf:
            self.measurement.measure(exposure, sources)

        # We always run astrometry; if you want to effectively turn it off, set
        # "forceKnownWcs=True" and "calibrateSip=False" in config.astrometry.
        # In that case, we'll still match to the reference catalog, but we won't update the Wcs.
        astromRet = self.astrometry.run(exposure, sources)
        matches = astromRet.matches
        matchMeta = astromRet.matchMeta

        # We always run star selection, but you can set 'starSelector.name = "catalog"'
        # to use stars from the astrometry.net catalog.
        # If we fail, we fall back to a default-constructed catalog star selector.
        psfCandidateList = self.starSelectors[filterName].selectStars(exposure, sources, matches=matches)
        if (len(psfCandidateList) < self.config.minPsfCandidates
            and filterConfig.starSelector.name != "catalog"):
            self.log.warn("'%s' PSF star selector found %d < %d candidates; trying catalog star selector" 
                          % (filterConfig.starSelector.name, len(psfCandidateList),
                             self.config.minPsfCandidates))
            self.metadata.add("StarSelectorStatus", "failed; had to fall back to catalog")
            selector = CatalogStarSelector()
            psfCandidateList = selector.selectStars(exposure, sources, matches=matches)
            self.log.log(self.log.INFO, "'catalog' PSF star selector found %d candidates" 
                         % len(psfCandidateList))
        else:
            self.log.log(self.log.INFO, "'%s' PSF star selector found %d candidates" 
                         % (filterConfig.starSelector.name, len(psfCandidateList)))
        if psfCandidateList and self.psfCandidateKey is not None:
            for cand in psfCandidateList:
                source = cand.getSource()
                source.set(self.psfCandidateKey, True)

        if self.config.doPsf:
            try:
                psf, cellSet = self.psfDeterminers[filterName].determinePsf(exposure, psfCandidateList,
                                                                            self.metadata,
                                                                            flagKey=self.psfUsedKey)
                self.log.log(self.log.INFO, "PSF determination using %d/%d stars." % 
                             (self.metadata.get("numGoodStars"), self.metadata.get("numAvailStars")))
                exposure.setPsf(psf)
            except Exception, err:
                self.log.warn("PSF determination failed; falling back to initial PSF: %s" % err)
                self.metadata.add("PsfDeterminerStatus", "failed; had to fall back to initial PSF")
                psf = initialPsf
        else:
            self.log.log(self.log.INFO, "Not running PSF determination")
            psf = initialPsf
        
        # If we aren't using the input PSF, we need to re-repair and re-measure before doing
        # aperture corrections and photometric calibration.
        if self.config.doPsf:
            self.repair.run(exposure, defects=defects, keepCRs=None)
            self.display('repair', exposure=exposure)
            self.measurement.run(exposure, sources)
            self.log.log(self.log.INFO, "Re-running astrometry after measurement with improved PSF.")
            astromRet = self.astrometry.run(exposure, sources)
            matches = astromRet.matches
            matchMeta = astromRet.matchMeta


        if self.config.doPhotoCal:
            photocalRet = self.photocal.run(exposure, matches)
            self.log.info("Photometric zero-point: %f" % photocalRet.calib.getMagnitude(1.0))
            exposure.getCalib().setFluxMag0(photocalRet.calib.getFluxMag0())
        self.display('calibrate', exposure=exposure, sources=sources, matches=matches)

        return pipeBase.Struct(
            exposure = exposure,
            backgrounds = backgrounds,
            psf = psf,
            sources = sources,
            matches = matches,
            matchMeta = matchMeta,
        )

    def makeCellSet(self, exposure, psfCandidateList):
        """
        Make a spatial cell set that's more-or-less like the one the a psfDeterminer
        would have created, if it had succeeded.

        TODO: this code should be refactored into meas_algorithms, probably as a way
        to construct aperture corrections without an input spatial cell set, or as
        an alternate method on PsfDeterminers that always succeeds or at least succeeds
        more often.
        """
        bbox = exposure.getBBox()
        filterName = exposure.getFilter().getName()
        filterConfig = getattr(self.config, filterName)
        # This is slightly dangerous: we assume PsfDeterminerConfig has sizeCellX/sizeCellY,
        # which may only be true for PcaPsfDeterminer, but that's the only one we're likely
        # to use now, and others will probably have at least those settings.
        psfDeterminerConfig = filterConfig.psfDeterminer.active
        cellSet = afwMath.SpatialCellSet(
            bbox, psfDeterminerConfig.sizeCellX, psfDeterminerConfig.sizeCellY
        )
        for i, psfCandidate in enumerate(psfCandidateList):
            try:
                cellSet.insertCandidate(psfCandidate)
            except Exception, e:
                self.log.logdebug("Skipping PSF candidate %d of %d: %s" % (i, len(psfCandidateList), e))
                continue
        return cellSet
