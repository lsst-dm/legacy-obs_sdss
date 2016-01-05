#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
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
import math

import lsst.daf.base as dafBase
import lsst.pipe.base as pipeBase
import lsst.pex.config as pexConfig
import lsst.meas.algorithms as measAlg
from lsst.meas.astrom.catalogStarSelector import CatalogStarSelector
import lsst.afw.table as afwTable
import lsst.afw.math as afwMath
from lsstDebug import getDebugFrame
from lsst.afw.display import getDisplay
from lsst.meas.base import BasePlugin, SingleFrameMeasurementTask, MeasureApCorrTask
from lsst.meas.astrom import AstrometryTask, displayAstrometry
from lsst.pipe.tasks.photoCal import PhotoCalTask
from lsst.pipe.tasks.calibrate import InitialPsfConfig, CalibrateTask
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
    doPsf = pexConfig.Field(
        dtype = bool,
        doc = "Perform PSF fitting?",
        default = True,
    )
    doMeasureApCorr = pexConfig.Field(
        dtype = bool,
        doc = "Compute aperture corrections?",
        default = True,
    )
    doPhotoCal = pexConfig.Field(
        dtype = bool,
        doc = "Compute photometric zeropoint?",
        default = True,
    )
    background = pexConfig.ConfigField(
        dtype = measAlg.estimateBackground.ConfigClass,
        doc = "Background estimation configuration"
        )
    repair = pexConfig.ConfigurableField(
        target = RepairTask,
        doc = "Interpolate over defects and cosmic rays",
    )
    detection = pexConfig.ConfigurableField(
        target = measAlg.SourceDetectionTask,
        doc = "Initial (high-threshold) detection phase for calibration",
    )
    initialMeasurement = pexConfig.ConfigurableField(
        target = SingleFrameMeasurementTask,
        doc = "Initial measurements used to feed PSF determination and aperture correction determination",
    )
    measurement = pexConfig.ConfigurableField(
        target = SingleFrameMeasurementTask,
        doc = "Post-PSF-determination measurements used to feed other calibrations",
    )
    measureApCorr = pexConfig.ConfigurableField(
        target = MeasureApCorrTask,
        doc = "subtask to measure aperture corrections"
    )
    astrometry = pexConfig.ConfigurableField(
        target = AstrometryTask,
        doc = "fit WCS of exposure",
    )
    photocal = pexConfig.ConfigurableField(
        target = PhotoCalTask,
        doc = "peform photometric calibration",
    )

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
    minPsfCandidates = pexConfig.Field(
        dtype=int, default=1,
        doc=("If the number of candidates returned by the star selector is less than this amount, "
             "retry with the catalog star selector")
    )
    
    def validate(self):
        pexConfig.Config.validate(self)
        if self.measurement.doApplyApCorr.startswith("yes") and not self.doMeasureApCorr:
            raise ValueError("Cannot set measurement.doApplyApCorr to 'yes...'"
                " unless doMeasureApCorr is True")

    def setDefaults(self):
        self.detection.includeThresholdMultiplier = 10.0
        self.initialMeasurement.algorithms.names = ["base_SdssCentroid", "base_PixelFlags",
                                                    "base_SdssShape", "base_PsfFlux",
                                                    "base_CircularApertureFlux"]
        self.initialMeasurement.algorithms["base_CircularApertureFlux"].radii = [7.0]
        self.initialMeasurement.slots.centroid = "base_SdssCentroid"
        self.initialMeasurement.slots.apFlux = "base_CircularApertureFlux_7_0"
        self.initialMeasurement.slots.modelFlux = None
        self.initialMeasurement.slots.instFlux = None
        self.initialMeasurement.slots.calibFlux = "base_CircularApertureFlux_7_0"
        self.initialMeasurement.doApplyApCorr = "no" # no aperture correction data yet
        self.measurement.doApplyApCorr = "yes"
        self.repair.doInterpolate = False
        self.repair.doCosmicRay = False
        # we rarely run PSF determination on SDSS data, so use the output of the star selector instead
        self.measureApCorr.inputFilterFlag = "calib_psfCandidate"

class SdssCalibrateTask(CalibrateTask):
    """SDSS-specific version of lsst.pipe.tasks.calibrate.CalibrateTask
    
    Changes from the default task include:
    - Add support for filter-specific configuration of star selection and PSF measurement
    - Always measure astrometry; it is not optional
    """
    ConfigClass = SdssCalibrateConfig

    def __init__(self, **kwargs):
        """!
        Create the calibration task

        \param **kwargs keyword arguments to be passed to lsst.pipe.base.task.Task.__init__
        """
        pipeBase.Task.__init__(self, **kwargs)

        # the calibrate Source Catalog is divided into two catalogs to allow measurement to be run twice
        # schema1 contains everything except what is added by the second measurement task.
        # Before the second measurement task is run, self.schemaMapper transforms the sources into
        # the final output schema, at the same time renaming the measurement fields to "initial_" 
        self.schema1 = afwTable.SourceTable.makeMinimalSchema()
        self.algMetadata = dafBase.PropertyList()
        self.makeSubtask("repair")
        self.makeSubtask("detection", schema=self.schema1)
        beginInitial = self.schema1.getFieldCount()
        self.makeSubtask("initialMeasurement", schema=self.schema1, algMetadata=self.algMetadata)
        endInitial = self.schema1.getFieldCount()

        # create subtasks that are run with schema1 (and possibly also the final schema)
        self.makeSubtask("astrometry")
        self.starSelectors = {}
        self.psfDeterminers = {}
        self.psfCandidateKey = self.schema1.addField(
            "calib_psfCandidate", type="Flag", 
            doc="Set if the source was selected by the star selector algorithm"
        )
        self.psfUsedKey = self.schema1.addField(
            "calib_psfUsed", type="Flag",
            doc="Set if the source was used in PSF determination"
        ) 

        for filterName in ("u", "g", "r", "i", "z"):
            # We don't pass a schema to the star selectors and PSF determiners (it's optional) because we
            # don't currently have a way to make them all share the same flag field.
            subConfig = getattr(self.config, filterName)
            self.starSelectors[filterName] = subConfig.starSelector.apply()
            self.psfDeterminers[filterName] = subConfig.psfDeterminer.apply()

        # create a schemaMapper to map schema1 into the final schema
        self.schemaMapper = afwTable.SchemaMapper(self.schema1)
        separator =  "_"
        count = 0
        for item in self.schema1:
            count = count + 1
            field = item.getField()
            name = field.getName()
            if count > beginInitial and count <= endInitial:
                name = "initial" + separator + name
            self.schemaMapper.addMapping(item.key, name)

        # create subtasks that are run only with the final schema
        schema = self.schemaMapper.editOutputSchema()
        self.makeSubtask("measurement", schema=schema, algMetadata=self.algMetadata)
        self.makeSubtask("measureApCorr", schema=schema)
        self.makeSubtask("photocal", schema=schema)

        # the final schema is the same as the schemaMapper output
        self.schema = self.schemaMapper.getOutputSchema()

    def getCalibKeys(self):
        """!
        Return a sequence of schema keys that represent fields that should be propagated from
        icSrc to src by ProcessCcdTask.
        """
        return (self.psfCandidateKey, self.psfUsedKey)

    @pipeBase.timeMethod
    def run(self, exposure, defects=None, idFactory=None, expId=0):
        """!Run the calibration task on an exposure

        \param[in,out]  exposure   Exposure to calibrate; measured PSF will be installed there as well
        \param[in]      defects    List of defects on exposure
        \param[in]      idFactory  afw.table.IdFactory to use for source catalog.
        \param[in]      expId      Exposure id used for random number generation. Note: Unused as psfs
                                   are loaded from disk.
        \return a pipeBase.Struct with fields:
        - exposure: Repaired exposure
        - backgrounds: A list of background models applied in the calibration phase
        - psf: Point spread function
        - sources: Sources used in calibration
        - matches: Astrometric matches
        - matchMeta: Metadata for astrometric matches
        - photocal: Output of photocal subtask

        It is moderately important to provide a decent initial guess for the seeing if you want to
        deal with cosmic rays.  If there's a PSF in the exposure it'll be used; failing that the
        CalibrateConfig.initialPsf is consulted (although the pixel scale will be taken from the
        WCS if available).

        If the exposure contains an lsst.afw.image.Calib object with the exposure time set, MAGZERO
        will be set in the task metadata.
        """

        psf = None
        matches = None
        matchMeta = None
        cellSet = None
        backgrounds = afwMath.BackgroundList()

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
        frame = getDebugFrame(self._display, "repair")
        if frame:
            getDisplay(frame).mtv(exposure)

        if self.config.doBackground:
            with self.timer("background"):
                bg, exposure = measAlg.estimateBackground(exposure, self.config.background, subtract=True)
                backgrounds.append(bg)
            frame = getDebugFrame(self._display, "background")
            if frame:
                getDisplay(frame).mtv(exposure)

        # Make both tables from the same detRet, since detection can only be run once
        table1 = afwTable.SourceTable.make(self.schema1, idFactory)
        table1.setMetadata(self.algMetadata)
        detRet = self.detection.makeSourceCatalog(table1, exposure)
        sources1 = detRet.sources
        if detRet.fpSets.background:
            backgrounds.append(detRet.fpSets.background)

        # do the initial measurement.  This is normally done for star selection, but do it 
        # even if the psf is not going to be calculated for consistency
        self.initialMeasurement.run(exposure, sources1, allowApCorr=False)
     
        # make a second table with which to do the second measurement
        # the schemaMapper will copy the footprints and ids, which is all we need.
        # Note that the old measurements fields are copied to "initial_" names
        table2 = afwTable.SourceTable.make(self.schema, idFactory)
        table2.setMetadata(self.algMetadata)
        sources = afwTable.SourceCatalog(table2)
        # transfer to a second table -- note that the slots do not have to be reset here
        # as long as measurement.run follows immediately
        sources.extend(sources1, self.schemaMapper)
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

        # Measurement gets run now if doPsf.  Otherwise it gets run conditionally later.
        # astrometry is always run, though it is run twice if doPsf, the last time on 
        # the results of measurement.
        if not self.config.doPsf:
            self.measurement.run(exposure, sources, allowApCorr=False)

        # We always run astrometry; if you want to effectively turn it off, set
        # "forceKnownWcs=True" in config.astrometry.
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
            frame = getDebugFrame(self._display, "repair")
            if frame:
                getDisplay(frame).mtv(exposure)
            self.measurement.run(exposure, sources, allowApCorr=False)
            self.log.log(self.log.INFO, "Re-running astrometry after measurement with improved PSF.")
            astromRet = self.astrometry.run(exposure, sources)
            matches = astromRet.matches
            matchMeta = astromRet.matchMeta

        if self.config.doMeasureApCorr:
            # Run measurement through all flux measurements (all have the same execution order),
            # then apply aperture corrections, then run the rest of the measurements
            self.measurement.run(exposure, sources, endOrder=BasePlugin.APCORR_ORDER)
            apCorrMap = self.measureApCorr.run(bbox=exposure.getBBox(), catalog=sources).apCorrMap
            exposure.getInfo().setApCorrMap(apCorrMap)
            self.measurement.run(exposure, sources, beginOrder=BasePlugin.APCORR_ORDER)
        else:
            self.measurement.run(exposure, sources, allowApCorr=False)

        if self.config.doPhotoCal:
            assert(matches is not None)
            try:
                photocalRet = self.photocal.run(exposure, matches)
            except Exception, e:
                self.log.warn("Failed to determine photometric zero-point: %s" % e)
                photocalRet = None
                self.metadata.set('MAGZERO', float("NaN"))

            if photocalRet:
                self.log.info("Photometric zero-point: %f" % photocalRet.calib.getMagnitude(1.0))
                exposure.getCalib().setFluxMag0(photocalRet.calib.getFluxMag0())
                metadata = exposure.getMetadata()
                # convert to (mag/sec/adu) for metadata
                try:
                    magZero = photocalRet.zp - 2.5 * math.log10(exposure.getCalib().getExptime() )
                    metadata.set('MAGZERO', magZero)
                except:
                    self.log.warn("Could not set normalized MAGZERO in header: no exposure time")
                metadata.set('MAGZERO_RMS', photocalRet.sigma)
                metadata.set('MAGZERO_NOBJ', photocalRet.ngood)
                metadata.set('COLORTERM1', 0.0)
                metadata.set('COLORTERM2', 0.0)
                metadata.set('COLORTERM3', 0.0)
        else:
            photocalRet = None

        frame = getDebugFrame(self._display, "calibrate")
        if frame:
            displayAstrometry(exposure=exposure, sourceCat=sources, matches=matches, frame=frame, pause=False)

        return pipeBase.Struct(
            exposure = exposure,
            backgrounds = backgrounds,
            psf = psf,
            sources = sources,
            matches = matches,
            matchMeta = matchMeta,
            photocal = photocalRet,
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
