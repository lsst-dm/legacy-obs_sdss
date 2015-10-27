#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
from .makeCamera import makeCamera
from lsst.pipe.tasks.processCcd import ProcessCcdTask

class ProcessCcdSdssConfig(ProcessCcdTask.ConfigClass):
    """Config for ProcessCcdSdss"""
    removePedestal = pexConfig.Field(dtype=bool, default=True, doc="Remove SDSS pedestal from fpC file")
    pedestalVal = pexConfig.Field(dtype=int, default=1000, doc="Number of counts in the SDSS pedestal")

    removeOverlap = pexConfig.Field(dtype=bool, default=True,
                                    doc="Remove SDSS field overlap from fpC file")
    overlapSize = pexConfig.Field(dtype=int, default=128,
                                  doc="Number of pixels to remove from top of the fpC file")
    loadSdssWcs = pexConfig.Field(dtype=bool, default=False,
                                  doc = ("Load WCS from asTrans; it can then be used as-is or updated by "
                                         "our own code, dependening on calibrate.astrometry parameters.")
    )

    def setDefaults(self):
        ProcessCcdTask.ConfigClass.setDefaults(self)
        self.doIsr = False

class ProcessCcdSdssTask(ProcessCcdTask):
    """Process a CCD for SDSS
    """
    ConfigClass = ProcessCcdSdssConfig
    _DefaultName = "processCcd"
    dataPrefix = ""

    def __init__(self, **kwargs):
        ProcessCcdTask.__init__(self, **kwargs)
        self.camera = makeCamera()

    @classmethod
    def _makeArgumentParser(cls):
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", datasetType="fpC", help="data ID, e.g. --id run=1 camcol=2 field=345")
        return parser

    @pipeBase.timeMethod
    def makeExp(self, sensorRef):
        image = sensorRef.get("fpC").convertF()
        if self.config.removePedestal:
            image -= self.config.pedestalVal
        mask = sensorRef.get("fpM")
        wcs = sensorRef.get("asTrans")
        calib, gain = sensorRef.get("tsField")
        var = afwImage.ImageF(image, True)
        var /= gain

        mi = afwImage.MaskedImageF(image, mask, var)

        if self.config.removeOverlap:
            bbox = mi.getBBox(afwImage.LOCAL)
            begin = bbox.getBegin()
            extent = bbox.getDimensions()
            extent -= afwGeom.Extent2I(0, self.config.overlapSize)
            tbbox = afwGeom.BoxI(begin, extent)
            mi = afwImage.MaskedImageF(mi, tbbox, True)

        exp = afwImage.ExposureF(mi, wcs)
        exp.setCalib(calib)

        det = sensorRef.get('camera')["%s%d"%(sensorRef.dataId["filter"], sensorRef.dataId["camcol"])]
        exp.setDetector(det)
        exp.setFilter(afwImage.Filter(sensorRef.dataId['filter']))

        # Install the SDSS PSF here; if we want to overwrite it later, we can.
        psf = sensorRef.get('psField')
        exp.setPsf(psf)

        return exp

    def setPostIsrExposure(self, sensorRef):
        """Load the post instrument signature removal image

        \param[in]  sensorRef        sensor-level butler data reference

        \return     postIsrExposure  exposure to be passed to processCcdExposure
        """
        if self.config.doCalibrate:
            exp = self.makeExp(sensorRef)
            if self.config.loadSdssWcs:
                self.log.info("Loading WCS from asTrans")
                wcs = sensorRef.get("asTrans")
                exp.setWcs(wcs)
        else:
            exp = None
        return exp

    @pipeBase.timeMethod
    def run(self, sensorRef):
        """Process a CCD: including source detection, photometry and WCS determination

        @param sensorRef: sensor-level butler data reference to SDSS fpC file
        @return pipe_base Struct containing these fields:
        - exposure: calibrated exposure (calexp): as computed if config.doCalibrate,
            else as upersisted and updated if config.doDetection, else None
        - calib: object returned by calibration process if config.doCalibrate, else None
        - sources: detected source if config.doPhotometry, else None
        """
        # delegate most of the work to ProcessCcdTask (which, in turn, delegates to ProcessImageTask)
        result = ProcessCcdTask.run(self, sensorRef)
        return result
