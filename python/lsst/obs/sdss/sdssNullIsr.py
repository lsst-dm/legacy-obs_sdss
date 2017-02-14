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
from lsst.pipe.tasks.processCcd import ProcessCcdTask


class SdssNullIsrConfig(ProcessCcdTask.ConfigClass):
    """Config for SdssNullIsrTask"""
    removePedestal = pexConfig.Field(
        dtype=bool,
        doc="Remove SDSS pedestal from fpC file?",
        default=True,
    )
    pedestalVal = pexConfig.Field(
        dtype=int,
        doc="Number of counts in the SDSS pedestal",
        default=1000,
    )
    removeOverlap = pexConfig.Field(
        dtype=bool,
        doc="Remove SDSS field overlap from fpC file?",
        default=True,
    )
    overlapSize = pexConfig.Field(
        dtype=int,
        doc="Number of pixels to remove from top of the fpC file",
        default=128,
    )
    doWrite = pexConfig.Field(
        dtype=bool,
        doc="Persist loaded data as a postISRCCD? The default is false, to avoid duplicating data.",
        default=False,
    )
    datasetType = pexConfig.Field(
        dtype=str,
        doc="Dataset type for input data; read by ProcessCcdTask; users will typically leave this alone",
        default="fpC",
    )


# \addtogroup LSST_task_documentation
# \{
# \page SdssNullIsrTask
# \ref SdssNullIsrTask_ "SdssNullIsrTask"
# \copybrief SdssNullIsrTask
# \}

class SdssNullIsrTask(pipeBase.Task):
    """!Load SDSS post-ISR data from fpC, fpM, asTrans, etc. files

    @anchor SdssNullIsrTask_

    @section pipe_tasks_sdssNullIsr_Contents  Contents

     - @ref pipe_tasks_sdssNullIsr_Purpose
     - @ref pipe_tasks_sdssNullIsr_Initialize
     - @ref pipe_tasks_sdssNullIsr_IO
     - @ref pipe_tasks_sdssNullIsr_Config

    @section pipe_tasks_sdssNullIsr_Purpose  Description

    Load "instcal" exposures from the community pipeline as a post-ISR exposure,
    and optionally persists it as a `postISRCCD`.

    This is used to retarget the `isr` subtask in `ProcessCcdTask` for SDSS
    because the LSST software stack is not capable of performing ISR on SDSS data.

    @section pipe_tasks_sdssNullIsr_Initialize  Task initialisation

    @copydoc \_\_init\_\_

    @section pipe_tasks_sdssNullIsr_IO  Invoking the Task

    The main method is `runDataRef`.

    @section pipe_tasks_sdssNullIsr_Config  Configuration parameters

    See @ref SdssNullIsrConfig
    """
    ConfigClass = SdssNullIsrConfig
    _DefaultName = "isr"

    @pipeBase.timeMethod
    def loadExposure(self, sensorRef):
        """Load SDSS data as a post-ISR exposure

        - Image is from fpC
        - Mask is from fpM
        - Wcs is from asTrans
        - Calib is from tsField
        - Psf is from psField
        """
        image = sensorRef.get("fpC").convertF()
        if self.config.removePedestal:
            image -= self.config.pedestalVal
        mask = sensorRef.get("fpM")
        wcs = sensorRef.get("asTrans")
        tsField = sensorRef.get("tsField")
        calib = tsField.calib
        gain = tsField.gain
        var = afwImage.ImageF(image, True)
        var /= gain

        mi = afwImage.MaskedImageF(image, mask, var)

        if self.config.removeOverlap:
            bbox = mi.getBBox()
            begin = bbox.getBegin()
            extent = bbox.getDimensions()
            extent -= afwGeom.Extent2I(0, self.config.overlapSize)
            tbbox = afwGeom.BoxI(begin, extent)
            mi = afwImage.MaskedImageF(mi, tbbox)

        exposure = afwImage.ExposureF(mi, wcs)
        expInfo = exposure.getInfo()
        expInfo.setCalib(calib)

        camera = sensorRef.get('camera')
        detector = camera["%(filter)s%(camcol)d" % sensorRef.dataId]
        expInfo.setDetector(detector)
        expInfo.setFilter(afwImage.Filter(sensorRef.dataId['filter']))

        visitInfo = afwImage.makeVisitInfo(
            exposureTime=tsField.exptime,
            date=tsField.dateAvg,
            boresightAirmass=tsField.airmass,
        )
        expInfo.setVisitInfo(visitInfo)

        # Install the SDSS PSF here; if we want to overwrite it later, we can.
        psf = sensorRef.get('psField')
        exposure.setPsf(psf)

        return exposure

    @pipeBase.timeMethod
    def runDataRef(self, sensorRef):
        """!Load SDSS data as post-ISR exposure and possibly persist it as a post-ISR CCD

        @param[in] sensorRef  butler data reference for post-ISR exposure
            (a daf.persistence.butlerSubset.ButlerDataRef)

        @return a pipeBase.Struct with fields:
        - exposure: the exposure after application of ISR
        """
        self.log.info("Loading SDSS asTrans file %s" % (sensorRef.dataId))

        exposure = self.loadExposure(sensorRef)

        if self.config.doWrite:
            sensorRef.put(exposure, "postISRCCD")

        return pipeBase.Struct(
            exposure=exposure,
        )
