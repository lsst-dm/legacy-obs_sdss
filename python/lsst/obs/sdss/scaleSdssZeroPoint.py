#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import numpy
import lsst.pex.config as pexConfig
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
from lsst.pipe.tasks.scaleZeroPoint import ScaleZeroPointTask
from .selectFluxMag0 import SelectSdssFluxMag0Task

__all__ = ["ScaleSdssZeroPointTask"]

class SdssImageScaler(object):
    """Multiplicative image scaler using interpolation over a grid of points.

    This version only interpolates in the X direction; it is designed for SDSS Stripe82 images
    which scan along the X direction.
    """
    def __init__(self, interpStyle, xList, yList, scaleList):
        """Construct an SdssImageScaler

        @warning: scaleErrList is presently not used

        @param[in] interpStyle: interpolation style (see lsst.afw.math.Interpolate for options)
        @param[in] xList: list of X pixel positions
        @param[in] yList: list of Y pixel positions
        @param[in] scaleList: list of multiplicative scales at (x,y)

        @raise RuntimeError if the lists have different lengths
        """
        if len(xList) != len(yList) or len(xList) != len(scaleList):
            raise RuntimeError(
                "len(xList)=%s len(yList)=%s, len(scaleList)=%s but all lists must have the same length" % \
                (len(xList), len(yList), len(scaleList)))

        self.interpStyle = getattr(afwMath.Interpolate, interpStyle)
        self._xList = xList
        self._yList = yList
        self._scaleList = scaleList

    def scaleMaskedImage(self, maskedImage):
        """Apply scale correction to the specified masked image

        @param[in,out] image to scale; scale is applied in place
        """
        scale = self.getInterpImage(maskedImage.getBBox())
        maskedImage *= scale

    def getInterpImage(self, bbox):
        """Return an image interpolated in R.A direction covering supplied bounding box

        @param[in] bbox: integer bounding box for image (afwGeom.Box2I)
        """

        npoints = len(self._xList)
        #sort by X coordinate
        if npoints < 1:
            raise RuntimeError("Cannot create scaling image. Found no fluxMag0s to interpolate")

        x, z = zip(*sorted(zip(self._xList, self._scaleList)))

        xvec = afwMath.vectorD(x)
        zvec = afwMath.vectorD(z)
        height = bbox.getHeight()
        width = bbox.getWidth()
        x0, y0 = bbox.getMin()

        interp = afwMath.makeInterpolate(xvec, zvec, self.interpStyle)
        interpValArr = numpy.zeros(width, dtype=numpy.float32)

        for i, xInd in enumerate(range(x0, x0 + width)):
            xPos = afwImage.indexToPosition(xInd)
            interpValArr[i] = interp.interpolate(xPos)

        # assume the maskedImage being scaled is MaskedImageF (which is usually true); see ticket #3070
        interpGrid = numpy.meshgrid(interpValArr, range(0, height))[0].astype(numpy.float32)
        image = afwImage.makeImageFromArray(interpGrid)
        image.setXY0(x0, y0)
        return image

class ScaleSdssZeroPointConfig(ScaleZeroPointTask.ConfigClass):
    """Config for ScaleSdssZeroPointTask
    """
    selectFluxMag0 = pexConfig.ConfigurableField(
        doc = "Task to select data to compute spatially varying photometric zeropoint",
        target = SelectSdssFluxMag0Task,
    )
    interpStyle = pexConfig.ChoiceField(
        dtype = str,
        doc = "Algorithm to interpolate the flux scalings;" \
              "Maps to an enum; see afw.math.Interpolate",
        default = "NATURAL_SPLINE",
        allowed={
             "CONSTANT" : "Use a single constant value",
             "LINEAR" : "Use linear interpolation",
             "CUBIC_SPLINE": "cubic spline",
             "NATURAL_SPLINE" : "cubic spline with zero second derivative at endpoints",
             "AKIMA_SPLINE": "higher-level nonlinear spline that is more robust to outliers",
             }
    )
    bufferWidth = pexConfig.Field(
        dtype = float,
        doc = "Buffer in the R.A. direction added to the region to be searched by selectFluxMag0" \
        "Units are multiples of SDSS field widths (1489pix). (e.g. if the exposure is 1000x1000pixels, " \
        "a bufferWidth of 2 results in a search region of 6956 x 1000, centered on the original position.",
        default = 3,
    )

class ScaleSdssZeroPointTask(ScaleZeroPointTask):
    """Compute spatially varying scale factor to scale exposures to a desired photometric zero point
    """
    ConfigClass = ScaleSdssZeroPointConfig
    _DefaultName = "scaleSdssZeroPoint"

    def __init__(self, *args, **kwargs):
        ScaleZeroPointTask.__init__(self, *args, **kwargs)
        self.makeSubtask("selectFluxMag0")
        self.FIELD_WIDTH = 1489.

    def computeImageScaler(self, exposure, dataRef):
        """Query a database for fluxMag0s and return a SdssImageScaler

        @param[in] exposure: exposure for which we want an image scaler
        @param[in] dataRef: dataRef of exposure

        First, triple the width (R.A. direction) of the patch bounding box. Query the database for
        overlapping fluxMag0s corresponding to the same run and filter.
        """
        wcs = exposure.getWcs()
        bbox = exposure.getBBox()
        buffer = int(self.config.bufferWidth * self.FIELD_WIDTH)
        biggerBbox = afwGeom.Box2I(afwGeom.Point2I(bbox.getBeginX()-buffer, bbox.getBeginY()),
                                   afwGeom.Extent2I(bbox.getWidth()+ buffer + buffer, bbox.getHeight()))
        cornerPosList = afwGeom.Box2D(biggerBbox).getCorners()
        coordList = [wcs.pixelToSky(pos) for pos in cornerPosList]
        fluxMagInfoList = self.selectFluxMag0.run(dataRef.dataId, coordList).fluxMagInfoList

        xList = []
        yList = []
        scaleList = []

        for fluxMagInfo in fluxMagInfoList:
            #find center of field in tract coordinates
            if not fluxMagInfo.coordList:
                raise RuntimeError("no x,y data for fluxMagInfo")
            ctr = afwGeom.Extent2D()
            for coord in fluxMagInfo.coordList:
                #accumulate x, y
                ctr += afwGeom.Extent2D(wcs.skyToPixel(coord))
            #and find average x, y as the center of the chip
            ctr = afwGeom.Point2D(ctr / len(fluxMagInfo.coordList))
            xList.append(ctr.getX())
            yList.append(ctr.getY())
            scaleList.append(self.scaleFromFluxMag0(fluxMagInfo.fluxMag0).scale)

        self.log.info("Found %d flux scales for interpolation: %s" % (len(scaleList),
                                                                      ["%0.4f"%(s) for s in scaleList]))
        return SdssImageScaler(
            interpStyle = self.config.interpStyle,
            xList = xList,
            yList = yList,
            scaleList = scaleList,

        )
