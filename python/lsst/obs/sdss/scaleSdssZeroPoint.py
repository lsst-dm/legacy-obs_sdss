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
import MySQLdb
import os

import numpy
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
from lsst.afw.coord import IcrsCoord
import lsst.afw.geom as afwGeom
from lsst.daf.persistence import DbAuth
import lsst.pipe.base as pipeBase
from lsst.pipe.tasks.selectImages import SelectImagesConfig, BaseExposureInfo
from lsst.coadd.utils import ImageScaler, ScaleZeroPointTask
from .selectFluxMag0 import SelectSdssfluxMag0Task

__all__ = ["ScaleSdssZeroPointTask"]

class SdssImageScaler(object):
    def __init__(self):
        """Construct a multiplicative scale factor.
        
        It consists of  a list of points in tract coordinates. Each point has an X, Y, and a scalefactor.
        """
        self.xList = []
        self.yList = []
        self.scaleList = []
        #self.scaleErrList = []

    def scaleMaskedImage(self, maskedImage):
        """Apply scale correction to the specified masked image
        
        @param[in,out] image to scale; scale is applied in place
        """
        scale = self.getInterpImage(maskedImage.getBBox(afwImage.PARENT))
        maskedImage *= scale

    def getInterpImage(self, bbox):
        """Return an image interpolated in R.A direction covering supplied bounding box
        
         Hard-coded to work with obs_sdss only
        """
        
        npoints = len(self.xList)
        xvec = afwMath.vectorF(self.xList)
        zvec = afwMath.vectorF(self.scaleList)      
        height = bbox.getHeight()
        width = bbox.getWidth()
        x0, y0 = bbox.getBegin()

        # I can't get makeInterpolate to work "Message: gsl_interp_init failed: invalid argument supplied by user [4]"
        # interp = afwMath.makeInterpolate(xvec, zvec, afwMath.Interpolate.LINEAR)
        # interp = afwMath.makeInterpolate(xvec, zvec) won't initialize spline. 
        # eval = interp.interpolate(range(y0, height))
        
        eval = numpy.interp(range(x0, x0 + width), xvec, zvec)
        
        evalGrid = numpy.meshgrid(eval.astype(numpy.float32),range(0, height))[0]
        image = afwImage.makeImageFromArray(evalGrid)
        image.setXY0(x0, y0)
        return image


class ScaleSdssZeroPointConfig(ScaleZeroPointTask.ConfigClass):
    selectFluxMag0 = pexConfig.ConfigurableField(
        doc = "Task to select data to compute spatially varying photometric zeropoint",
        target = SelectSdssfluxMag0Task,
    )


class ScaleSdssZeroPointTask(ScaleZeroPointTask):
    """Select SDSS images suitable for coaddition
    """
    ConfigClass = ScaleSdssZeroPointConfig
    _DefaultName = "scaleSdssZeroPoint"
    
    def __init__(self, *args, **kwargs):
        """Construct a ScaleZeroPointTask
        """
        pipeBase.Task.__init__(self, *args, **kwargs)
        self.makeSubtask("selectFluxMag0")
        

        fluxMag0 = 10**(0.4 * self.config.zeroPoint)
        self._calib = afwImage.Calib()
        self._calib.setFluxMag0(fluxMag0)

    def computeImageScaler(self, exposure, exposureId, wcs):
        """
        Query a database for fluxMag0s and return a SdssImageScaler

        First, triple the width (R.A. direction) of the patch bounding box. Query the database for
        overlapping fluxMag0s corresponding to the same run and filter.

        
        """
        imageScaler = SdssImageScaler()
        bbox = exposure.getBBox()
        buffer = 2*bbox.getWidth()
        biggerBbox = afwGeom.Box2I(afwGeom.Point2I(bbox.getBeginX()-buffer, bbox.getBeginY()),
                                   afwGeom.Extent2I(bbox.getWidth()+buffer, bbox.getHeight()))
        cornerPosList = afwGeom.Box2D(biggerBbox).getCorners()
        coordList = [wcs.pixelToSky(pos) for pos in cornerPosList]
        runArgDict = self.selectFluxMag0._runArgDictFromDataId(exposureId)
        
        fluxMagInfoList = self.selectFluxMag0.run(coordList, **runArgDict).fluxMagInfoList

        for fluxMagInfo in fluxMagInfoList:
            self.log.info("found %s, fluxMag0 %s"%(
                fluxMagInfo.dataId, self.scaleFromFluxMag0(fluxMagInfo.fluxMag0).scale))
            raCenter = (fluxMagInfo.coordList[0].getRa() +  fluxMagInfo.coordList[1].getRa() +
                        fluxMagInfo.coordList[2].getRa() +  fluxMagInfo.coordList[3].getRa())/ 4.
            decCenter = (fluxMagInfo.coordList[0].getDec() +  fluxMagInfo.coordList[1].getDec() +
                        fluxMagInfo.coordList[2].getDec() +  fluxMagInfo.coordList[3].getDec())/ 4.
            x, y = wcs.skyToPixel(raCenter,decCenter)
            imageScaler.xList.append(x)
            imageScaler.yList.append(y)          
            imageScaler.scaleList.append(self.scaleFromFluxMag0(fluxMagInfo.fluxMag0).scale)
            #self.imageScaler.scaleErrList.append(self.fluxMag0ToScale(fluxMagInfo.fluxMag0Sigma))

        return imageScaler
