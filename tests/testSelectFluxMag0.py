#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011, 2012, 2013 LSST Corporation.
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

"""Test lsst.obs.sdss.selectFluxMag0Task and integration with lsst.obs.sdss.scaleSdssZeroPointTask
"""
import unittest

import lsst.utils.tests as utilsTests
import lsst.daf.base
from lsst.daf.persistence import DbAuth
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.coord as afwCoord
from lsst.obs.sdss.scaleSdssZeroPoint import ScaleSdssZeroPointTask
from lsst.obs.sdss.selectFluxMag0 import SelectSdssFluxMag0Task


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
class WrapDataId():
    """A container for dataId that looks like dataRef to computeImageScaler()
    """
    def __init__(self, dataId):
        self.dataId = dataId

class ScaleSdssZeroPointTaskTestCase(unittest.TestCase):
    """A test case for ScaleSdssZeroPointTask
    """
    def makeTestExposure(self, xNumPix=2060, yNumPix=1967):
        """
        Create and return an exposure with wcs. Wcs is chosen such that the exposure is
        completely covered by the Science_Ccd_Exposure table in the database: test_select_sdss_images
        """
        metadata = lsst.daf.base.PropertySet()
        metadata.set("NAXIS", 2)
        metadata.set("RADECSYS", "ICRS")
        metadata.set("EQUINOX", 2000.)
        metadata.setDouble("CRVAL1", 315.)
        metadata.setDouble("CRVAL2", 0.)
        metadata.setDouble("CRPIX1", 68030.)
        metadata.setDouble("CRPIX2", 30.)
        metadata.set("CTYPE1", "RA---CEA")
        metadata.set("CTYPE2", "DEC--CEA")
        metadata.setDouble("CD1_1", -0.00011)
        metadata.setDouble("CD1_2", 0.000000)
        metadata.setDouble("CD2_2", 0.000110)
        metadata.setDouble("CD2_1", 0.000000)
        metadata.set("CUNIT1", "deg")
        metadata.set("CUNIT2", "deg")
        metadata.set("LTV1", -341970)
        metadata.set("LTV2", -11412)
        #exposure needs a wcs and a bbox
        wcs = afwImage.makeWcs(metadata)
        bbox = afwGeom.Box2I(afwGeom.Point2I(341970, 11412), afwGeom.Extent2I(xNumPix, yNumPix))
        exposure = afwImage.ExposureF(bbox, wcs)
        mi = exposure.getMaskedImage()
        mi.set(1.0)
        mi.getVariance().set(1.0)
        return exposure

    def testSelectFluxMag0(self):
        """Test SelectFluxMag0"""
        config = SelectSdssFluxMag0Task.ConfigClass()
        config.database = "test_select_sdss_images"
        task = SelectSdssFluxMag0Task(config=config)
        run = 4192
        filter = 'i'
        dataId = {'run': run, "filter": filter}
        coordList = [afwCoord.Coord(5.62839*afwGeom.radians, -5.66359e-05*afwGeom.radians),
                     afwCoord.Coord(5.62444*afwGeom.radians, -5.66359e-05*afwGeom.radians),
                     afwCoord.Coord(5.62444*afwGeom.radians, 0.00371974*afwGeom.radians),
                     afwCoord.Coord(5.62839*afwGeom.radians, 0.00371974*afwGeom.radians)]
        fmInfoStruct = task.run(dataId, coordList)
        fmInfoList = fmInfoStruct.fluxMagInfoList
        self.assertEqual(2, len(fmInfoList))


    def testScaleZeroPoint(self):
        ZEROPOINT = 27
        self.sctrl = afwMath.StatisticsControl()
        self.sctrl.setNanSafe(True)

        config = ScaleSdssZeroPointTask.ConfigClass()
        config.zeroPoint = ZEROPOINT
        config.interpStyle = "NATURAL_SPLINE"
        config.selectFluxMag0.database = "test_select_sdss_images"
        zpScaler = ScaleSdssZeroPointTask(config=config)


        outCalib = zpScaler.getCalib()
        self.assertAlmostEqual(outCalib.getMagnitude(1.0), ZEROPOINT)

        exposure = self.makeTestExposure()
        #create dataId for exposure. Visit is only field needed. Others ignored.
        exposureId = {'ignore_fake_key': 1234, 'run': 4192, 'filter': 'i'}

        #test methods: computeImageScale(), scaleMaskedImage(), getInterpImage()
        dataRef = WrapDataId(exposureId)
        imageScaler = zpScaler.computeImageScaler(exposure,dataRef)
        scaleFactorIm = imageScaler.getInterpImage(exposure.getBBox())

        predScale = 0.402867736 #image mean for "NATURAL_SPLINE"
        self.assertAlmostEqual(afwMath.makeStatistics(scaleFactorIm, afwMath.MEAN, self.sctrl).getValue(),
                               predScale)

        mi = exposure.getMaskedImage()
        imageScaler.scaleMaskedImage(mi)
        pixel11 = scaleFactorIm.getArray()[1,1]
        self.assertAlmostEqual(mi.get(1,1)[0], pixel11) #check image plane scaled
        self.assertAlmostEqual(mi.get(1,1)[2], pixel11**2) #check variance plane scaled

        exposure.setCalib(zpScaler.getCalib())
        self.assertAlmostEqual(exposure.getCalib().getFlux(ZEROPOINT), 1.0)


    def makeCalib(self, zeroPoint):
        calib = afwImage.Calib()
        fluxMag0 = 10**(0.4 * zeroPoint)
        calib.setFluxMag0(fluxMag0, 1.0)
        return calib


def suite():
    """Return a suite containing all the test cases in this module.
    """
    utilsTests.init()

    suites = [
        unittest.makeSuite(ScaleSdssZeroPointTaskTestCase),
        unittest.makeSuite(utilsTests.MemoryTestCase),
    ]

    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """Run the tests"""
    config = ScaleSdssZeroPointTask.ConfigClass()
    try:
        DbAuth.username(config.selectFluxMag0.host, str(config.selectFluxMag0.port)),
    except Exception, e:
        print "Warning: did not find host=%s, port=%s in your db-auth file; or %s " \
              "skipping unit tests" % \
            (config.selectFluxMag0.host, str(config.selectFluxMag0.port), e)
        return

    utilsTests.run(suite(), shouldExit)


if __name__ == "__main__":
    run(True)
