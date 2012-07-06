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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import unittest
import lsst.utils.tests as utilsTests

import os
from lsst.daf.persistence import DbAuth
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
from lsst.obs.sdss.selectSdssImages import SelectSdssImagesTask


def getCoordList(minRa, minDec, maxRa, maxDec):
    degList = (
        (minRa, minDec),
        (maxRa, minDec),
        (maxRa, maxDec),
        (minRa, maxDec),
    )
    return tuple(afwCoord.IcrsCoord(afwGeom.Point2D(d[0], d[1]), afwGeom.degrees) for d in degList)

class SdssMapperTestCase(unittest.TestCase):
    """A test case for SelectSdssImagesTask."""
    def testMaxFwhm(self):
        """Test config.maxFwhm
        """
        config = SelectSdssImagesTask.ConfigClass()
        for maxFwhm in (1.2, 1.3):
            config.maxFwhm = maxFwhm
            task = SelectSdssImagesTask(config=config)
            coordList = getCoordList(333.746,-0.63606,334.522,-0.41341)
            filter = "g"
            expInfoList = task.run(coordList, filter).exposureInfoList
            self.assertEqual(tuple(expInfo for expInfo in expInfoList if expInfo.fwhm > maxFwhm), ())

    def testMaxAirmass(self):
        """Test config.maxAirmass
        """
        config = SelectSdssImagesTask.ConfigClass()
        for maxAirmass in (1.2, 1.3):
            config.maxAirmass = maxAirmass
            task = SelectSdssImagesTask(config=config)
            coordList = getCoordList(333.746,-0.63606,334.522,-0.41341)
            filter = "g"
            expInfoList = task.run(coordList, filter).exposureInfoList
            self.assertEqual(tuple(expInfo for expInfo in expInfoList if expInfo.airmass > maxAirmass), ())

    def testMaxSky(self):
        """Test config.maxSky
        """
        config = SelectSdssImagesTask.ConfigClass()
        for maxSky in (5.0e-9, 1.0e-8):
            config.maxSky = maxSky
            task = SelectSdssImagesTask(config=config)
            coordList = getCoordList(333.746,-0.63606,334.522,-0.41341)
            filter = "g"
            expInfoList = task.run(coordList, filter).exposureInfoList
            self.assertEqual(tuple(expInfo for expInfo in expInfoList if expInfo.sky > maxSky), ())
    
    def testQuality(self):
        """Test config.quality
        """
        config = SelectSdssImagesTask.ConfigClass()
        for quality in (1, 2, 3):
            config.quality = quality
            task = SelectSdssImagesTask(config=config)
            coordList = getCoordList(333.746,-0.63606,334.522,-0.41341)
            filter = "g"
            expInfoList = task.run(coordList, filter).exposureInfoList
            self.assertEqual(tuple(expInfo for expInfo in expInfoList if expInfo.quality < quality), ())
    
    def testCullBlacklisted(self):
        """Test config.cullBlacklisted
        """
        config = SelectSdssImagesTask.ConfigClass()
        for cullBlacklisted in (False, True):
            config.quality = 1
            config.cullBlacklisted = cullBlacklisted
            task = SelectSdssImagesTask(config=config)
            coordList = getCoordList(333.746,-0.63606,334.522,-0.41341)
            filter = "g"
            expInfoList = task.run(coordList, filter).exposureInfoList
            blacklistedList = tuple(expInfo for expInfo in expInfoList if expInfo.isBlacklisted)
            if cullBlacklisted:
                self.assertEqual(blacklistedList, ())
            else:
                print "WARNING: test disabled because the blacklisted field is not yet populated"
                # self.assertGreater(len(blacklistedList), 0)
    
    def testCamcols(self):
        """Test config.camcols
        """
        config = SelectSdssImagesTask.ConfigClass()
        for camcols in ((1, 3, 4), (2,)):
            config.camcols = camcols
            task = SelectSdssImagesTask(config=config)
            coordList = getCoordList(333.746,-0.63606,334.522,-0.41341)
            filter = "g"
            expInfoList = task.run(coordList, filter).exposureInfoList
            self.assertEqual(tuple(expInfo for expInfo in expInfoList if expInfo.dataId["camcol"] not in camcols), ())

    def testStrip(self):
        """Test config.strip
        """
        config = SelectSdssImagesTask.ConfigClass()
        for strip in ("S", "N"):
            config.strip = strip
            task = SelectSdssImagesTask(config=config)
            coordList = getCoordList(333.746,-0.63606,334.522,-0.41341)
            filter = "g"
            expInfoList = task.run(coordList, filter).exposureInfoList
            self.assertEqual(tuple(expInfo for expInfo in expInfoList if expInfo.strip != strip), ())

    def testRejectWholeRuns(self):
        """Test config.rejectWholeRuns
        """
        config = SelectSdssImagesTask.ConfigClass()
        config.maxFwhm = 1.25 # make sure to cut out some partial runs
        config.rejectWholeRuns = True
        task = SelectSdssImagesTask(config=config)
        minRa = 333.746
        maxRa = 334.522
        coordList = getCoordList(minRa,-0.63606,maxRa,-0.41341)
        filter = "g"
        expInfoList = task.run(coordList, filter).exposureInfoList
        runExpInfoDict = dict()
        for expInfo in expInfoList:
            run = expInfo.dataId["run"]
            if run in runExpInfoDict:
                runExpInfoDict[run].append(expInfo)
            else:
                runExpInfoDict[run] = [expInfo]
            
        for run, expInfoList in runExpInfoDict.iteritems():
            raList = []
            for expInfo in expInfoList:
                raList += [coord.getLongitude().asDegrees() for coord in expInfo.coordList]
            raList.sort()
            self.assertGreaterEqual(minRa, raList[0])
            self.assertGreaterEqual(raList[-1], maxRa)
    
    def testMaxExposures(self):
        """Test config.maxExposures
        """
        config = SelectSdssImagesTask.ConfigClass()
        for maxExposures in (0, 6):
            config.maxExposures = maxExposures
            task = SelectSdssImagesTask(config=config)
            coordList = getCoordList(333.746,-0.63606,334.522,-0.41341)
            filter = "g"
            expInfoList = task.run(coordList, filter).exposureInfoList
            self.assertEqual(len(expInfoList), maxExposures)
        
    def testMaxRuns(self):                
        """Test config.maxRuns
        """
        config = SelectSdssImagesTask.ConfigClass()
        for maxRuns in (0, 2):
            config.maxRuns = maxRuns
            task = SelectSdssImagesTask(config=config)
            coordList = getCoordList(333.746,-0.63606,334.522,-0.41341)
            filter = "g"
            expInfoList = task.run(coordList, filter).exposureInfoList
            runSet = set(expInfo.dataId["run"] for expInfo in expInfoList)
            self.assertEqual(len(runSet), maxRuns)
    
    def testQScore(self):
        """Test QScore sorting
        """
        config = SelectSdssImagesTask.ConfigClass()
        config.quality = 1
        config.rejectWholeRuns = False
        task = SelectSdssImagesTask(config=config)
        coordList = getCoordList(333.746,-0.63606,334.522,-0.41341)
        filter = "g"
        expInfoList = task.run(coordList, filter).exposureInfoList
        qscoreList = list(expInfo.qscore for expInfo in expInfoList)
        self.assertEqual(qscoreList, sorted(qscoreList))
        bestExp = expInfoList[0]
        worstExp = expInfoList[-1]
        self.assertGreater(worstExp.fwhm, bestExp.fwhm)
        self.assertGreater(worstExp.sky, bestExp.sky)
        self.assertGreater(bestExp.quality, worstExp.quality)
        self.assertEqual(bestExp.quality, 3)
    
    def testConfigValidate(self):
        config = SelectSdssImagesTask.ConfigClass()
        for maxExposures in (None, 1):
            config.maxExposures = maxExposures
            for maxRuns in (None, 1):
                config.maxRuns = maxRuns
                if maxExposures and maxRuns:
                    self.assertRaises(Exception, config.validate)
                else:
                    config.validate() # should not raise an exception

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(SdssMapperTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    config = SelectSdssImagesTask.ConfigClass()
    try:
        user = DbAuth.username(config.host, str(config.port)),
    except Exception:
        print "Warning: did not find host=%s, port=%s in your db-auth file; skipping SelectSdssImagesTask unit tests" % \
            (config.host, str(config.port))
        return

    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
