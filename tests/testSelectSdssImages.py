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
from lsst.daf.persistence import DbAuth
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
from lsst.obs.sdss.selectSdssImages import SelectSdssImagesTask

Database = "test_select_sdss_images"


def getCoordList(minRa, minDec, maxRa, maxDec):
    degList = (
        (minRa, minDec),
        (maxRa, minDec),
        (maxRa, maxDec),
        (minRa, maxDec),
    )
    return tuple(afwCoord.IcrsCoord(afwGeom.Point2D(d[0], d[1]), afwGeom.degrees) for d in degList)

class SelectSdssImagesTestCase(unittest.TestCase):
    """A test case for SelectSdssImagesTask."""
    def testMaxFwhm(self):
        """Test config.maxFwhm
        """
        for maxFwhm in (1.2, 1.3):
            config = SelectSdssImagesTask.ConfigClass()
            config.database = Database
            config.maxFwhm = maxFwhm
            task = SelectSdssImagesTask(config=config)
            coordList = getCoordList(333.746,-0.63606,334.522,-0.41341)
            filter = "g"
            expInfoList = task.run(coordList=coordList, filter=filter).exposureInfoList
            self.assertEqual(tuple(expInfo for expInfo in expInfoList if expInfo.fwhm > maxFwhm), ())

    def testMaxAirmass(self):
        """Test config.maxAirmass
        """
        for maxAirmass in (1.2, 1.3):
            config = SelectSdssImagesTask.ConfigClass()
            config.database = Database
            config.maxAirmass = maxAirmass
            task = SelectSdssImagesTask(config=config)
            coordList = getCoordList(333.746,-0.63606,334.522,-0.41341)
            filter = "g"
            expInfoList = task.run(coordList=coordList, filter=filter).exposureInfoList
            self.assertEqual(tuple(expInfo for expInfo in expInfoList if expInfo.airmass > maxAirmass), ())

    def testMaxSky(self):
        """Test config.maxSky
        """
        for maxSky in (5.0e-9, 1.0e-8):
            config = SelectSdssImagesTask.ConfigClass()
            config.database = Database
            config.maxSky = maxSky
            task = SelectSdssImagesTask(config=config)
            coordList = getCoordList(333.746,-0.63606,334.522,-0.41341)
            filter = "g"
            expInfoList = task.run(coordList=coordList, filter=filter).exposureInfoList
            self.assertEqual(tuple(expInfo for expInfo in expInfoList if expInfo.sky > maxSky), ())

    def testQuality(self):
        """Test config.quality
        """
        for quality in (1, 2, 3):
            config = SelectSdssImagesTask.ConfigClass()
            config.database = Database
            config.quality = quality
            task = SelectSdssImagesTask(config=config)
            coordList = getCoordList(333.746,-0.63606,334.522,-0.41341)
            filter = "g"
            expInfoList = task.run(coordList=coordList, filter=filter).exposureInfoList
            self.assertEqual(tuple(expInfo for expInfo in expInfoList if expInfo.quality < quality), ())

    def testCullBlacklisted(self):
        """Test config.cullBlacklisted
        """
        for cullBlacklisted in (False, True):
            config = SelectSdssImagesTask.ConfigClass()
            config.database = Database
            config.quality = 1
            config.cullBlacklisted = cullBlacklisted
            task = SelectSdssImagesTask(config=config)
            coordList = getCoordList(300,-0.63606,302,-0.41341)
            filter = "g"
            expInfoList = task.run(coordList=coordList, filter=filter).exposureInfoList
            blacklistedList = tuple(expInfo for expInfo in expInfoList if expInfo.isBlacklisted)
            if cullBlacklisted:
                self.assertEqual(blacklistedList, ())
            else:
                self.assertGreater(len(blacklistedList), 0)

    def testCamcols(self):
        """Test config.camcols
        """
        for camcols in ((1, 3, 4), (2,)):
            config = SelectSdssImagesTask.ConfigClass()
            config.database = Database
            config.camcols = camcols
            task = SelectSdssImagesTask(config=config)
            coordList = getCoordList(333.746,-0.63606,334.522,-0.41341)
            filter = "g"
            expInfoList = task.run(coordList=coordList, filter=filter).exposureInfoList
            self.assertEqual(
                tuple(expInfo for expInfo in expInfoList if expInfo.dataId["camcol"] not in camcols), ())

    def testStrip(self):
        """Test config.strip
        """
        for strip in ("S", "N", "Auto", "Both"):
            config = SelectSdssImagesTask.ConfigClass()
            config.database = Database
            config.strip = strip
            task = SelectSdssImagesTask(config=config)
            coordList = getCoordList(333.746,-0.63606,334.522,-0.41341)
            dataId = {"filter": "g", "patch": "20,5"}
            runArgDict = task._runArgDictFromDataId(dataId)
            expInfoList = task.run(coordList=coordList, **runArgDict).exposureInfoList
            if strip in ("S","N"):
                self.assertEqual(tuple(expInfo for expInfo in expInfoList if expInfo.strip != strip), ())
            elif strip =="Auto":
                self.assertEqual(tuple(expInfo for expInfo in expInfoList if expInfo.strip != 'N'), ())
            #no assert for "Both"

    def testRejectWholeRuns(self):
        """Test config.rejectWholeRuns
        """
        config = SelectSdssImagesTask.ConfigClass()
        config.database = Database
        config.maxFwhm = 1.25 # make sure to cut out some partial runs due to bad exposures
        config.rejectWholeRuns = True
        task = SelectSdssImagesTask(config=config)
        minRa = 333.746
        maxRa = 334.522
        coordList = getCoordList(minRa,-0.63606,maxRa,-0.41341)
        filter = "g"
        expInfoList = task.run(coordList=coordList, filter=filter).exposureInfoList
        runExpInfoDict = dict()
        for expInfo in expInfoList:
            run = expInfo.dataId["run"]
            if run in runExpInfoDict:
                runExpInfoDict[run].append(expInfo)
            else:
                runExpInfoDict[run] = [expInfo]

        self.checkExpList(minRa, maxRa, runExpInfoDict)

    def testMaxExposures(self):
        """Test config.maxExposures
        """
        for maxExposures in (0, 6):
            config = SelectSdssImagesTask.ConfigClass()
            config.database = Database
            config.maxExposures = maxExposures
            task = SelectSdssImagesTask(config=config)
            coordList = getCoordList(333.746,-0.63606,334.522,-0.41341)
            filter = "g"
            expInfoList = task.run(coordList=coordList, filter=filter).exposureInfoList
            self.assertEqual(len(expInfoList), maxExposures)

    def testMaxRuns(self):
        """Test config.maxRuns
        """
        for maxRuns in (0, 2):
            config = SelectSdssImagesTask.ConfigClass()
            config.database = Database
            config.maxRuns = maxRuns
            task = SelectSdssImagesTask(config=config)
            coordList = getCoordList(333.746,-0.63606,334.522,-0.41341)
            filter = "g"
            expInfoList = task.run(coordList=coordList, filter=filter).exposureInfoList
            runSet = set(expInfo.dataId["run"] for expInfo in expInfoList)
            self.assertEqual(len(runSet), maxRuns)

    def testQScore(self):
        """Test QScore sorting
        """
        config = SelectSdssImagesTask.ConfigClass()
        config.database = Database
        config.quality = 1
        config.rejectWholeRuns = False
        task = SelectSdssImagesTask(config=config)
        coordList = getCoordList(333.746,-0.63606,334.522,-0.41341)
        filter = "g"
        expInfoList = task.run(coordList=coordList, filter=filter).exposureInfoList
        qscoreList = list(expInfo.qscore for expInfo in expInfoList)
        self.assertEqual(qscoreList, sorted(qscoreList))
        bestExp = expInfoList[0]
        worstExp = expInfoList[-1]
        self.assertGreater(worstExp.fwhm, bestExp.fwhm)
        self.assertGreater(worstExp.sky, bestExp.sky)
        self.assertGreater(bestExp.quality, worstExp.quality)
        self.assertEqual(bestExp.quality, 3)

    def testConfigValidate(self):
        """Test validation of config
        """
        for maxExposures in (None, 1):
            for maxRuns in (None, 1):
                config = SelectSdssImagesTask.ConfigClass()
                config.database = Database
                config.maxExposures = maxExposures
                config.maxRuns = maxRuns
                if maxExposures and maxRuns:
                    self.assertRaises(Exception, config.validate)
                else:
                    config.validate() # should not raise an exception

        config = SelectSdssImagesTask.ConfigClass()
        config.database = Database
        config.table = "invalid*name"
        self.assertRaises(Exception, config.validate)

    def testFilterValidation(self):
        """Test filter name validation
        """
        coordList = getCoordList(333.7,-0.6,333.71,-0.59)
        config = SelectSdssImagesTask.ConfigClass()
        config.database = Database
        task = SelectSdssImagesTask(config=config)
        for charVal in range(ord("a"), ord("z")+1):
            filter = chr(charVal)
            if filter in ("u", "g", "r", "i", "z"):
                task.run(coordList=coordList, filter=filter)
            else:
                self.assertRaises(Exception, task.run, coordList, filter)

    def testAcrossWrap(self):
        """Test rejectWholeRuns across the RA 0/360 boundary
        """
        config = SelectSdssImagesTask.ConfigClass()
        config.database = Database
        config.rejectWholeRuns = True
        task = SelectSdssImagesTask(config=config)
        minRa = 359
        maxRa = 1
        coordList = getCoordList(minRa,-0.63606,maxRa,-0.41341)
        filter = "g"
        expInfoList = task.run(coordList=coordList, filter=filter).exposureInfoList
        runExpInfoDict = dict()
        for expInfo in expInfoList:
            run = expInfo.dataId["run"]
            if run in runExpInfoDict:
                runExpInfoDict[run].append(expInfo)
            else:
                runExpInfoDict[run] = [expInfo]

        self.assertEqual(len(runExpInfoDict), 6)
        self.checkExpList(minRa, maxRa, runExpInfoDict)

    def checkExpList(self, minRa, maxRa, runExpInfoDict):
        """Check that all exposures runExpInfoDict are within the specified range

        @param[in] minRa: minimum RA (degrees)
        @param[in] maxRa: maxinum RA (degrees)
        @param[in] runExpInfoDict: a dict of run: list of ExposureInfo objects
        """
        minRaAngle = minRa * afwGeom.degrees
        maxRaAngle = maxRa * afwGeom.degrees
        minRaAngle.wrapNear(maxRaAngle)
        ctrRaAngle = (minRaAngle + maxRaAngle) * 0.5
        raDegList = []
        for expInfoList in runExpInfoDict.itervalues():
            for expInfo in expInfoList:
                raAngleList = [coord.getLongitude() for coord in expInfo.coordList]
                for raAngle in raAngleList:
                    raAngle.wrapNear(ctrRaAngle)
                raDegList += [raAngle.asDegrees() for raAngle in raAngleList]
            raDegList.sort()
        self.assertGreaterEqual(minRa, raDegList[0])
        self.assertGreaterEqual(raDegList[-1], maxRa)


    def testTable(self):
        """Test config.table
        """
        config = SelectSdssImagesTask.ConfigClass()
        config.table = "Bad_table_name_JutmgQEXm76O38VDtcNAICLrtQiSQ64y"
        task = SelectSdssImagesTask(config=config)
        for coordList in [None, getCoordList(333.746,-0.63606,334.522,-0.41341)]:
            filter = "g"
            self.assertRaises(Exception, task.run, coordList, filter)

    def testWholeSky(self):
        """Test whole-sky search (slow so don't do much)
        """
        config = SelectSdssImagesTask.ConfigClass()
        config.database = Database
        config.camcols = (2,)
        config.quality = 1
        config.rejectWholeRuns = False
        task = SelectSdssImagesTask(config=config)
        coordList = None
        filter = "g"
        expInfoList = task.run(coordList=coordList, filter=filter).exposureInfoList
        self.assertEqual(tuple(expInfo for expInfo in expInfoList if expInfo.quality < config.quality), ())
        print "found %s exposures" % (len(expInfoList),)
        self.assertEqual(tuple(expInfo for expInfo in expInfoList \
            if expInfo.dataId["camcol"] not in config.camcols), ())


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(SelectSdssImagesTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    config = SelectSdssImagesTask.ConfigClass()
    try:
        DbAuth.username(config.host, str(config.port)),
    except Exception:
        print "Warning: did not find host=%s, port=%s in your db-auth file; skipping SelectSdssImagesTask unit tests" % \
            (config.host, str(config.port))
        return

    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
