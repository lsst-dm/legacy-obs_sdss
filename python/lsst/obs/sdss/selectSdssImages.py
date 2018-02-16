from __future__ import absolute_import, division, print_function
from builtins import range
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
import os
import re

import MySQLdb
import numpy as np

import lsst.pex.config as pexConfig
from lsst.afw.coord import IcrsCoord
import lsst.afw.geom as afwGeom
from lsst.daf.persistence import DbAuth
import lsst.pipe.base as pipeBase
from lsst.pipe.tasks.selectImages import DatabaseSelectImagesConfig, BaseSelectImagesTask, BaseExposureInfo

__all__ = ["SelectSdssImagesTask"]


class SelectSdssImagesConfig(DatabaseSelectImagesConfig):
    """Config for SelectSdssImagesTask
    """
    table = pexConfig.Field(
        doc="Name of database table",
        dtype=str,
        default="SeasonFieldQuality_Test",
    )
    maxFwhm = pexConfig.Field(
        doc="maximum FWHM (arcsec)",
        dtype=float,
        optional=True,
    )
    maxSky = pexConfig.Field(
        doc="maximum sky level (maggies/arcsec^2)",
        dtype=float,
        optional=True,
    )
    maxAirmass = pexConfig.Field(
        doc="Maximum airmass",
        dtype=float,
        optional=True,
    )
    quality = pexConfig.ChoiceField(
        doc="SDSS quality flag",
        dtype=int,
        default=3,
        allowed={1: "All data", 2: "Flagged ACCEPTABLE or GOOD", 3: "Flagged GOOD"},
    )
    cullBlacklisted = pexConfig.Field(
        doc="Omit blacklisted images? (Some run/field combinations have been blacklisted even though "
            "the quality metrics may have missed them.)",
        dtype=bool,
        default=True,
    )
    camcols = pexConfig.ListField(
        doc="Which camcols to include? If None then include all",
        dtype=int,
        optional=True,
    )
    strip = pexConfig.Field(
        doc="Strip: N, S, Both for both, or Auto to pick the strip based on the dataId",
        dtype=str,
        default="Auto",
        optional=False,
    )
    rejectWholeRuns = pexConfig.Field(
        doc="If any exposure in the region is bad or the run does not cover the whole region, "
            "then reject the whole run?",
        dtype=bool,
        default=True,
    )
    maxRuns = pexConfig.Field(
        doc="maximum runs to select; ignored if None; cannot be specified with maxExposures",
        dtype=int,
        optional=True,
    )

    def setDefaults(self):
        BaseSelectImagesTask.ConfigClass.setDefaults(self)
        self.host = "lsst-db.ncsa.illinois.edu"
        self.port = 3306
        self.database = "krughoff_SDSS_quality_db"

    def validate(self):
        BaseSelectImagesTask.ConfigClass.validate(self)
        if (self.maxRuns is not None) and (self.maxExposures is not None):
            raise RuntimeError("maxRuns=%s or maxExposures=%s must be None" %
                               (self.maxRuns, self.maxExposures))
        if not re.match(r"^[a-zA-Z0-9_.]+$", self.table):
            raise RuntimeError("table=%r is an invalid name" % (self.table,))


class ExposureInfo(BaseExposureInfo):
    """Data about a selected exposure

    Data includes:
    - dataId: data ID of exposure (a dict)
    - coordList: a list of corner coordinates of the exposure (list of IcrsCoord)
    - fwhm: mean FWHM of exposure
    - quality: quality field from self.config.table
    """

    def __init__(self, result):
        """Set exposure information based on a query result from a db connection
        """
        self._ind = -1

        dataId = dict(
            run=result[self._nextInd],
            rerun=result[self._nextInd],
            camcol=result[self._nextInd],
            field=result[self._nextInd],
            filter=result[self._nextInd],
        )
        coordList = []
        for i in range(4):
            coordList.append(
                IcrsCoord(
                    afwGeom.Angle(result[self._nextInd], afwGeom.degrees),
                    afwGeom.Angle(result[self._nextInd], afwGeom.degrees),
                )
            )
        BaseExposureInfo.__init__(self, dataId, coordList)

        self.strip = result[self._nextInd]
        self.fwhm = result[self._nextInd]
        self.sky = result[self._nextInd]
        self.airmass = result[self._nextInd]
        self.quality = result[self._nextInd]
        self.isBlacklisted = result[self._nextInd]

        # compute RHL quality factors
        self.q = self.sky * (self.fwhm**2)
        self.qscore = None  # not known yet

    @property
    def _nextInd(self):
        """Return the next index

        This allows one to progress through returned database columns,
        repeatedly using _nextInd, e.g.:
            dataId = dict(
                raft = result[self._nextInd],
                visit = result[self._nextInd],
                sensor = result[self._nextInd],
                filter = result[self._nextInd],
            )
        """
        self._ind += 1
        return self._ind

    @staticmethod
    def getColumnNames():
        """Get database columns to retrieve, in a format useful to the database interface

        @return database column names as list of strings
        """
        return (
            "run rerun camcol field filter ra1 dec1 ra2 dec2 ra3 dec3 ra4 dec4".split() +
            "strip psfWidth sky airmass quality isblacklisted".split()
        )

    def __hash__(self):
        """
        Objects used in sets are required to be hashable, which requires (among other things) a
        ``__hash__`` method. In python 2, if no ``__hash`` method was explicitly set, the id is used.
        Python 3 does not have this behavior, so we explicitly return the id.
        """
        return id(self)


class SelectSdssImagesTask(BaseSelectImagesTask):
    """Select SDSS images suitable for coaddition
    """
    ConfigClass = SelectSdssImagesConfig

    @pipeBase.timeMethod
    def run(self, coordList, filter, strip=None):
        """Select SDSS images suitable for coaddition in a particular region

        @param[in] filter: filter for images (one of "u", "g", "r", "i" or "z")
        @param[in] coordList: list of coordinates defining region of interest

        @return a pipeBase Struct containing:
        - exposureInfoList: a list of ExposureInfo objects

        @raise RuntimeError if filter not one of "u", "g", "r", "i" or "z"
        """
        if filter not in set(("u", "g", "r", "i", "z")):
            raise RuntimeError("filter=%r is an invalid name" % (filter,))

        read_default_file = os.path.expanduser("~/.my.cnf")

        try:
            open(read_default_file)
            kwargs = dict(
                read_default_file=read_default_file,
            )
        except IOError:
            kwargs = dict(
                user=DbAuth.username(self.config.host, str(self.config.port)),
                passwd=DbAuth.password(self.config.host, str(self.config.port)),
            )

        db = MySQLdb.connect(
            host=self.config.host,
            port=self.config.port,
            db=self.config.database,
            **kwargs
        )
        cursor = db.cursor()

        columnNames = tuple(ExposureInfo.getColumnNames())
        if not columnNames:
            raise RuntimeError("Bug: no column names")
        queryStr = "select %s " % (", ".join(columnNames),)
        dataTuple = ()  # tuple(columnNames)

        if coordList is not None:
            # look for exposures that overlap the specified region

            # create table scisql.Region containing patch region
            coordStrList = ["%s, %s" % (c.getLongitude().asDegrees(),
                                        c.getLatitude().asDegrees()) for c in coordList]
            coordStr = ", ".join(coordStrList)
            coordCmd = "call scisql.scisql_s2CPolyRegion(scisql_s2CPolyToBin(%s), 10)" % (coordStr,)
            cursor.execute(coordCmd)
            cursor.nextset()  # ignore one-line result of coordCmd

            queryStr += """
from %s as ccdExp,
    (select distinct fieldid
    from SeasonFieldQuality_To_Htm10 as ccdHtm inner join scisql.Region
    on (ccdHtm.htmId10 between scisql.Region.htmMin and scisql.Region.htmMax)
    where ccdHtm.filter = %%s) as idList
where ccdExp.fieldid = idList.fieldid and """ % (self.config.table,)
            dataTuple += (filter,)
        else:
            # no region specified; look over the whole sky
            queryStr += """
from %s as ccdExp where """ % (self.config.table,)

        # compute where clauses as a list of (clause, data)
        whereDataList = [
            ("filter = %s", filter),
        ]

        if self.config.camcols is not None:
            whereDataList.append(_whereDataFromList("camcol", self.config.camcols))

        if strip is not None:  # None corresponds to query for both strips: no constraint added
            whereDataList.append(("strip = %s", strip))

        queryStr += " and ".join(wd[0] for wd in whereDataList)
        dataTuple += tuple(wd[1] for wd in whereDataList)

        self.log.info("queryStr=%r; dataTuple=%s" % (queryStr, dataTuple))

        cursor.execute(queryStr, dataTuple)
        exposureInfoList = [ExposureInfo(result) for result in cursor]

        runExpInfoSetDict = dict()

        for expInfo in exposureInfoList:
            run = expInfo.dataId["run"]
            expInfoSet = runExpInfoSetDict.get(run)
            if expInfoSet:
                expInfoSet.add(expInfo)
            else:
                runExpInfoSetDict[run] = set([expInfo])

        self.log.info("Before quality cuts found %d exposures in %d runs" %
                      (len(exposureInfoList), len(runExpInfoSetDict)))

        goodRunSet = set()
        goodExposureInfoList = []
        if self.config.rejectWholeRuns:
            # reject runs for which any exposure does not meet our quality criteria
            # or the run begins or ends in the region
            regionRaRange = None
            regionCtrRa = None
            if coordList is not None:
                regionRaRange = _computeRaRange(coordList)
                regionCtrRa = (regionRaRange[0] + regionRaRange[1]) * 0.5

            numRangeCuts = 0
            for run, expInfoSet in runExpInfoSetDict.items():
                runRaRange = None
                for expInfo in expInfoSet:
                    if self._isBadExposure(expInfo):
                        break

                    if regionRaRange is not None:
                        expRaRange = _computeRaRange(expInfo.coordList, ctrRa=regionCtrRa)
                        if runRaRange is None:
                            runRaRange = expRaRange
                        else:
                            runRaRange = (min(runRaRange[0], expRaRange[0]),
                                          max(runRaRange[1], expRaRange[1]))
                else:
                    # all exposures in this run are valid;
                    # if approriate, check that the run starts and ends outside the region
                    if regionRaRange is not None:
                        if (runRaRange[0] > regionRaRange[0]) or (runRaRange[1] < regionRaRange[1]):
                            numRangeCuts += 1
                            continue

                    goodExposureInfoList += list(expInfoSet)
                    goodRunSet.add(run)
            self.log.info("Rejected %d whole runs, including %d for incomplete range" %
                          (len(runExpInfoSetDict) - len(goodRunSet), numRangeCuts))
        else:
            # reject individual exposures which do not meet our quality criteria
            for expInfo in exposureInfoList:
                if not self._isBadExposure(expInfo):
                    goodExposureInfoList.append(expInfo)
                    goodRunSet.add(expInfo.dataId["run"])
            self.log.info("Rejected %d individual exposures" %
                          (len(exposureInfoList) - len(goodExposureInfoList),))

        exposureInfoList = goodExposureInfoList

        self.log.info("After quality cuts, found %d exposures in %d runs" %
                      (len(exposureInfoList), len(goodRunSet)))

        if exposureInfoList:
            # compute qscore according to RHL's formula and sort by it
            qArr = np.array([expInfo.q for expInfo in exposureInfoList])
            qMax = np.percentile(qArr, 95.0)
            for expInfo in exposureInfoList:
                expInfo.qscore = (expInfo.q / qMax) - expInfo.quality
            exposureInfoList.sort(key=lambda expInfo: expInfo.qscore)

            if self.config.maxExposures is not None:
                # select config.maxExposures exposures with the highest qscore
                exposureInfoList = exposureInfoList[0:self.config.maxExposures]
                self.log.info("After maxExposures cut, found %d exposures" %
                              (len(exposureInfoList),))

            elif self.config.maxRuns is not None:
                # select config.maxRuns runs with the highest median qscore
                # (of those exposures that overlap the patch)
                runQualListDict = dict()
                for expInfo in exposureInfoList:
                    run = expInfo.dataId["run"]
                    qualList = runQualListDict.get(run)
                    if qualList:
                        qualList.append(expInfo.qscore)
                    else:
                        runQualListDict[run] = [expInfo.qscore]

                if len(runQualListDict) > self.config.maxRuns:
                    qualRunList = []
                    for run, qualList in runQualListDict.items():
                        runQscore = np.median(qualList)
                        qualRunList.append((runQscore, run))
                    qualRunList.sort()
                    qualRunList = qualRunList[0:self.config.maxRuns]

                    goodRunSet = set(qr[1] for qr in qualRunList)
                    exposureInfoList = [ei for ei in exposureInfoList if ei.dataId["run"] in goodRunSet]

        return pipeBase.Struct(
            exposureInfoList=exposureInfoList,
        )

    def _isBadExposure(self, expInfo):
        """Return True of exposure does not meet quality criteria

        @param[in] expInfo: exposure info (an ExposureInfo)
        @return True if exposure does not meet quality criteria
        """
        return (expInfo.quality < self.config.quality) \
            or (self.config.cullBlacklisted and expInfo.isBlacklisted) \
            or ((self.config.maxFwhm is not None) and (expInfo.fwhm > self.config.maxFwhm)) \
            or ((self.config.maxSky is not None) and (expInfo.sky > self.config.maxSky)) \
            or ((self.config.maxAirmass is not None) and (expInfo.airmass > self.config.maxAirmass))

    def _runArgDictFromDataId(self, dataId):
        """Extract keyword arguments for run (other than coordList) from a data ID

        @param[in] dataId: a data ID dict
        @return keyword arguments for run (other than coordList), as a dict
        """
        patch = dataId["patch"]
        if self.config.strip.lower() == 'both':
            stripVal = None
        elif self.config.strip.lower() == 'auto':
            stripVal = 'S' if int(patch.split(",")[1]) % 2 == 0 else 'N'
        elif self.config.strip.lower() == 'n':
            stripVal = 'N'
        elif self.config.strip.lower() == 's':
            stripVal = 'S'
        else:
            raise RuntimeError("Invalid config.strip %s. Must be one of N, S, Auto or Both" %
                               (self.config.strip))
        return dict(
            filter=dataId["filter"],
            strip=stripVal
        )


def _formatList(valueList):
    """Format a value list as "v0,v1,v2...vlast"
    """
    return ",".join(str(v) for v in valueList)


def _whereDataFromList(name, valueList):
    """Return a where clause and associated value(s)

    For example if valueList has a single value then returns
        "name = %s", valueList[0]
    but if valueList has more values then returns
        "name in %s", valueList

    This function exists because MySQL requires multiple values for "in" clauses.

    @raise RuntimeError if valueList is None or has length 0
    """
    if not valueList:
        raise RuntimeError("%s valueList = %s; must be a list with at least one value" % (name, valueList))

    if len(valueList) == 1:
        return ("%s = %%s" % (name,), valueList[0])
    else:
        return ("%s in %%s" % (name,), tuple(valueList))


def _computeRaRange(coordList, ctrRa=None):
    """Compute RA range from a list of coords

    The angles are wrapped to be near ctrRa (or coordList[0] if refCoord is None).
    If refCoord is in the middle of your RA range then the coordinates in coordList can span
    essentially the entire range of RA with predictable results. But if refCoord is omitted
    or is not the center of your range, then it is safer to restrict your coordinates
    to span a range of ctrRa - pi < coord < ctrRa + pi

    @param[in] coordList: list of afwCoord.Coord
    @param[in] ctrRa: RA of center of range as an afwGeom.Angle; if None then
        coordList[0].toIcrs().getLongitude() is used after being wrapped into the range [-pi, pi)
    @return RA range as (minRA, maxRA), each an ICRS RA as an afwGeom.Angle
    """
    if not coordList:
        raise RuntimeError("coordList contains no elements")
    raList = [coord.toIcrs().getLongitude() for coord in coordList]
    if ctrRa is None:
        ctrRa = raList[0].wrapCtr()
    raList = [ra.wrapNear(ctrRa) for ra in raList]
    raRadList = np.array([ra.asRadians() for ra in raList])
    minAngRad = np.min(raRadList)
    maxAngRad = np.max(raRadList)
    return (minAngRad * afwGeom.radians, maxAngRad * afwGeom.radians)
