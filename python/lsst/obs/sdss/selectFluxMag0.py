from builtins import str
from builtins import range
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

from lsst.afw.coord import IcrsCoord
import lsst.afw.geom as afwGeom
from lsst.daf.persistence import DbAuth
import lsst.pipe.base as pipeBase
from lsst.pipe.tasks.selectImages import DatabaseSelectImagesConfig, BaseExposureInfo

__all__ = ["SelectSdssFluxMag0Task"]


class SelectSdssFluxMag0Config(DatabaseSelectImagesConfig):
    """Config for SelectSdssFluxMag0Task
    """

    def setDefaults(self):
        DatabaseSelectImagesConfig.setDefaults(self)
        self.database = "krughoff_early_prod_11032012"
        self.host = "lsst-db.ncsa.illinois.edu"
        self.port = 3306


class FluxMagInfo(BaseExposureInfo):
    """Data about a selected exposure

    Data includes:
    - dataId: data ID of exposure (a dict)
    - coordList: a list of corner coordinates of the exposure (list of IcrsCoord)
    - fluxMag0: float
    - fluxMag0Sigma: float
    """

    def __init__(self, result):
        """Set exposure information based on a query result from a db connection
        """
        result = list(result)
        dataId = dict(
            run=result.pop(0),
            camcol=result.pop(0),
            field=result.pop(0),
            filter=result.pop(0),
        )
        coordList = [IcrsCoord(afwGeom.Angle(result.pop(0), afwGeom.degrees),
                               afwGeom.Angle(result.pop(0), afwGeom.degrees)) for i in range(4)]

        BaseExposureInfo.__init__(self, dataId=dataId, coordList=coordList)
        self.fluxMag0 = result.pop(0)
        self.fluxMag0Sigma = result.pop(0)

    @staticmethod
    def getColumnNames():
        """Get database columns to retrieve, in a format useful to the database interface

        @return database column names as list of strings
        """
        return (
            "run camcol field  filterName".split() +
            "corner1Ra corner1Decl corner2Ra corner2Decl".split() +
            "corner3Ra corner3Decl corner4Ra corner4Decl".split() +
            "fluxMag0 fluxMag0Sigma".split()
        )


class SelectSdssFluxMag0Task(pipeBase.Task):
    """Select SDSS data suitable for computing fluxMag0
    """
    ConfigClass = SelectSdssFluxMag0Config
    _DefaultName = "selectFluxMag0"

    @pipeBase.timeMethod
    def run(self, dataId, coordList):
        """Select flugMag0's of SDSS images for a particular run

        @param[in] dataId: a dataId containing at least a run and filter
        @param[in] coordList: list of coordinates defining region of interest

        @return a pipeBase Struct containing:
        - fluxMagInfoList: a list of FluxMagInfo objects
        """
        argDict = self.runArgDictFromDataId(dataId)
        run = argDict["run"]
        filter = argDict["filter"]

        if filter not in set(("u", "g", "r", "i", "z")):
            raise RuntimeError("filter=%r is an invalid name" % (filter,))

        filterDict = {"u": 0,
                      "g": 1,
                      "r": 2,
                      "i": 3,
                      "z": 4}

        if self._display:
            self.log.info(self.config.database)

        db = MySQLdb.connect(
            host=self.config.host,
            port=self.config.port,
            db=self.config.database,
            user=DbAuth.username(self.config.host, str(self.config.port)),
            passwd=DbAuth.password(self.config.host, str(self.config.port)),
        )
        cursor = db.cursor()

        columnNames = tuple(FluxMagInfo.getColumnNames())

        queryStr = "select %s from Science_Ccd_Exposure where " % (", ".join(columnNames))
        dataTuple = ()

        if coordList is not None:
            # look for exposures that overlap the specified region
            for c in coordList:
                dataTuple += (c.getLongitude().asDegrees(), c.getLatitude().asDegrees())
            queryStr += " scisql_s2PtInCPoly(ra, decl"
            queryStr += ", %s, %s" * len(coordList)
            queryStr += ") = 1 and "

        # compute where clauses as a list of (clause, data)
        whereDataList = [
            ("filterId = %s", filterDict[filter]),
            ("run = %s", run),
        ]

        queryStr += " and ".join(wd[0] for wd in whereDataList)
        dataTuple += tuple(wd[1] for wd in whereDataList)

        queryStr += " order by field desc"
        if self._display:
            self.log.info("queryStr=%r; dataTuple=%s" % (queryStr, dataTuple))

        cursor.execute(queryStr, dataTuple)
        exposureInfoList = [FluxMagInfo(result) for result in cursor]
        if self._display:
            self.log.info("Found %d exposures" %
                          (len(exposureInfoList)))

        return pipeBase.Struct(
            fluxMagInfoList=exposureInfoList,
        )

    def runArgDictFromDataId(self, dataId):
        """Extract keyword arguments for run (other than coordList) from a data ID

        @param[in] dataId: a data ID dict
        @return keyword arguments (other than coordList), as a dict
        """
        return dict(
            filter=dataId["filter"],
            run=dataId["run"]
        )
