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

import lsst.pex.config as pexConfig
from lsst.afw.coord import IcrsCoord
import lsst.afw.geom as afwGeom
from lsst.daf.persistence import DbAuth
import lsst.pipe.base as pipeBase
from lsst.pipe.tasks.selectImages import SelectImagesConfig, BaseExposureInfo

__all__ = ["SelectSdssFluxMag0Task"]

class SelectSdssFluxMag0Config(SelectImagesConfig):
    table = pexConfig.Field(
        doc = "Name of database table",
        dtype = str,
        default = "Science_Ccd_Exposure",
    )

    def setDefaults(self):
        SelectImagesConfig.setDefaults(self)
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
        BaseExposureInfo.__init__(self)
        self.dataId = dict(
           run = result[self._nextInd],
           camcol = result[self._nextInd],
           field = result[self._nextInd],
           filter = result[self._nextInd],
        )
        self.coordList = []
        for i in range(4):
            self.coordList.append(
                IcrsCoord(
                    afwGeom.Angle(result[self._nextInd], afwGeom.degrees),
                    afwGeom.Angle(result[self._nextInd], afwGeom.degrees),
                )
            )
        self.fluxMag0 = result[self._nextInd]
        self.fluxMag0Sigma = result[self._nextInd]
                
    @staticmethod
    def getColumnNames():
        """Get database columns to retrieve, in a format useful to the database interface
        
        @return database column names as list of strings
        """
        return (
            "run camcol field  filterName".split() + \
            "corner1Ra corner1Decl corner2Ra corner2Decl corner3Ra corner3Decl corner4Ra corner4Decl".split() + \
            "fluxMag0 fluxMag0Sigma".split()
        )

class SelectSdssFluxMag0Task(pipeBase.Task):
    """Select SDSS data suitable for computing fluxMag0
    """
    ConfigClass = SelectSdssFluxMag0Config

    @pipeBase.timeMethod
    def run(self, filter, run, camcol):
        """Select flugMag0's of SDSS images for a particular run

        @param[in] filter: filter for images (one of "u", "g", "r", "i" or "z")
        @param[in] coordList: list of coordinates defining region of interest
        
        @return a pipeBase Struct containing:
        - fluxMagInfoList: a list of FluxMagInfo objects
        """
        if filter not in set(("u", "g", "r", "i", "z")):
            raise RuntimeError("filter=%r is an invalid name" % (filter,))

        filterDict = {"u": 0,
                      "g": 1,
                      "r": 2,
                      "i": 3,
                      "z": 4}
                      

        read_default_file=os.path.expanduser("~/.my.cnf")

        try:
            open(read_default_file)
            kwargs = dict(
                read_default_file=read_default_file,
                )
        except IOError:
            kwargs = dict(
                user = DbAuth.username(self.config.host, str(self.config.port)),
                passwd = DbAuth.password(self.config.host, str(self.config.port)),
                )

        self.log.info(self.config.table)    
        self.log.info(self.config.database)
        
        db = MySQLdb.connect(
            host = self.config.host,
            port = self.config.port,
            db = self.config.database,
            **kwargs
        )
        cursor = db.cursor()
        
        columnNames = tuple(FluxMagInfo.getColumnNames())
        if not columnNames:
            raise RuntimeError("Bug: no column names")
        queryStr = "select %s from %s where " % (", ".join(columnNames), self.config.table)
        dataTuple = () # tuple(columnNames)

        # compute where clauses as a list of (clause, data)
        whereDataList = [
            ("filterId = %s", filterDict[filter]),
            ("run = %s", run),
            ("camcol = %s", camcol),
        ]
        
        queryStr += " and ".join(wd[0] for wd in whereDataList)
        dataTuple += tuple(wd[1] for wd in whereDataList)
        
        queryStr += " order by field"

        self.log.info("queryStr=%r; dataTuple=%s" % (queryStr, dataTuple))
        
        cursor.execute(queryStr, dataTuple)
        exposureInfoList = [FluxMagInfo(result) for result in cursor]        
       
        self.log.info("Found %d exposures" % \
            (len(exposureInfoList)))
        
        return pipeBase.Struct(
            fluxMagInfoList = exposureInfoList,
        )
