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
import lsst.pex.config as pexConfig
from lsst.afw.coord import IcrsCoord
import lsst.afw.geom as afwGeom
import lsst.pipe.base as pipeBase
import lsst.daf.persistence as dafPersist
from lsst.pipe.tasks.selectImages import BaseSelectImagesTask, BaseExposureInfo

__all__ = ["SelectSdssImagesTask"]

class _BandSpecificConfig(pexConfig.Config):
    """Band-specific selection criteria
    """
    maxFwhm = pexConfig.Field(
        doc = "maximum FWHM (arcsec)",
        dtype = float,
    )
    maxSky = pexConfig.Field(
        doc = "maximum sky level (maggies/arcsec^2)",
        dtype = float,
        optional = True,
    )
    maxAirmass = pexConfig.Field(
        doc = "Maximum airmass",
        dtype = float,
        optional = True,
    )
        

class SelectSdssImagesConfig(BaseSelectImagesTask.ConfigClass):
    """Config for SelectSdssImagesTask
    """
    BAND_CONFIG_DICT = {}
    BAND_CONFIG_DICT['u'] = _BandSpecificConfig 
    BAND_CONFIG_DICT['g'] = _BandSpecificConfig 
    BAND_CONFIG_DICT['r'] = _BandSpecificConfig 
    BAND_CONFIG_DICT['i'] = _BandSpecificConfig 
    BAND_CONFIG_DICT['z'] = _BandSpecificConfig 
    band = pexConfig.ConfigChoiceField(
        doc = "Band-specific selection criteria",
        typemap = BAND_CONFIG_DICT,
    )
    table = pexConfig.Field(
        doc = "Name of database table",
        dtype = str,
        default = "SeasonFieldQuality",
    )
    quality = pexConfig.ChoiceField(
        doc = "SDSS quality flag",
        dtype = int,
        default = 3,
        allowed={1:"All data", 2:"Flagged ACCEPTABLE or GOOD", 3:"Flagged GOOD"},
    )
    cullBlacklisted = pexConfig.Field(
        doc = "Omit blacklisted images? (Some run/field combinations have been blacklisted even though the quality metrics may have missed them.)",
        dtype = bool,
        default = True,
    )
    camcols = pexConfig.ListField(
        doc = "Which camcols to include? If None then include all",
        dtype = int,
        optional = True,
    )

    def setDefaults(self):
        BaseSelectImagesTask.ConfigClass.setDefaults(self)
        self.host = "lsst-db.ncsa.illinois.edu"
        self.port = 3306
        self.database = "krughoff_SDSS_quality_db"
        # These defaults are the mean seeing for the GOOD (quality = 3) fields in stripe 82 (per band of course) 
        self.band['u'].maxFwhm = 1.52
        self.band['g'].maxFwhm = 1.43
        self.band['r'].maxFwhm = 1.31
        self.band['i'].maxFwhm = 1.25
        self.band['z'].maxFwhm = 1.29
        self.band.name = "g" # to allow instantiation; the code retrieves the data by band name


class ExposureInfo(BaseExposureInfo):
    """Data about a selected exposure
    
    Data includes:
    - dataId: data ID of exposure (a dict)
    - coordList: a list of corner coordinates of the exposure (list of IcrsCoord)
    - fwhm: mean FWHM of exposure
    - flags: flags field from Science_Ccd_Exposure table
    """
    def __init__(self, db, band):
        """Set exposure information based on a query result from a db connection
        """
        BaseExposureInfo.__init__(self)
        self.dataId = dict(
           run = db.getColumnByPosInt(self._nextInd),
           rerun = db.getColumnByPosInt(self._nextInd),
           camcol = db.getColumnByPosInt(self._nextInd),
           frame = db.getColumnByPosInt(self._nextInd),
           band = band,
        )
        self.coordList = []
        for i in range(4):
            self.coordList.append(
                IcrsCoord(
                    afwGeom.Angle(db.getColumnByPosDouble(self._nextInd), afwGeom.degrees),
                    afwGeom.Angle(db.getColumnByPosDouble(self._nextInd), afwGeom.degrees),
                )
            )
        self.fwhm = db.getColumnByPosFloat(self._nextInd)
        self.sky = db.getColumnByPosFloat(self._nextInd)
        self.airmass = db.getColumnByPosFloat(self._nextInd)
        self.quality = db.getColumnByPosInt(self._nextInd)
        self.isblacklisted = db.getColumnByPosChar(self._nextInd)

    @staticmethod
    def getColumnNames():
        """Get database columns to retrieve, in a format useful to the database interface
        
        @return database column names as list of strings
        """
        return "run rerun camcol field ra1 dec1 ra2 dec2 ra3 dec3 ra4 dec4".split() + \
            "psfWidth sky airmass quality isblacklisted".split()



class SelectSdssImagesTask(BaseSelectImagesTask):
    """Select SDSS images suitable for coaddition
    """
    ConfigClass = SelectSdssImagesConfig
    _DefaultName = "selectImages"
    
    @pipeBase.timeMethod
    def run(self, coordList, band):
        """Select SDSS images suitable for coaddition in a particular region
        
        @param[in] band: filter band for images (one of "u", "g", "r", "i" or "z")
        @param[in] coordList: list of coordinates defining region of interest
        
        @return a pipeBase Struct containing:
        - dataIdList: a list of data ID dicts
        """
        idList = []
        db = dafPersist.DbStorage()
        loc = dafPersist.LogicalLocation("mysql://lsst-db.ncsa.illinois.edu:3306/%s"%(self.config.database))

        db.setRetrieveLocation(loc) # was setPersistLocation, but K-T suggests setRetrieveLocation
        db.startTransaction()
        db.setTableForQuery(self.config.table)
        for colName in ExposureInfo.getColumnNames():
            db.outColumn(colName)
        wstr = self.getWhereString(coordList = coordList, band = band)
        db.setQueryWhere(wstr)
        db.query()
        exposureInfoList = []
        while db.next():
            exposureInfoList.append(ExposureInfo(db=db, band=band))

        return pipeBase.Struct(
            exposureInfoList = exposureInfoList,
        )

    def getWhereString(self, coordList, band):
        """Construct SQL query
        
        @return SQL query string
        """
        whereList = []
        if coordList is not None:
            skyBox = afwGeom.Box2D()
            for coord in coordList:
                skyBox.include(coord.getPosition(afwGeom.degrees))
            minx = skyBox.getMinX()
            miny = skyBox.getMinY()
            maxx = skyBox.getMaxX()
            maxy = skyBox.getMaxY()
            bboxstr = "POLYGON((%f %f,%f %f,%f %f,%f %f,%f %f))" % \
                (minx, miny, maxx, miny, maxx, maxy, minx, maxy, minx, miny)
            whereList.append("MBRIntersects(GeomFromText('%s'), bbox)" % bboxstr)
        
        whereList.append("filter = %r" % (band,))
        whereList.append("psfWidth < %f" % (self.config.band[band].maxFwhm,))
        
        if self.config.cullBlacklisted:
            whereList.append("isblacklisted is false")

        # It would be simpler to use str(tuple(range(...))), but that gives "(val,)" for a single value
        # and the final comma may not be compatible with SQL queries
        qualityFlagList = range(self.config.quality, 4)
        qualityFlagStr = ",".join(str(v) for v in qualityFlagList)
        whereList.append("quality in (%s)" % (qualityFlagStr,))
        
        if self.config.camcols is not None:
            camColStr = ",".join(str(v) for v in self.config.camcols)
            whereList.append("camcol in (%s)" % (camColStr,))

        if self.config.band[band].maxSky:
            whereList.append("sky < %f" % (self.config.band[band].maxSky,))

        if self.config.band[band].maxAirmass:
            whereList.append("airmass < %f" % (self.config.band[band].maxAirmass,))
        
        return " and ".join(whereList)

    def _runArgDictFromDataId(self, dataId):
        """Extract keyword arguments for run (other than coordList) from a data ID
        
        @return keyword arguments for run (other than coordList), as a dict
        """
        return dict(
            band = dataId["band"]
        )

def _formatList(valueList):
    """Format a value list as "v0,v1,v2...vlast"
    """
    return ",".join(str(v) for v in valueList)
    
