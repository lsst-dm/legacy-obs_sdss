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
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.pipe.base as pipeBase
import lsst.daf.persistence as dafPersist

__all__ = ["SelectSDSSImagesTask"]

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
        

class SelectSDSSImagesConfig(pexConfig.Config):
    """Config for SelectSDSSImagesTask
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
    database = pexConfig.Field(
        doc = "Name of database",
        dtype = str,
        default = "krughoff_SDSS_quality_db",
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
    def setDefaults(self):
        # These defaults are the mean seeing for the GOOD (quality = 3) fields in stripe 82 (per band of course) 
        self.band['u'].maxFwhm = 1.52
        self.band['g'].maxFwhm = 1.43
        self.band['r'].maxFwhm = 1.31
        self.band['i'].maxFwhm = 1.25
        self.band['z'].maxFwhm = 1.29



class SelectSDSSImagesTask(pipeBase.Task):
    """Select SDSS images suitable for coaddition
    """
    ConfigClass = SelectSDSSImagesConfig
    _DefaultName = "selectImages"
    
    @pipeBase.timeMethod
    def run(self, band, coordList):
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
        db.outColumn("run")
        db.outColumn("rerun")
        db.outColumn("camcol")
        db.outColumn("field")
        wstr = self.getWhereString(band, coordList)
        db.setQueryWhere(wstr)
        db.query()
        while db.next():
            run = db.getColumnByPosInt(0)
            rerun = db.getColumnByPosInt(1)
            camcol = db.getColumnByPosInt(2)
            field = db.getColumnByPosInt(3)
            dataId = {'run':run, 'rerun':rerun, 'camcol':camcol, 'frame':field, 'band':band}
            idList.append(dataId)

        return pipeBase.Struct(
            dataIdList = idList,
        )
    
    def runDataRef(self, dataRef, coordList):
        """Run based on a data reference
        
        @param[in] dataRef: data reference; must contain key "band"
        @param[in] coordList: list of coordinates defining region of interest
        @return a pipeBase Struct containing:
        - dataRefList: a list of data references
        """
        butler = dataRef.butlerSubset.butler
        band = dataRef["band"]
        dataIdList = self.run(band, coordList).dataIddList
        dataRefList = [butler.dataRef(dataId=dataId, level="sensor") for dataId in dataIdList]
        return pipeBase.Struct(
            dataRefList = dataRefList,
        )

    def getWhereString(self, band, coordList):
        """Construct SQL query
        
        @return SQL query string
        """
        skyBox = afwGeom.Box2D()
        for coord in coordList:
            skyBox.include(coord.getPosition(afwGeom.degrees))
        minx = skyBox.getMinX()
        miny = skyBox.getMinY()
        maxx = skyBox.getMaxX()
        maxy = skyBox.getMaxY()
        bboxstr = "POLYGON((%f %f,%f %f,%f %f,%f %f,%f %f))" % \
            (minx, miny, maxx, miny, maxx, maxy, minx, maxy, minx, miny)
        whereTemplate = "MBRIntersects(GeomFromText('%s'), bbox) and filter = '%s' and psfWidth < %f" % \
            (bboxstr, band, self.config.band[band].maxFwhm)
        
        if self.config.cullBlacklisted:
            whereTemplate += " and isblacklisted is false"

        if self.config.quality == 1:
            whereTemplate += " and quality in (1,2,3)"
        elif self.config.quality == 2:
            whereTemplate += " and quality in (2,3)"
        else:
            whereTemplate += " and quality = 3"

        if self.config.band[band].maxSky:
            whereTemplate += " and sky < %f"%(self.config.band[band].maxSky)

        if self.config.band[band].maxAirmass:
            whereTemplate += " and airmass < %f"%(self.config.band[band].maxAirmass)
        
        return whereTemplate
