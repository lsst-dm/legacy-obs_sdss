#!/usr/bin/env python

import lsst.afw.table as afwTable
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.dadf.persistence as dafPersist
from lsst.pipe.tasks.forcedPhot import ForcedPhotTask, ForcedPhotConfig
from lsst.pex.config import Field

class SdssForcedPhotConfig(ForcedPhotConfig):
    dbName = Field(dtype=str, doc="Name of database") # Note: no default, so must be set by an override
    dbUrl = Field(dtype=str, doc="URL for database (without the trailing database name)",
                  default="mysql://lsst10.ncsa.uiuc.edu:3390/")


class SdssForcedPhotTask(ForcedPhotTask):
    def getReferences(self, dataRef, exposure):
        """Get reference sources on (or close to) exposure"""
        coordList = self.getRaDecFromDatabase(dataRef)

        table = afwTable.SimpleTable.makeMinimalSchema()
        references = afwTable.SimpleCatalog(table)
        references.preallocate(len(dbRows))
        for coord in coordList:
            ref = table.makeRecord()
            ref.setCoord(coord)
            references.append(ref)

        return references

    def getRaDecFromDatabase(self, dataRef):
        """Get a list of RA, Dec from the database

        @param dataRef     Data reference, which includes the identifiers
        @return List of coordinates
        """
        dbFullUrl = self.config.dbUrl + self.config.dbName

        db = dafPersist.DbStorage()
        db.setPersistLocation(dafPersist.LogicalLocation(dbFullUrl))
        db.startTransaction()
        db.executeSql("""
           SELECT poly FROM Science_Ccd_Exposure
               WHERE scienceCcdExposureId = %d
               INTO @poly;""" % dataRef.get("ccdExposureId"))
        db.executeSql("CALL scisql.scisql_s2CPolyRegion(@poly, 20)")
        db.setTableListForQuery(["Object", "Region"])
        db.outColumn("ra")
        db.outColumn("dec")
        db.setQueryWhere("""
           Object.htmId20 BETWEEN Region.htmMin AND Region.htmMax
           AND scisql_s2PtInCPoly(Object.ra_PS, Object.decl_PS, @poly) = 1""")
        db.query()

        coordList = []
        while db.next():
            ra = db.getColumnByPosDouble(0) * afwGeom.degrees
            dec = db.getColumnByPosDouble(1) * afwGeom.degrees
            coordList.append(afwCoord.IcrsCoord(ra, dec))

        db.finishQuery()
        db.endTransaction()
        return coordList
