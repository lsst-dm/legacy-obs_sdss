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

        schema = afwTable.SimpleTable.makeMinimalSchema()
        references = afwTable.SimpleCatalog(schema)
        table = references.table
        references.preallocate(len(dbRows))
        for ident, coord in coordList:
            ref = table.makeRecord()
            ref.setId(ident)
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
        db.outColumn("objectId")
        db.outColumn("ra")
        db.outColumn("decl")
        db.setQueryWhere("""
           Object.htmId20 BETWEEN Region.htmMin AND Region.htmMax
           AND scisql_s2PtInCPoly(Object.ra_PS, Object.decl_PS, @poly) = 1""")
        db.query()

        coordList = []
        while db.next():
            ident = db.getColumnByPosInt64(0)
            ra = db.getColumnByPosDouble(1) * afwGeom.degrees
            dec = db.getColumnByPosDouble(2) * afwGeom.degrees
            row = (ident, afwCoord.IcrsCoord(ra, dec))
            coordList.append(row)

        db.finishQuery()
        db.endTransaction()
        return coordList


class TestSdssForcedPhotTask(SdssForcedPhotTask):
    """Get SDSS reference objects out of Mario Juric's
    custom database table.  Mario writes:

    {{{
    I just imported the DR7 Stripe 82 co-add catalog object table to table
    'Stripe82RefObject' in database 'juric_DR7_stripe82'. This is to support
    the work on forced photometry, until our coadd/detection code matures.

    It has the following columns:

        sdssObjectId,run,rerun,camcol,field,obj,mode,type,ra,decl,
        uMag,uErr,gMag,gErr,rMag,rErr,iMag,iErr,zMag,zErr
        htmId20,isStar

    where the magnitudes are SDSS model magnitude. The data has been
    downloaded from SDSS DR7 CAS website, Stripe82 database using:

        SELECT objID, run, rerun, camcol, field, obj, mode, type, ra, dec, u,
            g, r, i, z, err_u, err_g, err_r, err_i, err_z from PhotoObjAll where run
            in (106, 206)
    }}}

    Since we can, we will retrieve these by run,camcol,field instead of
    using a cone search.
    """
    
    def getRaDecFromDatabase(self, dataRef):
        """Get a list of RA, Dec from the database

        @param dataRef     Data reference, which includes the identifiers
        @return List of coordinates
        """
        dataId = dataRef.dataId

        dbFullUrl = self.config.dbUrl + self.config.dbName
        db = dafPersist.DbStorage()
        db.setPersistLocation(dafPersist.LogicalLocation(dbFullUrl))
        db.startTransaction()
        db.setTableListForQuery(["Stripe82RefObject"])
        db.outColumn("sdssObjectId")
        db.outColumn("ra")
        db.outColumn("decl")
        db.setQueryWhere("run = %(run)d AND camcol = %(camcol)d AND field = %(frame)d" % dataId)
        db.query()

        coordList = []
        while db.next():
            ident = db.getColumnByPosInt64(0)
            ra = db.getColumnByPosDouble(1) * afwGeom.degrees
            dec = db.getColumnByPosDouble(2) * afwGeom.degrees
            row = (ident, afwCoord.IcrsCoord(ra, dec))
            coordList.append(row)

        db.finishQuery()
        db.endTransaction()
        return coordList
