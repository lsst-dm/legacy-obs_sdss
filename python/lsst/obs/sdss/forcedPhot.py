#!/usr/bin/env python

import lsst.afw.table as afwTable
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.daf.persistence as dafPersist
from lsst.pipe.tasks.forcedPhot import ReferencesTask, ReferencesConfig
from lsst.pex.config import Field
import collections

RefSource = collections.namedtuple("RefSource", ["ident", "coord", "mag", "magErr"])


class SdssReferencesConfig(ReferencesConfig):
    dbName = Field(dtype=str, doc="Name of database") # Note: no default, so must be set by an override
    dbUrl = Field(dtype=str, doc="URL for database (without the trailing database name)",
                  default="mysql://lsst10.ncsa.uiuc.edu:3390/")


class SdssReferencesTask(ReferencesTask):
    ConfigClass = SdssReferencesConfig
    
    def getReferences(self, dataRef, exposure):
        """Get reference sources on (or close to) exposure

        @param dataRef     Data reference from butler
        @param exposure    Exposure that has been read
        @return Catalog (lsst.afw.table.SourceCatalog) of reference sources
        """
        sourceList = self.getRaDecFromDatabase(dataRef, exposure)

        schema = afwTable.SourceTable.makeMinimalSchema()
        magKey = schema.addField("refMag", float, "Magnitude from database", "mag")
        magErrKey = schema.addField("refMag.err", float, "Magnitude error from database", "mag")
        
        references = afwTable.SourceCatalog(schema)
        table = references.table
        references.preallocate(len(sourceList))
        for source in sourceList:
            ref = table.makeRecord()
            ref.setId(source.ident)
            ref.setCoord(source.coord)
            ref.set(magKey, source.mag)
            ref.set(magErrKey, source.magErr)
            references.append(ref)

        return references

    def getRaDecFromDatabase(self, dataRef, exposure):
        """Get a list of RA, Dec from the database

        @param dataRef     Data reference, which includes the identifiers
        @return List of RefSources
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
        # XXX Get magnitude, error
        db.setQueryWhere("""
           Object.htmId20 BETWEEN Region.htmMin AND Region.htmMax
           AND scisql_s2PtInCPoly(Object.ra_PS, Object.decl_PS, @poly) = 1""")
        db.query()

        sourceList = []
        while db.next():
            ident = db.getColumnByPosInt64(0)
            ra = db.getColumnByPosDouble(1) * afwGeom.degrees
            dec = db.getColumnByPosDouble(2) * afwGeom.degrees
            mag = 0.0
            magErr = 0.0
            sourceList.append(RefSource(ident, afwCoord.IcrsCoord(ra, dec), mag, magErr))

        db.finishQuery()
        db.endTransaction()
        return sourceList


class TestSdssReferencesTask(SdssReferencesTask):
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

    The Science_Ccd_Exposure table isn't populated, so we can't use the
    above query; we'll just do a simple cone search.
    """
    
    def getRaDecFromDatabase(self, dataRef, exposure):
        """Get a list of RA, Dec from the database

        @param dataRef     Data reference, which includes the identifiers
        @return List of RefSources
        """
        dataId = dataRef.dataId

        padding = 1.1 # Padding factor; could live in a Config except this is just a temporary test...

        filtName = exposure.getFilter().getName()

        wcs = exposure.getWcs()
        width, height = exposure.getWidth(), exposure.getHeight()
        center = wcs.pixelToSky(afwGeom.Point2D(width/2.0, height/2.0))
        radius = center.angularSeparation(wcs.pixelToSky(afwGeom.Point2D(0.0, 0.0)))

        dbFullUrl = self.config.dbUrl + self.config.dbName
        db = dafPersist.DbStorage()
        db.setPersistLocation(dafPersist.LogicalLocation(dbFullUrl))
        db.startTransaction()
        db.setTableListForQuery(["Stripe82RefObject"])
        db.outColumn("sdssObjectId")
        db.outColumn("ra")
        db.outColumn("decl")
        db.outColumn(filtName + "Mag")
        db.outColumn(filtName + "Err")
        db.setQueryWhere("scisql_s2PtInCircle(ra, decl, %f, %f, %f) = 1" %
                         (center.getLongitude().asDegrees(), center.getLatitude().asDegrees(),
                          radius.asDegrees() * padding))
        db.query()

        sourceList = []
        while db.next():
            ident = db.getColumnByPosInt64(0)
            ra = db.getColumnByPosDouble(1) * afwGeom.degrees
            dec = db.getColumnByPosDouble(2) * afwGeom.degrees
            mag = db.getColumnByPosFloat(3)
            magErr = db.getColumnByPosFloat(4)
            sourceList.append(RefSource(ident, afwCoord.IcrsCoord(ra, dec), mag, magErr))

        db.finishQuery()
        db.endTransaction()
        return sourceList
