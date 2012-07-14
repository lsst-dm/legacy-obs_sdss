#!/usr/bin/env python

import lsst.afw.table as afwTable
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.daf.persistence as dafPersist
from lsst.pipe.base import Struct
from lsst.pipe.tasks.forcedPhot import ReferencesTask, ReferencesConfig
from lsst.pex.config import Field
import collections

__all__ = ["SdssReferencesConfig", "SdssReferencesTask", "PolySdssReferencesTask", "TestSdssReferencesTask"]

RefSource = collections.namedtuple("RefSource", ["ident", "coord", "flux", "fluxErr"])

class SdssReferencesConfig(ReferencesConfig):
    dbName = Field(dtype=str, doc="Name of database") # Note: no default, so must be set by an override
    dbUrl = Field(dtype=str, doc="URL for database (without the trailing database name)",
                  default="mysql://lsst10.ncsa.uiuc.edu:3390/")
    padding = Field(dtype=float, doc="Padding factor for cone search", default=1.1,
                    check=lambda x: x >= 1.0)

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
        fluxKey = schema.addField("refFlux", float, "Flux from database", "counts")
        fluxErrKey = schema.addField("refFlux.err", float, "Flux error from database", "counts")
        
        references = afwTable.SourceCatalog(schema)
        table = references.table
        references.preallocate(len(sourceList))
        for source in sourceList:
            ref = table.makeRecord()
            ref.setId(source.ident)
            ref.setCoord(source.coord)
            ref.set(fluxKey, source.flux)
            ref.set(fluxErrKey, source.fluxErr)
            references.append(ref)

        return references

    def getTables(self):
        """Return a list of which tables to include in the query"""
        return ["Object"]

    def getColumns(self, dataRef, exposure):
        """Return which columns to get in the query

        @param dataRef     Data reference from butler
        @param exposure    Exposure that has been read
        @return List of column names
        """
        filterName = exposure.getFilter().getName()
        return ["objectId", "ra", "decl", filterName + "PsfFlux", filterName + "PsfFluxSigma"]

    def parseRow(self, db):
        """Parse a row returned from the query

        @param db Database handle
        @return RefSource
        """
        ident = db.getColumnByPosInt64(0)
        ra = db.getColumnByPosDouble(1) * afwGeom.degrees
        dec = db.getColumnByPosDouble(2) * afwGeom.degrees
        flux = db.getColumnByPosFloat(3)
        fluxErr = db.getColumnByPosFloat(4)
        return RefSource(ident, afwCoord.IcrsCoord(ra, dec), flux, fluxErr)

    def whereClause(self, base):
        """Update a basic spatial search where clause with other conditions.
        
        @param base (str)  SQL WHERE clause for spatial search
        @return (str)      Modified SQL WHERE clause"""

        # By default, do nothing
        return base

    def getRaDecFromDatabase(self, dataRef, exposure):
        """Get a list of RA, Dec from the database

        @param dataRef     Data reference, which includes the identifiers
        @param exposure    Exposure that has been read
        @return List of RefSources
        """
        dbFullUrl = self.config.dbUrl + self.config.dbName

        wcs = exposure.getWcs()
        width, height = exposure.getWidth(), exposure.getHeight()
        center = wcs.pixelToSky(afwGeom.Point2D(width/2.0, height/2.0))
        llc = wcs.pixelToSky(afwGeom.Point2D(0.0, 0.0))
        radius = center.angularSeparation(llc)

        db = dafPersist.DbStorage()
        db.setPersistLocation(dafPersist.LogicalLocation(dbFullUrl))
        db.startTransaction()
        db.setTableListForQuery(self.getTables())
        columns = self.getColumns(dataRef, exposure)
        for col in columns:
            db.outColumn(col)
        coneSearch = "scisql_s2PtInCircle(ra, decl, %f, %f, %f) = 1" % \
                (center.getLongitude().asDegrees(),
                        center.getLatitude().asDegrees(),
                        radius.asDegrees() * self.config.padding)
        db.setQueryWhere(self.whereClause(coneSearch))
        db.query()

        sourceList = []
        while db.next():
            source = self.parseRow(db)
            sourceList.append(source)

        db.finishQuery()
        db.endTransaction()
        return sourceList


class SdssCoaddReferencesConfig(SdssReferencesConfig):
    """Configuration for references from a coadd table"""
    coaddName = Field(dtype=str, default="goodSeeing", doc="Name of coadd reference")

class SdssCoaddReferencesTask(SdssReferencesTask):
    """Use a coadd table instead of the Object table."""

    ConfigClass = SdssCoaddReferencesConfig
    
    def getTables(self):
        table = self.config.coaddName + "Source"
        table = table[0].upper() + table[1:]
        return [table]

    def getColumns(self, dataRef, exposure):
        ident = self.config.coaddName + "SourceId"
        ident = ident[0].lower() + ident[1:]
        return [ident, "ra", "decl", "psfFlux", "psfFluxSigma"]

    def whereClause(self, base):
        # Use r filter sources only; cheat and don't join to Filter table.
        return base + " AND filterId = 2"

class SdssCoaddFileReferencesTask(SdssReferencesTask):
    """Forsake the database altogether, and use files available on disk.

    This can be useful for testing.
    """
    def getRaDecFromDatabase(self, dataRef, exposure):
        butler = dataRef.getButler()

        sourceList = []

        wcs = exposure.getWcs()
        width, height = exposure.getWidth(), exposure.getHeight()
        pointList = [(0,0), (width, 0), (width, height), (0, height)]
        coordList = [wcs.pixelToSky(afwGeom.Point2D(x, y)) for x, y in pointList]
        skymap = butler.get(self.config.coaddName + "Coadd_skyMap")
        tractPatchList = skymap.findTractPatchList(coordList)

        for tractInfo, patchInfoList in tractPatchList:
            for patchInfo in patchInfoList:
                catalog = butler.get(self.config.coaddName + "Coadd_src",
                                     filter=exposure.getFilter().getName(),
                                     tract=tractInfo.getId(), patch="%d,%d" % patchInfo.getIndex())
                for src in catalog:
                    sourceList.append(RefSource(src.getId(), src.getCoord(),
                                                src.getPsfFlux(), src.getPsfFluxErr()))

        return sourceList                                 


class PolySdssReferencesTask(SdssReferencesTask):
    """Use exposure polygons to get reference objects out of the database.
    
    This may be slightly more efficient since the polygons will describe the
    exposure on the sky much better than a simple cone search, saving us from
    doing more refinement.  However, I (PAP) haven't been able to get it to
    work yet.  The query was originally provided by KTL, but the schema has
    changed since then.
    """
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
        db.setTableListForQuery(["RefObject", "Region"])
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
    """Get SDSS reference objects out of Mario Juric's custom database table.

    This was originally provided for testing the forced photometry functionality,
    and serves as a useful example.  This one-off database schema is slightly
    different from the LSST database, so hence this override.

    Mario writes:

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
