#!/usr/bin/env python
import os

import MySQLdb

import lsst.afw.table as afwTable
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
from lsst.daf.persistence import DbAuth, DbStorage, LogicalLocation
from lsst.pipe.base import Struct
from lsst.pipe.tasks.forcedPhot import ReferencesTask, ReferencesConfig
from lsst.pex.config import Field
import collections

__all__ = ["SdssReferencesConfig", "SdssReferencesTask", "PolySdssReferencesTask", "TestSdssReferencesTask"]

RefSource = collections.namedtuple("RefSource", ["ident", "coord", "flux", "fluxErr"])

_FilterIdDict = dict(u = 0, g = 1, r = 2, i = 3, z = 4)

class SdssReferencesConfig(ReferencesConfig):
    host = Field(
        doc = "Database server host name",
        dtype = str,
        default = "lsst10.ncsa.uiuc.edu",
    )
    port = Field(
        doc = "Database server port",
        dtype = int,
        default = 3306,
    )
    dbName = Field(
        doc = "Name of database",
        dtype = str,
    )
    filterName = Field(
        dtype=str, 
        default=None, 
        doc="Name of the band to use for the references, if None use the band in the exposure",
    )

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

    def getTableName(self):
        """Return the table to include in the query"""
        return "Object"

    def getColumns(self, dataRef, exposure):
        """Return which columns to get in the query

        @param dataRef     Data reference from butler
        @param exposure    Exposure that has been read
        @return List of column names
        """
        filterName = exposure.getFilter().getName()
        return ["objectId", "ra", "decl", filterName + "PsfFlux", filterName + "PsfFluxSigma"]

    def parseRow(self, result):
        """Parse a result returned from the query

        @param result   result returned from MySQLDb query
        @return RefSource
        """
        ident = result[0]
        ra = afwGeom.Angle(result[1], afwGeom.degrees)
        dec = afwGeom.Angle(result[2], afwGeom.degrees)
        flux = result[3]
        fluxErr = result[4]
        return RefSource(ident, afwCoord.IcrsCoord(ra, dec), flux, fluxErr)

    def getRaDecFromDatabase(self, dataRef, exposure):
        """Get a list of RA, Dec from the database

        @param dataRef     Data reference, which includes the identifiers
        @param exposure    Exposure that has been read
        @return List of RefSources
        """
        # this ~/.my.cnf stuff is a Princetonism; LSST uses DbAuth
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

        db = MySQLdb.connect(
            host = self.config.host,
            port = self.config.port,
            db = self.config.dbName,
            **kwargs
        )
        cursor = db.cursor()

        wcs = exposure.getWcs()
        posBBox = afwGeom.Box2D(exposure.getBBox(afwImage.PARENT))
        coordList = [wcs.pixelToSky(pos) for pos in posBBox.getCorners()]
        coordStrList = ["%s, %s" % (c.getLongitude().asDegrees(),
                                    c.getLatitude().asDegrees()) for c in coordList]
        coordStr = ", ".join(coordStrList)
        coordCmd = "call scisql.scisql_s2CPolyRegion(scisql_s2CPolyToBin(%s), 20)" % (coordStr,)
        self.log.info("queryStr=%r" % (coordCmd,))
        cursor.execute(coordCmd)
        cursor.nextset() # ignore one-line result of coordCmd

        columnNames = self.getColumns(dataRef, exposure)
        colNameStr = ", ".join("s.%s" % (cn,) for cn in columnNames)

        tableName = self.getTableName()

        queryStr = """
SELECT
    %s
FROM
    %s AS s INNER JOIN
    scisql.Region AS r ON (s.htmId20 BETWEEN r.htmMin AND r.htmMax)
WHERE
    filterId = %%s
""" % (colNameStr, tableName)
        if not self.config.filterName:
            filterId = _FilterIdDict[exposure.getFilter().getName()]
        else:
            filterId = _FilterIdDict[self.config.filterName]
        dataTuple = (filterId,)

        self.log.info("queryStr=%r; dataTuple=%s" % (queryStr, dataTuple))

        cursor.execute(queryStr, dataTuple)

        sourceList = []
        for result in cursor:
            sourceList.append(self.parseRow(result))
        
        self.log.info("Found %s sources" % (len(sourceList),))

        return sourceList


class SdssCoaddReferencesConfig(SdssReferencesConfig):
    """Configuration for references from a coadd table"""
    coaddName = Field(dtype=str, default="goodSeeing", doc="Name of coadd reference")


class SdssCoaddReferencesTask(SdssReferencesTask):
    """Use a coadd table instead of the Object table."""

    ConfigClass = SdssCoaddReferencesConfig
    
    def getTableName(self):
        table = self.config.coaddName + "Source"
        table = table[0].upper() + table[1:]
        return table

    def getColumns(self, dataRef, exposure):
        ident = self.config.coaddName + "SourceId"
        ident = ident[0].lower() + ident[1:]
        return [ident, "ra", "decl", "psfFlux", "psfFluxSigma"]


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
