#!/usr/bin/env python

import lsst.afw.table as afwTable
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
from lsst.pipe.tasks.forcedPhot import ForcedPhotTask, ForcedPhotConfig

class SdssForcedPhotTask(ForcedPhotTask):
    def getReferences(self, dataRef, exposure):
        """Get reference sources on (or close to) exposure"""
        raise NotImplementedError("Don't know how to get reference sources for SDSS yet")

        table = afwTable.SimpleTable.makeMinimalSchema()
        references = afwTable.SimpleCatalog(table)

        dbRows = lookupDatabaseSomehow(exposure) # XXX

        references.preallocate(len(dbRows))
        for row in dbRows:
            ref = table.makeRecord()
            ref.setCoord(afwCoord.IcrsCoord(row.ra * afwGeom.degrees, row.dec * afwGeom.degrees))
            references.append(ref)

        return references
