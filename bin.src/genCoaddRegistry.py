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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import glob
from optparse import OptionParser
import os
import re
import shutil
try:
    import sqlite3
except ImportError:
    # try external pysqlite package; deprecated
    import sqlite as sqlite3
import sys
import lsst.afw.image as afwImage
import lsst.skypix as skypix


def process(dirList, inputRegistry, outputRegistry="registry.sqlite3"):
    if os.path.exists(outputRegistry):
        print >>sys.stderr, "Output registry exists; will not overwrite."
        sys.exit(1)
    if inputRegistry is not None:
        if not os.path.exists(inputRegistry):
            print >>sys.stderr, "Input registry does not exist."
            sys.exit(1)
        shutil.copy(inputRegistry, outputRegistry)

    conn = sqlite3.connect(outputRegistry)

    done = {}
    if inputRegistry is None:
        # Create tables in new output registry.
        cmd = """CREATE TABLE raw (id INTEGER PRIMARY KEY AUTOINCREMENT,
            run INT, filter TEXT, camcol INT, field INT)"""
        # cmd += ", unique(run, filter, camcol, field))"
        conn.execute(cmd)
        cmd = "CREATE TABLE raw_skyTile (id INTEGER, skyTile INTEGER)"
        # cmd += ", unique(id, skyTile), foreign key(id) references raw(id))"
        conn.execute(cmd)
    else:
        cmd = """SELECT run || '_B' || filter ||
            '_C' || camcol || '_F' || field FROM raw"""
        for row in conn.execute(cmd):
            done[row[0]] = True

    qsp = skypix.createQuadSpherePixelization()

    try:
        for dir in dirList:
            for filterDir in glob.iglob(os.path.join(dir, "*")):
                processBand(filterDir, conn, done, qsp)
    finally:
        print >>sys.stderr, "Cleaning up..."
        conn.execute("CREATE INDEX ix_skyTile_id ON raw_skyTile (id)")
        conn.execute("CREATE INDEX ix_skyTile_tile ON raw_skyTile (skyTile)")
        conn.commit()
        conn.close()


def processBand(filterDir, conn, done, qsp):
    nProcessed = 0
    nSkipped = 0
    nUnrecognized = 0
    print >>sys.stderr, filterDir, "... started"
    for fits in glob.iglob(
            os.path.join(filterDir, "fpC*_ts_coaddNorm_NN.fit.gz")):
        m = re.search(r'/([ugriz])/fpC-(\d{6})-\1(\d)-(\d{4})_ts_coaddNorm_NN.fit.gz', fits)
        if not m:
            print >>sys.stderr, "Warning: Unrecognized file:", fits
            nUnrecognized += 1
            continue

        (filter, run, camcol, field) = m.groups()
        camcol = int(camcol)
        run = int(run)
        field = int(field)
        key = "%d_B%s_C%d_F%d" % (run, filter, camcol, field)
        if done.has_key(key):
            nSkipped += 1
            continue

        md = afwImage.readMetadata(fits)
        conn.execute("""INSERT INTO raw VALUES
            (NULL, ?, ?, ?, ?)""", (run, filter, camcol, field))

        for row in conn.execute("SELECT last_insert_rowid()"):
            id = row[0]
            break

        wcs = afwImage.makeWcs(md)
        poly = skypix.imageToPolygon(wcs,
                                     md.get("NAXIS1"), md.get("NAXIS2"),
                                     padRad=0.000075)  # about 15 arcsec
        pix = qsp.intersect(poly)
        for skyTileId in pix:
            conn.execute("INSERT INTO raw_skyTile VALUES(?, ?)",
                         (id, skyTileId))

        nProcessed += 1
        if nProcessed % 100 == 0:
            conn.commit()

    conn.commit()
    print >>sys.stderr, filterDir, \
        "... %d processed, %d skipped, %d unrecognized" % \
        (nProcessed, nSkipped, nUnrecognized)

if __name__ == "__main__":
    parser = OptionParser(usage="""%prog [options] DIR ...

DIR should contain a directory per filter containing coadd pieces.""")
    parser.add_option("-i", dest="inputRegistry", help="input registry")
    parser.add_option("-o", dest="outputRegistry", default="registry.sqlite3",
                      help="output registry (default=registry.sqlite3)")
    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.error("Missing directory argument(s)")
    process(args, options.inputRegistry, options.outputRegistry)
