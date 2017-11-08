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

from __future__ import print_function
import glob
from optparse import OptionParser
import os
import re
import shutil
import sqlite3
import sys
import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.skypix as skypix


def process(dirList, inputRegistry, outputRegistry="registry.sqlite3"):
    if os.path.exists(outputRegistry):
        print("Output registry exists; will not overwrite.", file=sys.stderr)
        sys.exit(1)
    if inputRegistry is not None:
        if not os.path.exists(inputRegistry):
            print("Input registry does not exist.", file=sys.stderr)
            sys.exit(1)
        shutil.copy(inputRegistry, outputRegistry)

    conn = sqlite3.connect(outputRegistry)

    done = {}
    if inputRegistry is None:
        # Create tables in new output registry.
        cmd = """CREATE TABLE raw (id INTEGER PRIMARY KEY AUTOINCREMENT,
            run INT, rerun INT, filter TEXT, camcol INT, field INT,
            taiObs TEXT, strip TEXT)"""
        # cmd += ", unique(run, filter, camcol, field))"
        conn.execute(cmd)
        cmd = "CREATE TABLE raw_skyTile (id INTEGER, skyTile INTEGER)"
        # cmd += ", unique(id, skyTile), foreign key(id) references raw(id))"
        conn.execute(cmd)
    else:
        cmd = """SELECT run || '_R' || rerun || '_B' || filter ||
            '_C' || camcol || '_F' || field FROM raw"""
        for row in conn.execute(cmd):
            done[row[0]] = True

    qsp = skypix.createQuadSpherePixelization()

    try:
        for dir in dirList:
            if dir.endswith("runs"):
                for runDir in glob.iglob(os.path.join(dir, "*")):
                    processRun(runDir, conn, done, qsp)
            else:
                processRun(dir, conn, done, qsp)
    finally:
        print("Cleaning up...", file=sys.stderr)
        conn.execute("""CREATE UNIQUE INDEX uq_raw ON raw
                (run, filter, camcol, field)""")
        conn.execute("CREATE INDEX ix_skyTile_id ON raw_skyTile (id)")
        conn.execute("CREATE INDEX ix_skyTile_tile ON raw_skyTile (skyTile)")
        conn.commit()
        conn.close()


def processRun(runDir, conn, done, qsp):
    nProcessed = 0
    nSkipped = 0
    nUnrecognized = 0
    print(runDir, "... started", file=sys.stderr)
    for fits in glob.iglob(
            os.path.join(runDir, "*", "corr", "[1-6]", "fpC*.fit.gz")):
        m = re.search(r'(\d+)/corr/([1-6])/fpC-(\d{6})-([ugriz])\2-(\d{4}).fit.gz', fits)
        if not m:
            print("Warning: Unrecognized file:", fits, file=sys.stderr)
            nUnrecognized += 1
            continue

        (rerun, camcol, run, filter, field) = m.groups()
        rerun = int(rerun)
        camcol = int(camcol)
        run = int(run)
        field = int(field)
        key = "%d_R%d_B%s_C%d_F%d" % (run, rerun, filter, camcol, field)
        if key in done or rerun < 40:
            nSkipped += 1
            continue

        md = afwImage.readMetadata(fits)
        date = md.get("DATE-OBS")
        if date.find("-") != -1:
            (year, month, day) = md.get("DATE-OBS").split("-")
        else:
            (day, month, year) = md.get("DATE-OBS").split("/")
            year = 1900 + int(year)
        (hour, minute, second) = md.get("TAIHMS").split(":")
        seconds = float(second)
        second = int(seconds)
        taiObs = dafBase.DateTime(int(year), int(month), int(day), int(hour),
                                  int(minute), second, dafBase.DateTime.TAI)
        taiObs = dafBase.DateTime(taiObs.nsecs() +
                                  int((seconds - second) * 1000000000), dafBase.DateTime.TAI)
        taiObs = taiObs.toString(dafBase.DateTime.UTC)[:-1]
        strip = "%d%s" % (md.get('STRIPE'), md.get('STRIP'))
        conn.execute("""INSERT INTO raw VALUES
            (NULL, ?, ?, ?, ?, ?, ?, ?)""",
                     (run, rerun, filter, camcol, field, taiObs, strip))

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
    print(runDir,
          "... %d processed, %d skipped, %d unrecognized" %
          (nProcessed, nSkipped, nUnrecognized), file=sys.stderr)


if __name__ == "__main__":
    parser = OptionParser(usage="""%prog [options] DIR ...

DIR may be either a root directory containing a 'raw' subdirectory
or a visit subdirectory.""")
    parser.add_option("-i", dest="inputRegistry", help="input registry")
    parser.add_option("-o", dest="outputRegistry", default="registry.sqlite3",
                      help="output registry (default=registry.sqlite3)")
    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.error("Missing directory argument(s)")
    process(args, options.inputRegistry, options.outputRegistry)
