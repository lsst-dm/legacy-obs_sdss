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

import sys
import os
import re
from astropy.io import fits

import lsst.afw.image as afwImage
import lsst.geom as geom


class Span(object):

    def __init__(self, y, x1, x2):
        self.y = y
        self.x1 = x1
        self.x2 = x2


class Objmask(object):
    nperspan = 6

    def __init__(self, frow, cval, verbose=False):
        self.refcntr = frow[0]
        self.nspan = frow[1]
        self.row0 = frow[2]
        self.col0 = frow[3]
        self.rmin = frow[4]
        self.rmax = frow[5]
        self.cmin = frow[6]
        self.cmax = frow[7]
        self.npix = frow[8]
        self.span = frow[9]
        if len(self.span) == 0:
            self.nspan = 0  # some bogus fpM files

        self.spans = []
        npixcheck = 0
        for i in range(self.nspan):
            b1 = self.span[6*i + 0]
            b2 = self.span[6*i + 1]
            y = (b1 << 8) + b2
            b1 = self.span[6*i + 2]
            b2 = self.span[6*i + 3]
            x1 = (b1 << 8) + b2
            b1 = self.span[6*i + 4]
            b2 = self.span[6*i + 5]
            x2 = (b1 << 8) + b2
            self.spans.append(Span(y, x1, x2))
            for i in range(x1, x2+1):
                npixcheck += 1

        # Some fpM files are actually wrong and this test fails!
        # So warn, not assert
        # 5759/40/objcs/1/fpM-005759-r1-0011.fit
        # Plane S_MASK_NOTCHECKED
        if self.npix != npixcheck and verbose:
            print("WARNING: npix != npixcheck (%d != %d)" % (self.npix, npixcheck))

        self.cval = cval

    def setMask(self, mask):
        nrow = mask.getHeight()
        ncol = mask.getWidth()

        for i in range(self.nspan):
            y = int(self.spans[i].y - self.row0)

            if (y < 0 or y >= nrow):
                continue

            x1 = self.spans[i].x1 - self.col0
            x2 = self.spans[i].x2 - self.col0

            if (x1 < 0):
                x1 = 0

            if(x2 >= ncol):
                x2 = ncol - 1

            mask.array[y, x1: x2 + 1] |= self.cval


def convertfpM(infile, allPlanes=False):
    hdr = fits.open(infile)
    hdr[0].header['RUN']
    hdr[0].header['CAMCOL']
    hdr[0].header['FIELD']
    nRows = hdr[0].header['MASKROWS']
    nCols = hdr[0].header['MASKCOLS']
    hdr[0].header['NPLANE']

    names = hdr[-1].data.names
    if ("attributeName" not in names) or ("Value" not in names):
        raise LookupError("Missing data in fpM header")

    planes = hdr[-1].data.field("attributeName").tolist()
    mask = afwImage.Mask(geom.ExtentI(nCols, nRows))

    # Minimal sets of mask planes needed for LSST
    interpPlane = planes.index("S_MASK_INTERP") + 1
    satPlane = planes.index("S_MASK_SATUR") + 1
    crPlane = planes.index("S_MASK_CR") + 1

    interpBitMask = afwImage.Mask.getPlaneBitMask("INTRP")
    satBitMask = afwImage.Mask.getPlaneBitMask("SAT")
    crBitMask = afwImage.Mask.getPlaneBitMask("CR")

    listToSet = [(interpPlane, interpBitMask),
                 (satPlane, satBitMask),
                 (crPlane, crBitMask)]

    # Add the rest of the SDSS planes
    if allPlanes:
        for plane in ['S_MASK_NOTCHECKED', 'S_MASK_OBJECT', 'S_MASK_BRIGHTOBJECT',
                      'S_MASK_BINOBJECT', 'S_MASK_CATOBJECT', 'S_MASK_SUBTRACTED', 'S_MASK_GHOST']:
            idx = planes.index(plane) + 1
            planeName = re.sub("S_MASK_", "", plane)
            mask.addMaskPlane(planeName)
            planeBitMask = afwImage.Mask.getPlaneBitMask(planeName)
            listToSet.append((idx, planeBitMask))

    for plane, bitmask in listToSet:
        if len(hdr) < plane:
            continue

        if hdr[plane].data is None:
            continue

        nmask = len(hdr[plane].data)
        for i in range(nmask):
            frow = hdr[plane].data[i]
            Objmask(frow, bitmask).setMask(mask)

    return mask


if __name__ == '__main__':
    infile = sys.argv[1]
    outfile = sys.argv[2]

    if not os.path.isfile(infile):
        sys.exit(1)

    convertfpM(infile).writeFits(outfile)

    comparison = """  # noqa ignore the really long line
import lsst.afw.image as afwImage
import numpy as num
import os

readMask = "/astro/users/acbecker/LSST/lsst_devel/clue/Coadd/src/native/SDSS_PSFs/readAtlasImages_mod_by_Keith/src/read_mask"
cmd      = readMask + " %s INTERP /tmp/interp.fit"
os.system(cmd)
cmd      = readMask + " %s SATUR /tmp/satur.fit"
os.system(cmd)
cmd      = readMask + " %s CR /tmp/cr.fit"
os.system(cmd)

interp  = afwImage.ImageU("/tmp/interp.fit")
interp *= 2**2

sat     = afwImage.ImageU("/tmp/satur.fit")
sat    *= 2**1

cr      = afwImage.ImageU("/tmp/cr.fit")
cr     *= 2**3

interp += sat
interp += cr

lsst    = afwImage.ImageU("%s")
interp -= lsst
interp.writeFits("/tmp/mask_diff.fits")
print len(num.where(interp.getArray() != 0)[0])

""" % (infile, infile, infile, outfile)
    print(comparison)
