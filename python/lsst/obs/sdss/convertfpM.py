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

import sys, os
import pyfits
import numpy as num
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom

class Span(object):
    def __init__(self, y, x1, x2):
        self.y  = y
        self.x1 = x1
        self.x2 = x2

class Objmask(object):
    nperspan = 6

    def __init__(self, frow, cval):
        self.refcntr = frow[0]
        self.nspan   = frow[1]
        self.row0    = frow[2]
        self.col0    = frow[3]
        self.rmin    = frow[4]
        self.rmax    = frow[5]
        self.cmin    = frow[6]
        self.cmax    = frow[7]
        self.npix    = frow[8]
        self.span    = frow[9]

        self.spans   = []
        npixcheck    = 0
        for i in range(self.nspan):
            b1 = self.span[6*i + 0]
            b2 = self.span[6*i + 1]
            y  = (b1 << 8) + b2
            b1 = self.span[6*i + 2]
            b2 = self.span[6*i + 3]
            x1 = (b1 << 8) + b2
            b1 = self.span[6*i + 4]
            b2 = self.span[6*i + 5]
            x2 = (b1 << 8) + b2
            self.spans.append(Span(y, x1, x2))
            for i in range(x1, x2+1): npixcheck += 1

        assert(self.npix == npixcheck)

        self.cval    = cval

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
           
            x = int(x1)
            while x <= x2:
                mask.set(x, y, mask.get(x, y) | self.cval)
                x += 1


if __name__ == '__main__':
    infile  = sys.argv[1]
    outfile = sys.argv[2]
    
    if not os.path.isfile(infile):
        sys.exit(1)
    
    hdr    = pyfits.open(infile)
    run    = hdr[0].header['RUN']
    camcol = hdr[0].header['CAMCOL']
    field  = hdr[0].header['FIELD']
    nRows  = hdr[0].header['MASKROWS']
    nCols  = hdr[0].header['MASKCOLS']
    nPlane = hdr[0].header['NPLANE']

    planes = hdr[-1].data["attributeName"].tolist()
    values = hdr[-1].data["Value"].tolist()

    mask   = afwImage.MaskU(afwGeom.ExtentI(nCols, nRows))

    interpPlane = planes.index("S_MASK_INTERP") + 1
    satPlane    = planes.index("S_MASK_SATUR") + 1
    crPlane     = planes.index("S_MASK_CR") + 1

    interpBitMask = afwImage.MaskU_getPlaneBitMask("INTRP")
    satBitMask    = afwImage.MaskU_getPlaneBitMask("SAT")
    crBitMask     = afwImage.MaskU_getPlaneBitMask("CR")

    for plane, bitmask in [ (interpPlane, interpBitMask),
                            (satPlane, satBitMask),
                            (crPlane, crBitMask) ]:

        for frow in hdr[plane].data:
            Objmask(frow, bitmask).setMask(mask)
    
    mask.writeFits(outfile)





    comparison = """
import lsst.afw.image as afwImage

interp  = afwImage.ImageU("interp.fit")
interp *= 2**2

sat     = afwImage.ImageU("satur.fit")
sat    *= 2**1

cr      = afwImage.ImageU("cr.fit")
cr     *= 2**3

interp += sat
interp += cr
interp.writeFits("compare.fits")

"""
