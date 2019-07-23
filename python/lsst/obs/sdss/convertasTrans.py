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

from astropy.io import fits
import numpy as np

import lsst.afw.image as afwImage
from lsst.afw.geom import makeSkyWcs
import lsst.afw.table as afwTable
import lsst.geom as geom
import lsst.meas.astrom.sip as sip

deg2rad = np.pi / 180.
rad2deg = 180. / np.pi


class CoordinateMapper(object):
    # COMMENT mu nu are defined as:
    # COMMENT   r'-i' < riCut:
    # COMMENT       rowm = row+dRow0+dRow1*col+dRow2*(col^2)+dRow3*(col^3)+csRow*c
    # COMMENT       colm = col+dCol0+dCol1*col+dCol2*(col^2)+dCol3*(col^3)+csCol*c
    # COMMENT   r'-i' >= riCut
    # COMMENT       rowm = row+dRow0+dRow1*col+dRow2*(col^2)+dRow3*(col^3)+ccRow
    # COMMENT       colm = col+dCol0+dCol1*col+dCol2*(col^2)+dCol3*(col^3)+ccCol
    # COMMENT   mu = a + b * rowm + c * colm
    # COMMENT   nu = d + e * rowm + f * colm

    # We will do an "average" mapping for an r-i=0 object so that the
    # last cs/cc terms are not necessary

    def __init__(self, node_rad, incl_rad, dRow0, dRow1, dRow2, dRow3, dCol0, dCol1, dCol2, dCol3,
                 a, b, c, d, e, f, cOffset=+0.5):
        # Here cOffset reflects the differences between SDSS coords
        # (LLC = 0.5,0.5) and LSST coords (LLC = 0,0).  If SDSS
        # measures an object centered at (0.5,0.5) then LSST will
        # measure it at coordinate (0,0), but when using SDSS
        # astrometry we need to evaluate the equations at (x+0.5,
        # y+0.5)

        self.node_rad = node_rad
        self.incl_rad = incl_rad

        self.dRow0 = dRow0
        self.dRow1 = dRow1
        self.dRow2 = dRow2
        self.dRow3 = dRow3

        self.dCol0 = dCol0
        self.dCol1 = dCol1
        self.dCol2 = dCol2
        self.dCol3 = dCol3

        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.f = f

        self.cOff = cOffset

    def xyToMuNu(self, x, y):
        rowm = (y+self.cOff) + self.dRow0 + self.dRow1*(x+self.cOff) + self.dRow2*((x+self.cOff)**2) + \
            self.dRow3*((x+self.cOff)**3)
        colm = (x+self.cOff) + self.dCol0 + self.dCol1*(x+self.cOff) + self.dCol2*((x+self.cOff)**2) + \
            self.dCol3*((x+self.cOff)**3)

        mu_deg = self.a + self.b * rowm + self.c * colm
        nu_deg = self.d + self.e * rowm + self.f * colm
        mu_rad = mu_deg * deg2rad
        nu_rad = nu_deg * deg2rad

        return mu_rad, nu_rad

    def muNuToRaDec(self, mu_rad, nu_rad):
        x2 = np.cos(mu_rad - self.node_rad) * np.cos(nu_rad)
        y2 = np.sin(mu_rad - self.node_rad) * np.cos(nu_rad)
        z2 = np.sin(nu_rad)
        y1 = y2 * np.cos(self.incl_rad) - z2 * np.sin(self.incl_rad)
        z1 = y2 * np.sin(self.incl_rad) + z2 * np.cos(self.incl_rad)

        ra_rad = np.arctan2(y1, x2) + self.node_rad
        dec_rad = np.arcsin(z1)

        return ra_rad, dec_rad

    def xyToRaDec(self, x, y):
        mu_rad, nu_rad = self.xyToMuNu(x, y)
        return self.muNuToRaDec(mu_rad, nu_rad)


def createWcs(x, y, mapper, order=4, cOffset=1.0):
    # Here cOffset reflects the differences between FITS coords (LLC =
    # 1,1) and LSST coords (LLC = 0,0).  That is, when creating a Wcs
    # from scratch, we need to evaluate our WCS at coordinate 0,0 to
    # create CRVAL, but set CRPIX to 1,1.

    ra_rad, dec_rad = mapper.xyToRaDec(x, y)

    # Minimial table for sky coordinates
    catTable = afwTable.SimpleTable.make(afwTable.SimpleTable.makeMinimalSchema())

    # Minimial table + centroids for focal plane coordintes
    srcSchema = afwTable.SourceTable.makeMinimalSchema()
    centroidKey = afwTable.Point2DKey.addFields(srcSchema, "centroid", "centroid", "pixel")

    srcTable = afwTable.SourceTable.make(srcSchema)
    srcTable.defineCentroid("centroid")

    matches = []
    for i in range(len(x)):
        src = srcTable.makeRecord()
        src.set(centroidKey.getX(), x[i])
        src.set(centroidKey.getY(), y[i])

        cat = catTable.makeRecord()
        cat.set(catTable.getCoordKey().getRa(), geom.Angle(ra_rad[i], geom.radians))
        cat.set(catTable.getCoordKey().getDec(), geom.Angle(dec_rad[i], geom.radians))

        mat = afwTable.ReferenceMatch(cat, src, 0.0)
        matches.append(mat)

    # Need to make linear Wcs around which to expand solution

    # CRPIX1  = Column Pixel Coordinate of Ref. Pixel
    # CRPIX2  = Row Pixel Coordinate of Ref. Pixel
    crpix = geom.Point2D(x[0] + cOffset, y[0] + cOffset)

    # CRVAL1  = RA at Reference Pixel
    # CRVAL2  = DEC at Reference Pixel
    crval = geom.SpherePoint(ra_rad[0], dec_rad[0], geom.radians)

    # CD1_1   = RA  degrees per column pixel
    # CD1_2   = RA  degrees per row pixel
    # CD2_1   = DEC degrees per column pixel
    # CD2_2   = DEC degrees per row pixel
    LLl = mapper.xyToRaDec(0., 0.)
    ULl = mapper.xyToRaDec(0., 1.)
    LRl = mapper.xyToRaDec(1., 0.)

    LLc = geom.SpherePoint(LLl[0], LLl[1], geom.radians)
    ULc = geom.SpherePoint(ULl[0], ULl[1], geom.radians)
    LRc = geom.SpherePoint(LRl[0], LRl[1], geom.radians)

    cdN_1 = LLc.getTangentPlaneOffset(LRc)
    cdN_2 = LLc.getTangentPlaneOffset(ULc)
    cd1_1, cd2_1 = cdN_1[0].asDegrees(), cdN_1[1].asDegrees()
    cd1_2, cd2_2 = cdN_2[0].asDegrees(), cdN_2[1].asDegrees()

    cdMatrix = np.array([cd1_1, cd2_1, cd1_2, cd2_2], dtype=float)
    cdMatrix.shape = (2, 2)

    linearWcs = makeSkyWcs(crpix=crpix, crval=crval, cdMatrix=cdMatrix)
    wcs = sip.makeCreateWcsWithSip(matches, linearWcs, order).getNewWcs()

    return wcs


def validate(xs, ys, mapper, wcs):
    dists = []
    for i in range(len(xs)):
        tuple1 = mapper.xyToRaDec(xs[i], ys[i])
        coord1 = geom.SpherePoint(tuple1[0], tuple1[1], geom.radians)
        coord2 = wcs.pixelToSky(xs[i], ys[i])
        dist = coord1.separation(coord2).asArcseconds()
        dists.append(dist)

    print(np.mean(dists), np.std(dists))


def convertasTrans(infile, filt, camcol, field, stepSize=50, doValidate=False):
    with fits.open(infile) as hdulist:
        t0 = hdulist[0].header['ccdarray']
        if t0 != 'photo':
            raise RuntimeError('*** Cannot support ccdarray: %s' % (t0,))

        camcols = hdulist[0].header['camcols']
        filters = hdulist[0].header['filters']
        node_deg = hdulist[0].header['node']
        incl_deg = hdulist[0].header['incl']
        node_rad = node_deg * deg2rad
        incl_rad = incl_deg * deg2rad

        cList = [int(cc) for cc in camcols.split()]
        fList = filters.split()

        try:
            cIdx = cList.index(camcol)
        except Exception:
            print("Cannot extract data for camcol %s" % (camcol))
            return None

        try:
            fIdx = fList.index(filt)
        except Exception:
            print("Cannot extract data for filter %s" % (filt))
            return None

        ext = cIdx * len(fList) + (fIdx + 1)
        ehdr = hdulist[ext].header
        edat = hdulist[ext].data

    if ehdr['CAMCOL'] != camcol or ehdr['FILTER'] != filt:
        print("Extracted incorrect header; fix me")
        return None

    fields = edat.field('field').tolist()
    try:
        fIdx = fields.index(field)
    except Exception:
        print("Cannot extract data for field %d" % (field))
        return None

    dRow0 = edat.field('dRow0')[fIdx]
    dRow1 = edat.field('dRow1')[fIdx]
    dRow2 = edat.field('dRow2')[fIdx]
    dRow3 = edat.field('dRow3')[fIdx]

    dCol0 = edat.field('dCol0')[fIdx]
    dCol1 = edat.field('dCol1')[fIdx]
    dCol2 = edat.field('dCol2')[fIdx]
    dCol3 = edat.field('dCol3')[fIdx]

    a = edat.field('a')[fIdx]
    b = edat.field('b')[fIdx]
    c = edat.field('c')[fIdx]
    d = edat.field('d')[fIdx]
    e = edat.field('e')[fIdx]
    f = edat.field('f')[fIdx]

    # We need to fit for a TAN-SIP
    x = np.arange(0, 1489+stepSize, stepSize)
    y = np.arange(0, 2048+stepSize, stepSize)
    coords = np.meshgrid(x, y)
    xs = np.ravel(coords[0]).astype(np.float)
    ys = np.ravel(coords[1]).astype(np.float)
    mapper = CoordinateMapper(node_rad, incl_rad, dRow0, dRow1, dRow2, dRow3, dCol0, dCol1, dCol2, dCol3,
                              a, b, c, d, e, f)
    wcs = createWcs(xs, ys, mapper)

    if doValidate:
        validate(xs, ys, mapper, wcs)

    return wcs


if __name__ == '__main__':
    infile = sys.argv[1]
    filt = sys.argv[2]
    camcol = int(sys.argv[3])
    field = int(sys.argv[4])

    wcs = convertasTrans(infile, filt, camcol, field, doValidate=True)

    if len(sys.argv) > 5:
        fpC = sys.argv[5]
        image = afwImage.ImageF(fpC)
        mi = afwImage.MaskedImageF(image)
        exp = afwImage.ExposureF(mi, wcs)
        exp.writeFits("/tmp/exp.fits")
