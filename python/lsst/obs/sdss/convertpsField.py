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
import lsst.afw.math as afwMath
DEBUG = True

filtToHdu = {'u':1,'g':2,'r':3,'i':4,'z':5}

# Mapping from psField coefficient locations to PolynomialFunction2 locations
skMatrixPos2TriSeqPosT = [
	0,	2,	5,	9,	14,
	1,	4,	8,	13,	19,
	3,	7,	12,	18,	25,
	6,	11,	17,	24,	32,
	10,	16,	23,	31,	40
]

def convertpsField(infile, filt, trim = True, rcscale = 0.001, MAX_ORDER_B = 5, LSST_ORDER = 4):
    if not filt in filtToHdu.keys():
        print "INVALID FILTER", filt
        sys.exit(1)

    buff = open(infile, "rb")
    pstruct = pyfits.getdata(buff, ext = filtToHdu[filt])

    spaParList = [[]]*len(pstruct)
    kernelList = afwMath.KernelList()
    for i in range(len(pstruct)):
        nrow_b = pstruct[i][0]
        ncol_b = pstruct[i][1]
        cmat   = pstruct[i][2].reshape((MAX_ORDER_B, MAX_ORDER_B))
        e_val  = pstruct[i][3]
        krow   = pstruct[i][4] # RNROW
        kcol   = pstruct[i][5] # RNCOL
        karr   = pstruct[i][7].reshape((krow, kcol)).astype(num.float64)

        if trim:
            karr = karr[10:41, 10:41]
        kim  = afwImage.ImageD(karr)
        kern = afwMath.FixedKernel(kim)
        kernelList.push_back(kern)

        # NOTES:

        # Afw has the polynomial terms like:
        #
        # * f(x,y) =   c0                                                      (0th order)
        # *          + c1 x    + c2 y                                          (1st order)
        # *          + c3 x^2  + c4 x y    + c5 y^2                            (2nd order)
        # *          + c6 x^3  + c7 x^2 y  + c8 x y^2    + c9 y^3              (3rd order)
        # *          + c10 x^4 + c11 x^3 y + c12 x^2 y^2 + c13 x y^3 + c14 y^4 (4th order)
        #
        # So, ordered: x^0,y^0 x^1,y^0 x^0,y^1 x^2,y^0 x^1,y^1 x^0,y^2

        # SDSS has the terms ordered like, after reshape():
        #
        # x^0,y^0 x^0,y^1 x^0,y^2 
        # x^1,y^0 x^1,y^1 x^1,y^2
        # x^2,y^0 x^2,y^1 x^2,y^2
        #
        # So, it technically goes up to fourth order in LSST-speak.  OK, that is the trick.
        # 
        # Mapping:
        # cmat[0][0] = c0  x^0,y^0
        # cmat[1][0] = c1  x^1,y^0
        # cmat[2][0] = c3  x^2,y^0
        # cmat[0][1] = c2  x^0,y^1
        # cmat[1][1] = c4  x^1,y^1
        # cmat[2][1] = c7  x^2,y^1
        # cmat[0][2] = c5  x^0,y^2 
        # cmat[1][2] = c8  x^1,y^2
        # cmat[2][2] = c12 x^2,y^2
        #
        # This is quantified in skMatrixPos2TriSeqPosT

        spaParamsTri = num.zeros(MAX_ORDER_B * MAX_ORDER_B)
        for k in range(nrow_b * ncol_b):
            row         = k % nrow_b
            col         = k // nrow_b
            coeff       = cmat[row,col]
            scale       = pow(rcscale, row) * pow(rcscale, col)
            scaledCoeff = coeff * scale
            idx         = row * MAX_ORDER_B + col 

            #print row, col, cmat[row,col], idx, skMatrixPos2TriSeqPosT[idx], scaledCoeff

            spaParamsTri[skMatrixPos2TriSeqPosT[idx]] = scaledCoeff

        #print spaParamsTri
        nTerms = (LSST_ORDER + 1) * (LSST_ORDER + 2) // 2
        spaParamsTri  = spaParamsTri[:nTerms]
        spaParList[i] = spaParamsTri

    buff.close()
    spaFun = afwMath.PolynomialFunction2D(LSST_ORDER)
    spatialKernel = afwMath.LinearCombinationKernel(kernelList, spaFun)
    spatialKernel.setSpatialParameters(spaParList)
    return spatialKernel

def directCompare(infile, filt, x, y, soft_bias = 1000, amp = 30000, outfile = "/tmp/sdss_psf.fits"):
    if not filt in filtToHdu.keys():
        print "INVALID FILTER", filt
        sys.exit(1)
    
    # Make the kernel image from LSST
    kernel = convertpsField(infile, filt, trim = False)

    # Assumes you have built dervish and have read_PSF in your path
    cmd = "read_PSF %s %s %f %f %s" % (infile, filtToHdu[filt], x, y, outfile)
    os.system(cmd)
    if not os.path.isfile(outfile):
        print "Cannot find SDSS-derived kernel", outfile
        sys.exit(1)

    kImage1  = afwImage.ImageD(outfile)
    kImage1 -= soft_bias
    kImage1 /= (amp - soft_bias)
    maxVal   = afwMath.makeStatistics(kImage1, afwMath.MAX).getValue(afwMath.MAX)
    print "TEST 1", maxVal == 1.0
    kImage1.writeFits("/tmp/sdss_psf_scaled.fits")

    # 
    kImage2  = afwImage.ImageD(kernel.getDimensions())
    kernel.computeImage(kImage2, True, x, y)
    maxVal   = afwMath.makeStatistics(kImage2, afwMath.MAX).getValue(afwMath.MAX)
    kImage2 /= maxVal
    kImage2.writeFits("/tmp/kernel.fits")

    kImage2 -= kImage1
    kImage2.writeFits("/tmp/diff.fits")
    residSum = afwMath.makeStatistics(kImage2, afwMath.SUM).getValue(afwMath.SUM)
    print "TEST 2", residSum

if __name__ == '__main__':
    infile  = sys.argv[1]
    filt    = sys.argv[2]
    x       = float(sys.argv[3])
    y       = float(sys.argv[4])
    outfile = sys.argv[5]

    if not os.path.isfile(infile):
        sys.exit(1)
    
    if DEBUG:
        directCompare(infile, filt, x, y)
    else:
        kernel = convertpsField(infile, filt)
    
    # Persist the kernel at your own leisure
