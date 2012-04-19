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

import unittest
import lsst.utils.tests as utilsTests

import os
import lsst.daf.persistence as dafPersist
from lsst.obs.sdss import SdssMapper
import lsst.afw.image
import lsst.afw.detection

class SdssMapperTestCase(unittest.TestCase):
    """A test case for the SdssMapper."""

    def testGetDR7(self):
        butler = dafPersist.ButlerFactory(
                mapper=SdssMapper(root="/lsst7/stripe82/dr7/runs",
                    registry=os.path.join(os.environ['OBS_SDSS_DIR'],
                        "bin", "registry3.sqlite3")
                    )
                ).create()
        sub = butler.subset("fpC", run=5754, camcol=3, frame=280, band="r")
        self.assertEqual(len(sub), 1)
        for ref in sub:
            im = ref.get("fpC")
            w, h = im.getWidth(), im.getHeight()
            self.assertEqual(im.__class__, lsst.afw.image.ImageU)
            self.assertEqual(w, 2048)
            self.assertEqual(h, 1489)
        
            im_md = ref.get("fpC_md")
            self.assertEqual(im_md.get("RUN"), 5754)
            self.assertEqual(im_md.get("FRAME"), 280)
            self.assertEqual(im_md.get("STRIPE"), 82)
        
            msk = ref.get("fpM")
            w, h = msk.getWidth(), msk.getHeight()
            self.assertEqual(msk.__class__, lsst.afw.image.MaskU)
            self.assertEqual(w, 2048)
            self.assertEqual(h, 1489)

            psf = ref.get("psField")
            k = psf.getKernel()
            w, h = k.getWidth(), k.getHeight()
            self.assertEqual(psf.__class__,
                    lsst.afw.detection.detectionLib.KernelPsf)
            self.assertEqual(w, 31)
            self.assertEqual(h, 31)

    def testGetCoadd(self):
        butler = dafPersist.ButlerFactory(
                mapper=SdssMapper(root="/lsst7/stripe82/uw-coadds")).create()
        coadd = butler.get("coadd", run=6383, camcol=3, frame=280, band="r")
        w, h = coadd.getWidth(), coadd.getHeight()
        self.assertEqual(coadd.__class__, lsst.afw.image.ImageF)
        self.assertEqual(w, 1489)
        self.assertEqual(h, 2048)

        coadd_md = butler.get("coadd_md", run=6383, camcol=3, frame=280, band="r")
        self.assertAlmostEqual(coadd_md.get("CRPIX1"), -0.047789988341666, 9)
        self.assertAlmostEqual(coadd_md.get("CRPIX2"), 0.083988780414136, 9)

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(SdssMapperTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
