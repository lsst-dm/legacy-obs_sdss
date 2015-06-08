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
import os
import unittest

import lsst.utils
import lsst.utils.tests as utilsTests
import lsst.daf.base as dafBase
import lsst.daf.persistence as dafPersist
import lsst.afw.image
import lsst.afw.detection

class SdssMapperTestCase(unittest.TestCase):
    """A test case for the SdssMapper."""

    def testGetDR7(self):
        obsSdssDir = lsst.utils.getPackageDir('obs_sdss')
        butler = dafPersist.Butler(
            root=os.path.join(obsSdssDir, "tests", "data", "dr7", "runs"))
        sub = butler.subset("fpC", run=5754, camcol=3, field=280, filter="r")
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
            self.assertEqual(psf.__class__, lsst.meas.algorithms.PcaPsf)
            self.assertEqual(w, 31)
            self.assertEqual(h, 31)

            wcs = ref.get("asTrans")
            self.assertEqual(wcs.__class__, lsst.afw.image.TanWcs)
            self.assertFalse(wcs.isFlipped())
            self.assertAlmostEqual(wcs.getFitsMetadata().get("CRPIX1"), 1.0, 5)
            self.assertAlmostEqual(wcs.getFitsMetadata().get("CRPIX2"), 1.0, 5)

            calib, gain = ref.get("tsField")
            self.assertAlmostEqual(calib.getMidTime().mjd(),
                    53664.2260706 + 0.5 * 53.907456/3600/24, 7)
            self.assertAlmostEqual(calib.getExptime(), 53.907456, 6)
            self.assertAlmostEqual(gain, 4.72, 2)

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
