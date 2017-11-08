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
from __future__ import absolute_import, division, print_function
import os
import unittest

import lsst.utils
import lsst.utils.tests
import lsst.daf.persistence as dafPersist
from lsst.daf.base import DateTime
import lsst.afw.image
import lsst.afw.detection


class SdssMapperTestCase(lsst.utils.tests.TestCase):
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
            self.assertEqual(msk.__class__, lsst.afw.image.Mask[lsst.afw.image.MaskPixel])
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

            tsField = ref.get("tsField")
            self.assertAlmostEqual(tsField.gain, 4.72, 2)
            self.assertAlmostEqual(tsField.airmass, 1.2815132857671601)
            self.assertAlmostEqual(tsField.exptime, 53.907456)
            predDateStart = DateTime(53664.226070589997, DateTime.MJD, DateTime.TAI)
            predDateAvg = DateTime(predDateStart.nsecs(DateTime.TAI) + int(0.5e9*tsField.exptime),
                                   DateTime.TAI)
            self.assertAlmostEqual(tsField.dateAvg.get(), predDateAvg.get())


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
