#!/usr/bin/env python
"""
Generate the SDSS camera registry FITS files from the Yanny files in etc/.

Scons should have automatically run this when building obs_sdss. To produce
the same files that scons would have, run with no arguments.
"""
from __future__ import absolute_import, division, print_function

import os
import shutil
import argparse

import lsst.utils
from lsst.obs.sdss import makeCamera

if __name__ == "__main__":
    path = lsst.utils.getPackageDir('obs_sdss')
    outDir = os.path.join(path, 'description', 'camera')

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--clobber", action="store_true", dest="clobber", default=False,
                        help=("remove and re-create the output directory if it already exists?"))
    args = parser.parse_args()

    def makeDir(dirPath, doClobber=False):
        """Make a directory; if it exists then clobber or fail, depending on doClobber

        @param[in] dirPath: path of directory to create
        @param[in] doClobber: what to do if dirPath already exists:
            if True and dirPath is a dir, then delete it and recreate it, else raise an exception
        @throw RuntimeError if dirPath exists and doClobber False
        """
        if os.path.exists(dirPath):
            if doClobber and os.path.isdir(dirPath):
                print("Clobbering directory %r" % (dirPath,))
                shutil.rmtree(dirPath)
            else:
                raise RuntimeError("Directory %r exists. Will not overwrite." % (dirPath,))
        print("Creating directory %r" % (dirPath,))
        os.makedirs(dirPath)

    makeDir(dirPath=outDir, doClobber=args.clobber)

    makeCamera.makeCamera(outputDir=outDir)
