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
import sys, os
import numpy as num

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from convertfpM import convertfpM
from convertasTrans import convertasTrans

try:
    import pymssql 
except:
    print "You need pymssql installed to access the DB"
    sys.exit(1)

rootdir = "/astro/net/pogo1/stripe82/imaging"
db      = pymssql.connect(user="clue-1", password="wlZH2xWy", host="fatboy.npl.washington.edu", database="clue", port=1433)
cursor  = db.cursor()

class FieldMatch(object):
    def __init__(self, args):
        self.run    = args[0]
        self.rerun  = args[1]
        self.filt   = args[2]
        self.camcol = args[3]
        self.field  = args[4]
        self.strip  = args[5]

        self.fpC    = None
        self.fpM    = None
        self.wcs    = None

    def loadfpC(self):
        self.fpC = getfpC(self.run, self.rerun, self.filt, self.camcol, self.field)

    def loadfpM(self):
        self.fpM = getfpM(self.run, self.rerun, self.filt, self.camcol, self.field)

    def loadWcs(self):
        asTrans = getasTrans(self.run, self.rerun)
        if asTrans:
            self.wcs = convertasTrans(asTrans, self.filt, self.camcol, self.field)

class MatchBackgroundsConfig(pexConfig.Config):
    kernelName = pexConfig.Field(
        dtype = str,
        doc = """Type of kernel for remapping""",
        default = "lanczos3"
    )


class MatchBackgrounds(pipeBase.Task):
    ConfigClass = MatchBackgroundsConfig
    def __init__(self, refrun, rerun, camcol, filt, *args, **kwargs):
        pipeBase.Task.__init__(self, *args, **kwargs)
        self.refrun  = refrun
        self.rerun   = rerun
        self.camcol  = camcol
        self.filt    = filt
        self.asTrans = getasTrans(self.refrun, self.rerun)

        self.refExp      = {} # per field
        self.stitchedExp = {} # per run
        self.warpedExp   = {} # per run

        self.warper = afwMath.Warper(self.config.kernelName)
        self.writeFits = True

    @pipeBase.timeMethod
    def run(self, fields, **kwargs):

        for field in fields:
            fpCRef = getfpC(self.refrun, self.rerun, self.filt, self.camcol, field)
            if not fpCRef:
                continue
            
            fpMRef = getfpM(self.refrun, self.rerun, self.filt, self.camcol, field)
            if not fpMRef:
                continue
            
            wcsRef = convertasTrans(self.asTrans, self.filt, self.camcol, field)
            if not wcsRef:
                continue
        
            # Assemble an exposure out of this
            # Unknown gain for now
            varRef = afwImage.ImageF(fpCRef, True)
            miRef  = afwImage.MaskedImageF(fpCRef, fpMRef, varRef)
            self.refExp[field] = afwImage.ExposureF(miRef, wcsRef)
            
            matches = self.queryClue(fpCRef.getBBox(), wcsRef, self.filt)
            self.processMatches(matches, field)

    @pipeBase.timeMethod
    def queryClue(self, bbox, wcs, filt):
        LLC = wcs.pixelToSky(bbox.getMinX(), bbox.getMinY())
        ULC = wcs.pixelToSky(bbox.getMinX(), bbox.getMaxY())
        URC = wcs.pixelToSky(bbox.getMaxX(), bbox.getMaxY())
        LRC = wcs.pixelToSky(bbox.getMaxX(), bbox.getMinY())
    
        sql  = "select run,rerun,filter,camCol,field,strip from clue.dbo.SeasonStripColor_bboxNoView WITH(INDEX(idx_bbox))"
        sql += " where bbox.STIntersects(geography::STGeomFromText('POLYGON (("
        sql += " %f %f," % (LLC[0].asDegrees(), LLC[1].asDegrees())
        sql += " %f %f," % (ULC[0].asDegrees(), ULC[1].asDegrees())
        sql += " %f %f," % (URC[0].asDegrees(), URC[1].asDegrees())
        sql += " %f %f," % (LRC[0].asDegrees(), LRC[1].asDegrees())
        sql += " %f %f ))', 4326))=1" % (LLC[0].asDegrees(), LLC[1].asDegrees())
        sql += " and filter='%s'" % (filt)
    
    
        sql += " order by run asc, field asc;"   
        
        print sql
        cursor.execute(sql)
        results = cursor.fetchall()
    
        # Note: lets just add each strip up and we'll use the N/S overlap
        # to match coadd backgroudns
    
        amatches = []
        for result in results:
            amatches.append(FieldMatch(result))
    
        strips = num.array([x.strip for x in amatches])
        runs   = num.array([x.run for x in amatches])
        idx    = num.where(runs == self.refrun)[0]
        strip  = list(set(strips[idx]))
        if len(strip) != 1:
            print "ERROR in strips"
            sys.exit(1)
    
        idxs     = num.where(strips == strip)[0]
        smatches = []
        for idx in idxs:
            smatches.append(amatches[idx])
        return smatches

    @pipeBase.timeMethod
    def processMatches(self, matches, field, gain = 1.0, overlap = 128, testme = 0):
        runs  = num.array([x.run for x in matches])
        uruns = list(set(runs))
        uruns.sort()

        if self.writeFits:
            self.refExp[field].writeFits("/tmp/exp-%06d-%s%d-%04d.fits" % (self.refrun, self.filt, self.camcol, field))
    
        for run in uruns:
            print "RUNNING", run, "vs.", self.refrun

            if run == self.refrun:
                continue

            if run != 1755:
                continue
            
            idxs = num.where(runs == run)[0]
            if len(idxs) == 0:
                continue

            runMatches = []
            for idx in idxs:
                runMatches.append(matches[idx])

            # We need to pad the ends, it appears.  This is a bit
            # blunt, but for lack of time just do it for now...
            if False:
                firstMatch = FieldMatch((runMatches[0].run, runMatches[0].rerun, 
                                         runMatches[0].filt, runMatches[0].camcol, 
                                         runMatches[0].field-1, runMatches[0].strip))
                lastMatch  = FieldMatch((runMatches[-1].run, runMatches[-1].rerun, 
                                         runMatches[-1].filt, runMatches[-1].camcol, 
                                         runMatches[-1].field+1, runMatches[-1].strip))
                runMatches.insert(0, firstMatch)
                runMatches.append(lastMatch)

    
            nloaded = 0
            for match in runMatches:
                match.loadfpC()
                match.loadfpM()
                match.loadWcs()
                if match.fpC and match.fpM and match.wcs:
                    nloaded += 1
    
            if nloaded != len(runMatches):
                print "Not able to load all images, skipping to next run"
                continue
    
            width  = runMatches[0].fpC.getWidth()
            height = runMatches[0].fpC.getHeight() * nloaded - overlap * (nloaded - 1) + testme * nloaded
            stitch = afwImage.MaskedImageF(width, height)
    
            for i in range(len(runMatches)):
                match  = runMatches[i]
    
                symin  = (i + 0) * match.fpC.getHeight() + (i * testme)
                symax  = (i + 1) * match.fpC.getHeight() + (i * testme)
    
                iymin  = 0
                iymax  = match.fpC.getHeight()
                
                if i > 0:
                    iymin   = overlap
                    symin  -= (i - 1) * overlap
                    symax  -= (i - 0) * overlap
    
                # Note transpose of getArray()
                try:
                    stitch.getImage().getArray()[symin:symax,:] = match.fpC.getArray()[iymin:iymax,:]
                    stitch.getMask().getArray()[symin:symax,:] = match.fpM.getArray()[iymin:iymax,:]
                except:
                    import pdb; pdb.set_trace()
    
    
            # Variance from gain
            var  = stitch.getVariance()
            var  = afwImage.ImageF(stitch.getImage(), True)
            var /= gain
    
            # Keep Wcs of first image
            exp = afwImage.ExposureF(stitch, runMatches[0].wcs)
            self.stitchedExp[run] = exp
            self.warpedExp[run]   = self.warper.warpExposure(self.refExp[field].getWcs(), 
                                                             exp, 
                                                             destBBox = self.refExp[field].getBBox(afwImage.PARENT))

            if self.writeFits:
                self.stitchedExp[run].writeFits("/tmp/match-%06d-%s%d-%04d-r%06d.fits" % (self.refrun, self.filt, self.camcol, field, run))
                self.warpedExp[run].writeFits("/tmp/warp-%06d-%s%d-%04d-r%06d.fits" % (self.refrun, self.filt, self.camcol, field, run))
    
######
######
######
######

def getfpC(run, rerun, filt, camcol, field):
    fname = os.path.join(rootdir, str(run), str(rerun), "corr", str(camcol), "fpC-%06d-%s%d-%04d.fit.gz" % (run, filt, camcol, field))
    print fname
    if os.path.isfile(fname):
        im  = afwImage.ImageF(fname)
        im -= 1000 # damn pedestal
        return im
    return None

def getfpM(run, rerun, filt, camcol, field):
    fname1 = os.path.join(rootdir, str(run), str(rerun), "objcs", str(camcol), "fpM-%06d-%s%d-%04d.fit" % (run, filt, camcol, field))
    fname2 = os.path.join(rootdir, str(run), str(rerun), "objcs", str(camcol), "fpM-%06d-%s%d-%04d.fit.gz" % (run, filt, camcol, field))
    for fname in (fname1, fname2):
        print fname
        if os.path.isfile(fname):
            try:
                return convertfpM(fname, allPlanes = True)
            except:
                return None
    return None

def getasTrans(run, rerun):
    fname = os.path.join(rootdir, str(run), str(rerun), "astrom", "asTrans-%06d.fit" % (run))
    print fname
    if os.path.isfile(fname):
        return fname
    return None


if __name__ == '__main__':
    refrun  = int(sys.argv[1])
    camcol  = int(sys.argv[2])
    filt    = sys.argv[3]

    matcher = MatchBackgrounds(refrun, 40, camcol, filt)
    fields  = range(1, 1000)
    fields  = [11,]
    matcher.run(fields)
