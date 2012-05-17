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
import sys, os, re
import numpy as num

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDetect
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.ip.diffim as ipDiffim
import lsst.coadd.utils as coaddUtils

from lsst.pipe.tasks.coadd import CoaddTask
from convertfpM import convertfpM
from convertasTrans import convertasTrans
from convertpsField import convertpsField
from converttsField import converttsField
from scipy.interpolate import SmoothBivariateSpline
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
        self.psf    = None
        self.gain   = None
        self.calib  = None

    def loadfpC(self):
        self.fpC = getfpC(self.run, self.rerun, self.filt, self.camcol, self.field)

    def loadfpM(self):
        self.fpM = getfpM(self.run, self.rerun, self.filt, self.camcol, self.field)

    def loadWcs(self):
        asTrans = getasTrans(self.run, self.rerun)
        if asTrans:
            self.wcs = convertasTrans(asTrans, self.filt, self.camcol, self.field)

    def loadPsf(self):
        self.psf = getpsField(self.run, self.rerun, self.filt, self.camcol, self.field)

    def loadCalib(self):
        self.calib, self.gain = gettsField(self.run, self.rerun, self.filt, self.camcol, self.field)

    def createExp(self):
        var  = afwImage.ImageF(self.fpC, True)
        var /= self.gain
        mi   = afwImage.MaskedImageF(self.fpC, self.fpM, var)
        exp  = afwImage.ExposureF(mi, self.wcs)
        exp.setPsf(self.psf)
        exp.setCalib(self.calib)
        return exp

class SigmaClippedCoaddConfig(CoaddTask.ConfigClass):
    sigmaClip = pexConfig.Field(
        dtype = float,
        doc = "sigma for outlier rejection",
        default = 5.0,
        optional = None,
    )
    clipIter = pexConfig.Field(
        dtype = int,
        doc = "number of iterations of outlier rejection",
        default = 2,
        optional = False,
    )

class MatchBackgroundsConfig(pexConfig.Config):
    warpingKernelName = pexConfig.Field(
        dtype = str,
        doc = """Type of kernel for remapping""",
        default = "lanczos3"
    )
    backgroundOrder = pexConfig.Field(
        dtype = int,
        doc = """Order of background Chebyshev""",
        default = 2
    )
    backgroundBinsize = pexConfig.Field(
        dtype = int,
        doc = """Bin size for background matching""",
        default = 256
    )
    writeFits = pexConfig.Field(
        dtype = bool,
        doc = """Write output fits files""",
        default = False
    )
    outputPath = pexConfig.Field(
        dtype = str,
        doc = """Location of output files""",
        default = "/tmp"
    )

    refPsfSize = pexConfig.Field(
        dtype = int,
        doc = """Size of reference Psf matrix; must be same size as SDSS Psfs""",
        default = 31
    )
    refPsfSigma = pexConfig.Field(
        dtype = float,
        doc = """Gaussian sigma for Psf FWHM (pixels)""",
        default = 4.0
    )

    # With linear background model, this should fail
    # /astro/net/pogo1/stripe82/imaging/6447/40/corr/1/fpC-006447-r1-0718.fit.gz
    maxBgRms = pexConfig.Field(
        dtype = float,
        doc = """Maximum RMS of matched background differences, in counts""",
        default = 5.0
    )

    # Clouds
    # /astro/net/pogo1/stripe82/imaging/7071/40/corr/1/fpC-007071-r1-0190.fit.gz
    minFluxMag0 = pexConfig.Field(
        dtype = float,
        doc = """Minimum flux for star of mag 0""",
        default = 1.0e+10
    )

class SigmaClippedCoaddTask(CoaddTask):
    ConfigClass = SigmaClippedCoaddConfig
    def __init__(self, *args, **kwargs):
        CoaddTask.__init__(self, *args, **kwargs)

        # Stats object for sigma clipping
        self.statsCtrl = afwMath.StatisticsControl()
        self.statsCtrl.setNumSigmaClip(self.config.sigmaClip)
        self.statsCtrl.setNumIter(self.config.clipIter)
        self.statsCtrl.setAndMask(afwImage.MaskU.getPlaneBitMask(self.config.coadd.badMaskPlanes))

    @pipeBase.timeMethod   
    def normalizeForCoadd(self, exp):
        # WARNING: IF YOU HAVE BACKGROUND MATCHED IMAGES THIS WILL
        # DESTROY THEIR MATCHING.  RUN THIS BEFORE BACKGROUND MATCHING
        calib = exp.getCalib()
        scaleFac = 1.0 / calib.getFlux(self.config.coadd.coaddZeroPoint)
        mi  = afwImage.MaskedImageF(exp.getMaskedImage(), True)
        mi *= scaleFac
        self.log.log(self.log.INFO, "Normalized using scaleFac=%0.3g" % (scaleFac))

        # Copy Psf, Wcs, Calib (bad idea?)
        normExp = afwImage.ExposureF(exp, True)
        normExp.setMaskedImage(mi)
        return normExp

    @pipeBase.timeMethod   
    def weightForCoadd(self, exp, weightFactor = 1.0):
        statObj = afwMath.makeStatistics(exp.getMaskedImage().getVariance(), exp.getMaskedImage().getMask(),
                                         afwMath.MEANCLIP, self.statsCtrl)
        meanVar, meanVarErr = statObj.getResult(afwMath.MEANCLIP)
        weight = weightFactor / float(meanVar)
        return weight

    @pipeBase.timeMethod   
    def run(self, refExp, expList):

        # Calib object for coadd (i.e. zeropoint)
        coaddCalib  = coaddUtils.makeCalib(self.config.coadd.coaddZeroPoint)

        # Destination for the coadd
        refMi    = refExp.getMaskedImage()
        coaddMi  = refMi.Factory(refMi.getBBox(afwImage.PARENT))

        # Vectors for images to coadd and their weights
        maskedImageList = afwImage.vectorMaskedImageF()
        weightList = []
        
        # Add reference image first
        maskedImageList.append(refExp.getMaskedImage())
        weightList.append(self.weightForCoadd(refExp))
        
        # All the other matched images
        for exp in expList:
            maskedImageList.append(exp.getMaskedImage())
            weightList.append(self.weightForCoadd(exp))

        try:
            coaddMi = afwMath.statisticsStack(maskedImageList, afwMath.MEANCLIP, self.statsCtrl, weightList)
        except Exception, e:
            self.log.log(self.log.ERR, "Outlier rejected coadd failed: %s" % (e,))
            sys.exit(1)

        # Post processing
        coaddUtils.setCoaddEdgeBits(coaddMi.getMask(), coaddMi.getVariance())

        coaddExp = afwImage.ExposureF(coaddMi, refExp.getWcs())
        coaddExp.setCalib(coaddCalib)
        self.log.log(self.log.INFO, "COADD with %d images" % (len(expList) + 1))
        self.log.log(self.log.INFO, "")
        self.log.log(self.log.INFO, self.metadata.toString())
        return coaddExp

            
class MatchBackgrounds(pipeBase.Task):
    ConfigClass = MatchBackgroundsConfig
    def __init__(self, refrun, rerun, camcol, filt, field, *args, **kwargs):
        pipeBase.Task.__init__(self, *args, **kwargs)
        self.refrun  = refrun
        self.rerun   = rerun
        self.camcol  = camcol
        self.filt    = filt
        self.field   = field
        self.asTrans = getasTrans(self.refrun, self.rerun)

        self.warper  = afwMath.Warper(self.config.warpingKernelName)
        self.refPsf  = afwDetect.createPsf("DoubleGaussian", self.config.refPsfSize, self.config.refPsfSize, self.config.refPsfSigma)
        
        config          = ipDiffim.ModelPsfMatchTask.ConfigClass()
        config.kernel.active.kernelSize = self.config.refPsfSize // 2
        self.psfMatcher = ipDiffim.ModelPsfMatchTask(config=config)

        self.coadder    = SigmaClippedCoaddTask()

    @pipeBase.timeMethod
    def run(self, nMax = 10, **kwargs):
        fpCRef = getfpC(self.refrun, self.rerun, self.filt, self.camcol, self.field)
        if not fpCRef:
            return
        
        fpMRef = getfpM(self.refrun, self.rerun, self.filt, self.camcol, self.field)
        if not fpMRef:
            return
        
        wcsRef = convertasTrans(self.asTrans, self.filt, self.camcol, self.field)
        if not wcsRef:
            return

        psfRef = getpsField(self.refrun, self.rerun, self.filt, self.camcol, self.field)
        if not psfRef:
            return

        calib, gain = gettsField(self.refrun, self.rerun, self.filt, self.camcol, self.field)
        if (not calib) or (not gain):
            return

        # Assemble an exposure out of this
        varRef  = afwImage.ImageF(fpCRef, True)
        varRef /= gain
        miRef   = afwImage.MaskedImageF(fpCRef, fpMRef, varRef)
        exp     = afwImage.ExposureF(miRef, wcsRef)
        exp.setPsf(psfRef)
        exp.setCalib(calib)

        # Ref exposure for this field
        self.refExp          = self.psfMatcher.run(exp, self.refPsf).psfMatchedExposure
        
        # Matches per run
        matches = self.queryClue(fpCRef.getBBox(), wcsRef, self.filt)
        self.warpedExp       = {}
        self.bgMatchedExp    = {}
        
        self.loadMatches(matches, nMax = nMax)
        self.matchBackgrounds()
        self.createCoadd()
        self.log.log(self.log.INFO, "")
        self.log.log(self.log.INFO, self.metadata.toString())

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
        
        self.log.log(self.log.INFO, "SQL: %s" % (sql))
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
    def loadMatches(self, matches, nMax = None):
        """Finds matches on disk; Psf matches them to common Psf;
        stitches them together; remaps to fiducial field; stores in
        self.warpedExp"""

        runs  = num.array([x.run for x in matches])
        uruns = list(set(runs))
        uruns.sort()

        nProc = 0
        for run in uruns:
            if nMax and nProc >= nMax:
                break

            if run == self.refrun:
                continue

            self.log.log(self.log.INFO, "") # spacer
            self.log.log(self.log.INFO, "RUNNING %d vs. %d" % (run, self.refrun))

            #if run != 4136:
            #    continue
            
            idxs = num.where(runs == run)[0]
            if len(idxs) == 0:
                continue

            runMatches = []
            for idx in idxs:
                runMatches.append(matches[idx])

            nloaded = 0
            for match in runMatches:
                match.loadfpC()
                match.loadfpM()
                match.loadWcs()
                match.loadPsf()
                match.loadCalib()
                if (match.fpC and match.fpM and 
                    match.wcs and match.psf and 
                    match.gain and match.calib and 
                    match.calib.getFluxMag0() > self.config.minFluxMag0):
                    nloaded += 1
    
            if nloaded != len(runMatches):
                self.log.log(self.log.INFO, "Not able to load all images, skipping to next run")
                continue
            else:
                self.log.log(self.log.INFO, "OK")

            stitch = self.stitchMatches(runMatches)

            # Keep Wcs and Calib of first image
            exp = afwImage.ExposureF(stitch, runMatches[0].wcs)
            exp.setPsf(self.refPsf)
            exp.setCalib(runMatches[0].calib)

            if self.config.writeFits:
                exp.writeFits(os.path.join(self.config.outputPath, 
                                           "psfmatch-%06d-%s%d-%04d-r%06d.fits" % 
                                           (self.refrun, self.filt, self.camcol, self.field, run)))

            # Do need to keep this
            self.warpedExp[run] = self.warper.warpExposure(self.refExp.getWcs(), 
                                                           exp, 
                                                           destBBox = self.refExp.getBBox(afwImage.PARENT))

            # Do after warping, since it loses it in warping
            self.warpedExp[run].setPsf(self.refPsf)

            if self.config.writeFits:
                self.warpedExp[run].writeFits(os.path.join(self.config.outputPath, 
                                                           "warp-%06d-%s%d-%04d-r%06d.fits" % 
                                                           (self.refrun, self.filt, self.camcol, self.field, run)))

            nProc += 1


    @pipeBase.timeMethod
    def stitchMatches(self, runMatches, overlap = 128, testme = 0):
        # Stitching together neighboring images from the matching run
        nloaded = len(runMatches)
        width   = runMatches[0].fpC.getWidth()
        height  = runMatches[0].fpC.getHeight() * nloaded - overlap * (nloaded - 1) + testme * nloaded
        stitch  = afwImage.MaskedImageF(width, height)

        for i in range(len(runMatches)):
            match  = runMatches[i]
            matchExp = match.createExp()
            # Psf match before stitching!
            psfmatchedExp = self.psfMatcher.run(matchExp, self.refPsf).psfMatchedExposure
            
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
                stitch.getImage().getArray()[symin:symax,:]    = psfmatchedExp.getMaskedImage().getImage().getArray()[iymin:iymax,:]
                stitch.getMask().getArray()[symin:symax,:]     = psfmatchedExp.getMaskedImage().getMask().getArray()[iymin:iymax,:]
                stitch.getVariance().getArray()[symin:symax,:] = psfmatchedExp.getMaskedImage().getVariance().getArray()[iymin:iymax,:]
            except:
                import pdb; pdb.set_trace()
                
            # Clear up memory
            match.fpC = None
            match.fpM = None
        
        return stitch

    @pipeBase.timeMethod   
    def matchBackgrounds(self):
        """Puts images on a common zeropoint; then background matches
        them; saves results in self.bgMatchedExp after checking some
        quality flags"""

        # IMPORTANT DETAIL : match zeropoints before matching backgrounds!
        self.refExp = self.coadder.normalizeForCoadd(self.refExp)
        for run in self.warpedExp.keys():
            self.warpedExp[run] = self.coadder.normalizeForCoadd(self.warpedExp[run])
        # IMPORTANT DETAIL : match zeropoints before matching b0ackgrounds!

        if self.config.writeFits:
            self.refExp.writeFits(os.path.join(self.config.outputPath, 
                                               "exp-%06d-%s%d-%04d.fits" % 
                                               (self.refrun, self.filt, self.camcol, self.field)))

        refMask   = self.refExp.getMaskedImage().getMask().getArray()
        refArr    = self.refExp.getMaskedImage().getImage().getArray()

        runsToMatch = self.warpedExp.keys()
        expsToMatch = self.warpedExp.values()
        
        skyMask  = num.sum(num.array([x.getMaskedImage().getMask().getArray() for x in expsToMatch]), 0)
        skyArr   = num.array([x.getMaskedImage().getImage().getArray() for x in expsToMatch])
        Nim      = len(skyArr)

        # Find all unmasked (sky) pixels
        idx = num.where((refMask + skyMask) == 0)

        width  = self.refExp.getMaskedImage().getWidth()
        height = self.refExp.getMaskedImage().getHeight()
        nbinx  = width  // self.config.backgroundBinsize
        nbiny  = height // self.config.backgroundBinsize

        bgX  = num.zeros((nbinx*nbiny, Nim)) # coord
        bgY  = num.zeros((nbinx*nbiny, Nim)) # coord
        bgZ  = num.zeros((nbinx*nbiny, Nim)) # value
        bgdZ = num.zeros((nbinx*nbiny, Nim)) # unc

        for biny in range(nbiny):
            ymin = biny * self.config.backgroundBinsize
            ymax = min((biny + 1) * self.config.backgroundBinsize, height)
            idxy = num.where( (idx[0] >= ymin) & (idx[0] < ymax) )[0]

            for binx in range(nbinx):
                xmin   = binx * self.config.backgroundBinsize
                xmax   = min((binx + 1) * self.config.backgroundBinsize, width)
                idxx   = num.where( (idx[1] >= xmin) & (idx[1] < xmax) )[0]
                inreg  = num.intersect1d(idxx, idxy)

                Aij    = num.zeros((Nim, Nim))
                Eij    = num.ones((Nim, Nim))

                area0 = refArr[idx[0][inreg],idx[1][inreg]]
                for i in range(Nim):
                    areai     = skyArr[i][idx[0][inreg],idx[1][inreg]]
                    area      = area0 - areai
                    
                    bgX [binx + biny * nbinx, i] = 0.5 * (xmin + xmax)
                    bgY [binx + biny * nbinx, i] = 0.5 * (ymin + ymax) 
                    bgZ [binx + biny * nbinx, i] = num.mean(area) 
                    bgdZ[binx + biny * nbinx, i] = num.std(area) 

        # Function for each image
        for i in range(Nim):

            # Function for this image
            bbox  = afwGeom.Box2D(self.refExp.getMaskedImage().getBBox())
            poly  = afwMath.Chebyshev1Function2D(self.config.backgroundOrder, bbox)          
            terms = list(poly.getParameters())

            Nall  = nbiny * nbinx
            Ncell = num.sum(num.isfinite(bgZ[:,i]))
            Nterm = len(terms)

            # Mx = b; solve for x
            m  = num.zeros((Ncell, Nterm))
            b  = num.zeros((Ncell))
            iv = num.zeros((Ncell))

            # One constraint for each cell
            nc = 0
            for na in range(Nall):
                if not num.isfinite(bgZ[:,i][na]):
                    continue

                for nt in range(Nterm):
                    terms[nt] = 1.0
                    poly.setParameters(terms)
                    m[nc, nt] = poly(bgX[:,i][na], bgY[:,i][na])
                    terms[nt] = 0.0
                b[nc]  = bgZ[:,i][na]
                iv[nc] = 1.0 / (bgdZ[:,i][na])**2
                nc += 1
            #import pdb; pdb.set_trace()

            M    = num.dot(num.dot(m.T, num.diag(iv)), m)
            B    = num.dot(num.dot(m.T, num.diag(iv)), b)
            Minv = num.linalg.inv(M)
            Soln = num.dot(Minv, B)
            poly.setParameters(Soln)

            run = runsToMatch[i]
            exp = expsToMatch[i]
            im  = exp.getMaskedImage()
            im += poly


            if self.config.writeFits:
                exp.writeFits(os.path.join(self.config.outputPath, 
                                           "match-%06d-%s%d-%04d-r%06d.fits" % 
                                           (self.refrun, self.filt, self.camcol, self.field, run)))

            # DEBUGGING INFO
            tmp  = afwImage.MaskedImageF(im, True)
            tmp -= self.refExp.getMaskedImage()

            if self.config.writeFits:
                tmp.writeFits(os.path.join(self.config.outputPath, 
                                           "diff-%06d-%s%d-%04d-r%06d.fits" % 
                                           (self.refrun, self.filt, self.camcol, self.field, run)))
            
            # Lets see some stats!
            area = tmp.getImage().getArray()[idx]
            self.log.log(self.log.INFO, "Diff BG %06d: mean=%0.3g med=%0.3g std=%0.3g npts=%d" % (
                    run, num.mean(area), num.median(area), num.std(area), len(area))
            )

            if num.std(area) < self.config.maxBgRms:
                self.bgMatchedExp[run] = exp

            self.warpedExp[run] = None # Clear memory


    @pipeBase.timeMethod
    def createCoadd(self):
        coaddExp = self.coadder.run(self.refExp, self.bgMatchedExp.values())
        coaddExp.setPsf(self.refPsf)
        coaddExp.writeFits(os.path.join(self.config.outputPath, 
                                        "coadd-%06d-%s%d-%04d.fits" % 
                                        (self.refrun, self.filt, self.camcol, self.field)))
        

                    #Aij[0][i] = num.mean(area)
                    #Eij[0][i] = num.std(area)
                    #Aij[i][0] = -1 * Aij[0][i]
                    #Eij[i][0] = +1 * Eij[0][i]

# NN2 BELOW, NOT YET WORKING            
        
#                    print binx, biny, i, num.mean(area0), num.mean(areai), num.mean(area), num.mean(area0)-num.mean(areai)
#                    bgs1X.append(num.mean(area0) - num.mean(areai))
#
#                    for j in range(i+1, N):
#                        areaj     = skyArr[j-1][idx[0][inreg],idx[1][inreg]]
#                        area      = areai - areaj
#                        Aij[i][j] = num.mean(area)
#                        Eij[i][j] = num.std(area)
#                        Aij[j][i] = -1 * Aij[i][j]
#                        Eij[j][i] = +1 * Eij[i][j]
#                        
#                Eij2     = Eij**2
#                InvEhat2 = 2.0 / (N * (N - 1.0)) / num.sum(num.triu(Eij2)) # NOTE TRIU DOES INCLUDE DIAGONALS!!!
#                InvEhat2 = 1.0
#
#                Cij   = -1. / Eij2
#                for i in range(N): Cij[i][i] = 0.0
#                Cij  += num.diag(1.0 / num.sum(Eij2, 0))
#                Cij  += InvEhat2
#                Cinv  = num.linalg.inv(Cij)
#                AEij  = Aij / Eij2
#                for i in range(N): AEij[i][i] = 0.0
#                bgs   = num.sum(num.dot(Cinv, AEij), 0)
#                unc   = num.sqrt(num.diagonal(Cinv))
#                
#                import pdb; pdb.set_trace()
#
#                print num.sum(bgs)
#                bgs  -= bgs[0]
#                print bgs
#                print
#                # Note that these numbers are differences between
#                # image reference in a lightcurve sense.  So negative
#                # numbers means the background is fainter.  So you
#                # need to subtract the numbers from the images to get
#                # them to match.
#                bgs *= -1
#
#                bgs1.append(bgs1X)
#                bgs2.append(bgs)
#
#        bgs1     = num.array(bgs1)
#        bgs2     = num.array(bgs2)
#        offsets1 = num.mean(bgs1, 0)
#        offsets2 = num.mean(bgs2, 0)
#
#        # debugging!
#        self.refExp[field].writeFits("/tmp/refexp.fits")
#        for i in range(1, len(offsets1)):
#            comp1  = afwImage.MaskedImageF(self.warpedExp[field].values()[i-1].getMaskedImage(), True)
#            comp1 += offsets1[i]
#            comp1.writeFits("/tmp/comp_avg_%d.fits" % (i))
#
#            comp2  = afwImage.MaskedImageF(self.warpedExp[field].values()[i-1].getMaskedImage(), True)
#            comp2 += offsets2[i]
#            comp2.writeFits("/tmp/comp_nn2_%d.fits" % (i))
#
#        import pdb; pdb.set_trace()
        

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

def getpsField(run, rerun, filt, camcol, field):
    fname = os.path.join(rootdir, str(run), str(rerun), "objcs", str(camcol), "psField-%06d-%d-%04d.fit" % (run, camcol, field))
    print fname
    if os.path.isfile(fname):
        return convertpsField(fname, filt)
    return None

def gettsField(run, rerun, filt, camcol, field):
    fname = os.path.join(rootdir, str(run), str(rerun), "calibChunks", str(camcol), "tsField-%06d-%d-%d-%04d.fit" % (run, camcol, rerun, field))
    print fname
    if os.path.isfile(fname):
        return converttsField(fname, filt)
    return None, None


if __name__ == '__main__':
    refrun  = int(sys.argv[1])
    camcol  = int(sys.argv[2])
    filt    = sys.argv[3]
    field   = int(sys.argv[4])
    nMax    = int(sys.argv[5])

    matcher = MatchBackgrounds(refrun, 40, camcol, filt, field)
    if True:
        matcher.run(nMax = nMax)
        sys.exit(1)

    else:
        # If you have stuff in the output dir, use it...
        matcher.warpedExp[fields[0]] = {}
        matcher.bgMatchedExp[fields[0]] = {}
        tmpdir  = "/tmp/"
        nmax    = 100
        nfound  = 0
        for f in os.listdir(tmpdir):
            if f.startswith("exp-%06d" % (refrun)):
                refExp = afwImage.ExposureF(os.path.join(tmpdir, f))
                matcher.refExp[fields[0]] = refExp
            elif f.startswith("warp-%06d" % (refrun)) and nfound < nmax:
                skyExp = afwImage.ExposureF(os.path.join(tmpdir, f))
                skyrun  = int(re.sub("r", "", f.split("-")[4].split(".")[0]))
                matcher.warpedExp[fields[0]][skyrun] = skyExp
                nfound += 1
        matcher.matchBackgrounds(fields[0])
    
