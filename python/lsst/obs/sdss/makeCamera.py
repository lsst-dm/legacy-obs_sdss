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

import lsst.afw.geom as afwGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
from lsst.afw.cameraGeom import makeCameraFromCatalogs, CameraConfig, DetectorConfig, SCIENCE, PIXELS, PUPIL, FOCAL_PLANE
import lsst.afw.table as afwTable
from lsst.obs.sdss.convertOpECalib import SdssCameraState

#
# Make an Amp
#
def addAmp(ampCatalog, i, eparams):
    """ Add an amplifier to an AmpInfoCatalog

    @param i which amplifier? (i == 0 ? left : right)
    @param eparams Electronic parameters for this amplifier
    """
    #
    # Layout of active and overclock pixels in the as-readout data  The layout is:
    #     Amp0 || extended | overclock | data || data | overclock | extended || Amp1
    # for each row; all rows are identical in drift-scan data
    #
    height = 1361      # number of rows in a frame
    width = 1024       # number of data pixels read out through one amplifier
    nExtended = 8             # number of pixels in the extended register
    nOverclock = 32           # number of (horizontal) overclock pixels 
    #
    # Construct the needed bounding boxes given that geometrical information.
    #
    # Note that all the offsets are relative to the origin of this amp, not to its eventual
    # position in the CCD
    #
    record = ampCatalog.addNew()
    xtot = width + nExtended + nOverclock
    allPixels = afwGeom.BoxI(afwGeom.PointI(0, 0), afwGeom.ExtentI(xtot, height))
    biasSec = afwGeom.BoxI(afwGeom.PointI(nExtended if i == 0 else width, 0),
                           afwGeom.ExtentI(nOverclock, height)) 
    dataSec = afwGeom.BoxI(afwGeom.PointI(nExtended + nOverclock if i == 0 else 0, 0),
                           afwGeom.ExtentI(width, height))
    emptyBox = afwGeom.BoxI()
    bbox = afwGeom.BoxI(afwGeom.PointI(0, 0), afwGeom.ExtentI(width, height))
    bbox.shift(afwGeom.Extent2I(width*i, 0))

    shiftp = afwGeom.Extent2I(xtot*i, 0)
    allPixels.shift(shiftp)
    biasSec.shift(shiftp)
    dataSec.shift(shiftp)

    record.setBBox(bbox)
    record.setRawXYOffset(afwGeom.ExtentI(0,0))
    record.setName('left' if i==0 else 'right')
    record.setReadoutCorner(afwTable.LL if i==0 else afwTable.LR)    
    record.setGain(eparams['gain'])
    record.setReadNoise(eparams['readNoise'])
    record.setSaturation(int(eparams['fullWell']))
    record.setLinearityType('Proportional')
    record.setLinearityCoeffs([1.,])
    record.setHasRawInfo(True)
    record.setRawFlipX(False)
    record.setRawFlipY(False)
    record.setRawBBox(allPixels)
    record.setRawDataBBox(dataSec)
    record.setRawHorizontalOverscanBBox(biasSec)
    record.setRawVerticalOverscanBBox(emptyBox)
    record.setRawPrescanBBox(emptyBox)

#
# Make a Ccd out of 2 Amps
#
def makeCcd(ccdName, ccdId, offsetPoint):
    opDir = os.path.join(os.environ['OBS_SDSS_DIR'], "etc")
    sc = SdssCameraState(opDir, "opConfig-50000.par", "opECalib-50000.par")
    eparams = sc.getEParams(ccdName)
    width = 1024*2
    height = 1361

    pixelSize = 24e-3                   # pixel size in mm
    schema = afwTable.AmpInfoTable.makeMinimalSchema()
    ampCatalog = afwTable.AmpInfoCatalog(schema)
    for i in range(2):
        addAmp(ampCatalog, i, eparams[i][1])
    detConfig = DetectorConfig()
    detConfig.name = ccdName
    detConfig.id = ccdId
    detConfig.bbox_x0 = 0
    detConfig.bbox_y0 = 0
    detConfig.bbox_x1 = width - 1
    detConfig.bbox_y1 = height - 1
    detConfig.serial = ccdName
    detConfig.detectorType = SCIENCE
    detConfig.offset_x = offsetPoint.getX()
    detConfig.offset_y = offsetPoint.getY()
    detConfig.refpos_x = (width-1)/2.
    detConfig.refpos_y = (height-1)/2.
    detConfig.yawDeg = 0.
    detConfig.pitchDeg = 0.
    detConfig.rollDeg = 0.
    detConfig.pixelSize_x = pixelSize
    detConfig.pixelSize_y = pixelSize
    detConfig.transposeDetector = False
    detConfig.transformDict.nativeSys = PIXELS.getSysName()
    return {'ccdConfig':detConfig, 'ampInfo':ampCatalog}

#
# Make a Camera out of 6 dewars and 5 chips per dewar
#
def makeCamera(name="SDSS", outputDir=None):
    camConfig = CameraConfig()
    camConfig.name = name
    camConfig.detectorList = {}
    camConfig.plateScale = 16.5 # arcsec/mm
    pScaleRad = afwGeom.arcsecToRad(camConfig.plateScale)
    radialDistortCoeffs = [0.0, 1.0/pScaleRad]
    tConfig = afwGeom.TransformConfig()
    tConfig.transform.name = 'inverted'
    radialClass = afwGeom.xyTransformRegistry['radial']
    tConfig.transform.active.transform.retarget(radialClass)
    tConfig.transform.active.transform.coeffs = radialDistortCoeffs
    tmc = afwGeom.TransformMapConfig()
    tmc.nativeSys = FOCAL_PLANE.getSysName()
    tmc.transforms = {PUPIL.getSysName():tConfig}
    camConfig.transformDict = tmc


    ccdId = 0
    ampInfoCatDict = {}
    for i in range(6):
        dewarName = str(i+1)
        filters = "riuzg"
        for j, c in enumerate(reversed(filters)):
            ccdName = "%s%s" % (c, dewarName)
            offsetPoint = afwGeom.Point2D(25.4*2.5*(2.5-i), 25.4*2.1*(2.0 - j))
            ccdInfo = makeCcd(ccdName, ccdId, offsetPoint)
            ampInfoCatDict[ccdName] = ccdInfo['ampInfo']
            camConfig.detectorList[ccdId] = ccdInfo['ccdConfig']
            ccdId += 1
    if outputDir is not None:
        camConfig.save(os.path.join(outputDir, 'camera.py'))
        for k in ampInfoCatDict:
            ampInfoCatDict[k].writeFits(os.path.join(outputDir, "%s.fits"%(k)))
    return makeCameraFromCatalogs(camConfig, ampInfoCatDict)

#************************************************************************************************************
#
# Print a Ccd
#
def printCcd(title, ccd, trimmed=True, indent=""):
    print indent, title, "CCD: ", ccd.getName()
    if trimmed:
        allPixels = ccd.getBBox()
    else:
        allPixels = cameraGeomUtils.calcRawCcdBBox(ccd)
    print indent, "Total size: %dx%d" % (allPixels.getWidth(), allPixels.getHeight())
    for i, amp in enumerate(ccd):
        biasSec = amp.getRawHorizontalOverscanBBox() 
        dataSec = amp.getRawDataBBox() 

        print indent, "   Amp: %s gain: %g" % (amp.getName(),
                                               amp.getGain())

        print indent,"   bias sec: %dx%d+%d+%d" % (biasSec.getWidth(), biasSec.getHeight(),
                                                    biasSec.getMinX(), biasSec.getMinY())

        print indent, "   data sec: %dx%d+%d+%d" % (dataSec.getWidth(), dataSec.getHeight(),
                                                     dataSec.getMinX(), dataSec.getMinY())
        if i == 0:
            print

#
# Print a Camera
#
def printCamera(title, camera):
    print title, "Camera:", camera.getName()

    for det in camera:
        print "%s %dx%d centre (mm): %s" % \
            (det.getName(),
             det.getBBox().getWidth(), det.getBBox().getHeight(),
             det.getCenter(FOCAL_PLANE).getPoint())

#************************************************************************************************************

def main():
    camera = makeCamera("SDSS") 

    print
    printCamera("", camera) 

    ccd = camera["r1"]

    printCcd("Raw ", ccd, trimmed=False) 

    print
    printCcd("Trimmed ", ccd, trimmed=True) 

if __name__ == "__main__":
    main()

