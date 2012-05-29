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
import lsst.afw.geom as afwGeom
import lsst.afw.cameraGeom as cameraGeom
from lsst.obs.sdss.convertOpECalib import SdssCameraState

#
# Make an Amp
#
def makeAmp(i, eparams):
    """ Make an amplifier

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
    allPixels = afwGeom.BoxI(afwGeom.PointI(0, 0), afwGeom.ExtentI(width + nExtended + nOverclock, height))
    biasSec = afwGeom.BoxI(afwGeom.PointI(nExtended if i == 0 else width, 0),
                           afwGeom.ExtentI(nOverclock, height)) 
    dataSec = afwGeom.BoxI(afwGeom.PointI(nExtended + nOverclock if i == 0 else 0, 0),
                           afwGeom.ExtentI(width, height))
    
    return cameraGeom.Amp(cameraGeom.Id(i), allPixels, biasSec, dataSec, eparams)

#
# Make a Ccd out of 2 Amps
#
def makeCcd(ccdName):
    sc = SdssCameraState("/lsst7/stripe82/dr7/opfiles", "opConfig-50000.par", "opECalib-50000.par")
    eparams = sc.getEParams(ccdName)

    pixelSize = 24e-3                   # pixel size in mm
    ccd = cameraGeom.Ccd(cameraGeom.Id(ccdName), pixelSize)

    for i in range(2):
        ccd.addAmp(i, 0, makeAmp(i, eparams[i][1]))

    return ccd 

#
# Make a Raft (== SDSS dewar) out of 5 Ccds
#
def makeRaft(raftName):
    dewar = cameraGeom.Raft(cameraGeom.Id(str(raftName)), 1, 5)

    filters = "riuzg"
    for i, c in enumerate(reversed(filters)):
        ccdName = "%s%d" % (c, raftName)
        dewar.addDetector(afwGeom.PointI(0, i), cameraGeom.FpPoint(0.0, 25.4*2.1*(2.0 - i)),
                          cameraGeom.Orientation(0), makeCcd(ccdName))

    return dewar 

#
# Make a Camera out of 6 Rafts (==dewars)
#
def makeCamera(name="SDSS"):
    camera = cameraGeom.Camera(cameraGeom.Id(name), 6, 1)

    for i in range(6):
        dewarName = (i + 1)
        camera.addDetector(afwGeom.PointI(i, 0), cameraGeom.FpPoint(25.4*2.5*(2.5 - i), 0.0),
                           cameraGeom.Orientation(0), makeRaft(dewarName))

    return camera 

#************************************************************************************************************
#
# Print a Ccd
#
def printCcd(title, ccd, indent=""):
    print indent, title, "CCD: ", ccd.getId().getName()
    allPixels = ccd.getAllPixels() 
    print indent, "Total size: %dx%d" % (allPixels.getWidth(), allPixels.getHeight())
    for i, amp in enumerate(ccd):
        biasSec = amp.getBiasSec() 
        dataSec = amp.getDataSec() 

        print indent, "   Amp: %s gain: %g" % (amp.getId().getSerial(),
                                               amp.getElectronicParams().getGain())

        print indent,"   bias sec: %dx%d+%d+%d" % (biasSec.getWidth(), biasSec.getHeight(),
                                                    biasSec.getMinX(), biasSec.getMinY())

        print indent, "   data sec: %dx%d+%d+%d" % (dataSec.getWidth(), dataSec.getHeight(),
                                                     dataSec.getMinX(), dataSec.getMinY())
        if i == 0:
            print

#
# Print a Dewar
#
def printDewar(title, dewar, indent=""):
    print indent, title, "Dewar: ", dewar.getId().getName()

    for det in dewar:
        print indent, "%s %dx%d centre (mm): %s" % \
            (det.getId().getName(),
             det.getAllPixels(True).getWidth(), det.getAllPixels(True).getHeight(),
             det.getCenter().getMm())

#
# Print a Camera
#
def printCamera(title, camera):
    print title, "Camera:", camera.getId().getName()

    for raft in camera:
        printDewar("\n", cameraGeom.cast_Raft(raft), "    ") 

#************************************************************************************************************

def main():
    ccd = makeCcd("r1") 

    printCcd("Raw ", ccd) 

    ccd.setTrimmed(True) 

    print
    printCcd("Trimmed ", ccd) 
    #
    # The SDSS camera has 6 independent dewars, each with 5 CCDs mounted on a piece of invar
    #
    dewar = makeRaft(2)                 # dewar 2

    print
    printDewar("Single ", dewar) 
    #
    # On to the camera
    #
    camera = makeCamera("SDSS") 

    print
    printCamera("", camera) 

if __name__ == "__main__":
    main()

