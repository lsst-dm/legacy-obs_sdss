"""
Utilities to use SDSS data with the cameraGeom utilities.  For example:

import lsst.daf.persistence as dafPersist
import lsst.afw.cameraGeom.utils as cgUtils
import lsst.obs.sdss.makeCamera as makeCamera
from lsst.obs.sdss import SdssMapper;

butler = dafPersist.ButlerFactory(mapper=SdssMapper(root="/lsst7/stripe82/dr7/runs")).create()

camera = makeCamera.makeCamera()
cgUtils.showCamera(camera, imageSource=SdssCcdImage(butler=butler, run=94,
field=101), field=1)

"""
import lsst.afw.cameraGeom.utils as cgUtils
import lsst.afw.geom as afwGeom

class SdssCcdImage(cgUtils.GetCcdImage):
    raise NotImplementedError("This will be re-implemented in Summer 2014")
    """A class to return an Image of a given SDSS Ccd by using the butler"""
    ''' 
    def __init__(self, butler, run, field, *args):
        """Initialise"""
        super(SdssCcdImage, self).__init__(*args)
        self.butler = butler
        self.run = run
        self.field = field

    def getImage(self, ccd, amp=None, imageFactory=None):
        """Return the image of the chip with cameraGeom.Id == id; if provided only read the given amp"""

        filter, camCol = list(ccd.getName()) # list splits the name
        camCol = 3                                 # XXX
        fpC = self.butler.get("fpC", dict(run=self.run, camcol=camCol,
            filter=filter, field=self.field))

        if amp:
            if amp.getId().getSerial() == 0:
                origin = afwGeom.PointI(0, 0)
            else:
                origin = afwGeom.PointI(1024, 0)

            fpC = fpC.Factory(fpC, afwGeom.BoxI(origin, afwGeom.ExtentI(1024, 1361)))

        return fpC
    '''
