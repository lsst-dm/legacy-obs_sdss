import numpy as np
import os
from lsst.obs.sdss.yanny import yanny as Yanny

class SdssCameraState(Yanny):
    _filters = dict(u = 1, g = 2, r = 3, i = 4, z = 5)

    def __init__(self, opDir, opConfig, opECalib):
        self._ECalib = Yanny(os.path.join(opDir, opECalib))["ECALIB"]
        self._CcdConfig = Yanny(os.path.join(opDir, opConfig))["CCDCONFIG"]

    def _splitCcd(self, ccdName):
        filter, camCol = list(ccdName)

        return filter, int(camCol)

    def _getCamRow(self, filter):
        return SdssCameraState._filters[filter]

    def getCcdIndex(self, ECALIB, ccdName):
        """Return the index for the given ccd (e.g. g1) into the arrays returned by a Yanny object"""
        filter, camCol = self._splitCcd(ccdName)
        camRow = self._getCamRow(filter)

        me = np.where(np.logical_and(np.equal(ECALIB["camCol"], camCol), np.equal(ECALIB["camRow"], camRow)))
        if len(me) != 1:
            raise RuntimeError("Unable to lookup index for ccd %s" % ccdName)

        return me[0]

    def getEParams(self, ccdName):
        """Return a pair of ampId dict of electronic params for both amps of a named CCD (e.g. z4)"""
        ECALIB = self._ECalib
        me = self.getCcdIndex(ECALIB, ccdName)

        eparams = []
        for i in range(4):
            if int(self._CcdConfig["amp%d" % i][me]):
                gain = ECALIB["gain%d" % i][me]
                readNoise = ECALIB["readNoiseDN%d" % i][me]
                fullWell = ECALIB["fullWellDN%d" % i][me]

                eparams.append((i, {'gain':gain, 'readNoise':readNoise, 'fullWell':fullWell}))

        if len(eparams) == 1:
            eparams.append((1, eparams[0][1]))

        return eparams

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

if __name__ == "__main__":
    sc = SdssCameraState("/lsst7/stripe82/dr7/opfiles", "opConfig-50000.par", "opECalib-50000.par")
    print [(i, ep.getGain(), ep.getReadNoise(), ep.getSaturationLevel()) for i, ep in sc.getEParams("g2")]
