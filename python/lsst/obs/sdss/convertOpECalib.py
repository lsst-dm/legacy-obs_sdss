import numpy as np
import os
from lsst.obs.sdss.yanny import yanny as Yanny

class SdssConfig(Yanny):
    _bands = dict(u = 0, g = 1, r = 2, i = 3, z = 4)

    def __init__(self, opDir, opConfig, opECalib):
        self._ECalib = Yanny(os.path.join(opDir, opECalib))["ECALIB"]
        self._CcdConfig = Yanny(os.path.join(opDir, opConfig))["CCDCONFIG"]

    def _splitCcd(self, ccdName):
        band, camCol = list(ccdName)

        return band, int(camCol)

    def _getCamRow(self, band):
        return SdssConfig._bands[band]

    def getCcdIndex(self, ECALIB, ccdName):
        """Return the index for the given ccd (e.g. g1) into the arrays returned by a Yanny object"""
        band, camCol = self._splitCcd(ccdName)
        camRow = self._getCamRow(band)

        me = np.where(np.logical_and(np.equal(ECALIB["camCol"], camCol), np.equal(ECALIB["camRow"], camRow)))
        if len(me) != 1:
            raise RuntimeError("Unable to lookup index for ccd %s" % ccdName)

        return me[0]

    def getGain(self, ccdName):
        """Return an array of the (amp, gain) values for the given ccd (e.g. r2)"""
        ECALIB = self._ECalib
        me = self.getCcdIndex(ECALIB, ccdName)

        gains = []
        for i in range(4):
            if int(self._CcdConfig["amp%d" % i][me]):
                gains.append((i, ECALIB["gain%d" % i][me]))

        return gains
        
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

if __name__ == "__main__":
    sc = SdssConfig("/lsst7/stripe82/dr7/opfiles", "opConfig-50000.par", "opECalib-50000.par")
    print sc.getGain("g2")
