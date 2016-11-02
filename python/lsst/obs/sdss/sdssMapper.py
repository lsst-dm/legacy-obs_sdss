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

import re

import lsst.pex.policy as pexPolicy
from lsst.obs.base import CameraMapper, exposureFromImage
from lsst.obs.sdss.convertfpM import convertfpM
from lsst.obs.sdss.convertpsField import convertpsField
from lsst.obs.sdss.convertasTrans import convertasTrans
from lsst.obs.sdss.converttsField import converttsField
import lsst.afw.image.utils as afwImageUtils

# Solely to get boost serialization registrations for Measurement subclasses
import lsst.meas.algorithms as measAlgo  # flake8: noqa


class SdssMapper(CameraMapper):
    packageName = 'obs_sdss'

    def __init__(self, inputPolicy=None, **kwargs):
        policyFile = pexPolicy.DefaultPolicyFile(self.packageName, "SdssMapper.paf", "policy")
        policy = pexPolicy.Policy(policyFile)

        self.doFootprints = False
        if inputPolicy is not None:
            for kw in inputPolicy.paramNames(True):
                if kw == "doFootprints":
                    self.doFootprints = True
                else:
                    kwargs[kw] = inputPolicy.get(kw)

        super(SdssMapper, self).__init__(policy, policyFile.getRepositoryPath(), **kwargs)
        # define filters?
        self.filterIdMap = dict(u=0, g=1, r=2, i=3, z=4)

        afwImageUtils.defineFilter('u', lambdaEff=380)
        afwImageUtils.defineFilter('g', lambdaEff=450)
        afwImageUtils.defineFilter('r', lambdaEff=600)
        afwImageUtils.defineFilter('i', lambdaEff=770)
        afwImageUtils.defineFilter('z', lambdaEff=900)

    def _computeCcdExposureId(self, dataId):
        """Compute the 64-bit (long) identifier for a CCD exposure.

        @param dataId (dict) Data identifier with run, rerun, filter, camcol, field
        """
        return ((long(dataId['run'])
                 * 10 + self.filterIdMap[dataId['filter']])
                * 10 + dataId['camcol']) \
            * 10000 + dataId['field']

    def _computeCoaddExposureId(self, dataId, singleFilter):
        """Compute the 64-bit (long) identifier for a coadd.

        @param dataId (dict)       Data identifier with tract and patch.
        @param singleFilter (bool) True means the desired ID is for a single-
                                   filter coadd, in which case dataId
                                   must contain filter.
        """
        tract = long(dataId['tract'])
        if tract < 0 or tract >= 128:
            raise RuntimeError('tract not in range [0,128)')
        patchX, patchY = map(int, dataId['patch'].split(','))
        for p in (patchX, patchY):
            if p < 0 or p >= 2**13:
                raise RuntimeError('patch component not in range [0, 8192)')
        id = (tract * 2**13 + patchX) * 2**13 + patchY
        if singleFilter:
            return id * 8 + self.filterIdMap[dataId['filter']]
        return id

    def _setCcdExposureId(self, propertyList, dataId):
        propertyList.set("Computed_ccdExposureId", self._computeCcdExposureId(dataId))
        return propertyList

    def _standardizeExposure(self, mapping, item, dataId, filter=True,
                             trimmed=True):
        """Default standardization function for images.
        @param mapping (lsst.obs.base.Mapping)
        @param[in,out] item (lsst.afw.image.Exposure)
        @param dataId (dict) Dataset identifier
        @param filter (bool) Set filter?
        @param trimmed (bool) Should detector be marked as trimmed?
        @return (lsst.afw.image.Exposure) the standardized Exposure"""

        if (re.search(r'Exposure', mapping.python) and re.search(r'Image', mapping.persistable)):
            item = exposureFromImage(item)
        return item

###############################################################################

    def bypass_fpM(self, datasetType, pythonType, location, dataId):
        return convertfpM(location.getLocations()[0])

    def bypass_psField(self, datasetType, pythonType, location, dataId):
        return convertpsField(location.getLocations()[0], dataId['filter'])

    def bypass_asTrans(self, datasetType, pythonType, location, dataId):
        return convertasTrans(location.getLocations()[0], dataId['filter'],
                              dataId['camcol'], dataId['field'])

    def bypass_tsField(self, datasetType, pythonType, location, dataId):
        return converttsField(location.getLocations()[0], dataId['filter'])

    def bypass_ccdExposureId(self, datasetType, pythonType, location, dataId):
        return self._computeCcdExposureId(dataId)

    def bypass_ccdExposureId_bits(self, datasetType, pythonType, location, dataId):
        return 38

    def bypass_goodSeeingCoaddId(self, datasetType, pythonType, location, dataId):
        return self._computeCoaddExposureId(dataId, True)

    def bypass_goodSeeingCoaddId_bits(self, datasetType, pythonType, location, dataId):
        return 1 + 7 + 13*2 + 3

    # Deep coadds use tract, patch, and filter just like good-seeing coadds
    bypass_deepCoaddId = bypass_goodSeeingCoaddId
    bypass_deepCoaddId_bits = bypass_goodSeeingCoaddId_bits

    def bypass_chiSquaredCoaddId(self, datasetType, pythonType, location, dataId):
        return self._computeCoaddExposureId(dataId, False)

    def bypass_chiSquaredCoaddId_bits(self, datasetType, pythonType, location, dataId):
        return 1 + 7 + 13*2

    # Keith coadds use run, camcol, field, filter just like CCD exposures
    bypass_keithCoaddId = bypass_ccdExposureId
    bypass_keithCoaddId_bits = bypass_ccdExposureId_bits


###############################################################################


for dsType in ("fpC", "fpM", "calexp"):
    setattr(SdssMapper, "std_" + dsType + "_md",
            lambda self, item, dataId: self._setCcdExposureId(item, dataId))
