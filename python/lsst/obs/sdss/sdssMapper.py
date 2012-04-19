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
from lsst.daf.butlerUtils import CameraMapper
from lsst.obs.sdss.convertfpM import convertfpM
from lsst.obs.sdss.convertpsField import convertpsField
from lsst.obs.sdss.convertasTrans import convertasTrans

# Solely to get boost serialization registrations for Measurement subclasses
import lsst.meas.algorithms as measAlgo

class SdssMapper(CameraMapper):
    def __init__(self, inputPolicy=None, **kwargs):
        policyFile = pexPolicy.DefaultPolicyFile("obs_sdss", "SdssMapper.paf", "policy")
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

    def _computeCcdExposureId(self, dataId):
        """Compute the 64-bit (long) identifier for a CCD exposure.

        @param dataId (dict) Data identifier with run, rerun, band, camcol, frame
        """

        return ((long(run) \
                * 10 + self.filterIdMap[band]) \
                * 10 + camcol) \
                * 10000 + frame

    def _setCcdExposureId(self, propertyList, dataId):
        propertyList.set("Computed_ccdExposureId", self._computeCcdExposureId(dataId))
        return propertyList

    def _standardizeExposure(self, mapping, item, dataId, filter=True,
            trimmed=True):

        """Default standardization function for images.
        @param mapping (lsst.daf.butlerUtils.Mapping)
        @param[in,out] item (lsst.afw.image.Exposure)
        @param dataId (dict) Dataset identifier
        @param filter (bool) Set filter?
        @param trimmed (bool) Should detector be marked as trimmed?
        @return (lsst.afw.image.Exposure) the standardized Exposure"""

        if (re.search(r'Exposure', mapping.python) and re.search(r'Image',mapping.persistable)):
            item = exposureFromImage(item)
        return item

###############################################################################

    def bypass_fpM(self, datasetType, pythonType, location, dataId):
        return convertfpM(location.getLocations()[0])

    def bypass_psField(self, datasetType, pythonType, location, dataId):
        return convertpsField(location.getLocations()[0], dataId['band'])

    def bypass_asTrans(self, datasetType, pythonType, location, dataId):
        return convertasTrans(location.getLocations()[0], dataId['band'],
                dataId['camcol'], dataid['frame'])

###############################################################################

    def _addSources(self, dataId):
        """Generic 'add' function to add ampExposureId and filterId"""
        # Note that sources are identified by what is called an ampExposureId,
        # but in this case all we have is a CCD.
        ampExposureId = self._computeCcdExposureId(dataId)
        filterId = self.filterIdMap[pathId['band']]
        ad = dict(ampExposureId=ampExposureId, filterId=filterId)
        if self.doFootprints:
            ad["doFootprints"] = True
        return ad

    def _addSkytile(self, dataId):
        """Generic 'add' function to add skyTileId"""
        return {"skyTileId": dataId['skyTile']}

for dsType in ("icSrc", "src"):
    setattr(SdssMapper, "add_" + dsType, SdssMapper._addSources)
for dsType in ("source", "badSource", "invalidSource", "object", "badObject"):
    setattr(SdssMapper, "add_" + dsType, SdssMapper._addSkytile)

###############################################################################


for dsType in ("fpC", "fpM", "calexp"):
    setattr(SdssMapper, "std_" + dsType + "_md",
            lambda self, item, dataId: self._setCcdExposureId(item))
