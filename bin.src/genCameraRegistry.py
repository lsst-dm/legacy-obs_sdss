from lsst.obs.sdss import makeCamera
import lsst.utils
import os

path = lsst.utils.getPackageDir('obs_sdss')
repoPath = os.path.join(path, 'description', 'camera')
if os.path.exists(repoPath):
    raise RuntimeError("Path, $s, exists.  Will not overwrite"%(repoPath))
else:
    os.makedirs(repoPath)

makeCamera.makeCamera(outputDir=repoPath)
