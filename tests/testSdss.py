import lsst.daf.persistence as dafPersist
from lsst.obs.sdss import SdssMapper

butler = dafPersist.ButlerFactory(
        mapper=SdssMapper(root="/lsst7/stripe82/dr7/runs")).create()
exp = butler.get("corr", run=5754, camcol=3, frame=280, band="r", rerun=40)
print exp.getWidth(), exp.getHeight()
msk = butler.get("mask", run=5754, camcol=3, frame=280, band="r", rerun=40)
print msk.getWidth(), msk.getHeight()
# psf = butler.get("origPsf", run=5754, camcol=3, frame=280, band="r", rerun=40)
# do something with psf

butler2 = dafPersist.ButlerFactory(
        mapper=SdssMapper(root="/lsst7/stripe82/uw-coadds")).create()
exp = butler2.get("coadd", run=6383, camcol=3, frame=280, band="r")
print exp.getWidth(), exp.getHeight()
