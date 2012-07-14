# Currently set up to use Mario's test database
from lsst.obs.sdss.forcedPhot import SdssReferencesTask, SdssCoaddReferencesTask, SdssCoaddFileReferencesTask
#root.references.retarget(SdssReferencesTask)
root.references.retarget(SdssCoaddReferencesTask)
root.references.dbUrl = "mysql://lsst10.ncsa.uiuc.edu:3306/"

# Copy database columns over
for col in ("refFlux", "refFlux.err"):
    root.copyColumns[col] = col

root.measurement.algorithms.names += ("centroid.record",)
