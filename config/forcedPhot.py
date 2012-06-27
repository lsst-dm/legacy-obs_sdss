# Currently set up to use Mario's test database
from lsst.obs.sdss import SdssReferencesTask
root.references.retarget(SdssReferencesTask)
root.references.dbUrl = "mysql://lsst10.ncsa.uiuc.edu:3306/"

# Copy database columns over
for col in ("refMag", "refMag.err"):
    root.copyColumns[col] = col

root.measurement.algorithms.names += ("centroid.record",)
