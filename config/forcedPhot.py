# Currently set up to use Mario's test database
from lsst.obs.sdss import TestSdssReferencesTask
root.references.retarget(TestSdssReferencesTask)
root.references.dbUrl = "mysql://lsst10.ncsa.uiuc.edu:3306/"
root.references.dbName = "juric_DR7_stripe82"

# Copy database columns over
for col in ("refMag", "refMag.err"):
    root.copyColumns[col] = col
