# Currently set up to use Mario's test database
from lsst.obs.sdss.forcedPhot import SdssCoaddReferencesTask
root.references.retarget(SdssCoaddReferencesTask)

# Copy database columns over
for col in ("refFlux", "refFlux.err"):
    root.copyColumns[col] = col
root.measurement.algorithms.names.add("centroid.record")
