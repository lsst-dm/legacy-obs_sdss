import lsst.afw.cameraGeom.cameraConfig
assert type(config) == lsst.afw.cameraGeom.cameraConfig.CameraConfig, 'config is of type %s.%s instead of lsst.afw.cameraGeom.cameraConfig.CameraConfig' % (
    type(config).__module__, type(config).__name__)
# Plate scale of the camera in arcsec/mm
config.plateScale = 16.5

# Name of native coordinate system
config.transformDict.nativeSys = 'FocalPlane'

config.transformDict.transforms = {}
config.transformDict.transforms['FieldAngle'] = lsst.afw.geom.transformConfig.TransformConfig()
config.transformDict.transforms['FieldAngle'].transform['multi'].transformDict = None
# x, y translation vector
config.transformDict.transforms['FieldAngle'].transform['affine'].translation = [0.0, 0.0]

# 2x2 linear matrix in the usual numpy order;
#             to rotate a vector by theta use: cos(theta), sin(theta), -sin(theta), cos(theta)
config.transformDict.transforms['FieldAngle'].transform['affine'].linear = [1.0, 0.0, 0.0, 1.0]

# Coefficients for the radial polynomial; coeff[0] must be 0
config.transformDict.transforms['FieldAngle'].transform['radial'].coeffs = None

config.transformDict.transforms['FieldAngle'].transform['inverted'].transform.retarget(target=lsst.afw.geom.transformRegistry['radial'])
# Coefficients for the radial polynomial; coeff[0] must be 0
config.transformDict.transforms['FieldAngle'].transform['inverted'].transform.coeffs = [0.0, 12500.89734830887]

config.transformDict.transforms['FieldAngle'].transform.name = 'inverted'
config.detectorList = {}
config.detectorList[0] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[0].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[0].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[0].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[0].bbox_x0 = 0

# Name of detector slot
config.detectorList[0].name = 'g1'

# Pixel size in the x dimension in mm
config.detectorList[0].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[0].transformDict.nativeSys = 'Pixels'

config.detectorList[0].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[0].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[0].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[0].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[0].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[0].offset_x = 158.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[0].offset_y = 106.67999999999999

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[0].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[0].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[0].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[0].serial = 'g1'

# pitch (rotation about y) of the detector in degrees
config.detectorList[0].pitchDeg = 0.0

# ID of detector slot
config.detectorList[0].id = 0

config.detectorList[1] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[1].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[1].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[1].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[1].bbox_x0 = 0

# Name of detector slot
config.detectorList[1].name = 'z1'

# Pixel size in the x dimension in mm
config.detectorList[1].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[1].transformDict.nativeSys = 'Pixels'

config.detectorList[1].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[1].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[1].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[1].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[1].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[1].offset_x = 158.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[1].offset_y = 53.339999999999996

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[1].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[1].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[1].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[1].serial = 'z1'

# pitch (rotation about y) of the detector in degrees
config.detectorList[1].pitchDeg = 0.0

# ID of detector slot
config.detectorList[1].id = 1

config.detectorList[2] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[2].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[2].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[2].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[2].bbox_x0 = 0

# Name of detector slot
config.detectorList[2].name = 'u1'

# Pixel size in the x dimension in mm
config.detectorList[2].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[2].transformDict.nativeSys = 'Pixels'

config.detectorList[2].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[2].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[2].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[2].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[2].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[2].offset_x = 158.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[2].offset_y = 0.0

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[2].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[2].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[2].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[2].serial = 'u1'

# pitch (rotation about y) of the detector in degrees
config.detectorList[2].pitchDeg = 0.0

# ID of detector slot
config.detectorList[2].id = 2

config.detectorList[3] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[3].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[3].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[3].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[3].bbox_x0 = 0

# Name of detector slot
config.detectorList[3].name = 'i1'

# Pixel size in the x dimension in mm
config.detectorList[3].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[3].transformDict.nativeSys = 'Pixels'

config.detectorList[3].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[3].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[3].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[3].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[3].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[3].offset_x = 158.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[3].offset_y = -53.339999999999996

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[3].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[3].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[3].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[3].serial = 'i1'

# pitch (rotation about y) of the detector in degrees
config.detectorList[3].pitchDeg = 0.0

# ID of detector slot
config.detectorList[3].id = 3

config.detectorList[4] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[4].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[4].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[4].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[4].bbox_x0 = 0

# Name of detector slot
config.detectorList[4].name = 'r1'

# Pixel size in the x dimension in mm
config.detectorList[4].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[4].transformDict.nativeSys = 'Pixels'

config.detectorList[4].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[4].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[4].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[4].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[4].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[4].offset_x = 158.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[4].offset_y = -106.67999999999999

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[4].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[4].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[4].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[4].serial = 'r1'

# pitch (rotation about y) of the detector in degrees
config.detectorList[4].pitchDeg = 0.0

# ID of detector slot
config.detectorList[4].id = 4

config.detectorList[5] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[5].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[5].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[5].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[5].bbox_x0 = 0

# Name of detector slot
config.detectorList[5].name = 'g2'

# Pixel size in the x dimension in mm
config.detectorList[5].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[5].transformDict.nativeSys = 'Pixels'

config.detectorList[5].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[5].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[5].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[5].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[5].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[5].offset_x = 95.25

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[5].offset_y = 106.67999999999999

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[5].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[5].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[5].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[5].serial = 'g2'

# pitch (rotation about y) of the detector in degrees
config.detectorList[5].pitchDeg = 0.0

# ID of detector slot
config.detectorList[5].id = 5

config.detectorList[6] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[6].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[6].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[6].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[6].bbox_x0 = 0

# Name of detector slot
config.detectorList[6].name = 'z2'

# Pixel size in the x dimension in mm
config.detectorList[6].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[6].transformDict.nativeSys = 'Pixels'

config.detectorList[6].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[6].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[6].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[6].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[6].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[6].offset_x = 95.25

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[6].offset_y = 53.339999999999996

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[6].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[6].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[6].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[6].serial = 'z2'

# pitch (rotation about y) of the detector in degrees
config.detectorList[6].pitchDeg = 0.0

# ID of detector slot
config.detectorList[6].id = 6

config.detectorList[7] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[7].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[7].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[7].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[7].bbox_x0 = 0

# Name of detector slot
config.detectorList[7].name = 'u2'

# Pixel size in the x dimension in mm
config.detectorList[7].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[7].transformDict.nativeSys = 'Pixels'

config.detectorList[7].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[7].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[7].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[7].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[7].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[7].offset_x = 95.25

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[7].offset_y = 0.0

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[7].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[7].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[7].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[7].serial = 'u2'

# pitch (rotation about y) of the detector in degrees
config.detectorList[7].pitchDeg = 0.0

# ID of detector slot
config.detectorList[7].id = 7

config.detectorList[8] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[8].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[8].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[8].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[8].bbox_x0 = 0

# Name of detector slot
config.detectorList[8].name = 'i2'

# Pixel size in the x dimension in mm
config.detectorList[8].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[8].transformDict.nativeSys = 'Pixels'

config.detectorList[8].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[8].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[8].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[8].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[8].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[8].offset_x = 95.25

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[8].offset_y = -53.339999999999996

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[8].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[8].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[8].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[8].serial = 'i2'

# pitch (rotation about y) of the detector in degrees
config.detectorList[8].pitchDeg = 0.0

# ID of detector slot
config.detectorList[8].id = 8

config.detectorList[9] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[9].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[9].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[9].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[9].bbox_x0 = 0

# Name of detector slot
config.detectorList[9].name = 'r2'

# Pixel size in the x dimension in mm
config.detectorList[9].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[9].transformDict.nativeSys = 'Pixels'

config.detectorList[9].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[9].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[9].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[9].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[9].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[9].offset_x = 95.25

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[9].offset_y = -106.67999999999999

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[9].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[9].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[9].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[9].serial = 'r2'

# pitch (rotation about y) of the detector in degrees
config.detectorList[9].pitchDeg = 0.0

# ID of detector slot
config.detectorList[9].id = 9

config.detectorList[10] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[10].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[10].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[10].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[10].bbox_x0 = 0

# Name of detector slot
config.detectorList[10].name = 'g3'

# Pixel size in the x dimension in mm
config.detectorList[10].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[10].transformDict.nativeSys = 'Pixels'

config.detectorList[10].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[10].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[10].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[10].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[10].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[10].offset_x = 31.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[10].offset_y = 106.67999999999999

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[10].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[10].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[10].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[10].serial = 'g3'

# pitch (rotation about y) of the detector in degrees
config.detectorList[10].pitchDeg = 0.0

# ID of detector slot
config.detectorList[10].id = 10

config.detectorList[11] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[11].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[11].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[11].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[11].bbox_x0 = 0

# Name of detector slot
config.detectorList[11].name = 'z3'

# Pixel size in the x dimension in mm
config.detectorList[11].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[11].transformDict.nativeSys = 'Pixels'

config.detectorList[11].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[11].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[11].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[11].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[11].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[11].offset_x = 31.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[11].offset_y = 53.339999999999996

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[11].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[11].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[11].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[11].serial = 'z3'

# pitch (rotation about y) of the detector in degrees
config.detectorList[11].pitchDeg = 0.0

# ID of detector slot
config.detectorList[11].id = 11

config.detectorList[12] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[12].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[12].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[12].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[12].bbox_x0 = 0

# Name of detector slot
config.detectorList[12].name = 'u3'

# Pixel size in the x dimension in mm
config.detectorList[12].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[12].transformDict.nativeSys = 'Pixels'

config.detectorList[12].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[12].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[12].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[12].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[12].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[12].offset_x = 31.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[12].offset_y = 0.0

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[12].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[12].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[12].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[12].serial = 'u3'

# pitch (rotation about y) of the detector in degrees
config.detectorList[12].pitchDeg = 0.0

# ID of detector slot
config.detectorList[12].id = 12

config.detectorList[13] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[13].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[13].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[13].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[13].bbox_x0 = 0

# Name of detector slot
config.detectorList[13].name = 'i3'

# Pixel size in the x dimension in mm
config.detectorList[13].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[13].transformDict.nativeSys = 'Pixels'

config.detectorList[13].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[13].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[13].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[13].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[13].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[13].offset_x = 31.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[13].offset_y = -53.339999999999996

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[13].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[13].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[13].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[13].serial = 'i3'

# pitch (rotation about y) of the detector in degrees
config.detectorList[13].pitchDeg = 0.0

# ID of detector slot
config.detectorList[13].id = 13

config.detectorList[14] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[14].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[14].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[14].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[14].bbox_x0 = 0

# Name of detector slot
config.detectorList[14].name = 'r3'

# Pixel size in the x dimension in mm
config.detectorList[14].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[14].transformDict.nativeSys = 'Pixels'

config.detectorList[14].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[14].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[14].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[14].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[14].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[14].offset_x = 31.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[14].offset_y = -106.67999999999999

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[14].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[14].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[14].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[14].serial = 'r3'

# pitch (rotation about y) of the detector in degrees
config.detectorList[14].pitchDeg = 0.0

# ID of detector slot
config.detectorList[14].id = 14

config.detectorList[15] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[15].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[15].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[15].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[15].bbox_x0 = 0

# Name of detector slot
config.detectorList[15].name = 'g4'

# Pixel size in the x dimension in mm
config.detectorList[15].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[15].transformDict.nativeSys = 'Pixels'

config.detectorList[15].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[15].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[15].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[15].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[15].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[15].offset_x = -31.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[15].offset_y = 106.67999999999999

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[15].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[15].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[15].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[15].serial = 'g4'

# pitch (rotation about y) of the detector in degrees
config.detectorList[15].pitchDeg = 0.0

# ID of detector slot
config.detectorList[15].id = 15

config.detectorList[16] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[16].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[16].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[16].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[16].bbox_x0 = 0

# Name of detector slot
config.detectorList[16].name = 'z4'

# Pixel size in the x dimension in mm
config.detectorList[16].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[16].transformDict.nativeSys = 'Pixels'

config.detectorList[16].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[16].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[16].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[16].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[16].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[16].offset_x = -31.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[16].offset_y = 53.339999999999996

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[16].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[16].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[16].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[16].serial = 'z4'

# pitch (rotation about y) of the detector in degrees
config.detectorList[16].pitchDeg = 0.0

# ID of detector slot
config.detectorList[16].id = 16

config.detectorList[17] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[17].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[17].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[17].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[17].bbox_x0 = 0

# Name of detector slot
config.detectorList[17].name = 'u4'

# Pixel size in the x dimension in mm
config.detectorList[17].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[17].transformDict.nativeSys = 'Pixels'

config.detectorList[17].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[17].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[17].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[17].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[17].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[17].offset_x = -31.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[17].offset_y = 0.0

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[17].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[17].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[17].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[17].serial = 'u4'

# pitch (rotation about y) of the detector in degrees
config.detectorList[17].pitchDeg = 0.0

# ID of detector slot
config.detectorList[17].id = 17

config.detectorList[18] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[18].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[18].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[18].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[18].bbox_x0 = 0

# Name of detector slot
config.detectorList[18].name = 'i4'

# Pixel size in the x dimension in mm
config.detectorList[18].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[18].transformDict.nativeSys = 'Pixels'

config.detectorList[18].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[18].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[18].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[18].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[18].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[18].offset_x = -31.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[18].offset_y = -53.339999999999996

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[18].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[18].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[18].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[18].serial = 'i4'

# pitch (rotation about y) of the detector in degrees
config.detectorList[18].pitchDeg = 0.0

# ID of detector slot
config.detectorList[18].id = 18

config.detectorList[19] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[19].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[19].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[19].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[19].bbox_x0 = 0

# Name of detector slot
config.detectorList[19].name = 'r4'

# Pixel size in the x dimension in mm
config.detectorList[19].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[19].transformDict.nativeSys = 'Pixels'

config.detectorList[19].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[19].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[19].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[19].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[19].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[19].offset_x = -31.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[19].offset_y = -106.67999999999999

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[19].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[19].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[19].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[19].serial = 'r4'

# pitch (rotation about y) of the detector in degrees
config.detectorList[19].pitchDeg = 0.0

# ID of detector slot
config.detectorList[19].id = 19

config.detectorList[20] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[20].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[20].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[20].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[20].bbox_x0 = 0

# Name of detector slot
config.detectorList[20].name = 'g5'

# Pixel size in the x dimension in mm
config.detectorList[20].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[20].transformDict.nativeSys = 'Pixels'

config.detectorList[20].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[20].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[20].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[20].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[20].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[20].offset_x = -95.25

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[20].offset_y = 106.67999999999999

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[20].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[20].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[20].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[20].serial = 'g5'

# pitch (rotation about y) of the detector in degrees
config.detectorList[20].pitchDeg = 0.0

# ID of detector slot
config.detectorList[20].id = 20

config.detectorList[21] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[21].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[21].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[21].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[21].bbox_x0 = 0

# Name of detector slot
config.detectorList[21].name = 'z5'

# Pixel size in the x dimension in mm
config.detectorList[21].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[21].transformDict.nativeSys = 'Pixels'

config.detectorList[21].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[21].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[21].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[21].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[21].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[21].offset_x = -95.25

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[21].offset_y = 53.339999999999996

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[21].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[21].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[21].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[21].serial = 'z5'

# pitch (rotation about y) of the detector in degrees
config.detectorList[21].pitchDeg = 0.0

# ID of detector slot
config.detectorList[21].id = 21

config.detectorList[22] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[22].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[22].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[22].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[22].bbox_x0 = 0

# Name of detector slot
config.detectorList[22].name = 'u5'

# Pixel size in the x dimension in mm
config.detectorList[22].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[22].transformDict.nativeSys = 'Pixels'

config.detectorList[22].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[22].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[22].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[22].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[22].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[22].offset_x = -95.25

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[22].offset_y = 0.0

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[22].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[22].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[22].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[22].serial = 'u5'

# pitch (rotation about y) of the detector in degrees
config.detectorList[22].pitchDeg = 0.0

# ID of detector slot
config.detectorList[22].id = 22

config.detectorList[23] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[23].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[23].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[23].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[23].bbox_x0 = 0

# Name of detector slot
config.detectorList[23].name = 'i5'

# Pixel size in the x dimension in mm
config.detectorList[23].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[23].transformDict.nativeSys = 'Pixels'

config.detectorList[23].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[23].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[23].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[23].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[23].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[23].offset_x = -95.25

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[23].offset_y = -53.339999999999996

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[23].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[23].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[23].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[23].serial = 'i5'

# pitch (rotation about y) of the detector in degrees
config.detectorList[23].pitchDeg = 0.0

# ID of detector slot
config.detectorList[23].id = 23

config.detectorList[24] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[24].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[24].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[24].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[24].bbox_x0 = 0

# Name of detector slot
config.detectorList[24].name = 'r5'

# Pixel size in the x dimension in mm
config.detectorList[24].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[24].transformDict.nativeSys = 'Pixels'

config.detectorList[24].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[24].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[24].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[24].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[24].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[24].offset_x = -95.25

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[24].offset_y = -106.67999999999999

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[24].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[24].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[24].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[24].serial = 'r5'

# pitch (rotation about y) of the detector in degrees
config.detectorList[24].pitchDeg = 0.0

# ID of detector slot
config.detectorList[24].id = 24

config.detectorList[25] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[25].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[25].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[25].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[25].bbox_x0 = 0

# Name of detector slot
config.detectorList[25].name = 'g6'

# Pixel size in the x dimension in mm
config.detectorList[25].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[25].transformDict.nativeSys = 'Pixels'

config.detectorList[25].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[25].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[25].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[25].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[25].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[25].offset_x = -158.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[25].offset_y = 106.67999999999999

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[25].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[25].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[25].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[25].serial = 'g6'

# pitch (rotation about y) of the detector in degrees
config.detectorList[25].pitchDeg = 0.0

# ID of detector slot
config.detectorList[25].id = 25

config.detectorList[26] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[26].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[26].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[26].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[26].bbox_x0 = 0

# Name of detector slot
config.detectorList[26].name = 'z6'

# Pixel size in the x dimension in mm
config.detectorList[26].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[26].transformDict.nativeSys = 'Pixels'

config.detectorList[26].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[26].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[26].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[26].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[26].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[26].offset_x = -158.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[26].offset_y = 53.339999999999996

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[26].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[26].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[26].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[26].serial = 'z6'

# pitch (rotation about y) of the detector in degrees
config.detectorList[26].pitchDeg = 0.0

# ID of detector slot
config.detectorList[26].id = 26

config.detectorList[27] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[27].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[27].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[27].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[27].bbox_x0 = 0

# Name of detector slot
config.detectorList[27].name = 'u6'

# Pixel size in the x dimension in mm
config.detectorList[27].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[27].transformDict.nativeSys = 'Pixels'

config.detectorList[27].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[27].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[27].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[27].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[27].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[27].offset_x = -158.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[27].offset_y = 0.0

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[27].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[27].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[27].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[27].serial = 'u6'

# pitch (rotation about y) of the detector in degrees
config.detectorList[27].pitchDeg = 0.0

# ID of detector slot
config.detectorList[27].id = 27

config.detectorList[28] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[28].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[28].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[28].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[28].bbox_x0 = 0

# Name of detector slot
config.detectorList[28].name = 'i6'

# Pixel size in the x dimension in mm
config.detectorList[28].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[28].transformDict.nativeSys = 'Pixels'

config.detectorList[28].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[28].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[28].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[28].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[28].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[28].offset_x = -158.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[28].offset_y = -53.339999999999996

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[28].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[28].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[28].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[28].serial = 'i6'

# pitch (rotation about y) of the detector in degrees
config.detectorList[28].pitchDeg = 0.0

# ID of detector slot
config.detectorList[28].id = 28

config.detectorList[29] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
# y0 of pixel bounding box
config.detectorList[29].bbox_y0 = 0

# y1 of pixel bounding box
config.detectorList[29].bbox_y1 = 1360

# x1 of pixel bounding box
config.detectorList[29].bbox_x1 = 2047

# x0 of pixel bounding box
config.detectorList[29].bbox_x0 = 0

# Name of detector slot
config.detectorList[29].name = 'r6'

# Pixel size in the x dimension in mm
config.detectorList[29].pixelSize_x = 0.024

# Name of native coordinate system
config.detectorList[29].transformDict.nativeSys = 'Pixels'

config.detectorList[29].transformDict.transforms = None
# x position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[29].refpos_x = 1023.5

# y position of the reference point in the detector in pixels in transposed coordinates.
config.detectorList[29].refpos_y = 680.0

# Pixel size in the y dimension in mm
config.detectorList[29].pixelSize_y = 0.024

# Detector type: SCIENCE=0, FOCUS=1, GUIDER=2, WAVEFRONT=3
config.detectorList[29].detectorType = 0

# x offset from the origin of the camera in mm in the transposed system.
config.detectorList[29].offset_x = -158.75

# y offset from the origin of the camera in mm in the transposed system.
config.detectorList[29].offset_y = -106.67999999999999

# Transpose the pixel grid before orienting in focal plane?
config.detectorList[29].transposeDetector = False

# yaw (rotation about z) of the detector in degrees. This includes any
# necessary rotation to go from detector coordinates to camera coordinates
# after optional transposition.
config.detectorList[29].yawDeg = 0.0

# roll (rotation about x) of the detector in degrees
config.detectorList[29].rollDeg = 0.0

# Serial string associated with this specific detector
config.detectorList[29].serial = 'r6'

# pitch (rotation about y) of the detector in degrees
config.detectorList[29].pitchDeg = 0.0

# ID of detector slot
config.detectorList[29].id = 29

# Coefficients for radial distortion
config.radialCoeffs = None

# Name of this camera
config.name = 'SDSS'
