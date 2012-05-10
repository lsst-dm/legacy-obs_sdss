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
import os, sys
import pyfits
import cPickle
import numpy as num
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

from scipy.interpolate import UnivariateSpline
import pylab
from matplotlib.ticker import FormatStrFormatter

from convertfpM import convertfpM
from convertasTrans import convertasTrans

import lsst.afw.display.ds9 as ds9

fFormatter = FormatStrFormatter('%d')
rootdir = "/astro/net/pogo1/stripe82/imaging"

class RunBackgroundBinner(object):
    def __init__(self, run, rerun, filt, camcol):
        self.run     = run
        self.rerun   = rerun
        self.filt    = filt
        self.camcol  = camcol

        self.asTrans = getasTrans(self.run, self.rerun)
        
        self.ras     = []
        self.decs    = []
        self.medians = []
        self.iqrs    = []
        self.fields  = []
 
    def convertRa(self, ra):
        if ra > 180:
            return ra - 360
        return ra

    def bin(self, image, mask, field, nbin = 4):
        ny   = image.getHeight() // 4
        for n in range(nbin):
            ymin = n * ny
            ymax = min((n + 1) * ny, image.getHeight())

            idx  = num.where(mask.getArray()[ymin:ymax,:] == 0)
            data = image.getArray()[ymin:ymax,:][idx]
            stat = afwMath.makeStatistics(data.tolist(), afwMath.MEDIAN | afwMath.IQRANGE)
            med  = stat.getValue(afwMath.MEDIAN)
            iqr  = stat.getValue(afwMath.IQRANGE)

            wcs       = convertasTrans(self.asTrans, self.filt, self.camcol, field)
            center    = wcs.pixelToSky(image.getWidth()//2, 0.5 * (ymin + ymax))
            

            self.ras.append(self.convertRa(center[0].asDegrees()))
            self.decs.append(center[1].asDegrees())
            self.medians.append(med)
            self.iqrs.append(iqr)
            self.fields.append(field)

    def process(self):
        if len(self.ras) == 0:
            return
        self.ras     = num.array(self.ras)
        self.decs    = num.array(self.decs)
        self.medians = num.array(self.medians)
        self.iqrs    = num.array(self.iqrs)
        self.spline  = self.doSpline(self.ras, self.medians, 0.741 * self.iqrs)

    def doSpline(self, x, y, dy):
        return UnivariateSpline(x, y, w = 1./dy)

    def doPlot(self, show = True):
        if len(self.ras) == 0:
            return

        fig = pylab.figure()
        ax1 = fig.add_subplot(111)
        ax1.errorbar(self.ras, self.medians, yerr = 0.741 * self.iqrs, fmt='ro', ms=2)

        ax2 = ax1.twiny()
        ax2.set_xlabel("Field")
        ax2.set_xlim(min(self.fields), max(self.fields))
        ax2.xaxis.set_major_formatter( fFormatter )

        xspl = num.arange(num.min(self.ras), num.max(self.ras), 0.05)
        yspl = self.spline(xspl)
        ax1.plot(xspl, yspl, "k-")
        fig.suptitle("%06d-%s%d   Dec = %.5f" % (self.run, self.filt, self.camcol, 
                                                 afwMath.makeStatistics(self.decs, afwMath.MEDIAN).getValue(afwMath.MEDIAN)))
        ax1.set_xlabel("RA")
        ax1.set_ylabel("Sky")
        if show:
            pylab.show()
                    

def getfpC(run, rerun, filt, camcol, field):
    fname = os.path.join(rootdir, str(run), str(rerun), "corr", str(camcol), "fpC-%06d-%s%d-%04d.fit.gz" % (run, filt, camcol, field))
    print fname
    if os.path.isfile(fname):
        im  = afwImage.ImageF(fname)
        im -= 1000 # damn pedestal
        return im
    return None

def getfpM(run, rerun, filt, camcol, field):
    fname = os.path.join(rootdir, str(run), str(rerun), "objcs", str(camcol), "fpM-%06d-%s%d-%04d.fit" % (run, filt, camcol, field))
    print fname
    if os.path.isfile(fname):
        try:
            return convertfpM(fname, allPlanes = True)
        except:
            return None
    return None

def getasTrans(run, rerun):
    fname = os.path.join(rootdir, str(run), str(rerun), "astrom", "asTrans-%06d.fit" % (run))
    print fname
    if os.path.isfile(fname):
        return fname
    return None
 
if __name__ == '__main__':
    cc    = int(sys.argv[1]) # parallelize
    rerun = 40
    runs  = [94, 125, 1033, 1056, 1752, 1755, 1894, 2385, 2570, 2578, 2579,
             2583, 2585, 2589, 2649, 2650, 2659, 2662, 2677, 2700, 2708, 2709,
             2728, 2738, 2768, 2820, 2855, 2861, 2873, 2886, 2960, 2968, 3325,
             3355, 3360, 3362, 3384, 3388, 3427, 3430, 3434, 3437, 3438, 3460,
             3461, 3465, 4128, 4136, 4145, 4153, 4157, 4184, 4187, 4188, 4192,
             4198, 4203, 4207, 4247, 4253, 4263, 4288, 4797, 4849, 4858, 4868,
             4874, 4894, 4895, 4899, 4905, 4917, 4927, 4930, 4933, 4948, 5042,
             5052, 5566, 5582, 5590, 5597, 5603, 5607, 5610, 5619, 5622, 5628,
             5633, 5637, 5642, 5646, 5654, 5658, 5665, 5666, 5670, 5675, 5681,
             5698, 5702, 5709, 5713, 5719, 5729, 5730, 5731, 5732, 5743, 5744,
             5745, 5754, 5759, 5760, 5763, 5765, 5770, 5771, 5776, 5777, 5781,
             5782, 5786, 5792, 5797, 5800, 5807, 5808, 5813, 5820, 5823, 5836,
             5842, 5847, 5853, 5865, 5866, 5870, 5871, 5872, 5878, 5882, 5889,
             5895, 5898, 5902, 5905, 5909, 5915, 5918, 5924, 6281, 6283, 6287,
             6293, 6313, 6314, 6330, 6348, 6349, 6353, 6355, 6360, 6362, 6363,
             6367, 6370, 6373, 6374, 6377, 6383, 6391, 6400, 6401, 6402, 6404,
             6408, 6409, 6412, 6414, 6417, 6418, 6421, 6422, 6425, 6430, 6433,
             6435, 6441, 6444, 6447, 6448, 6450, 6453, 6458, 6461, 6464, 6468,
             6471, 6474, 6476, 6479, 6480, 6484, 6488, 6494, 6501, 6504, 6508,
             6513, 6518, 6522, 6524, 6525, 6528, 6530, 6533, 6534, 6537, 6542,
             6545, 6548, 6552, 6555, 6556, 6559, 6564, 6565, 6568, 6571, 6577,
             6580, 6584, 6590, 6592, 6596, 6600, 6604, 6609, 6615, 6618, 6920,
             6921, 6930, 6933, 6934, 6947, 6951, 6955, 6958, 6961, 6962, 6963,
             6964, 6976, 6981, 6982, 6985, 7003, 7006, 7013, 7016, 7018, 7024,
             7033, 7034, 7037, 7038, 7043, 7047, 7051, 7054, 7057, 7060, 7069,
             7071, 7074, 7076, 7077, 7080, 7081, 7084, 7092, 7095, 7096, 7101,
             7106, 7110, 7111, 7112, 7117, 7121, 7124, 7127, 7130, 7133, 7136,
             7140, 7142, 7145, 7150, 7151, 7152, 7155, 7158, 7161, 7164, 7167,
             7170, 7173, 7176, 7177, 7182, 7183, 7188, 7195, 7199, 7202]
    
    filts   = ["u", "g", "r", "i", "z"]
    camcols = [cc,]
    fields  = range(1, 1000)

    for run in runs:
        for filt in filts:
            for camcol in camcols:
                for run in runs:
                    outfile = "bg-%06d-%s%d.pickle" % (run, filt, camcol)
                    if os.path.isfile(outfile):
                        continue

                    binner = RunBackgroundBinner(run, rerun, filt, camcol)
            
                    for field in fields:
                        image   = getfpC(run, rerun, filt, camcol, field)
                        mask    = getfpM(run, rerun, filt, camcol, field)
                        if (not image) or (not mask):
                            continue
                        try:
                            binner.bin(image, mask, field)
                        except:
                            print "FAIL", 

                    binner.process()

                    buff = open(outfile, "wb")
                    cPickle.dump(binner, buff)
                    buff.close()
                    
