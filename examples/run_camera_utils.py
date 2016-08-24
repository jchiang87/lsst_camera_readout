from __future__ import print_function, absolute_import, division
import astropy.io.fits as fits
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.obs.lsstSim as lsstSim
from desc.lsst_camera_readout import ImageSource, set_itl_bboxes

eimage_file = '../data/lsst_e_921297_f2_R22_S11_E000_transpose.fits.gz'
seg_file = '../data/segmentation_itl.txt'

# Create an ImageSource object directly from the eimage file.
#image_source = ImageSource.create_from_eimage(eimage_file, seg_file=seg_file)

# Create an ImageSource object from a numpy array.
imarr = fits.open(eimage_file)[0].data
exptime = 30.
sensor_id = 'R22_S11'
image_source = ImageSource(imarr, exptime, sensor_id, seg_file=seg_file)

output = fits.HDUList()
output.append(fits.PrimaryHDU())
for col in '10':
    for row in '01234567':
        amp_name = 'R22_S11_C' + col + row
        outfile = 'C%s%s_image.fits' % (col, row)
        image_source.write_amplifier_image(amp_name, outfile)
        output.append(fits.open(outfile)[0])
output.writeto('mef.fits', clobber=True)

# Flip the eimage in the x-direction so that we are viewing it through L3.
eimage = fits.open(eimage_file)
eimage[0].data = eimage[0].data[:, ::-1]
eimage.writeto('eimage_xflip.fits', clobber=True)

# Use cameraGeomUtils to display a segment.
mapper = lsstSim.LsstSimMapper()
camera = mapper.camera
raft = 'R:2,2'
ccd = 'S:1,1'
sensor = camera[' '.join((raft, ccd))]

amp = sensor['0,3']
amp = set_itl_bboxes(amp)
cameraGeomUtils.showAmp(amp, imageSource=image_source,
                        display=ds9.getDisplay(frame=0),
                        imageFactory=afwImage.ImageI)
