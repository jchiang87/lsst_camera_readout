from __future__ import print_function, absolute_import, division
import astropy.io.fits as fits
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.obs.lsstSim as lsstSim
from desc.lsst_camera_readout import ImageSource, set_itl_bboxes

mapper = lsstSim.LsstSimMapper()
camera = mapper.camera

raft = 'R:2,2'
ccd = 'S:1,1'
sensor = camera[' '.join((raft, ccd))]

eimage_file = '../data/lsst_e_200_f2_R22_S11_E000.fits.gz'
image_source = ImageSource(eimage_file)
gain = 1  # use a common gain for all segments

output = fits.HDUList()
output.append(fits.PrimaryHDU())
for col in '10':
    for row in '01234567':
        amp_id = '%s,%s' % (col, row)
        print('processing', amp_id)
        amp = sensor[amp_id]
#        amp = set_itl_bboxes(amp)
        amp.setGain(gain)
        outfile = 'C%s%s_image.fits' % (col, row)
        image_source.write_ampliflier_image(amp, outfile)
        output.append(fits.open(outfile)[0])
output.writeto('mef.fits', clobber=True)

# Create a transposed eimage file for blinking in ds9 against the
# mosaicked MEF.
eimage = fits.open(eimage_file)
eimage[0].data = eimage[0].data.transpose()
eimage.writeto('eimage_transpose.fits', clobber=True)

# Use cameraGeomUtils to display one segment with different image factories.
amp = sensor['0,0']
amp = set_itl_bboxes(amp)
cameraGeomUtils.showAmp(amp, imageSource=image_source,
                        display=ds9.getDisplay(frame=0),
                        imageFactory=afwImage.ImageI)

cameraGeomUtils.showAmp(amp, imageSource=image_source,
                        display=ds9.getDisplay(frame=1),
                        imageFactory=afwImage.ImageF)

cameraGeomUtils.showAmp(amp, imageSource=image_source,
                        display=ds9.getDisplay(frame=2),
                        imageFactory=afwImage.ImageD)

cameraGeomUtils.showAmp(amp, imageSource=image_source,
                        display=ds9.getDisplay(frame=3),
                        imageFactory=afwImage.ImageU)
