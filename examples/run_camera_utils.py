from __future__ import print_function, absolute_import, division
import astropy.io.fits as fits
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.obs.lsstSim as lsstSim
from desc.lsst_camera_readout import ImageSource, set_itl_bboxes,\
    set_phosim_bboxes

#eimage_file = '../data/lsst_e_200_f2_R22_S11_E000.fits.gz'
#eimage_file = '../work/phosim-phosim_release-cc088644cb36/output/lsst_e_99999999_f2_R22_S11_E000.fits.gz'
eimage_file = '../work/phosim-phosim_release-cc088644cb36/work/lsst_e_99999999_f2_R22_S11_E000.fits.gz'
image_source = ImageSource(eimage_file, seg_file='../bin/segmentation.txt_new')

output = fits.HDUList()
output.append(fits.PrimaryHDU())
for col in '10':
    for row in '01234567':
        amp_name = 'R22_S11_C' + col + row
        outfile = 'C%s%s_image.fits' % (col, row)
        image_source.write_amplifier_image(amp_name, outfile)
        output.append(fits.open(outfile)[0])
output.writeto('mef.fits', clobber=True)

# Create a transposed eimage file for blinking in ds9 against the
# mosaicked MEF.
eimage = fits.open(eimage_file)
#eimage[0].data = eimage[0].data.transpose()
eimage.writeto('eimage_transpose.fits', clobber=True)

# Use cameraGeomUtils to display one segment with different image factories.
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

for col in '10':
    for row in '01234567':
        print(image_source.fp_props.get_amp('R22_S11_C'+col+row).mosaic_section)

#cameraGeomUtils.showAmp(amp, imageSource=image_source,
#                        display=ds9.getDisplay(frame=1),
#                        imageFactory=afwImage.ImageF)
#
#cameraGeomUtils.showAmp(amp, imageSource=image_source,
#                        display=ds9.getDisplay(frame=2),
#                        imageFactory=afwImage.ImageD)
#
#cameraGeomUtils.showAmp(amp, imageSource=image_source,
#                        display=ds9.getDisplay(frame=3),
#                        imageFactory=afwImage.ImageU)
