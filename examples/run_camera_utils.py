import lsst.afw.cameraGeom.utils as cameraGeomUtils
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.obs.lsstSim as lsstSim
from desc.lsst_camera_readout import ImageSource, set_raw_bboxes

mapper = lsstSim.LsstSimMapper()
camera = mapper.camera
display = ds9.getDisplay()

raft = 'R:2,2'
ccd = 'S:1,1'
sensor = camera[' '.join((raft, ccd))]

image_source = ImageSource('../data/lsst_e_200_f2_R22_S11_E000.fits.gz')

amp_id = '0,0'
amp = sensor[amp_id]
#amp = set_raw_bboxes(amp)

cameraGeomUtils.showAmp(amp, imageSource=image_source, display=display,
                        imageFactory=afwImage.ImageI)
