"""
Classes to convert an eimage to individual sensor segments, applying
electronics readout effects.

Steps:
  Read in single sensor eimage
  Extract geometry from amplifier record
  Copy imaging segment pixels
  Add defects (bright defects, dark defects, traps)
  Apply CTE
  Apply crosstalk
  Add bias
  Add dark current
  Apply gain
  Write FITS file for each amplifier
"""
from __future__ import print_function, absolute_import, division
import numpy as np
import astropy.io.fits as fits
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage

__all__ = ['ImageSource', 'set_raw_bboxes']

class ImageSource(object):
    '''
    Class to create single segment images based on the pixel geometry
    described by a Camera object from lsst.afw.cameraGeom.
    '''
    def __init__(self, eimage_file):
        "Constructor"
        self.eimage = fits.open(eimage_file)
        self.eimage_data = self.eimage[0].data
        self._amp_images = {}

    def getAmpImage(self, amp, imageFactory=afwImage.ImageI):
        if not self._amp_images.has_key(amp.getName()):
            self._amp_images[amp.getName()] \
                = self._makeAmpImage(amp, imageFactory)
        return self._amp_images[amp.getName()]

    def _makeAmpImage(self, amp, imageFactory):
        "Return the segment image for the amplier geometry specified in amp."
        bbox = amp.getBBox()
        full_segment = afwImage.ImageF(amp.getRawBBox())

        # Get the imaging segment (i.e., excluding pre- and overscan
        # regions), and fill with data from the eimage.
        imaging_segment = full_segment.Factory(full_segment,
                                               amp.getRawDataBBox())
        data = self.eimage_data[bbox.getMinY():bbox.getMaxY()+1,
                                bbox.getMinX():bbox.getMaxX()+1].copy()

        # Apply flips in x and y relative to assembled eimage.
        if amp.getRawFlipX():
            data = data[:, ::-1]
        if amp.getRawFlipY():
            data = data[::-1, :]

        imaging_segment.getArray()[:] = data

        # Add read noise.
        full_arr = full_segment.getArray()
        full_arr += np.random.normal(scale=amp.getReadNoise(),
                                     size=full_arr.shape)
        # Convert to ADU.
        full_arr /= amp.getGain()

        # Return image as ints.
        output_image = imageFactory(amp.getRawBBox())
        output_image.getArray()[:] = np.array(full_arr, dtype=np.int)

        return output_image

    def write_ampliflier_image(self, amp, outfile, clobber=True):
        output = fits.HDUList()
        output.append(self.eimage[0])
        amp_image = self.getAmpImage(amp)
        output[0].data = amp_image.getArray()
        output[0].header['DATASEC'] \
            = self._noao_section_keyword(amp.getRawDataBBox())
        output[0].header['DETSEC'] \
            = self._noao_section_keyword(amp.getBBox(),
                                         flipx=amp.getRawFlipX(),
                                         flipy=amp.getRawFlipY())
        output[0].header['BIASSEC'] \
            = self._noao_section_keyword(amp.getRawHorizontalOverscanBBox())
        output[0].header['GAIN'] = amp.getGain()
        output.writeto(outfile, clobber=clobber)

    @staticmethod
    def _noao_section_keyword(bbox, flipx=False, flipy=False):
        xmin, xmax = bbox.getMinX()+1, bbox.getMaxX()+1
        ymin, ymax = bbox.getMinY()+1, bbox.getMaxY()+1
        if flipx:
            xmin, xmax = xmax, xmin
        if flipy:
            ymin, ymax = ymax, ymin
        return '[%i:%i,%i:%i]' % (xmin, xmax, ymin, ymax)

def set_raw_bboxes(amp):
    "Apply realistic pixel geometry for ITL sensors."
    amp.setRawBBox(afwGeom.Box2I(afwGeom.Point2I(0, 0),
                                 afwGeom.Extent2I(532, 2020)))
    amp.setRawDataBBox(afwGeom.Box2I(afwGeom.Point2I(3, 0),
                                     afwGeom.Extent2I(509, 2000)))
    amp.setRawHorizontalOverscanBBox(afwGeom.Box2I(afwGeom.Point2I(512, 0),
                                                   afwGeom.Extent2I(20, 2000)))
    amp.setRawVerticalOverscanBBox(afwGeom.Box2I(afwGeom.Point2I(0, 2000),
                                                 afwGeom.Extent2I(532, 20)))
    amp.setRawPrescanBBox(afwGeom.Box2I(afwGeom.Point2I(0, 0),
                                        afwGeom.Extent2I(3, 2000)))
    return amp
