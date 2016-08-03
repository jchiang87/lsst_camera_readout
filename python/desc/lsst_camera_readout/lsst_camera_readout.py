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
import copy
import numpy as np
import astropy.io.fits as fits
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage

__all__ = ['ImageSource', 'set_itl_bboxes', 'set_e2v_bboxes']

class ImageSource(object):
    '''
    Class to create single segment images based on the pixel geometry
    described by a Camera object from lsst.afw.cameraGeom.

    Parameters
    ----------
    eimage_file : str
        Filename of the eimage FITS file from which the amplifier images
        will be extracted.

    Attributes
    ----------
    eimage : astropy.io.fits.HDUList
        The input eimage data.
    _amp_images : dict
        Dictionary of amplifier images to serve as a cache so that each
        amplifier image is constructed only once.
    '''
    def __init__(self, eimage_file):
        """
        Class constructor.
        """
        self.eimage = fits.open(eimage_file)
        # The eimage data from phosim seems to have x- and y-directions
        # swapped, so transpose it.
        self.eimage_data = self.eimage[0].data.transpose()
        self._amp_images = {}

    def getAmpImage(self, amp, imageFactory=afwImage.ImageI):
        """
        Return an amplifier afwImage.Image object with electronics readout
        applied.

        Parameters
        ----------
        amp : lsst.afw.table.tableLib.AmpInfoRecord
            Data structure containing the amplifier information such as
            pixel geometry, gain, noise, etc..
        imageFactory : lsst.afw.image.Image[FIU], optional
            Image factory to be used for creating the return value.

        Returns
        -------
        lsst.afw.Image[FIU]
            The image object containing the pixel data.
        """
        self._check_amp_geometry(amp)
        if not self._amp_images.has_key((amp.getName(), imageFactory)):
            self._amp_images[(amp.getName(), imageFactory)] \
                = self._make_amp_image(amp, imageFactory)
        return self._amp_images[(amp.getName(), imageFactory)]

    def _make_amp_image(self, amp, imageFactory):
        """
        Create the segment image for the amplier geometry specified in amp.

        Parameters
        ----------
        amp : lsst.afw.table.tableLib.AmpInfoRecord
            Data structure containing the amplifier information such as
            pixel geometry, gain, noise, etc..
        imageFactory : lsst.afw.image.Image[FIU], optional
            Image factory to be used for creating the return value.

        Returns
        -------
        lsst.afw.Image[FIU]
            The image object containing the pixel data.
        """
        bbox = amp.getBBox()
        full_segment = afwImage.ImageF(amp.getRawBBox())

        # Get the imaging segment (i.e., excluding prescan and
        # overscan regions), and fill with data from the eimage.
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

        # Return image as type given by imageFactory.
        output_image = imageFactory(amp.getRawBBox())
        output_image.getArray()[:] = full_arr

        return output_image

    def _check_amp_geometry(self, amp):
        """
        Check that the imaging section of amp is consistent with the
        eimage geometry.

        Parameters
        ----------
        amp : lsst.afw.table.tableLib.AmpInfoRecord
            Data structure containing the amplifier information such as
            pixel geometry, gain, noise, etc..

        Raises
        ------
        RuntimeError
            If the amplifier geometry is inconsistent with the eimage
            geometry.
        """
        bbox = amp.getBBox()
        xsize = bbox.getMaxX() - bbox.getMinX() + 1
        ysize = bbox.getMaxY() - bbox.getMinY() + 1
        if (8*xsize != self.eimage[0].header['NAXIS2'] or
            2*ysize != self.eimage[0].header['NAXIS1']):
            raise RuntimeError('amplifier geometry is inconsistent with eimage')

    def write_ampliflier_image(self, amp, outfile, clobber=True):
        """
        Write the pixel data for the specified amplifier as FITS image.

        Parameters
        ----------
        amp : lsst.afw.table.tableLib.AmpInfoRecord
            Data structure containing the amplifier information such as
            pixel geometry, gain, noise, etc..
        outfile : str
            Filename of the FITS file to be written.
        clobber : bool, optional
            Flag whether to overwrite an existing output file.
        """
        output = fits.HDUList()
        output.append(copy.deepcopy(self.eimage[0]))
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
        """
        Convert bounding boxes into NOAO section keywords.

        Parameters
        ----------
        bbox : lsst.afw.geom.Box2I
            Bounding box.
        flipx : bool
            Flag to indicate that data should be flipped in the x-direction.
        flipy : bool
            Flag to indicate that data should be flipped in the y-direction.
        """
        xmin, xmax = bbox.getMinX()+1, bbox.getMaxX()+1
        ymin, ymax = bbox.getMinY()+1, bbox.getMaxY()+1
        if flipx:
            xmin, xmax = xmax, xmin
        if flipy:
            ymin, ymax = ymax, ymin
        return '[%i:%i,%i:%i]' % (xmin, xmax, ymin, ymax)

def set_itl_bboxes(amp):
    """
    Function to apply realistic pixel geometry for ITL sensors.

    Parameters
    ----------
    amp : lsst.afw.table.tableLib.AmpInfoRecord
        Data structure containing the amplifier information such as
        pixel geometry, gain, noise, etc..

    Returns
    -------
    lsst.afw.table.tableLib.AmpInfoRecord
        The updated AmpInfoRecord.
    """
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

def set_e2v_bboxes(amp):
    """
    Function to apply realistic pixel geometry for e2v sensors.

    Parameters
    ----------
    amp : lsst.afw.table.tableLib.AmpInfoRecord
        Data structure containing the amplifier information such as
        pixel geometry, gain, noise, etc..

    Returns
    -------
    lsst.afw.table.tableLib.AmpInfoRecord
        The updated AmpInfoRecord.
    """
    amp.setRawBBox(afwGeom.Box2I(afwGeom.Point2I(0, 0),
                                 afwGeom.Extent2I(542, 2022)))
    amp.setRawDataBBox(afwGeom.Box2I(afwGeom.Point2I(10, 0),
                                     afwGeom.Extent2I(522, 2002)))
    amp.setRawHorizontalOverscanBBox(afwGeom.Box2I(afwGeom.Point2I(522, 0),
                                                   afwGeom.Extent2I(20, 2002)))
    amp.setRawVerticalOverscanBBox(afwGeom.Box2I(afwGeom.Point2I(0, 2002),
                                                 afwGeom.Extent2I(542, 20)))
    amp.setRawPrescanBBox(afwGeom.Box2I(afwGeom.Point2I(0, 0),
                                        afwGeom.Extent2I(10, 2002)))
    return amp
