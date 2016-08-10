"""
Code to convert an eimage to individual sensor segments, applying
electronics readout effects.

 * Read in single sensor eimage
 * Extract geometry from amplifier record
 * Copy imaging segment pixels
 * Add bias
 * Add dark current
 * Add defects (bright defects, dark defects, traps)
 * Apply CTE
 * Apply crosstalk
 * Apply gain
 * Write FITS file for each amplifier

"""
from __future__ import print_function, absolute_import, division
import os
import copy
import numpy as np
import astropy.io.fits as fits
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.utils as lsstUtils
from .focalplane_readout import FocalPlaneReadout

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
    seg_file : str, optional
        Full path of segmentation.txt file, the PhoSim-formatted file
        that describes the properties of the sensors in the focal
        plane.  If None, then the version in obs_lsstSim/description
        will be used.
    raft_sensor : str, optional
        The raft and sensor identifier, e.g., 'R:2,2_S:1,1'.  If None,
        then it will be extracted from the eimage_file name.
    add_read_noise : bool, optional
        Flag to add read noise.

    Attributes
    ----------
    eimage : astropy.io.fits.HDUList
        The input eimage data.
    _amp_images : dict
        Dictionary of amplifier images to serve as a cache so that each
        amplifier image is constructed only once.
    fp_props : FocalPlaneReadout object
        Object containing the readout properties of the sensors in the
        focal plane, extracted from the segmentation.txt file.
    '''
    def __init__(self, eimage_file, seg_file=None, raft_sensor=None,
                 add_read_noise=True):
        """
        Class constructor.
        """
        self.eimage = fits.open(eimage_file)
        # The eimage data from phosim seems to have x- and y-directions
        # swapped, so transpose it.
        self.eimage_data = self.eimage[0].data.transpose()

        if seg_file is None:
            seg_file = os.path.join(lsstUtils.getPackageDir('obs_lsstSim'),
                                    'description', 'segmentation.txt')
        self.fp_props = FocalPlaneReadout.read_phosim_seg_file(seg_file)

        if raft_sensor is None:
            self.raft, self.sensor = self.extract_sensor_id(eimage_file)
        else:
            self.raft, self.sensor = raft_sensor.split('_')

        self._make_amp_images(add_read_noise)

    @staticmethod
    def extract_sensor_id(eimage_file):
        """
        Extract the raft and sensor ids from the eimage filename.

        Parameters
        ----------
        eimage_file : str
            Filename of the eimage FITS file from which the amplifier
            images will be extracted.

        Returns
        -------
        tuple
            A tuple containing (raft, sensor).
        """
        tokens = os.path.basename(eimage_file).split('_')
        raft = 'R:%s,%s' % (tokens[4][1], tokens[4][2])
        sensor = 'S:%s,%s' % (tokens[5][1], tokens[5][2])
        return raft, sensor

    def _exptime(self):
        """
        The exposure time of the frame in seconds.

        Returns
        -------
        float
            The exposure time of the frame in seconds from the eimage_file.
        """
        return self.eimage[0].header['EXPTIME']

    def getAmpImage(self, amp, imageFactory=afwImage.ImageI):
        """
        Return an amplifier afwImage.Image object with electronics
        readout effects applied.

        Parameters
        ----------
        amp : lsst.afw.table.tableLib.AmpInfoRecord
            Data structure containing the amplifier information such as
            pixel geometry, gain, noise, etc..
        imageFactory : lsst.afw.image.Image[DFIU], optional
            Image factory to be used for creating the return value.

        Returns
        -------
        lsst.afw.Image[DFIU]
            The image object containing the pixel data.
        """
        float_image = self._amp_images[amp.getName()]
        if imageFactory == afwImage.ImageF:
            return float_image
        # Return image as the type given by imageFactory.
        output_image = imageFactory(amp.getRawBBox())
        output_image.getArray()[:] = float_image.getArray()
        return output_image

    def _make_amp_images(self, add_read_noise):
        sensor_props = self.fp_props.get_sensor(self.raft, self.sensor)
        for amp_name in sensor_props.amp_names:
            self._make_amp_image(amp_name, add_read_noise)
        self._apply_crosstalk()

    def _make_amp_image(self, amp_name, add_read_noise):
        """
        Create the segment image for the amplier geometry specified in amp.

        Parameters
        ----------
        amp_name : str
            The amplifier name, e.g., "R22_S11_C00".
        add_read_noise : bool
            Flag to add read noise.

        """
        amp_props = self.fp_props.get_amp(amp_name)
        bbox = amp_props.mosiac_section
        full_segment = afwImage.ImageF(amp_props.full_segment)

        # Get the imaging segment (i.e., excluding prescan and
        # overscan regions), and fill with data from the eimage.
        imaging_segment = full_segment.Factory(full_segment, amp_props.imaging)
        data = self.eimage_data[bbox.getMinY():bbox.getMaxY()+1,
                                bbox.getMinX():bbox.getMaxX()+1].copy()

        # Apply flips in x and y relative to assembled eimage in order
        # to have the pixels in readout order.
        if amp_props.flipx:
            data = data[:, ::-1]
        if amp_props.flipy:
            data = data[::-1, :]

        imaging_segment.getArray()[:] = data
        full_arr = full_segment.getArray()

        # Add dark current.
        full_arr += np.random.poisson(amp_props.dark_current*self._exptime(),
                                      size=full_arr.shape)

        # Add defects.

        # Apply CTE.

        # Convert to ADU.
        full_arr /= amp_props.gain

        # Add read noise.
        if add_read_noise:
            full_arr += np.random.normal(scale=amp_props.read_noise,
                                         size=full_arr.shape)
        # Add bias level.
        full_arr += amp_props.bias_level

        self._amp_images[amp_name] = full_segment

    def _apply_crosstalk(self):
        """
        Apply inter-amplifier crosstalk using the cross-talk matrix
        from segmentation.txt.  This should be run only once and
        only after ._make_amp_image has been run for each amplifier.
        """
        sensor_props = self.fp_props.get_sensor(self.raft, self.sensor)
        imarrs = np.array([self._amp_images[amp_name].getArray()
                           for amp_name in sensor_props.amp_names])
        for amp_name in sensor_props.amp_names:
            amp_props = self.fp_props.get_amplifier(self.raft, self.sensor,
                                                    channel(amp_name))
            self.amp_images[amp_name].getArray()[:, :] \
                = sum(imarrs*amp_props.crosstalk)

    def check_amp_geometry(self, amp):
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

def channel(amp_name):
    """
    Extract the channel id from the full amp_name used by Phosim.

    Parameters
    ----------
    amp_name : str
        Amplifier name used by PhoSim, e.g., "R22_S11_C00"

    Returns
    -------
    str
        The channel name formatted DM-style, e.g., "C:0,0"
    """
    tokens = amp_name.split('_')
    return "C:%s,%s" % (tokens[1], tokens[2])
