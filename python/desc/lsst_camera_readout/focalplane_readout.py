"""
Code to access sensor properties, including pixel layout and electronics
readout information.
"""
from __future__ import print_function, absolute_import, division
import numpy as np
import lsst.afw.geom as afwGeom

__all__ = ['FocalPlaneReadout']

class FocalPlaneReadout(object):
    """
    Class to serve up electronics readout properties of sensors based
    on a focal plane description, such as the segmentation.txt file
    used by PhoSim.

    Attributes
    ----------
    sensors : dict
        A dictionary of sensor properties, keyed by sensor id in the
        LSST focal plane, e.g., "R22_S11".
    amps : dict
        A diction of amplifier properties, keyed by amplifier id, e.g.,
        "R22_S11_C00".
    """
    def __init__(self):
        self.sensors = {}
        self.amps = {}

    def get_sensor(self, sensor_id):
        """
        Access to the specified SensorProperties object.

        Parameters
        ----------
        sensor_id : str
            Sensor ID of the form "Rrr_Sss", e.g., "R22_S11".

        Returns
        -------
        SensorProperties object
            The object containing the sensor-wide properties.
        """
        return self.sensors[sensor_id]

    def get_amp(self, amp_id):
        """
        Access to the specified AmplifierProperties object.

        Parameters
        ----------
        amp_id : str
            Amplifier ID of the form "Rrr_Sss_Ccc", e.g., "R22_S11_C00"

        Returns
        -------
        AmplifierProperties object
            The object containing the amplifier properties.
        """
        return self.amps[amp_id]

    @staticmethod
    def sensor_id(raft, ccd):
        """
        Convert the lsst.afw.cameraGeom specifiers to the name used
        by PhoSim.

        Parameters
        ----------
        raft : str
            Raft id using lsst.afw.cameraGeom syntax, e.g., "R:2,2".
        ccd : str
            Sensor id using lsst.afw.cameraGeom syntax, e.g., "S:1,1".

        Returns
        -------
        str
            e.g., "R22_S11"
        """
        return 'R%s%s_S%s%s' % (raft[2], raft[4], ccd[2], ccd[4])

    @staticmethod
    def amp_id(raft, ccd, chan):
        """
        Convert the lsst.afw.cameraGeom specifiers to the name used
        by PhoSim.

        Parameters
        ----------
        raft : str
            Raft id using lsst.afw.cameraGeom syntax, e.g., "R:2,2".
        ccd : str
            Sensor id using lsst.afw.cameraGeom syntax, e.g., "S:1,1".
        chan : str
            Amplifier channel id using lsst.afw.cameraGeom syntax, e.g.,
            "C:0,0".

        Returns
        -------
        str
            e.g., "R22_S11_C00"
        """
        return 'R%s%s_S%s%s_C%s%s' % (raft[2], raft[4], ccd[2], ccd[4],
                                      chan[2], chan[4])

    @staticmethod
    def read_phosim_seg_file(seg_file):
        """
        Factory method to create a FocalPlaneReadout object which has
        been filled with the data from a PhoSim segmentation.txt file.

        Parameters
        ----------
        seg_file : str
            The PhoSim formatted segmentation.txt file.

        Returns
        -------
        FocalPlaneReadout object
            The filled FocalPlaneReadout object.
        """
        my_self = FocalPlaneReadout()
        with open(seg_file, 'r') as f:
            lines = [line for line in f.readlines() if not line.startswith('#')]
        i = -1
        while True:
            try:
                i += 1
                sensor_props = SensorProperties(lines[i])
                my_self.sensors[sensor_props.name] = sensor_props
                for j in range(sensor_props.num_amps):
                    i += 1
                    amp_props = AmplifierProperties(lines[i])
                    my_self.amps[amp_props.name] = amp_props
                    sensor_props.append_amp(amp_props)
            except IndexError:
                break
        return my_self

class SensorProperties(object):
    """
    Class to contain the properties of a sensor.

    Attributes
    ----------
    name : str
        The sensor name, e.g., "R22_S11".
    num_amps : int
        The number of amplifiers in this sensor.
    height : int
        The number of physical sensor pixels in the parallel direction.
    width : int
        The number of physical sensor pixels in the serial direction.
    amp_names : tuple
        The amplifier names in the order in which they were added
        to self.  For data read in from segmentation.txt, this ordering
        is used for the crosstalk matrix column ordering.
    """
    def __init__(self, line):
        tokens = line.strip().split()
        self.name = tokens[0]
        self.num_amps = int(tokens[1])
        self.height = int(tokens[2])
        self.width = int(tokens[3])
        self._amp_names = []

    @property
    def amp_names(self):
        return tuple(self._amp_names)

    def append_amp(self, amp_props):
        """
        Append an amplifier to the ._amp_order list.

        Parameters
        ----------
        amp_props : AmplifierProperties object
            The object containing the amplifier properties.
        """
        self._amp_names.append(amp_props.name)

class AmplifierProperties(object):
    """
    Class to contain the properties of an amplifier.

    Attributes
    ----------
    name : str
        The amplifier name, e.g., "R22_S11_C00".
    imaging : lsst.afw.geom.Box2I
        The imaging region bounding box.
    full_segment : lsst.afw.geom.Box2I
        The bounding box for the full segment.
    gain : float
        The amplifier gain in units of e-/ADU.
    bias_level : float
        The bias level in e-/pixel.
    dark_current : float
        The dark current in units of e-/pixel/s
    crosstalk : numpy.array
        The row of the intrasensor crosstalk matrix for this amplifier.
    flip_x : bool
        Flag to indicate that pixel ordering in x-direction should be
        reversed relative to mosaicked image.
    flip_y : bool
        Flag to indicate that pixel ordering in y-direction should be
        reversed relative to mosaicked image.
    """
    def __init__(self, line):
        tokens = line.strip().split()
        self.name = tokens[0]
        ymin, ymax, xmin, xmax = (int(x) for x in tokens[1:5])
        xsize = xmax - xmin + 1
        ysize = ymax - ymin + 1
        self.mosaic_section = afwGeom.Box2I(afwGeom.Point2I(xmin, ymin),
                                            afwGeom.Extent2I(xsize, ysize))
        parallel_prescan = int(tokens[15])
        serial_overscan = int(tokens[16])
        serial_prescan = int(tokens[17])
        parallel_overscan = int(tokens[18])
        self.imaging = afwGeom.Box2I(afwGeom.Point2I(serial_prescan,
                                                     parallel_prescan),
                                     afwGeom.Extent2I(xsize, ysize))
        self.full_segment \
            = afwGeom.Box2I(afwGeom.Point2I(0, 0),
                            afwGeom.Extent2I(xsize + serial_prescan
                                             + serial_overscan,
                                             ysize + parallel_prescan
                                             + parallel_overscan))
        self.prescan = afwGeom.Box2I(afwGeom.Point2I(serial_prescan, 0),
                                     afwGeom.Extent2I(xsize, parallel_prescan))
        self.serial_overscan = afwGeom.Box2I(afwGeom.Point2I(0, parallel_prescan),
                                             afwGeom.Extent2I(serial_prescan, ysize))
        self.parallel_overscan = afwGeom.Box2I(afwGeom.Point2I(0, 0),
                                               afwGeom.Extent2I(0, 0))
        self.gain = float(tokens[7])
        self.bias_level = float(tokens[9])
        self.read_noise = float(tokens[11])
        self.dark_current = float(tokens[13])
        self.crosstalk = np.array([float(x) for x in tokens[21:]])
        self.flip_x = (tokens[5] == '1')
        self.flip_y = (tokens[6] == '-1')
