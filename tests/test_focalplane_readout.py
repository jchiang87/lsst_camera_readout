import os
import unittest
import lsst.obs.lsstSim as lsstSim
import lsst.utils as lsstUtils
from desc.lsst_camera_readout import FocalPlaneReadout

class FocalPlaneReadoutTestCase(unittest.TestCase):
    def setUp(self):
        self.camera = lsstSim.LsstSimMapper().camera

    def tearDown(self):
        del self.camera

    def test_read_segmentation_txt(self):
        seg_file = os.path.join(lsstUtils.getPackageDir('obs_lsstSim'),
                                'description', 'segmentation.txt')
        readout_props = FocalPlaneReadout.read_phosim_seg_file(seg_file)
        raft = 'R:2,2'
        ccd = 'S:1,1'
        sensor = self.camera[' '.join((raft, ccd))]
        sensor_id = readout_props.sensor_id(raft, ccd)
        sensor_props = readout_props.get_sensor(sensor_id)
        self.assertEqual(sensor.getBBox().getWidth(), sensor_props.width)
        self.assertEqual(sensor.getBBox().getHeight(), sensor_props.height)
        self.assertEqual(len(sensor), sensor_props.num_amps)
        amp_num = 0
        for col in '01':
            for row in '01234567':
                amp_num += 1
                chan = 'C:%s,%s' % (col, row)
                amp_id = readout_props.amp_id(raft, ccd, chan)
                amp_props = readout_props.get_amp(amp_id)
                self.assertEqual(amp_props.gain, 1.7)
                self.assertEqual(amp_props.bias_level, 1000.)
                self.assertEqual(amp_props.read_noise, 7.)
                self.assertEqual(amp_props.dark_current, 0.02)
                crosstalk = [0]*len(sensor)
                crosstalk[amp_num - 1] = 1.0
                self.assertListEqual(amp_props.crosstalk.tolist(), crosstalk)

                # Check pixel geometry against cameraGeom.
                amp = sensor[chan[2:]]
                self.assertEqual(amp.getRawDataBBox(), amp_props.imaging)
#                self.assertEqual(amp.getRawBBox(), amp_props.full_segment)
#                self.assertEqual(amp.getRawPrescanBBox(), amp_props.prescan)
#                self.assertEqual(amp.getRawHorizontalOverscanBBox(),
#                                 amp_props.serial_overscan)
#                self.assertEqual(amp.getRawVerticalOverscanBBox(),
#                                 amp_props.parallel_overscan)
                self.assertEqual(amp.getRawFlipX(), amp_props.flip_x)
                self.assertEqual(amp.getRawFlipY(), amp_props.flip_y)

if __name__ == '__main__':
    unittest.main()
