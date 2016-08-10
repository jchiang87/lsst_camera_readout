"""
Example unit tests for lsst_camera_readout package
"""
import unittest
from desc.lsst_camera_readout import ImageSource

class ImageSourceTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_extract_sensor_id(self):
        eimage_file = '../data/lsst_e_200_f2_R23_S12_E000.fits.gz'
        raft, sensor = ImageSource.extract_sensor_id(eimage_file)
        self.assertEqual(raft, 'R:2,3')
        self.assertEqual(sensor, 'S:1,2')

if __name__ == '__main__':
    unittest.main()
