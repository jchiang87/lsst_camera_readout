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
        sensor_id = ImageSource.extract_sensor_id(eimage_file)
        self.assertEqual(sensor_id, 'R23_S12')

if __name__ == '__main__':
    unittest.main()
