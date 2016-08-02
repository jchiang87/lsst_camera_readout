"""
Example unit tests for lsst_camera_readout package
"""
import unittest
import desc.lsst_camera_readout

class lsst_camera_readoutTestCase(unittest.TestCase):
    def setUp(self):
        self.message = 'Hello, world'
        
    def tearDown(self):
        pass

    def test_run(self):
        foo = desc.lsst_camera_readout.lsst_camera_readout(self.message)
        self.assertEquals(foo.run(), self.message)

    def test_failure(self):
        self.assertRaises(TypeError, desc.lsst_camera_readout.lsst_camera_readout)
        foo = desc.lsst_camera_readout.lsst_camera_readout(self.message)
        self.assertRaises(RuntimeError, foo.run, True)

if __name__ == '__main__':
    unittest.main()
