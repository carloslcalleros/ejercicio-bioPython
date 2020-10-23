import unittest
import scrip

class Prueba(unittest.TestCase):
	def test_summarize_contents(self):
		s = summarize_contents("data/AF323668.gbk")
		self.assertEqual("\nfile: AF323668.gbk\npath: /mnt/c/Users/carlo/desktop/biopython-notebook/notebooks/data\nnum_records: 1\nrecord(s):\n----------------------------------------------------------\n- id:AF323668.1\nname: AF323668\ndescription: Bacteriophage bIL285, complete genome",s)
