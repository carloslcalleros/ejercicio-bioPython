from Bio.Seq import Seq
from Bio.SeqFeauture import SeqFeature,FeatureLocation
from Bio import SeqIO 
from Bio.SeqRecord import SeqRecord
def summarize_contents(filename):
        prueba = SeqIO.read(filename,"genbank")
        print(record)
summarize_contents(filename)
