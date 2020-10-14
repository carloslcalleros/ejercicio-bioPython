from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
#filename es la direccion del archivo o lista de archivos que la funcion tomará para completar su función
filename = "/mnt/c/Users/carlo/Desktop/Bionfo/biopython-notebook/notebooks/data/ls_orchid.gbk"

def summarize_contents(filename):
        all_records=[]
        records = list(SeqIO.parse(filename, "genbank"))        
        print ("Path: ", os.path.dirname(filename))
        print("num_records = %i records" % len(records))
        print("\n")
        for seq_record in SeqIO.parse(filename, "genbank"):
                all_records.append(seq_record.name)
                print("Name: ", seq_record.name)
                print("ID:",seq_record.id)
                print("Location:")
                for seq_feature in seq_record.features :
                        print('Start: %d, Stop: %d'%(int(seq_feature.location.start),int(seq_feature.location.end)))
                print("\n\n")

summarize_contents(filename)
