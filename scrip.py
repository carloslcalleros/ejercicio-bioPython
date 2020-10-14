from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

#filename es la direccion de la lista de archivos que la funcion tomará para completar su función
filename= "/mnt/c/Users/carlo/Desktop/BioInf/ls_orchid.gbk"

def summarize_contents(filename):
        all_records=[]
        records = list(SeqIO.parse(filename, "genbank"))
        print ("Path: ", os.path.dirname(filename))
        print("num_records = %i records" % len(records))
        print("\n")
        for seq_record in SeqIO.parse(filename, "genbank"):
                all_records.append(seq_record.name)
                print(counter, ".-")
                print("Name: ", seq_record.name)
                print("ID:",seq_record.id)
                print("Description: ", seq_record.description)
                
summarize_contents(filename)


