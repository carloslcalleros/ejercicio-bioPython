import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

#filename es la direccion de la lista de archivos que la funcion tomará para completar su función
filename = "/mnt/c/Users/carlo/Desktop/ejercicio-biopython/ls_orchid.gbk"
lista = []

def summarize_contents(filename):
        lista = os.path.split(filename)
        print("file:", all_records[1])
        all_records = []
        record = list(SeqIO.parse(filename,"genbank"))
        print("path: ", os.path.dirname(filename))
        print("num_records = %i records" % len(record))
        print("record:")
        for seq_record in SeqIO.parse(filename, "genbank"):
                all_records.append(seq_record.name)
                print("ID:",seq_record.id)
                print ("Name: ", seq_record.name)
                print("Description: ", seq_record.description)
                print("\n\n")
summarize_contents(filename)


