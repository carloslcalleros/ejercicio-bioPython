from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

# Archivo .gbk de mi escritorio, cambiar a la dirección del nuevo archivo .gbk a leer
filename =  "/mnt/c/Users/carlo/desktop/biopython-notebook/notebooks/data/ls_orchid.gbk"

# Definición de la función summarize_contents
def summarize_contents(filename):
	lista = []
	lista = os.path.split(filename)
	cadena = " "
	cadena = ("file: "+ lista[1] + "\npath: " + lista[0])
	all_records=[]
	records = list(SeqIO.parse(filename, "genbank"))
	cadena += ("\nnum_records: " + str(len(records)))
	cadena += ("\nrecords:")
	for seq_record in SeqIO.parse(filename, "genbank"):
		all_records.append(seq_record.name)
		cadena += ("\n- id:" + str(seq_record.id))
		cadena += ("\nname: " + seq_record.name)
		cadena += ("\ndescription: " + str(seq_record.description))
		cadena += ("\n")
	return cadena
resultado = summarize_contents(filename)
print(resultado)
