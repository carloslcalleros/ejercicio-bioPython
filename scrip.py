from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

# Archivo de mi escritorio, se debe de cambiar a la dirección del nuevo archivo .gbk a leer
filename =  "/mnt/c/Users/carlo/desktop/ejercicio-biopython/data/opuntia.fasta"

def summarize_contents(filename):
	lista = []
	lista = os.path.split(filename)
	cadena = " "
	cadena = ("\nfile: "+ lista[1] + "\npath: " + lista[0])
	all_records=[]
	records = list(SeqIO.parse(filename, "genbank"))
	cadena += ("\nnum_records: " + str(len(records)))
	cadena += ("\nrecord(s): ")
	cadena += ("\n----------------------------------------------------------")
	for seq_record in SeqIO.parse(filename, "genbank"):
		all_records.append(seq_record.name)
		cadena += ("\n- id: " + str(seq_record.id))
		cadena += ("\nname: " + seq_record.name)
		cadena += ("\ndescription: " + str(seq_record.description))
	return cadena
if __name__=="__main__":
	resultado = summarize_contents(filename)
	print(resultado)
