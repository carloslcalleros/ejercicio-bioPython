from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

# Archivo de mi escritorio, se debe de cambiar a la dirección del nuevo archivo .gbk o .fasta a leer
filename = "/mnt/c/users/carlo/desktop/ejercicio-biopython/data/NC_002703.gbk" 

def summarize_contents2(filename):
	File_List = []
	File_Extension = []

	File_List = os.path.split(filename)
	File_Extension = os.path.splitext(filename)

	if(File_Extension[1] == ".gbk"):
		type_file= "genbank"
	else:
		type_file= "fasta"

	record = list(SeqIO.parse(filename, type_file))
#Diccionario
	d = {}
	d['File:'] = File_List[1]
	d['Path:'] = File_List[0]
	d['Num_records:'] = len(record)
	
	d['Names:'] = []
	d['IDs:'] = []
	d['Descriptions'] = []
	for seq_rcd in SeqIO.parse(filename,type_file):
		d['Names:'].append(seq_rcd.name)
		d['IDs:'].append(seq_rcd.id)
		d['Descriptions'].append(seq_rcd.description)
	return d

if __name__ == "__main__":
	resultados = summarize_contents2(filename)
	print(resultados)
