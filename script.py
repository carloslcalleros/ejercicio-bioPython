from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

# NC_002703.gbk solo es un archivo dentro de la carpeta data, se puede cambiar por otros archivos que esten dentro de esta carpeta 
filename = os.path.abspath("data/NC_002703.gbk") 

#Secuencias ejemplo para probar la función concatenate_and_get_reverse_of_complement
seq1="ATCGACTC"
seq2="CAaCATacgagaaatAGG"

#----------------------------FUNCIÓN 1-----------------------------------
def summarize_contents(filename):
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
	dictionary = {}
	dictionary['File:'] = File_List[1]
	dictionary['Path:'] = File_List[0]
	dictionary['Num_records:'] = len(record)
	
	dictionary['Names:'] = []
	dictionary['IDs:'] = []
	dictionary['Descriptions'] = []
	for seq_rcd in SeqIO.parse(filename,type_file):
		dictionary['Names:'].append(seq_rcd.name)
		dictionary['IDs:'].append(seq_rcd.id)
		dictionary['Descriptions'].append(seq_rcd.description)
	return dictionary

if __name__ == "__main__":
	resultados = summarize_contents(filename)
	print(resultados)


#----------------------------FUNCIÓN 2-----------------------------------

def concatenate_and_get_reverse_of_complement(secuencia1,secuencia2):
	concatenate = Seq(secuencia1 + secuencia2)
	reverse = concatenate.reverse_complement()
	return reverse.upper()

if __name__ == "__main__":
	resultado = concatenate_and_get_reverse_of_complement(seq1,seq2)
	print(resultado)
