from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
from Bio.Data import CodonTable

# NC_002703.gbk solo es un archivo dentro de la carpeta data, se puede cambiar por otros archivos que esten dentro de esta carpeta 
filename = os.path.abspath("data/NC_002703.gbk") 

#Secuencias ejemplo para probar la función concatenate_and_get_reverse_of_complement
seq1="TTTTTAAAAUUUGGGKSKSKSKn"
seq2="uATG"
DNA = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
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
#	print(resultados)


#----------------------------FUNCIÓN 2-----------------------------------

def concatenate_and_get_reverse_of_complement(secuencia1,secuencia2):
	concatenate = Seq(secuencia1+secuencia2)
	reverse = concatenate.reverse_complement()
	return reverse.upper()

if __name__ == "__main__":
	resultado = concatenate_and_get_reverse_of_complement(seq1,seq2)
	print(resultado)

#--------------------------FUNCIÓN 3------------------------------------
def print_protein_and_codons_using_standard_table(sec):
	# Convertir secuencia en objeto Seq
	secuencia = Seq(sec)	
	# Crea diccionario
	diccionario = {}
	# mRNA, proteínas y codones de paro
	diccionario['mRNA'] = secuencia.transcribe()
	diccionario['proteins'] = []
	diccionario['stop_codons'] = []
		
	# Proceso de búsqueda de proteínas
	aminoacidos = secuencia.translate(table = 1, stop_symbol = "-")
	posibles_proteinas = aminoacidos.split('-')
	lista_proteinas = []
	for i in range(len(posibles_proteinas)):
		empieza = posibles_proteinas[i].find('M')
		if empieza != -1:
			lista_proteinas.append(posibles_proteinas[i][empieza:])	
	diccionario['proteins'] = lista_proteinas
	if (diccionario['proteins'] == []):
		diccionario['proteins'] = "Not found proteins"
		
	# Proceso de búsqueda de stop codons
	# Exportando codones de inicio y de parada
	codons_table = CodonTable.unambiguous_dna_by_id[1]
	start_codons_list = codons_table.start_codons
	stop_codons_list = codons_table.stop_codons
	start_codon, stop_codon = False, False
	i = 0
	pos_inicio = None
	pos_final = None
	while i < len(aminoacidos):
		if secuencia[i*3:i*3+3].upper() in start_codons_list:
			start_codon = True
			pos_inicio = i
			if i+1 == len(aminoacidos):
				diccionario['proteins'].append(aminoacidos[i:])
				break
			j = i+1
			while j < len(aminoacidos):
				if secuencia[j*3:j*3+3].upper() in stop_codons_list:
					stop_codon = True
					pos_final = j
					diccionario['proteins'].append(aminoacidos[i:j])
					diccionario['stop_codons'].append(secuencia[j*3:j*3+3])
					start_codon, stop_codon = False, False
					i = j
					break
				j += 1
		if(start_codon == True and stop_codon == False):
			diccionario['proteins'].append(aminoacidos[i:])
			break
		i += 1

	if diccionario['stop_codons'] == []:
		diccionario['stop_codons'] = "Not found stop codons"
	return diccionario

if __name__ == "__main__":
	resultado = print_protein_and_codons_using_standard_table(DNA)
#	print(resultado)
	
