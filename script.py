from Bio import SeqIO
from Bio.Seq import Seq,reverse_complement
from Bio.SeqRecord import SeqRecord
import os
from Bio.Data import CodonTable



# NC_002703.gbk solo es un archivo dentro de la carpeta data, se puede cambiar por otros archivos que esten dentro de esta carpeta 
filename = os.path.abspath("data/NC_002703.gbk") 

#Secuencias ejemplo para probar la función concatenate_and_get_reverse_of_complement
seq1="GATCA"
seq2="GATCA"
DNA = "ATGGCCTAGTTGAGCAgTATATAA"
h="ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

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
	concatenate = Seq(secuencia1+secuencia2)
	reverse = concatenate.reverse_complement()
	return reverse.upper()

if __name__ == "__main__":
	resultados = concatenate_and_get_reverse_of_complement(seq1,seq2)
	print(resultados)

#--------------------------FUNCIÓN 3------------------------------------

def print_protein_and_codons_using_standard_table(seq):
	sequence = Seq(seq)	
	dictionary = {'mRNA': sequence.upper().transcribe(),'proteins':[], 'stop_codons': []}
	aminoacids_seq = sequence.translate(table = 1)	
# Exportando codones de inicio y de parada de la tabla estandar
	Table = CodonTable.unambiguous_dna_by_name["Standard"]
	start_codon = False
	stop_codon = False
	begin = None
	end = None
	c = 0
	while c < len(aminoacids_seq):
		if (sequence[c*3:c*3+3].upper() == "TTG") or (sequence[c*3:c*3+3].upper() == "CTG") or (sequence[c*3:c*3+3].upper() == "ATG"):
			start_codon = True
			begin = c
			if c+1 == len(aminoacids_seq):
				dictionary['proteins'].append(aminoacids_seq[c:])
				break
			k = c+1
			while k < len(aminoacids_seq):
				if (sequence[k*3:k*3+3].upper() == "TAA") or (sequence[k*3:k*3+3].upper() == "TGA") or (sequence[k*3:k*3+3].upper() == "TAG"):
					stop_codon = True
					end = k
					dictionary['proteins'].append(aminoacids_seq[c:k])
					dictionary['stop_codons'].append(sequence[k*3:k*3+3])
					start_codon, stop_codon = False, False
					c = k
					break
				k += 1
		if(start_codon == True and stop_codon == False):
			dictionary['proteins'].append(aminoacids_seq[c:])
			break
		c += 1

	if (dictionary['proteins'] == []):
		dictionary['proteins'] = "No proteins were found in this sequence"
	if dictionary['stop_codons'] == []:
		dictionary['stop_codons'] = "No stop codons were found un this sequence"
	return dictionary

if __name__ == "__main__":
	resultado = print_protein_and_codons_using_standard_table(h)
	print(resultado)

#--------------------------- FUNCIÓN 4 -----------------------------------

def print_protein_and_codons_using_mitocondrial_yeast_table(seq):
	sequence = Seq(seq)	
	dictionary = {'mRNA': sequence.upper().transcribe(),'proteins':[], 'stop_codons': []}
	aminoacids_seq = sequence.translate(table = 3)	
# Exportando codones de inicio y de parada de la tabla de mictocondria de levadura
	Table = CodonTable.unambiguous_dna_by_name["Yeast Mitochondrial"]
	start_codon = False
	stop_codon = False
	begin = None
	end = None
	c = 0
	while c < len(aminoacids_seq):
		if (sequence[c*3:c*3+3].upper() == "ATA") or (sequence[c*3:c*3+3].upper() == "ATG") or (sequence[c*3:c*3+3].upper() == "GTG"):
			start_codon = True
			begin = c
			if c+1 == len(aminoacids_seq):
				dictionary['proteins'].append(aminoacids_seq[c:])
				break
			k = c+1
			while k < len(aminoacids_seq):
				if (sequence[k*3:k*3+3].upper() == "TAA") or (sequence[k*3:k*3+3].upper() == "TAG"):
					stop_codon = True
					end = k
					dictionary['proteins'].append(aminoacids_seq[c:k])
					dictionary['stop_codons'].append(sequence[k*3:k*3+3])
					start_codon, stop_codon = False, False
					c = k
					break
				k += 1
		if(start_codon == True and stop_codon == False):
			dictionary['proteins'].append(aminoacids_seq[c:])
			break
		c += 1

	if (dictionary['proteins'] == []):
		dictionary['proteins'] = "No proteins were found in this sequence"
	if dictionary['stop_codons'] == []:
		dictionary['stop_codons'] = "No stop codons were found un this sequence"
	return dictionary

if __name__ == "__main__":
	resultado = print_protein_and_codons_using_mitocondrial_yeast_table(h)
	print(resultado)

#--------------------FUNCIÓN 5---------------------------------
def extract_sequences(file,formato): 

	Direction = os.path.abspath(file)
	File_Extension = os.path.splitext(file)
	if (formato != "genbank"):
		return "Error el formato de salida debe de ser genbank"
	else:
		if (File_Extension[1] != ".fasta"):
			return "Error la extencion del archivo debe de ser .fasta"
		else:
			GBK_File= "auxiliar.gbk"
			SeqIO.convert(file, "fasta" , GBK_File , formato, molecule_type= "DNA")
			Direction = os.path.abspath(GBK_File)
			GBK_Records = list(SeqIO.parse( Direction , formato))
			c=0
			while c < len(GBK_Records):
				name_file = "sequence" + str(c+1) + ".gbk"
				file = open(name_file, "w")
				file.write(str(GBK_Records[c].format("genbank")))
				file.close()
				c+=1
			os.remove (GBK_File)
if __name__ == "__main__":
	extract_sequences("data/sequences.fasta","genbank")
	
#-----------------------FUNCION 6----------------------------
def extract_sequences_revcomp(file):

	Extention = os.path.splitext(file)
	if (Extention[1] != ".fasta"):
		return "Error el formato del archivo debe de ser .fasta"
	else:
		Direction = os.path.abspath(file)
		File_records = list(SeqIO.parse(Direction, "fasta"))
	
		name_file = "sequence_rev_comp" + ".fasta"
		file = open(name_file, "w")
	
		for i in range(len(File_records)):
			Nucleotides = Seq(str(File_records[i].seq))
			file.write(">" + File_records[i].id)
			R = Nucleotides.reverse_complement()
			file.write(str( os.linesep + R + os.linesep))
		file.close()

if __name__ == "__main__":
	extract_sequences_revcomp("data/sequences.fasta")