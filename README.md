# ejercicio-biopython
Descripcion del funcionamiento del scrip.py: importamos las bibliotecas necesarias para utilizar biopython, se define una función llamada summarize_contents que imprime un resumen 
del contenido de un archivo en formato .gbk (genbank), el output de esta fución es el siguiente: file: [nombre de archivo] path: [ruta al archivo] num_records: [numero de registros] records: - id: [id del primer registro] name: [nombre] description: [descripción] - id: [id del segundo registro] name: [nombre] description: [descripción] - ...

Descripción de la función summarize_contents: recibe el nombre de un archivo que se asume que está en formato genbank, se utiliza una lista, que guarda la ruta y el nombre del archivo. Después se imprime el número de records, y cada uno de ellos con un ciclo "for" que itera hasta el último registro. Imprimiendo  en el formato antes dicho.

