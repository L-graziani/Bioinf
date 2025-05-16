import xml.etree.ElementTree as ET
import os   
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

def obtener_mejores_alineaciones(xml_file, num_mejores=10):
    # Parsear el archivo XML
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Inicializar lista para almacenar las alineaciones
    alineaciones = []
    query_sequence = None

    # Recorrer todos los hits en el archivo XML
    for iteration in root.findall(".//Hit"):
        hit_id = iteration.find("Hit_id").text
        hit_def = iteration.find("Hit_def").text
        hit_evalue = float(iteration.find("Hit_hsps/Hsp/Hsp_evalue").text)
        
        # Recorrer los HSPs dentro de cada Hit
        for hsp in iteration.findall(".//Hsp"):
            # Extraer la secuencia alineada de la consulta (qseq) y del hit (hseq)
            secuencia_qseq = hsp.find("Hsp_qseq").text
            secuencia_hseq = hsp.find("Hsp_hseq").text
            
            # Guardar la secuencia de consulta si es la primera vez que la encontramos
            if query_sequence is None:
                query_sequence = secuencia_qseq
            
            # Guardar los datos de la alineación junto con las secuencias
            alineaciones.append((hit_id, secuencia_qseq, secuencia_hseq, hit_evalue))

    # Filtrar las alineaciones para excluir aquellas cuyo hit es igual a la secuencia de consulta
    alineaciones = [alineacion for alineacion in alineaciones if alineacion[1] != query_sequence]

    # Ordenar las alineaciones por el valor E (de menor a mayor)
    alineaciones_ordenadas = sorted(alineaciones, key=lambda x: x[3])

    # Obtener las mejores alineaciones (10 mejores)
    mejores_alineaciones = alineaciones_ordenadas[:num_mejores]

    return query_sequence, mejores_alineaciones

def guardar_en_fasta(query_sequence, alineaciones, output_fasta):
    # Guardar las secuencias en formato FASTA
    with open(output_fasta, "w") as fasta_file:
        # Escribir la secuencia de consulta
        fasta_file.write(f">consulta\n")
        fasta_file.write(f"{query_sequence}\n")

        # Guardar las mejores secuencias de los hits
        for idx, alineacion in enumerate(alineaciones, 1):
            hit_id, _, secuencia_hseq, _ = alineacion
            fasta_file.write(f">hit_{hit_id}\n")
            fasta_file.write(f"{secuencia_hseq}\n")
    
    print(f"Las secuencias alineadas han sido guardadas en '{output_fasta}'")
    return output_fasta
def run_clustalw(fasta_path: str,
                  clustalw_exe: str = 'clustalw2',
                  output_format: str = 'clustal',
                  outfile: str = None) -> AlignIO.MultipleSeqAlignment:
    """
    Realiza un alineamiento múltiple con ClustalW sobre un archivo FASTA.

    Parámetros:
    ----------
    fasta_path : str
        Ruta al archivo FASTA de entrada (con la secuencia query y los 10 mejores resultados de BLAST).
    clustalw_exe : str, opcional
        Nombre o ruta al ejecutable de ClustalW (por defecto 'clustalw2').
    output_format : str, opcional
        Formato de salida para el alineamiento (p. ej. 'clustal', 'fasta', 'phylip').
    outfile : str, opcional
        Ruta del archivo de salida. Si no se especifica, ClustalW generará un archivo .aln en el mismo directorio.

    Retorna:
    -------
    MultipleSeqAlignment
        Objeto con el alineamiento resultante.

    Ejemplo de uso:
    --------------
    >>> alignment = run_clustalw('resultados_blast.fasta')
    >>> print(alignment)
    """
    # Si no se proporciona outfile, ClustalW generará one basado en el nombre del FASTA (.aln)
    if outfile is None:
        base = fasta_path.rsplit('.', 1)[0]+"Alineado"
        outfile = f"{base}.aln"

    # Construir la línea de comando
    cline = ClustalwCommandline(clustalw_exe, infile=fasta_path, outfile=outfile, output=output_format)

    # Ejecutar ClustalW
    stdout, stderr = cline()
    if stderr:
        raise RuntimeError(f"ClustalW error: {stderr}")

    # Leer el alineamiento resultante
    alignment = AlignIO.read(outfile, output_format)
    return alignment


# Ejemplo de uso
script_dir=os.path.dirname(os.path.abspath(__file__))
archivo_xml = os.path.join(script_dir, "NM_001302688.2_frame+2_74-1105.xml")
query_sequence, mejores = obtener_mejores_alineaciones(archivo_xml, num_mejores=10)
clusdir=os.path.join(script_dir, "clustalw2","bin", "clustalw2")
# Guardar las mejores secuencias en un archivo FASTA
fasta=guardar_en_fasta(query_sequence, mejores, "Blast.fasta")
run_clustalw(fasta,clusdir)
