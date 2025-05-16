import xml.etree.ElementTree as ET
import os   
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

def obtener_mejores_alineaciones(xml_file, num_mejores=10):
    # Parsear el archivo XML
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Intentar extraer la secuencia completa de la consulta (si está en el XML)
    # Ajusta la ruta si en tu XML el tag es distinto
    query_sequence = root.findtext(".//Iteration/Iteration_query-sequence")
    
    alineaciones = []

    # Recorrer todos los hits
    for hit in root.findall(".//Hit"):
        hit_id    = hit.findtext("Hit_id")
        evalue    = float(hit.findtext("Hit_hsps/Hsp/Hsp_evalue"))

        for hsp in hit.findall(".//Hsp"):
            qseq = hsp.findtext("Hsp_qseq")
            hseq = hsp.findtext("Hsp_hseq")

            # Si no había query_sequence, la guardo a partir de este primer HSP
            if query_sequence is None:
                query_sequence = qseq

            alineaciones.append((hit_id, qseq, hseq, evalue))

    # Excluir auto‑alineaciones (donde qseq == query_sequence)
    alineaciones = [a for a in alineaciones if a[1] != query_sequence]

    # Ordenar por E‑value y quedarnos con las mejores
    mejores = sorted(alineaciones, key=lambda x: x[3])[:num_mejores]

    return query_sequence, mejores

def guardar_en_fasta(query_sequence, alineaciones, output_fasta):
    with open(output_fasta, "w") as out:
        # Escribo la consulta
        out.write(">consulta\n")
        out.write(query_sequence + "\n")

        # Escribo cada hit
        for idx, (hit_id, qseq, hseq, evalue) in enumerate(alineaciones, start=1):
            # Si quieres también la parte de consulta alineada:

            out.write(f">hit_{idx}|{hit_id}_h e={evalue}\n")
            out.write(hseq + "\n")

    print(f"Guardado en {output_fasta}")
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
