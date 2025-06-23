import xml.etree.ElementTree as ET
import os
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO


def obtener_mejores_alineaciones(xml_file: str, num_mejores: int = 10):
    """
    Parsea un archivo XML de resultados de BLAST y extrae las mejores alineaciones.

    Parámetros:
    ----------
    xml_file : str
        Ruta al archivo XML de entrada.
    num_mejores : int, opcional
        Número de mejores alineaciones a devolver (por defecto 10).

    Retorna:
    --------
    tuple
        Secuencia de consulta (str) y lista de alineaciones:
        [(hit_id, qseq, hseq, evalue), ...].
    """
    # Parseamos el XML y obtenemos la raíz
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Intentamos extraer la secuencia de consulta original
    query_sequence = root.findtext(".//Iteration/Iteration_query-sequence")
    alineaciones = []

    # Recorremos todos los hits del BLAST
    for hit in root.findall(".//Hit"):
        hit_id = hit.findtext("Hit_id")
        evalue = float(hit.findtext("Hit_hsps/Hsp/Hsp_evalue"))

        # Por cada HSP dentro del hit
        for hsp in hit.findall(".//Hsp"):
            qseq = hsp.findtext("Hsp_qseq")
            hseq = hsp.findtext("Hsp_hseq")

            # Si no se extrajo la query originalmente, la tomamos de aquí
            if query_sequence is None:
                query_sequence = qseq

            alineaciones.append((hit_id, qseq, hseq, evalue))

    # Excluimos auto-alineaciones (qseq == query_sequence)
    alineaciones = [a for a in alineaciones if a[1] != query_sequence]

    # Ordenamos por e-value ascendente y seleccionamos los mejores
    mejores = sorted(alineaciones, key=lambda x: x[3])[:num_mejores]

    return query_sequence, mejores


def guardar_en_fasta(query_sequence: str, alineaciones: list, output_fasta: str) -> str:
    """
    Escribe la secuencia de consulta y las mejores alineaciones en un archivo FASTA.

    Parámetros:
    ----------
    query_sequence : str
        Secuencia de consulta.
    alineaciones : list
        Lista de tuplas (hit_id, qseq, hseq, evalue).
    output_fasta : str
        Nombre o ruta del archivo FASTA de salida.

    Retorna:
    --------
    str
        Ruta del archivo FASTA generado.
    """
    with open(output_fasta, "w") as out:
        # Secuencia de consulta
        out.write(">consulta\n")
        out.write(query_sequence + "\n")

        # Cada hit alineado
        for idx, (hit_id, qseq, hseq, evalue) in enumerate(alineaciones, start=1):
            header = f">hit_{idx}|{hit_id} e={evalue}"
            out.write(header + "\n")
            out.write(hseq + "\n")

    print(f"Guardado en {output_fasta}")
    return output_fasta


def run_clustalw(
    fasta_path: str,
    clustalw_exe: str = 'clustalw2',
    output_format: str = 'clustal',
    outfile: str = None
) -> AlignIO.MultipleSeqAlignment:
    """
    Realiza un alineamiento múltiple con ClustalW sobre un archivo FASTA.

    Parámetros:
    ----------
    fasta_path : str
        Ruta al archivo FASTA de entrada (query + mejores hits).
    clustalw_exe : str, opcional
        Ejecutable de ClustalW (por defecto 'clustalw2').
    output_format : str, opcional
        Formato de salida ('clustal', 'fasta', 'phylip', etc.).
    outfile : str, opcional
        Ruta del archivo de salida. Si no se especifica, se genera un .aln.

    Retorna:
    -------
    AlignIO.MultipleSeqAlignment
        Objeto con el alineamiento resultante.
    """
    # Definimos el archivo de salida si no se proporcionó
    if outfile is None:
        base = fasta_path.rsplit('.', 1)[0] + "_alineado"
        outfile = f"{base}.aln"

    # Construimos y ejecutamos la línea de comando de ClustalW
    cline = ClustalwCommandline(
        clustalw_exe,
        infile=fasta_path,
        outfile=outfile,
        output=output_format
    )
    stdout, stderr = cline()

    if stderr:
        raise RuntimeError(f"ClustalW error: {stderr}")

    # Leemos y devolvemos el alineamiento generado
    alignment = AlignIO.read(outfile, output_format)
    return alignment


# Ejecución principal del script
script_dir = os.path.dirname(os.path.abspath(__file__))
archivo_xml = os.path.join(
    script_dir,
    "NM_001302688.2_frame+2_74-1105.xml"
)

# Obtenemos la secuencia y las mejores alineaciones
query_sequence, mejores = obtener_mejores_alineaciones(
    archivo_xml,
    num_mejores=10
)

# Ruta al ejecutable de ClustalW
clusdir = os.path.join(
    script_dir,
    "clustalw2",
    "bin",
    "clustalw2"
)

# Guardamos en FASTA y ejecutamos ClustalW
fasta = guardar_en_fasta(query_sequence, mejores, "Blast.fasta")
run_clustalw(fasta, clusdir)
