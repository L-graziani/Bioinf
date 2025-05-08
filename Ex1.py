from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

Entrez.email = "lgraziani@itba.edu.ar"

def descargar_mrnas_apoe(name,email: str, organism: str="Homo sapiens"):
    Entrez.email = email
    term = f'{name}[symbol] AND "{organism}"[Organism]'
    print(f"🔍 Buscando Gene APOE con: {term}")
    
    h = Entrez.esearch(db="gene", term=term, retmax=1)
    rec = Entrez.read(h); h.close()
    if not rec["IdList"]:
        print("❌ No se encontró el gen.")
        return None
    gene_id = rec["IdList"][0]
    print(f"✔ Gene ID correcto: {gene_id}")
    
    sum_h = Entrez.esummary(db="gene", id=gene_id)
    summary = Entrez.read(sum_h); sum_h.close()
    symbol = summary['DocumentSummarySet']['DocumentSummary'][0]['NomenclatureSymbol']
    if symbol.upper() != name:
        print(f"❌ El Gene ID {gene_id} no corresponde a {name} (sino a {symbol}).")
        return None
    
    link_h = Entrez.elink(dbfrom="gene",
                         db="nucleotide",
                         id=gene_id,
                         linkname="gene_nuccore_refseqrna")
    links = Entrez.read(link_h); link_h.close()
    try:
        mrna_ids = [l["Id"] for l in links[0]["LinkSetDb"][0]["Link"]]
    except (IndexError, KeyError):
        print("❌ No hay mRNAs RefSeq asociados a este Gene ID.")
        return None

    print(f"✔ IDs de mRNA RefSeq encontrados: {mrna_ids}")

    fetch_h = Entrez.efetch(db="nucleotide",
                            id=",".join(mrna_ids),
                            rettype="gb",
                            retmode="text")
    filename = f"{name}_mrnas_refseq.gb"
    with open(filename, "w") as out:
        out.write(fetch_h.read())
    fetch_h.close()
    print(f"💾 Guardado en '{filename}'")
    return filename

def convert_gb_to_fasta(gb_filename: str) -> None:

    # Determinar ruta absoluta y nombre de salida
    base, _ = os.path.splitext(gb_filename)
    fasta_filename = f"{base}.fasta"

    try:
        # Parsear todas las secuencias desde el archivo GenBank
        records = list(SeqIO.parse(gb_filename, "genbank"))
        if not records:
            print(f"No se encontraron registros en '{gb_filename}'")
            return

        # Escribir todas las secuencias en formato FASTA
        SeqIO.write(records, fasta_filename, "fasta")
        print(f"Se han guardado {len(records)} secuencias en '{fasta_filename}'")

    except Exception as e:
        print(f"Error al procesar '{gb_filename}': {e}")
    return fasta_filename

def filtrar_fasta(fasta_path):
    """
    Lee un archivo .fasta, filtra las secuencias que tienen "transcript variant 1"
    en el encabezado y escribe un nuevo archivo con dichas secuencias.

    Parámetros:
        fasta_path (str): Ruta al archivo .fasta de entrada.

    Crea:
        Un archivo .fasta nuevo con el sufijo _variant1.fasta.
    """
    salida_path = fasta_path.rsplit('.', 1)[0] + '_variant1.fasta'

    with open(fasta_path, 'r') as f, open(salida_path, 'w') as salida:
        guardar = False
        buffer = ""
        for linea in f:
            if linea.startswith(">"):
                if guardar and buffer:
                    salida.write(buffer)
                    buffer = ""
                if "transcript variant 1" in linea:
                    guardar = True
                    buffer = linea
                else:
                    guardar = False
            elif guardar:
                buffer += linea
        if guardar and buffer:
            salida.write(buffer)

    print(f"Archivo filtrado creado: {salida_path}")
    return salida_path

def traducir_6_frames(input_fasta, output_fasta):
    records_orf = []
    stop_codons = {"TAA", "TAG", "TGA"}

    for record in SeqIO.parse(input_fasta, "fasta"):
        seq = record.seq.upper()

        # Inicializar diccionario para guardar la mejor ORF por cada uno de los seis marcos
        mejores_orfs = {}
        for i in range(3):
            mejores_orfs[f'+{i+1}'] = (None, 0)
            mejores_orfs[f'-{i+1}'] = (None, 0)

        # Procesar los tres marcos en sentido directo
        for frame in range(3):
            seq_frame = seq[frame:]
            mejor_prot, max_len = None, 0
            for i in range(0, len(seq_frame) - 2, 3):
                if str(seq_frame[i:i+3]) == "ATG":
                    for j in range(i+3, len(seq_frame) - 2, 3):
                        codon2 = str(seq_frame[j:j+3])
                        if codon2 in stop_codons:
                            orf_dna = seq_frame[i:j+3]
                            orf_prot = orf_dna.translate(to_stop=False)
                            if len(orf_prot) > max_len:
                                max_len = len(orf_prot)
                                mejor_prot = (orf_prot, frame, i, j+3)
                            break
            if mejor_prot:
                prot, frm, start, end = mejor_prot
                id_orf = f"{record.id}_frame+{frm+1}_{frm+start+1}-{frm+end}"
                mejores_orfs[f'+{frm+1}'] = (Seq(str(prot)), id_orf)

        # Procesar los tres marcos en reverso complementario
        seq_rc = seq.reverse_complement()
        for frame in range(3):
            seq_frame = seq_rc[frame:]
            mejor_prot, max_len = None, 0
            for i in range(0, len(seq_frame) - 2, 3):
                if str(seq_frame[i:i+3]) == "ATG":
                    for j in range(i+3, len(seq_frame) - 2, 3):
                        codon2 = str(seq_frame[j:j+3])
                        if codon2 in stop_codons:
                            orf_dna = seq_frame[i:j+3]
                            orf_prot = orf_dna.translate(to_stop=False)
                            if len(orf_prot) > max_len:
                                max_len = len(orf_prot)
                                mejor_prot = (orf_prot, frame, i, j+3)
                            break
            if mejor_prot:
                prot, frm, start, end = mejor_prot
                id_orf = f"{record.id}_frame-{frm+1}_{frm+start+1}-{frm+end}"
                mejores_orfs[f'-{frm+1}'] = (Seq(str(prot)), id_orf)

        # Agregar las mejores ORFs de los seis marcos
        for key, (prot_seq, orf_id) in mejores_orfs.items():
            if prot_seq:
                records_orf.append(SeqRecord(prot_seq, id=orf_id, description=""))

    # Guardar las ORFs en un nuevo archivo FASTA
    SeqIO.write(records_orf, output_fasta, "fasta")
    print(f"Extracción de ORFs más largas completa. Archivo generado: {output_fasta}")


# Parámetros: modifica aquí según tu búsqueda
term = 'APOE'
email = "lgraziani@itba.edu.ar"
organism="Homo sapiens"

asd=descargar_mrnas_apoe(term,email,organism)
asd=convert_gb_to_fasta(asd)
asd=filtrar_fasta(asd)
traducir_6_frames(asd, "A_"+asd)

