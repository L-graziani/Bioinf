# six_frame_translation.py

import sys
import os
import tempfile
import subprocess
from Bio import SeqIO


def six_frame_translation(input_fasta, output_fasta, emboss_dir=None, table=1, clean=True):
    """
    Traduce la secuencia de nucleótidos en las seis fases de lectura usando EMBOSS Transeq vía subprocess.

    :param input_fasta: Ruta al fasta de entrada (nucleótidos)
    :param output_fasta: Ruta al fasta de salida (proteínas)
    :param emboss_dir: Directorio donde está instalado EMBOSS (si no está en PATH)
    :param table: Código de tabla de traducción (default 1: estándar)
    :param clean: Si True, elimina archivos temporales.
    """

    # Crear un archivo temporal para la salida de Transeq
    tmp_fd, tmp_path = tempfile.mkstemp(suffix=".faa")
    os.close(tmp_fd)

    # Determinar ruta al ejecutable transeq
    transeq_exec = os.path.join(emboss_dir, 'transeq')

    # Construir el comando de Transeq
    cmd = [
        transeq_exec,
        '-sequence', input_fasta,
        '-outseq', tmp_path,
        '-table', str(table),
        '-frame', '6'
    ]

    # Ejecutar Transeq
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.stderr.write(f"Error al ejecutar transeq: {result.stderr}\n")
        if not clean:
            sys.stderr.write(f"Archivo temporal conservado en: {tmp_path}\n")
        sys.exit(result.returncode)

    # Leer el resultado y escribir en el fasta de salida
    proteins = list(SeqIO.parse(tmp_path, "fasta"))
    if not proteins:
        sys.stderr.write("No se encontraron secuencias traducidas en la salida de Transeq.\n")
        sys.exit(1)

    SeqIO.write(proteins, output_fasta, "fasta")
    print(f"Traducción completada. Se han generado {len(proteins)} secuencias proteicas.")

    # Limpiar temporal si corresponde
    if clean:
        os.remove(tmp_path)



def contar_secuencias_fasta(fasta_path):
    """Cuenta el número de secuencias en un archivo FASTA"""
    count = 0
    secuencias_info = []
    
    with open(fasta_path, 'r') as f:
        current_header = None
        current_seq_length = 0
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    secuencias_info.append((current_header, current_seq_length))
                current_header = line[1:]  # Remove '>'
                current_seq_length = 0
                count += 1
            elif line and not line.startswith('>'):
                current_seq_length += len(line)
        
        # Add the last sequence
        if current_header:
            secuencias_info.append((current_header, current_seq_length))
    
    return count, secuencias_info

def extraer_secuencias_individuales(fasta_path):
    """Extrae cada secuencia del FASTA como archivos temporales individuales"""
    secuencias = []
    current_header = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header and current_seq:
                    secuencias.append((current_header, ''.join(current_seq)))
                current_header = line
                current_seq = []
            elif line:
                current_seq.append(line)
        
        # Add the last sequence
        if current_header and current_seq:
            secuencias.append((current_header, ''.join(current_seq)))
    
    return secuencias
    """Extrae cada secuencia del FASTA como archivos temporales individuales"""
    secuencias = []
    current_header = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header and current_seq:
                    secuencias.append((current_header, ''.join(current_seq)))
                current_header = line
                current_seq = []
            elif line:
                current_seq.append(line)
        
        # Add the last sequence
        if current_header and current_seq:
            secuencias.append((current_header, ''.join(current_seq)))
    
    return secuencias

def analizar_secuencias_individuales(fasta_path, prosite_dat_path, emboss_bin_path):
    """Analiza cada secuencia del FASTA por separado para asegurar que todas se procesen"""
    # First try the normal approach
    try:
        resultado_normal = analizar_motifs_dominios(fasta_path, prosite_dat_path, emboss_bin_path)
        
        # Check if results contain multiple sequences by counting sequence headers in output
        sequence_count_in_results = resultado_normal.count('Sequence:') + resultado_normal.count('# Sequence')
        fasta_sequence_count, _ = contar_secuencias_fasta(fasta_path)
        
        if sequence_count_in_results >= fasta_sequence_count:
            print(f"✓ Método normal procesó todas las secuencias ({sequence_count_in_results}/{fasta_sequence_count})")
            return resultado_normal
        else:
            print(f"⚠ Método normal solo procesó {sequence_count_in_results}/{fasta_sequence_count} secuencias")
            print("Intentando método individual...")
    except Exception as e:
        print(f"Método normal falló: {e}")
        print("Intentando método individual...")
    
    # Individual sequence approach
    secuencias = extraer_secuencias_individuales(fasta_path)
    resultados_combinados = []
    
    patmatmotifs_path = os.path.join(emboss_bin_path, "patmatmotifs")
    
    for i, (header, seq) in enumerate(secuencias, 1):
        print(f"Procesando secuencia {i}/{len(secuencias)}: {header[:50]}...")
        
        # Create temporary file for this sequence
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".fasta") as temp_fasta:
            temp_fasta.write(f"{header}\n")
            # Write sequence in lines of 80 characters
            for j in range(0, len(seq), 80):
                temp_fasta.write(seq[j:j+80] + "\n")
            temp_fasta_path = temp_fasta.name
        
        with tempfile.NamedTemporaryFile(delete=False, suffix=".out") as temp_out:
            output_path = temp_out.name
        
        try:
            subprocess.run([
                patmatmotifs_path,
                "-sequence", temp_fasta_path,
                "-outfile", output_path,
                "-full", "Y",
                "-auto", "Y"
            ], check=True)
            
            with open(output_path, "r") as result_file:
                resultado_individual = result_file.read()
                resultados_combinados.append(f"\n{'='*60}\n{header}\n{'='*60}\n{resultado_individual}")
        
        except subprocess.CalledProcessError as e:
            resultados_combinados.append(f"\n{'='*60}\n{header}\n{'='*60}\nError procesando esta secuencia: {e}\n")
        
        finally:
            # Clean up temporary files
            if os.path.exists(temp_fasta_path):
                os.remove(temp_fasta_path)
            if os.path.exists(output_path):
                os.remove(output_path)
    
    return "\n".join(resultados_combinados)

def analizar_motifs_dominios(fasta_path, prosite_dat_path, emboss_bin_path):
    if not os.path.isfile(fasta_path):
        raise FileNotFoundError(f"Archivo FASTA no encontrado: {fasta_path}")
    if not os.path.isfile(prosite_dat_path):
        raise FileNotFoundError(f"Archivo prosite.dat no encontrado: {prosite_dat_path}")
    if not os.path.isdir(emboss_bin_path):
        raise NotADirectoryError(f"Directorio EMBOSS bin no válido: {emboss_bin_path}")
    
    # Get the directory containing prosite.dat
    prosite_dir = os.path.dirname(prosite_dat_path)
    
    # Set EMBOSS_DATA environment variable
    os.environ["EMBOSS_DATA"] = prosite_dir
    
    # Check if prosextract exists
    prosextract_path = os.path.join(emboss_bin_path, "prosextract")
    if not os.path.isfile(prosextract_path):
        raise FileNotFoundError(f"No se encontró 'prosextract' en: {prosextract_path}")
    
    # Create PROSITE subdirectory if it doesn't exist
    prosite_output_dir = os.path.join(prosite_dir, "PROSITE")
    if not os.path.exists(prosite_output_dir):
        print(f"Creando directorio PROSITE: {prosite_output_dir}")
        os.makedirs(prosite_output_dir, exist_ok=True)
    
    # Run prosextract with the directory path, not the file path
    try:
        print(f"Ejecutando prosextract en directorio: {prosite_dir}")
        subprocess.run([
            prosextract_path,
            "-prositedir", prosite_dir
        ], check=True, cwd=prosite_dir)
        print("prosextract ejecutado exitosamente")
    except subprocess.CalledProcessError as e:
        print(f"Error ejecutando prosextract: {e}")
        # Try alternative approach without explicit -prositedir flag
        try:
            print("Intentando enfoque alternativo...")
            subprocess.run([prosextract_path], check=True, cwd=prosite_dir)
            print("prosextract ejecutado exitosamente (enfoque alternativo)")
        except subprocess.CalledProcessError as e2:
            # Try one more approach with explicit output directory
            try:
                print("Intentando con directorio de salida explícito...")
                subprocess.run([
                    prosextract_path,
                    "-prositedir", prosite_dir,
                    "-outdir", prosite_output_dir
                ], check=True, cwd=prosite_dir)
                print("prosextract ejecutado exitosamente (con outdir)")
            except subprocess.CalledProcessError as e3:
                raise RuntimeError(f"Falló prosextract con todos los enfoques: {e}, {e2}, {e3}")
    
    # Check if patmatmotifs exists
    patmatmotifs_path = os.path.join(emboss_bin_path, "patmatmotifs")
    if not os.path.isfile(patmatmotifs_path):
        raise FileNotFoundError(f"No se encontró 'patmatmotifs' en: {patmatmotifs_path}")
    
    # Create temporary output file
    with tempfile.NamedTemporaryFile(delete=False, suffix=".out") as temp_out:
        output_path = temp_out.name
    
    try:
        print(f"Ejecutando patmatmotifs...")
        # Force analysis of all sequences with -sformat and -sbegin/-send parameters
        result = subprocess.run([
            patmatmotifs_path,
            "-sequence", fasta_path,
            "-outfile", output_path,
            "-full", "Y",
            "-auto", "Y",
            "-sformat", "fasta",  # Explicitly specify FASTA format
            "-sbegin", "1",       # Start from first position
            "-send", "0"          # Process until end (0 means end)
        ], check=True, capture_output=True, text=True)
        print("patmatmotifs ejecutado exitosamente")
        
        # Print any warnings or messages from patmatmotifs
        if result.stderr:
            print(f"Mensajes de patmatmotifs: {result.stderr}")
        
        # Read results
        with open(output_path, "r") as result_file:
            resultado = result_file.read()
            print(f"Tamaño del archivo de resultados: {len(resultado)} caracteres")
            
    except subprocess.CalledProcessError as e:
        print(f"Error ejecutando patmatmotifs: {e}")
        if e.stderr:
            print(f"Error detallado: {e.stderr}")
        raise
    finally:
        # Clean up temporary file
        if os.path.exists(output_path):
            os.remove(output_path)
    
    return resultado




if __name__ == "__main__":
    script_dir=os.path.dirname(os.path.abspath(__file__))
    fasta_file = os.path.join(script_dir, "APOE_mrnas_refseq_variant1.fasta")
    output_fasta = os.path.join(script_dir, "APOE_translated.fasta")
    emboss_dir = os.path.join(script_dir, 'emboss', 'bin')
    PROSITE = os.path.join(script_dir, "prosite.dat")
    table = 1
    clean = True

    six_frame_translation(
        fasta_file,
        output_fasta,
        emboss_dir,  # si es None, dentro de la función caerá en DEFAULT_EMBOSS_DIR
        table=table,
        clean=clean)
    
    
    
    try:
    	# Use the enhanced analysis method that ensures all sequences are processed
    	resultado = analizar_secuencias_individuales(output_fasta, PROSITE, emboss_dir)

    
    	# Save results to a text file
    	output_file = os.path.join(script_dir, "motifs_analysis_results.txt")
    	with open(output_file, "w", encoding="utf-8") as f:
        	f.write("=== ANÁLISIS DE MOTIFS Y DOMINIOS PROSITE ===\n")
        	f.write(f"Archivo FASTA analizado: {output_fasta}\n")
        	f.write(f"Base de datos PROSITE: {PROSITE}\n")
        	f.write(f"Fecha y hora del análisis: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        
        	# Add sequence information
        	if os.path.exists(output_fasta):
            		num_seq, info_seq = contar_secuencias_fasta(output_fasta)
            		f.write(f"Número de secuencias analizadas: {num_seq}\n")
            		for i, (header, length) in enumerate(info_seq, 1):
                		f.write(f"  Secuencia {i}: {header} (longitud: {length} aa)\n")
        
        	f.write("=" * 70 + "\n\n")
        	f.write(resultado)
    
    	print(f"\nResultados guardados en: {output_file}")
    except Exception as e:
    	print(f"Error en el análisis: {e}")
    
