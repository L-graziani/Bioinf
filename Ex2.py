# Importación de librerías necesarias
import os                      # Para operaciones del sistema de archivos
import subprocess              # Para ejecutar comandos externos (BLASTP en este caso)
import tempfile                # Para crear archivos temporales
import time                    # Para manejar esperas entre reintentos
from pathlib import Path       # Para manejar rutas de forma más elegante
from Bio import SeqIO          # Para leer y escribir archivos FASTA (de Biopython)

def blastp_per_sequence(fasta_path,
                        blastp_exe, 
                        db_path,
                        output_dir,
                        evalue=1e-5):
    """
    Realiza BLASTP para cada secuencia individual del archivo FASTA.
    Crea archivos temporales por cada secuencia y los elimina después.
    Devuelve una lista de rutas de archivos XML con los resultados.
    """
    
    # Si no se especifica output_dir, usar el directorio actual
    if output_dir is None:
        output_dir = os.getcwd()
    os.makedirs(output_dir, exist_ok=True)  # Crear el directorio si no existe

    xml_files = []  # Lista para guardar rutas a archivos de salida (XML)
    
    # Iterar por cada secuencia del archivo FASTA
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_id = record.id  # ID de la secuencia
        
        # Crear archivo temporal para escribir la secuencia individual
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_fasta:
            SeqIO.write(record, tmp_fasta, "fasta")  # Guardar secuencia en el archivo temporal
            tmp_fasta_path = tmp_fasta.name          # Guardar ruta del archivo temporal

        # Definir nombre del archivo de salida (formato XML, outfmt 5)
        xml_out = Path(output_dir) / f"{seq_id}.xml"
        
        # Comando para ejecutar BLASTP
        cmd = [
            blastp_exe,
            '-query', tmp_fasta_path,
            '-db', db_path,
            '-out', str(xml_out),
            '-evalue', str(evalue),
            '-outfmt', '5',  # Formato XML
        ]

        try:
            # Ejecutar el comando BLASTP, cerrando file descriptors
            subprocess.run(cmd, check=True, close_fds=True)
            xml_files.append(str(xml_out))  # Guardar salida si el comando fue exitoso
        except subprocess.CalledProcessError as e:
            print(f"Error ejecutando BLASTP para {seq_id}: {e}")
        finally:
            # Eliminar el archivo temporal
            os.remove(tmp_fasta_path)

    return xml_files  # Devolver lista de archivos XML generados

def borrar_xml_sin_hits(directorio, retries=5, delay=0.2):
    """
    Borra los archivos XML que contengan la cadena "No hits found".
    Reintenta en caso de errores de permisos (PermissionError).
    """
    for archivo in os.listdir(directorio):
        if not archivo.endswith(".xml"):  # Ignorar archivos que no sean XML
            continue
        
        ruta = os.path.join(directorio, archivo)

        try:
            # Abrir archivo en modo lectura
            with open(ruta, 'r', encoding='utf-8') as f:
                # Buscar la cadena "No hits found" en alguna línea
                if not any("No hits found" in line for line in f):
                    continue  # Si no está, no se borra el archivo
        except Exception as e:
            print(f"Error al abrir {archivo}: {e}")
            continue

        # Intentar borrar el archivo con reintentos si es necesario
        for intento in range(1, retries+1):
            try:
                os.remove(ruta)
                print(f"Archivo eliminado: {archivo}")
                break
            except PermissionError:
                if intento < retries:
                    time.sleep(delay)  # Esperar antes de reintentar
                else:
                    print(f"No se pudo eliminar {archivo} tras {retries} intentos")
            except Exception as e:
                print(f"Error borrando {archivo}: {e}")
                break

# ==== Parámetros de entrada (ajusta según tu entorno) ====
fasta_file = r"C:\Users\grazi\projects\helloworld\A_APOE_mrnas_refseq_variant1.fasta"
output_directory = r"C:\Users\grazi\projects\helloworld"
blastp_path = r"C:\Program Files\NCBI\blast-2.16.0+\bin\blastp.exe"
database_path = r"C:\PRUEBABLAST\swissprot_local"

# ==== Ejecutar BLASTP por secuencia ====
xmls = blastp_per_sequence(
    fasta_path=fasta_file,
    blastp_exe=blastp_path,
    db_path=database_path,
    output_dir=output_directory,
    evalue=1e-5
)

# Imprimir las rutas de los archivos XML generados
print("Archivos XML generados:")
for x in xmls:
    print(x)

# Eliminar resultados sin coincidencias
borrar_xml_sin_hits(output_directory)
