#!/usr/bin/env python3
"""
Módulo para traducción de secuencias de nucleótidos en seis marcos de lectura
y análisis de motifs y dominios usando herramientas EMBOSS.
"""

import sys
import os
import tempfile
import subprocess
import datetime
from pathlib import Path
from typing import List, Tuple, Optional
from Bio import SeqIO


class SixFrameTranslator:
    """Clase para manejar la traducción en seis marcos de lectura."""
    
    def __init__(self, emboss_dir: Optional[str] = None):
        """
        Inicializa el traductor.
        
        Args:
            emboss_dir: Directorio donde está instalado EMBOSS (si no está en PATH)
        """
        self.emboss_dir = emboss_dir
        
    def translate(self, input_fasta: str, output_fasta: str, 
                 table: int = 1, clean: bool = True) -> int:
        """
        Traduce la secuencia de nucleótidos en las seis fases de lectura usando EMBOSS Transeq.

        Args:
            input_fasta: Ruta al fasta de entrada (nucleótidos)
            output_fasta: Ruta al fasta de salida (proteínas)
            table: Código de tabla de traducción (default 1: estándar)
            clean: Si True, elimina archivos temporales

        Returns:
            Número de secuencias proteicas generadas
            
        Raises:
            FileNotFoundError: Si no se encuentra el archivo de entrada o transeq
            subprocess.CalledProcessError: Si falla la ejecución de transeq
        """
        if not os.path.isfile(input_fasta):
            raise FileNotFoundError(f"Archivo FASTA de entrada no encontrado: {input_fasta}")

        # Crear archivo temporal para la salida de Transeq
        tmp_fd, tmp_path = tempfile.mkstemp(suffix=".faa")
        os.close(tmp_fd)

        try:
            # Determinar ruta al ejecutable transeq
            transeq_exec = self._get_transeq_path()

            # Construir y ejecutar el comando de Transeq
            cmd = [
                transeq_exec,
                '-sequence', input_fasta,
                '-outseq', tmp_path,
                '-table', str(table),
                '-frame', '6'
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Leer el resultado y escribir en el fasta de salida
            proteins = list(SeqIO.parse(tmp_path, "fasta"))
            if not proteins:
                raise ValueError("No se encontraron secuencias traducidas en la salida de Transeq")

            SeqIO.write(proteins, output_fasta, "fasta")
            print(f"Traducción completada. Se han generado {len(proteins)} secuencias proteicas.")
            return len(proteins)

        except subprocess.CalledProcessError as e:
            error_msg = f"Error al ejecutar transeq: {e.stderr if e.stderr else str(e)}"
            if not clean:
                error_msg += f"\nArchivo temporal conservado en: {tmp_path}"
            raise subprocess.CalledProcessError(e.returncode, e.cmd, error_msg)
        
        finally:
            # Limpiar temporal si corresponde
            if clean and os.path.exists(tmp_path):
                os.remove(tmp_path)

    def _get_transeq_path(self) -> str:
        """Obtiene la ruta al ejecutable transeq."""
        if self.emboss_dir:
            transeq_path = os.path.join(self.emboss_dir, 'transeq')
            if not os.path.isfile(transeq_path):
                raise FileNotFoundError(f"No se encontró transeq en: {transeq_path}")
            return transeq_path
        return 'transeq'  # Asume que está en PATH


class SequenceAnalyzer:
    """Clase para analizar motifs y dominios en secuencias proteicas."""
    
    def __init__(self, emboss_bin_path: str):
        """
        Inicializa el analizador.
        
        Args:
            emboss_bin_path: Ruta al directorio bin de EMBOSS
        """
        self.emboss_bin_path = emboss_bin_path
        self._validate_emboss_tools()
    
    def _validate_emboss_tools(self):
        """Valida que las herramientas EMBOSS necesarias estén disponibles."""
        tools = ['prosextract', 'patmatmotifs']
        for tool in tools:
            tool_path = os.path.join(self.emboss_bin_path, tool)
            if not os.path.isfile(tool_path):
                raise FileNotFoundError(f"No se encontró '{tool}' en: {tool_path}")

    def count_sequences(self, fasta_path: str) -> Tuple[int, List[Tuple[str, int]]]:
        """
        Cuenta el número de secuencias en un archivo FASTA.
        
        Args:
            fasta_path: Ruta al archivo FASTA
            
        Returns:
            Tupla con (número_de_secuencias, lista_de_info_secuencias)
        """
        sequences_info = []
        
        try:
            with open(fasta_path, 'r') as f:
                current_header = None
                current_seq_length = 0
                
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_header:
                            sequences_info.append((current_header, current_seq_length))
                        current_header = line[1:]  # Remove '>'
                        current_seq_length = 0
                    elif line and not line.startswith('>'):
                        current_seq_length += len(line)
                
                # Agregar la última secuencia
                if current_header:
                    sequences_info.append((current_header, current_seq_length))
        
        except FileNotFoundError:
            raise FileNotFoundError(f"Archivo FASTA no encontrado: {fasta_path}")
        
        return len(sequences_info), sequences_info

    def extract_individual_sequences(self, fasta_path: str) -> List[Tuple[str, str]]:
        """
        Extrae cada secuencia del FASTA como tuplas (header, secuencia).
        
        Args:
            fasta_path: Ruta al archivo FASTA
            
        Returns:
            Lista de tuplas (header, secuencia)
        """
        sequences = []
        current_header = None
        current_seq = []
        
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_header and current_seq:
                        sequences.append((current_header, ''.join(current_seq)))
                    current_header = line
                    current_seq = []
                elif line:
                    current_seq.append(line)
            
            # Agregar la última secuencia
            if current_header and current_seq:
                sequences.append((current_header, ''.join(current_seq)))
        
        return sequences

    def _run_prosextract(self, prosite_dir: str):
        """Ejecuta prosextract para preparar la base de datos PROSITE."""
        prosextract_path = os.path.join(self.emboss_bin_path, "prosextract")
        
        # Crear directorio PROSITE si no existe
        prosite_output_dir = os.path.join(prosite_dir, "PROSITE")
        os.makedirs(prosite_output_dir, exist_ok=True)
        
        # Intentar diferentes enfoques para ejecutar prosextract
        commands_to_try = [
            [prosextract_path, "-prositedir", prosite_dir],
            [prosextract_path],
            [prosextract_path, "-prositedir", prosite_dir, "-outdir", prosite_output_dir]
        ]
        
        for i, cmd in enumerate(commands_to_try, 1):
            try:
                print(f"Ejecutando prosextract (intento {i})...")
                subprocess.run(cmd, check=True, cwd=prosite_dir, 
                             capture_output=True, text=True)
                print("prosextract ejecutado exitosamente")
                return
            except subprocess.CalledProcessError as e:
                if i == len(commands_to_try):
                    raise RuntimeError(f"Falló prosextract con todos los enfoques: {e}")
                continue

    def analyze_individual_sequences(self, fasta_path: str, prosite_dat_path: str) -> str:
        """
        Analiza cada secuencia del FASTA por separado para asegurar procesamiento completo.
        
        Args:
            fasta_path: Ruta al archivo FASTA
            prosite_dat_path: Ruta al archivo prosite.dat
            
        Returns:
            Resultados combinados del análisis
        """
        print("Procesando secuencias individualmente...")
        return self._analyze_sequences_individually(fasta_path, prosite_dat_path)

    def _analyze_sequences_individually(self, fasta_path: str, prosite_dat_path: str) -> str:
        """Analiza cada secuencia individualmente."""
        sequences = self.extract_individual_sequences(fasta_path)
        combined_results = []
        
        patmatmotifs_path = os.path.join(self.emboss_bin_path, "patmatmotifs")
        
        # Preparar base de datos PROSITE
        prosite_dir = os.path.dirname(prosite_dat_path)
        os.environ["EMBOSS_DATA"] = prosite_dir
        self._run_prosextract(prosite_dir)
        
        for i, (header, seq) in enumerate(sequences, 1):
            print(f"Procesando secuencia {i}/{len(sequences)}: {header[:50]}...")
            
            # Crear archivo temporal para esta secuencia
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".fasta") as temp_fasta:
                temp_fasta.write(f"{header}\n")
                # Escribir secuencia en líneas de 80 caracteres
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
                ], check=True, capture_output=True, text=True)
                
                with open(output_path, "r") as result_file:
                    resultado_individual = result_file.read()
                    combined_results.append(
                        f"\n{'='*60}\n{header}\n{'='*60}\n{resultado_individual}"
                    )
            
            except subprocess.CalledProcessError as e:
                combined_results.append(
                    f"\n{'='*60}\n{header}\n{'='*60}\n"
                    f"Error procesando esta secuencia: {e}\n"
                )
            
            finally:
                # Limpiar archivos temporales
                for temp_file in [temp_fasta_path, output_path]:
                    if os.path.exists(temp_file):
                        os.remove(temp_file)
        
        return "\n".join(combined_results)


def save_analysis_results(results: str, output_fasta: str, prosite_path: str, 
                         script_dir: str, analyzer: SequenceAnalyzer) -> str:
    """
    Guarda los resultados del análisis en un archivo de texto.
    
    Args:
        results: Resultados del análisis
        output_fasta: Ruta al archivo FASTA analizado
        prosite_path: Ruta a la base de datos PROSITE
        script_dir: Directorio del script
        analyzer: Instancia del analizador para obtener información de secuencias
        
    Returns:
        Ruta del archivo de resultados guardado
    """
    output_file = os.path.join(script_dir, "motifs_analysis_results.txt")
    
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("=== ANÁLISIS DE MOTIFS Y DOMINIOS PROSITE ===\n")
        f.write(f"Archivo FASTA analizado: {output_fasta}\n")
        f.write(f"Base de datos PROSITE: {prosite_path}\n")
        f.write(f"Fecha y hora del análisis: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        
        # Agregar información de secuencias
        if os.path.exists(output_fasta):
            try:
                num_seq, info_seq = analyzer.count_sequences(output_fasta)
                f.write(f"Número de secuencias analizadas: {num_seq}\n")
                for i, (header, length) in enumerate(info_seq, 1):
                    f.write(f"  Secuencia {i}: {header} (longitud: {length} aa)\n")
            except Exception as e:
                f.write(f"Error al obtener información de secuencias: {e}\n")
        
        f.write("=" * 70 + "\n\n")
        f.write(results)
    
    return output_file


def main():
    """Función principal del script."""
    # Configuración de rutas
    script_dir = os.path.dirname(os.path.abspath(__file__))
    fasta_file = os.path.join(script_dir, "APOE_mrnas_refseq_variant1.fasta")
    output_fasta = os.path.join(script_dir, "APOE_translated.fasta")
    emboss_dir = os.path.join(script_dir, 'emboss', 'bin')
    prosite_path = os.path.join(script_dir, "prosite.dat")
    
    # Configuración de parámetros
    table = 1
    clean = True

    try:
        # Traducción en seis marcos de lectura
        print("=== INICIANDO TRADUCCIÓN EN SEIS MARCOS DE LECTURA ===")
        translator = SixFrameTranslator(emboss_dir)
        num_proteins = translator.translate(fasta_file, output_fasta, table, clean)
        print(f"Traducción completada: {num_proteins} secuencias proteicas generadas")

        # Análisis de motifs y dominios
        print("\n=== INICIANDO ANÁLISIS DE MOTIFS Y DOMINIOS ===")
        analyzer = SequenceAnalyzer(emboss_dir)
        results = analyzer.analyze_individual_sequences(output_fasta, prosite_path)

        # Guardar resultados
        output_file = save_analysis_results(
            results, output_fasta, prosite_path, script_dir, analyzer
        )
        print(f"\nResultados guardados en: {output_file}")
        
    except Exception as e:
        print(f"Error en la ejecución: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
