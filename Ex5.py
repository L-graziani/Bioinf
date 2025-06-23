#!/usr/bin/env python3
"""
Diseñador de Primers para qPCR
Script para diseñar primers específicos a partir de secuencias de ARNm en formato FASTA
"""

import json
import os
import math
from typing import List, Dict, Tuple
from dataclasses import dataclass
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

def modify_fasta_position(input_fasta: str,
                          output_fasta: str,
                          position: int = 538,
                          new_base: str = "C") -> None:
    """
    Lee un archivo FASTA, modifica el nucleótido/aminoácido en la posición indicada
    (0-based) y escribe el resultado en un archivo nuevo.

    Parámetros:
    -----------
    input_fasta : str
        Ruta al archivo FASTA de entrada.
    output_fasta : str
        Ruta al archivo FASTA de salida donde se guardará la secuencia modificada.
    position : int, opcional
        Índice 0-based de la posición a modificar (por defecto 538).
    new_base : str, opcional
        Nuevo carácter a colocar en esa posición (por defecto "C").
    """
    modified_records = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_str = str(record.seq)
        if position < 0 or position >= len(seq_str):
            raise IndexError(f"La posición {position} está fuera del rango de la secuencia (longitud {len(seq_str)}).")

        # Convertimos a lista para poder mutar el carácter
        seq_list = list(seq_str)
        seq_list[position] = new_base
        # Reconstruimos Seq y asignamos al record
        record.seq = Seq("".join(seq_list))
        modified_records.append(record)
    print("HOLA QUE HACE")
    # Escribimos todas las secuencias modificadas al archivo de salida
    SeqIO.write(modified_records, output_fasta, "fasta")
@dataclass
class PrimerConfig:
    """Configuración para el diseño de primers"""
    min_length: int
    max_length: int
    min_gc_content: float
    max_gc_content: float
    max_tm: float
    avoid_gc_clamp: bool
    num_primers: int


@dataclass
class Primer:
    """Clase para representar un primer"""
    sequence: str
    start_pos: int
    end_pos: int
    length: int
    gc_content: float
    tm: float
    
    def __str__(self):
        return (f"Secuencia: {self.sequence}\n"
                f"Posición: {self.start_pos}-{self.end_pos}\n"
                f"Longitud: {self.length} pb\n"
                f"Contenido GC: {self.gc_content:.1f}%\n"
                f"Tm: {self.tm:.1f}°C")


class PrimerDesigner:
    """Clase principal para el diseño de primers"""
    
    def __init__(self, config: PrimerConfig):
        self.config = config
    
    def read_fasta(self, fasta_file: str) -> str:
        """Lee el archivo FASTA y extrae la secuencia"""
        sequence = ""
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line.startswith('>'):
                    sequence += line.upper()
        
        if not sequence:
            raise ValueError("No se encontró secuencia válida en el archivo FASTA")
        
        return sequence
    
    def calculate_gc_content(self, sequence: str) -> float:
        """Calcula el contenido de GC de una secuencia"""
        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count / len(sequence)) * 100
    
    def calculate_tm(self, sequence: str) -> float:
        """
        Calcula la temperatura de melting usando la fórmula de Wallace
        Para primers cortos (< 14 nt): Tm = (A+T) × 2 + (G+C) × 4
        Para primers largos: Tm = 64.9 + 41 × (G+C-16.4)/N
        """
        length = len(sequence)
        gc_count = sequence.count('G') + sequence.count('C')
        at_count = sequence.count('A') + sequence.count('T')
        
        if length < 14:
            tm = at_count * 2 + gc_count * 4
        else:
            tm = 64.9 + 41 * (gc_count - 16.4) / length
        
        return tm
    
    def has_gc_clamp(self, sequence: str) -> bool:
        """Verifica si el primer tiene GC clamp en los extremos (últimos 5 nucleótidos)"""
        if len(sequence) < 5:
            return False
        
        # Verificar extremo 3' (últimos 5 nucleótidos)
        last_5 = sequence[-5:]
        gc_count_end = last_5.count('G') + last_5.count('C')
        
        # Verificar extremo 5' (primeros 5 nucleótidos)
        first_5 = sequence[:5]
        gc_count_start = first_5.count('G') + first_5.count('C')
        
        # Evitar más de 3 GC en cualquier extremo
        return gc_count_end > 3 or gc_count_start > 3
    
    def has_secondary_structures(self, sequence: str) -> bool:
        """Verifica estructuras secundarias simples (repeticiones y palíndromos)"""
        # Verificar repeticiones de 4 o más nucleótidos
        for i in range(len(sequence) - 3):
            for j in range(i + 4, len(sequence)):
                if sequence[i:i+4] == sequence[j:j+4]:
                    return True
        
        # Verificar palíndromos (complemento reverso)
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        rev_comp = ''.join([complement.get(base, base) for base in reversed(sequence)])
        
        # Buscar complementariedad de 4 o más bases
        for i in range(len(sequence) - 3):
            if sequence[i:i+4] in rev_comp:
                return True
        
        return False
    
    def validate_primer(self, sequence: str) -> bool:
        """Valida si un primer cumple con todos los criterios"""
        # Verificar longitud
        if not (self.config.min_length <= len(sequence) <= self.config.max_length):
            return False
        
        # Verificar contenido GC
        gc_content = self.calculate_gc_content(sequence)
        if not (self.config.min_gc_content <= gc_content <= self.config.max_gc_content):
            return False
        
        # Verificar temperatura de melting
        tm = self.calculate_tm(sequence)
        if tm > self.config.max_tm:
            return False
        
        # Verificar GC clamp si está habilitado
        if self.config.avoid_gc_clamp and self.has_gc_clamp(sequence):
            return False
        
        # Verificar estructuras secundarias
        if self.has_secondary_structures(sequence):
            return False
        
        return True
    
    def generate_primers(self, sequence: str) -> List[Primer]:
        """Genera primers válidos a partir de la secuencia"""
        primers = []
        
        # Generar primers de diferentes longitudes y posiciones
        for length in range(self.config.min_length, self.config.max_length + 1):
            for start in range(0, len(sequence) - length + 1, 1):  # Paso de 5 para eficiencia
                primer_seq = sequence[start:start + length]
                
                if self.validate_primer(primer_seq):
                    primer = Primer(
                        sequence=primer_seq,
                        start_pos=start + 1,  # Posición 1-indexed
                        end_pos=start + length,
                        length=length,
                        gc_content=self.calculate_gc_content(primer_seq),
                        tm=self.calculate_tm(primer_seq)
                    )
                    primers.append(primer)
        
        return primers
    
    def select_best_primers(self, primers: List[Primer]) -> List[Primer]:
        """Selecciona los mejores primers basado en criterios de calidad"""
        if len(primers) <= self.config.num_primers:
            return primers
        
        # Puntuar primers basado en múltiples criterios
        scored_primers = []
        
        for primer in primers:
            score = 0
            
            # Favorecer Tm cercano a 60°C
            tm_score = 10 - abs(primer.tm - 60)
            score += max(0, tm_score)
            
            # Favorecer GC content cercano al 55%
            gc_score = 10 - abs(primer.gc_content - 55)
            score += max(0, gc_score)
            
            # Favorecer longitud intermedia (20-22 pb)
            length_score = 10 - abs(primer.length - 21)
            score += max(0, length_score)
            
            scored_primers.append((score, primer))
        
        # Ordenar por puntuación y tomar los mejores
        scored_primers.sort(key=lambda x: x[0], reverse=True)
        
        # Seleccionar primers no solapantes
        selected_primers = []
        for score, primer in scored_primers:
            # Verificar que no solape significativamente con primers ya seleccionados
            overlap = False
            for selected in selected_primers:
                if (abs(primer.start_pos - selected.start_pos) < 10):  # Mínimo 10 pb de separación
                    overlap = True
                    break
            
            if not overlap:
                selected_primers.append(primer)
                if len(selected_primers) >= self.config.num_primers:
                    break
        
        return selected_primers


def load_config(config_file: str) -> PrimerConfig:
    """Carga la configuración desde un archivo JSON"""
    try:
        with open(config_file, 'r') as f:
            config_data = json.load(f)
        
        return PrimerConfig(
            min_length=config_data.get('min_length', 18),
            max_length=config_data.get('max_length', 24),
            min_gc_content=config_data.get('min_gc_content', 50.0),
            max_gc_content=config_data.get('max_gc_content', 60.0),
            max_tm=config_data.get('max_tm', 67.0),
            avoid_gc_clamp=config_data.get('avoid_gc_clamp', True),
            num_primers=config_data.get('num_primers', 5)
        )
    except FileNotFoundError:
        print(f"No se encontró el archivo de configuración {config_file}")
        print("Usando configuración por defecto...")
        return PrimerConfig(
            min_length=18,
            max_length=24,
            min_gc_content=50.0,
            max_gc_content=60.0,
            max_tm=67.0,
            avoid_gc_clamp=True,
            num_primers=5
        )


def create_default_config(filename: str = "config.json"):
    """Crea un archivo de configuración por defecto"""
    default_config = {
        "min_length": 18,
        "max_length": 24,
        "min_gc_content": 50.0,
        "max_gc_content": 60.0,
        "max_tm": 67.0,
        "avoid_gc_clamp": True,
        "num_primers": 5
    }
    
    with open(filename, 'w') as f:
        json.dump(default_config, f, indent=4)
    
    print(f"Archivo de configuración creado: {filename}")


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    fasta_path = os.path.join(script_dir, 'APOE_mrnas_refseq_variant1.fasta')
    config_path = os.path.join(script_dir, 'config.json')
    output_path = os.path.join(script_dir, 'resultados.json')
    modify_fasta_position(
        fasta_path,
        fasta_path,
        position=538,
        new_base="C"
    )
    # Cambia esto a True si quieres generar un archivo de configuración por defecto
    create_config = False

    if create_config:
        create_default_config(config_path)
        print(f"Archivo de configuración creado en: {config_path}")
        return

    try:
        # Cargar configuración
        config = load_config(config_path)

        # Crear diseñador de primers
        designer = PrimerDesigner(config)

        # Leer secuencia FASTA
        print(f"Leyendo secuencia desde {fasta_path}...")
        sequence = designer.read_fasta(fasta_path)
        print(f"Secuencia cargada: {len(sequence)} nucleótidos")

        # Generar primers
        print("Generando primers candidatos...")
        all_primers = designer.generate_primers(sequence)
        print(f"Se encontraron {len(all_primers)} primers candidatos")

        if not all_primers:
            print("No se encontraron primers que cumplan con los criterios especificados.")
            print("Considere ajustar los parámetros en el archivo de configuración.")
            return

        # Seleccionar mejores primers
        best_primers = designer.select_best_primers(all_primers)

        print(f"\n=== MEJORES {len(best_primers)} PRIMERS DISEÑADOS ===\n")

        # Mostrar resultados
        results = []
        for i, primer in enumerate(best_primers, 1):
            print(f"PRIMER {i}:")
            print(primer)
            print("-" * 40)

            results.append({
                'id': f'Primer_{i}',
                'sequence': primer.sequence,
                'position': f'{primer.start_pos}-{primer.end_pos}',
                'length': primer.length,
                'gc_content': round(primer.gc_content, 1),
                'tm': round(primer.tm, 1)
            })

        # Guardar resultados
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=4)
        print(f"\nResultados guardados en: {output_path}")

    except Exception as e:
        print(f"Error: {str(e)}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
