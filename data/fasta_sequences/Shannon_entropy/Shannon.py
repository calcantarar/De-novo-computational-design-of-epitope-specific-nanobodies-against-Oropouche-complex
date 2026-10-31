import math
import csv
from Bio import AlignIO
from collections import Counter

def shannon_entropy(column):
    """Calcula la entropía de Shannon para una columna del alineamiento."""
    # Filtramos los gaps (-) para no inflar la entropía por falta de datos
    residues = [res for res in column if res != "-"]
    if not residues:
        return 0.0
    
    n = len(residues)
    counts = Counter(residues)
    entropy = 0.0
    
    for res in counts:
        p = counts[res] / n
        entropy -= p * math.log2(p)
    
    return entropy

def export_entropy_to_csv(fasta_msa, output_csv):
    try:
        # Cargamos el MSA (formato FASTA alineado)
        alignment = AlignIO.read(fasta_msa, "fasta")
        alignment_len = alignment.get_alignment_length()
        
        with open(output_csv, mode='w', newline='') as file:
            writer = csv.writer(file)
            # Cabecera del CSV
            writer.writerow(['Posicion', 'Entropia', 'Residuo_Consenso', 'Conservacion_Porcentaje'])
            
            for i in range(alignment_len):
                col = alignment[:, i]
                entropy = shannon_entropy(col)
                
                # Cálculo de consenso y porcentaje
                counts = Counter(col)
                most_common_res, count = counts.most_common(1)[0]
                conservation_pct = (count / len(col)) * 100
                
                writer.writerow([i + 1, round(entropy, 4), most_common_res, round(conservation_pct, 2)])
        
        print(f"Análisis completado. Datos guardados en: {output_csv}")

    except Exception as e:
        print(f"Error al procesar el archivo: {e}")

# Instrucciones de uso:
# 1. Alinea tus secuencias (puedes usar MAFFT o ClustalOmega online).
# 2. Guarda el resultado como 'msa_oropouche.fasta'.
# 3. Corre esta función:
export_entropy_to_csv("msa_oropouche.fasta", "resultados_shannon.csv")
