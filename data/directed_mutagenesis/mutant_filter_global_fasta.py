#!/usr/bin/env python3
"""
🧬 SCRIPT AUTOMÁTICO: Evaluador Multi-Cluster de Candidatos VHH (OROV-N)
Versión 2.0: Corrección de detección de epítopos y Modo Rescate
"""

import os, glob, re
import pandas as pd
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import warnings
warnings.filterwarnings('ignore')

# 🎯 UMBRALES ESTRICTOS (Para Epítopo C)
MIN_PLDDT = 0.75   
MIN_PTM   = 0.70   
MAX_PAE   = 10.0   
MAX_RMSD  = 2.5    

# 🎯 UMBRALES RELAJADOS (Para Epítopo B)
RELAXED_THRESHOLDS = {
    'MIN_PLDDT': 0.65,  # Relajado aún más
    'MIN_PTM': 0.55,    # Relajado aún más
    'MAX_PAE': 15.0,    # Relajado aún más
    'MAX_RMSD': 5.0     # Relajado aún más
}

# 📍 RUTA BASE EN WSL
BASE_DIR = "/mnt/c/Users/craff/OneDrive/Desktop/inmunosero/MUTAGENESIS"

def parse_fasta_header(header):
    header = header.lstrip('>').strip()
    metrics = {}
    design_match = re.search(r'design:(\d+)', header)
    if design_match: metrics['design'] = int(design_match.group(1))
    n_match = re.search(r'n:(\d+)', header)
    if n_match: metrics['n'] = int(n_match.group(1))
    
    parts = header.split('|')
    for part in parts:
        if ':' in part:
            key, value = part.split(':', 1)
            key = key.strip().lower()
            try: metrics[key] = float(value)
            except ValueError: metrics[key] = value
    return metrics

def parse_fasta_file(fasta_path):
    records = []
    with open(fasta_path, 'r') as f:
        current_header = None
        current_seq = []
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                if current_header is not None:
                    metrics = parse_fasta_header(current_header)
                    metrics['seq'] = ''.join(current_seq)
                    records.append(metrics)
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            metrics = parse_fasta_header(current_header)
            metrics['seq'] = ''.join(current_seq)
            records.append(metrics)
    return pd.DataFrame(records) if records else pd.DataFrame()

def compute_vhh_metrics(row):
    raw_seq = str(row.get('seq', ''))
    vhh_seq = raw_seq.split('/')[-1] if '/' in raw_seq else raw_seq
    vhh_seq = re.sub(r'[^A-Z]', '', vhh_seq.upper())
    if not vhh_seq:
        return pd.Series({'vhh_seq_clean': '', 'cys_count': 0, 'cdr3_seq': '', 'cdr3_length': 0, 'cdr3_hydrophobicity': 0.0, 'pI': 0.0, 'mw_kDa': 0.0, 'composite_score': 0.0})
    
    cys_count = vhh_seq.count('C')
    analisis = ProteinAnalysis(vhh_seq)
    pi = round(analisis.isoelectric_point(), 2)
    mw = round(analisis.molecular_weight() / 1000, 2)
    cdr3_match = re.search(r'C(.{8,18})[WF]', vhh_seq)
    cdr3 = cdr3_match.group(1) if cdr3_match else ""
    kd = {'I':4.5,'V':4.2,'L':3.8,'F':2.8,'C':2.5,'M':1.9,'A':1.8,'G':-0.4,'T':-0.7,'S':-0.8,'W':-0.9,'Y':-1.3,'P':-1.6,'H':-3.2,'E':-3.5,'Q':-3.5,'D':-3.5,'N':-3.5,'K':-3.9,'R':-4.5}
    cdr3_hydro = np.mean([kd.get(aa, 0) for aa in cdr3]) if cdr3 else 0.0
    
    plddt_norm = row.get('plddt', 0) / 100 if row.get('plddt', 0) > 1 else row.get('plddt', 0)
    score = round(row.get('ptm', 0) * 0.5 + plddt_norm * 0.5, 3)
    
    return pd.Series({'vhh_seq_clean': vhh_seq, 'cys_count': cys_count, 'cdr3_seq': cdr3, 'cdr3_length': len(cdr3), 'cdr3_hydrophobicity': round(cdr3_hydro, 2), 'pI': pi, 'mw_kDa': mw, 'composite_score': score})

def apply_filters(df, thresholds):
    plddt_cutoff = thresholds['MIN_PLDDT'] * 100 if df['plddt'].max() > 1.0 else thresholds['MIN_PLDDT']
    mask = (
        (df['plddt'] >= plddt_cutoff) &
        (df['ptm'] >= thresholds['MIN_PTM']) &
        (df['pae'] <= thresholds['MAX_PAE']) &
        (df['rmsd'] <= thresholds['MAX_RMSD'])
    )
    return df[mask].copy()

# 🔍 1. BÚSQUEDA AUTOMÁTICA
print("🔎 Buscando archivos FASTA de diseños...")
patterns = [
    os.path.join(BASE_DIR, "**", "design.fasta"),
    os.path.join(BASE_DIR, "**", "*.fasta"),
]
fasta_files = []
for pattern in patterns:
    fasta_files.extend(glob.glob(pattern, recursive=True))
fasta_files = list(dict.fromkeys(fasta_files))

if not fasta_files:
    raise FileNotFoundError(f"❌ No se encontraron archivos FASTA en:\n{BASE_DIR}")

print(f"✅ Se encontraron {len(fasta_files)} archivos FASTA.\n")

all_passed_designs = []
diagnostics_summary = []

# 🧪 2. PROCESAMIENTO
for fasta_path in fasta_files:
    path_normalized = fasta_path.replace("\\", "/")
    path_parts = path_normalized.split("/")
    exp_name = next((part for part in path_parts if ".result" in part), path_parts[-2])
    
    # ✅ CORRECCIÓN CLAVE: Detectar epítopo por el sufijo _B o _C antes de .result
    is_epitope_b = bool(re.search(r'_b\d*(?=\.result|$)', exp_name.lower()))
    
    if is_epitope_b:
        thresholds = RELAXED_THRESHOLDS
        threshold_label = "RELAJADOS (Epítopo B)"
    else:
        thresholds = {'MIN_PLDDT': MIN_PLDDT, 'MIN_PTM': MIN_PTM, 'MAX_PAE': MAX_PAE, 'MAX_RMSD': MAX_RMSD}
        threshold_label = "ESTRICTOS (Epítopo C)"
    
    df_raw = parse_fasta_file(fasta_path)
    if df_raw.empty: continue
    
    total_inicial = len(df_raw)
    df_passed = apply_filters(df_raw, thresholds)
    pass_count = len(df_passed)
    
    print(f"📂 Experimento: {exp_name}")
    print(f"   ├─ Epítopo detectado: {'B' if is_epitope_b else 'C'} | Umbrales: {threshold_label}")
    print(f"   ├─ Diseños evaluados: {total_inicial}")
    print(f"   └─ Diseños válidos: {pass_count}/{total_inicial}")
    
    # 💡 MODO RESCATE: Si es Epítopo B y da 0, mostrar los "menos malos"
    if pass_count == 0 and is_epitope_b and total_inicial > 0:
        plddt_cutoff = thresholds['MIN_PLDDT'] * 100 if df_raw['plddt'].max() > 1.0 else thresholds['MIN_PLDDT']
        print(f"   ⚠️  DIAGNÓSTICO: Ninguno pasó. Mejores valores en este lote:")
        print(f"      pLDDT máx: {df_raw['plddt'].max():.3f} (corte: {plddt_cutoff})")
        print(f"      pTM máx:   {df_raw['ptm'].max():.3f} (corte: {thresholds['MIN_PTM']})")
        print(f"      PAE mín:   {df_raw['pae'].min():.2f} (corte: {thresholds['MAX_PAE']})")
        print(f"      RMSD mín:  {df_raw['rmsd'].min():.2f} (corte: {thresholds['MAX_RMSD']})")
        
        # Calcular score temporal para rankear los mejores de los peores
        plddt_norm = df_raw['plddt'] / 100 if df_raw['plddt'].max() > 1.0 else df_raw['plddt']
        df_raw['temp_score'] = df_raw['ptm'] * 0.4 + plddt_norm * 0.4 - (df_raw['rmsd'] / 10) * 0.2
        top3 = df_raw.sort_values('temp_score', ascending=False).head(3)
        print(f"   💡 TOP 3 'Menos malos' (Candidatos de rescate para MD):")
        for _, row in top3.iterrows():
            print(f"      - Design {row.get('design', '?')} | pTM:{row['ptm']:.2f} | pLDDT:{row['plddt']:.2f} | PAE:{row['pae']:.1f} | RMSD:{row['rmsd']:.1f}")
            
    diagnostics_summary.append({
        'experimento': exp_name, 'tipo_epitopo': 'B' if is_epitope_b else 'C',
        'umbrales': threshold_label, 'total_evaluados': total_inicial,
        'validos': pass_count, 'ruta_fasta': fasta_path
    })
    
    if pass_count > 0:
        df_passed['experimento'] = exp_name
        df_passed['tipo_epitopo'] = 'B' if is_epitope_b else 'C'
        metrics = df_passed.apply(compute_vhh_metrics, axis=1)
        df_passed = pd.concat([df_passed, metrics], axis=1)
        all_passed_designs.append(df_passed)

print("\n" + "="*60)

# 📊 3. EXPORTAR
df_diag = pd.DataFrame(diagnostics_summary)
df_diag.to_csv("global_diagnostic.csv", index=False)

if all_passed_designs:
    df_global = pd.concat(all_passed_designs, ignore_index=True)
    df_global = df_global.sort_values('composite_score', ascending=False).reset_index(drop=True)
    df_global['rank_global'] = df_global.index + 1
    cols_export = ['rank_global', 'experimento', 'tipo_epitopo', 'design', 'n', 'composite_score', 'plddt', 'ptm', 'pae', 'rmsd', 'pI', 'mw_kDa', 'cdr3_seq', 'vhh_seq_clean']
    cols_finales = [c for c in cols_export if c in df_global.columns]
    df_global[cols_finales].to_excel('fasta_candidates.xlsx', index=False)
    df_global[cols_finales].to_csv('fasta_candidates.csv', index=False)
    
    print(f"\n🏆 ¡SE ENCONTRARON {len(df_global)} CANDIDATOS VÁLIDOS EN TOTAL!")
    print(f"   • Epítopo B: {len(df_global[df_global['tipo_epitopo']=='B'])}")
    print(f"   • Epítopo C: {len(df_global[df_global['tipo_epitopo']=='C'])}")
else:
    print("\n❌ NINGÚN candidato superó los filtros.")
