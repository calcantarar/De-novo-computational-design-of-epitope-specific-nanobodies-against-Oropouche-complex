# MD_figures.py
import os
import numpy as np
import matplotlib.pyplot as plt

OUT = '/root/figuras_poster'          # ext4: evita bloqueo WSL/OneDrive
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    'font.size': 14, 'axes.titlesize': 18, 'axes.labelsize': 15,
    'xtick.labelsize': 12, 'ytick.labelsize': 12, 'legend.fontsize': 11,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'axes.linewidth': 1.2, 'xtick.major.width': 1.2, 'ytick.major.width': 1.2,
})

C_B, C_C = '#E74C3C', '#3498DB'
CDR_COLS = {'CDR1': '#27AE60', 'CDR2': '#F39C12', 'CDR3': '#9B59B6'}  # ← CAMBIO 3

db = '/mnt/c/Users/craff/OneDrive/Desktop/inmunosero/MM-PBSA/9ezu-b-7'
dc = '/mnt/c/Users/craff/OneDrive/Desktop/inmunosero/MM-PBSA/9ezu-c-0'

def load_xvg(fn):
    d = []
    with open(fn) as f:
        for line in f:
            if line.startswith(('#', '@')): continue
            p = line.split()
            if len(p) >= 2:
                try: d.append([float(p[0]), float(p[1])])
                except: continue
    return np.array(d) if d else None

def rolling(y, w=250): return np.convolve(y, np.ones(w)/w, mode='same')

# ---------------- Carga (idéntica a tu v2) ----------------
rmsd_b = load_xvg(f"{db}/rmsd_interfaz.xvg"); rmsd_b[:,1] *= 10
rmsd_c = load_xvg(f"{dc}/rmsd_interfaz.xvg"); rmsd_c[:,1] *= 10
hb_b = load_xvg(f"{db}/hbonds_interfaz.xvg"); hb_b[:,0] /= 1000.0
hb_c = load_xvg(f"{dc}/hbonds_interfaz.xvg"); hb_c[:,0] /= 1000.0
vhh_b = load_xvg(f"{db}/rmsf_vhh.xvg");  vhh_c = load_xvg(f"{dc}/rmsf_vhh.xvg")
ag_b  = load_xvg(f"{db}/rmsf_orovn.xvg"); ag_c  = load_xvg(f"{dc}/rmsf_orovn.xvg")

for n,a in [("hb_B",hb_b),("hb_C",hb_c)]:
    assert len(a) > 4000 and 90 < a[-1,0] <= 101, f"❌ {n}: frames/tiempo inválidos"
print(f"✅ B: {len(hb_b)} frames, avg H-bonds {hb_b[:,1].mean():.2f}")
print(f"✅ C: {len(hb_c)} frames, avg H-bonds {hb_c[:,1].mean():.2f}")

# ================= FIG 7 =================
fig,(ax1,ax2) = plt.subplots(1,2, figsize=(14,5))

# --- Panel A: CAMBIO 1 (leyenda 2 entradas + promedios al borde) ---
ax1.plot(rmsd_b[:,0], rmsd_b[:,1], color=C_B, alpha=0.5, lw=0.5, label='Epitope B')
ax1.plot(rmsd_c[:,0], rmsd_c[:,1], color=C_C, alpha=0.5, lw=0.5, label='Epitope C')
mb, mc = rmsd_b[:,1].mean(), rmsd_c[:,1].mean()
ax1.axhline(mb, color=C_B, ls='--', lw=1.8)                      # sin label
ax1.axhline(mc, color=C_C, ls='--', lw=1.8)                      # sin label
ax1.annotate(f'{mb:.2f} Å', xy=(100,mb), xytext=(101.5,mb-0.5), color=C_B,
             fontsize=11, fontweight='bold', va='center', annotation_clip=False)
ax1.annotate(f'{mc:.2f} Å', xy=(100,mc), xytext=(101.5,mc+0.5), color=C_C,
             fontsize=11, fontweight='bold', va='center', annotation_clip=False)
ax1.set_xlabel('Time (ns)'); ax1.set_ylabel('Interface RMSD (Å)')
ax1.set_title('A) Structural stability', fontweight='bold')
ax1.set_xlim(0,100); ax1.grid(alpha=0.3)
ax1.legend(loc='upper left', framealpha=0.9)                     # solo 2 entradas
ax1.text(0.98,0.02,'dashed = trajectory mean', transform=ax1.transAxes,
         ha='right', va='bottom', fontsize=8.5, style='italic', color='gray')

# --- Panel B: CAMBIO 2 (leyenda 2 entradas + promedios al borde) ---
ax2.plot(hb_b[:,0], hb_b[:,1], color=C_B, alpha=0.22, lw=0.4)    # raw sin label
ax2.plot(hb_c[:,0], hb_c[:,1], color=C_C, alpha=0.22, lw=0.4)
ax2.plot(hb_b[:,0], rolling(hb_b[:,1]), color=C_B, lw=2, label='Epitope B')
ax2.plot(hb_c[:,0], rolling(hb_c[:,1]), color=C_C, lw=2, label='Epitope C')
mb2, mc2 = hb_b[:,1].mean(), hb_c[:,1].mean()
ax2.axhline(mb2, color=C_B, ls='--', lw=1.5)                     # sin label
ax2.axhline(mc2, color=C_C, ls='--', lw=1.5)                     # sin label
ax2.annotate(f'{mb2:.2f}', xy=(100,mb2), xytext=(101.5,mb2-0.7), color=C_B,
             fontsize=11, fontweight='bold', va='center', annotation_clip=False)
ax2.annotate(f'{mc2:.2f}', xy=(100,mc2), xytext=(101.5,mc2+0.7), color=C_C,
             fontsize=11, fontweight='bold', va='center', annotation_clip=False)
ax2.set_xlabel('Time (ns)'); ax2.set_ylabel('Number of H-bonds')
ax2.set_title('B) Interfacial H-bonds', fontweight='bold')
ax2.set_xlim(0,100); ax2.set_ylim(0, max(hb_b[:,1].max(), hb_c[:,1].max())+2)
ax2.grid(alpha=0.3)
ax2.legend(loc='upper right', framealpha=0.9)                    # solo 2 entradas
ax2.text(0.02,0.98,'faint: raw per-frame  •  solid: 5-ns rolling  •  dashed: mean',
         transform=ax2.transAxes, va='top', fontsize=8.5, style='italic', color='gray')

plt.tight_layout()
plt.savefig(f'{OUT}/FIG7_rmsd_hbonds.png'); plt.savefig(f'{OUT}/FIG7_rmsd_hbonds.svg')
plt.close()
print('✅ FIG7 generada (leyendas de 2 entradas)')

# ================= RMSF =================
fig, axes = plt.subplots(2,1, figsize=(12,8), sharey=False)

# --- Panel A: SIN CAMBIOS (tu versión validada) ---
axes[0].plot(ag_b[:,0], ag_b[:,1]*10, color=C_B, alpha=0.8, lw=1, label='Epitope B')
axes[0].plot(ag_c[:,0], ag_c[:,1]*10, color=C_C, alpha=0.8, lw=1, label='Epitope C')
axes[0].axvspan(75, 85, color=C_B, alpha=0.15, label='Ep B region')
axes[0].axvspan(171, 186, color=C_C, alpha=0.15, label='Ep C region')
axes[0].set_xlabel('Residue number (OROV-N)'); axes[0].set_ylabel('RMSF (Å)')
axes[0].set_title('Antigen (OROV-N) per-residue flexibility', fontweight='bold')
axes[0].legend(loc='upper right', ncol=2); axes[0].grid(alpha=0.3)

# --- Panel B: CAMBIO 3 (CDRs con colores distintos) ---
axes[1].plot(vhh_b[:,0], vhh_b[:,1]*10, color=C_B, alpha=0.8, lw=1, label='Anti-Ep B')
axes[1].plot(vhh_c[:,0], vhh_c[:,1]*10, color=C_C, alpha=0.8, lw=1, label='Anti-Ep C')
axes[1].axvspan(24, 32,  color=CDR_COLS['CDR1'], alpha=0.18, label='CDR1')
axes[1].axvspan(50, 57,  color=CDR_COLS['CDR2'], alpha=0.18, label='CDR2')
axes[1].axvspan(98, 112, color=CDR_COLS['CDR3'], alpha=0.18, label='CDR3')
axes[1].set_xlabel('Residue number (9EZU VHH)'); axes[1].set_ylabel('RMSF (Å)')
axes[1].set_title('Nanobody (9EZU scaffold) per-residue flexibility', fontweight='bold')
axes[1].legend(loc='upper right', ncol=3); axes[1].grid(alpha=0.3)
# Opcional (legibilidad a distancia en póster): descomenta para nombres sobre bandas
# y1 = axes[1].get_ylim()[1]
# for s,e,nm in [(24,32,'CDR1'),(50,57,'CDR2'),(98,112,'CDR3')]:
#     axes[1].text((s+e)/2, y1*0.94, nm, color=CDR_COLS[nm], ha='center',
#                  fontweight='bold', fontsize=11)

plt.tight_layout()
plt.savefig(f'{OUT}/FIG_rmsf.png'); plt.savefig(f'{OUT}/FIG_rmsf.svg')
plt.close()
print('✅ FIG_rmsf generada (CDR1 verde / CDR2 naranja / CDR3 morado)')
print(f'\n📁 Figuras guardadas en: {OUT}')
