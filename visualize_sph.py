#!/usr/bin/env python3
"""
visualize_sph.py  --  Visualizacion de la simulacion SPH

Uso:
    python3 visualize_sph.py                   # produce sph_animation.mp4
    python3 visualize_sph.py --frames          # produce frames PNG en ./frames/
    python3 visualize_sph.py --step 50         # muestra un unico paso en pantalla
    python3 visualize_sph.py --color density   # colorea por densidad (default: speed)
    python3 visualize_sph.py --nfluid 1600     # numero de particulas de fluido (default: 1600)
    python3 visualize_sph.py --skip 5          # usa 1 de cada 5 frames (animacion mas rapida)

Columnas del archivo de estado:
    id  x  y  vx  vy  ax  ay  rho  mass  p  c  u
    0   1  2   3   4   5   6   7    8    9  10 11
"""

import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# --------------------------------------------------------------------------
# Parseo de argumentos
# --------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Visualiza la simulacion SPH")
parser.add_argument("--frames",  action="store_true",
                    help="Guarda frames PNG en ./frames/ en lugar de un MP4")
parser.add_argument("--step",    type=int, default=None,
                    help="Muestra un unico paso (numero entero) y sale")
parser.add_argument("--color",   choices=["speed", "density", "pressure"],
                    default="speed",
                    help="Variable usada para colorear las particulas")
parser.add_argument("--nfluid",  type=int, default=1600,
                    help="Numero de particulas de fluido (las primeras N del archivo)")
parser.add_argument("--skip",    type=int, default=1,
                    help="Saltar pasos: 1 = todos, 5 = uno de cada cinco")
parser.add_argument("--outdir",  default="./output",
                    help="Directorio con los archivos state_XXXX (default: ./output)")
parser.add_argument("--quiver",  action="store_true",
                    help="Superpone vectores de velocidad sobre el fluido")
args = parser.parse_args()

# --------------------------------------------------------------------------
# Localizacion de los archivos de estado
# --------------------------------------------------------------------------
pattern = os.path.join(args.outdir, "state_*")
all_files = sorted(glob.glob(pattern))

if not all_files:
    sys.exit(f"[ERROR] No se encontraron archivos en '{args.outdir}/state_*'.\n"
             f"Ejecuta la simulacion primero: ./sph <pasos>")

# Aplica el salto si se pidio
files = all_files[::args.skip]
print(f"Archivos encontrados: {len(all_files)}  |  A procesar: {len(files)}")

# --------------------------------------------------------------------------
# Lectura de un archivo de estado
# Devuelve un dict con arrays numpy por columna
# --------------------------------------------------------------------------
def read_state(path):
    """Lee un archivo state_XXXX y devuelve un diccionario con arrays."""
    # np.loadtxt es robusto para columnas numericas separadas por espacios
    data = np.loadtxt(path)
    return {
        "id":      data[:, 0].astype(int),
        "x":       data[:, 1],   # posicion en x [m]
        "y":       data[:, 2],   # posicion en y [m]
        "vx":      data[:, 3],   # velocidad en x [m/s]
        "vy":      data[:, 4],   # velocidad en y [m/s]
        "ax":      data[:, 5],   # aceleracion en x
        "ay":      data[:, 6],   # aceleracion en y
        "rho":     data[:, 7],   # densidad [kg/m^3]
        "mass":    data[:, 8],   # masa [kg]
        "p":       data[:, 9],   # presion [Pa]
        "c":       data[:, 10],  # velocidad del sonido
        "u":       data[:, 11],  # energia interna
    }

# --------------------------------------------------------------------------
# Seleccion de la variable de color y su etiqueta
# --------------------------------------------------------------------------
def get_color_data(state, color_var):
    """Devuelve el array de valores y la etiqueta para la barra de color."""
    if color_var == "speed":
        # modulo de la velocidad: sqrt(vx^2 + vy^2)
        values = np.sqrt(state["vx"]**2 + state["vy"]**2)
        label  = "Velocidad [m/s]"
    elif color_var == "density":
        values = state["rho"]
        label  = "Densidad [kg/m³]"
    elif color_var == "pressure":
        values = state["p"]
        label  = "Presion [Pa]"
    return values, label

# --------------------------------------------------------------------------
# Calculo del rango de color a partir del PRIMER y ULTIMO estado
# para que la barra no salte entre frames
# --------------------------------------------------------------------------
print("Calculando rango de color (primer y ultimo estado)...")
first_state = read_state(files[0])
last_state  = read_state(files[-1])

c0, _ = get_color_data(first_state, args.color)
cN, _ = get_color_data(last_state,  args.color)

# Se usa solo el fluido para fijar la escala (las fronteras distorsionan)
fluid_mask_first = first_state["id"] < args.nfluid
fluid_mask_last  = last_state["id"]  < args.nfluid

vmin = min(c0[fluid_mask_first].min(), cN[fluid_mask_last].min())
vmax = max(c0[fluid_mask_first].max(), cN[fluid_mask_last].max())

# Evita que vmin == vmax (puede ocurrir en el estado inicial con vel=0)
if np.isclose(vmin, vmax):
    vmax = vmin + 1e-12

print(f"Rango de color ({args.color}): [{vmin:.4e}, {vmax:.4e}]")

# --------------------------------------------------------------------------
# Configuracion de la figura
# --------------------------------------------------------------------------
Lx = 1e-3   # dominio en x [m]  (debe coincidir con el valor en el codigo C)
Ly = 1e-3   # dominio en y [m]

fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlim(-0.05*Lx*1e3, 1.05*Lx*1e3)   # en mm
ax.set_ylim(-0.05*Ly*1e3, 1.05*Ly*1e3)
ax.set_xlabel("x [mm]")
ax.set_ylabel("y [mm]")
ax.set_aspect("equal")
ax.set_facecolor("#111111")
fig.patch.set_facecolor("#1a1a1a")
ax.tick_params(colors="white")
ax.xaxis.label.set_color("white")
ax.yaxis.label.set_color("white")
for spine in ax.spines.values():
    spine.set_edgecolor("#555555")

# Barra de color
norm    = Normalize(vmin=vmin, vmax=vmax)
sm      = ScalarMappable(cmap="plasma", norm=norm)
sm.set_array([])
cbar    = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
cbar.ax.tick_params(colors="white")
_, clabel = get_color_data(first_state, args.color)
cbar.set_label(clabel, color="white")

# Scatter inicial: se divide en fluido y frontera para usar marcadores distintos
s0 = read_state(files[0])
fluid_mask = s0["id"] < args.nfluid
bound_mask = ~fluid_mask

# Fronteras: puntos pequenos grises
sc_bound = ax.scatter(
    s0["x"][bound_mask] * 1e3,   # conversion m -> mm
    s0["y"][bound_mask] * 1e3,
    s=4, c="#888888", zorder=2, label="Frontera"
)

# Fluido: puntos medianos coloreados
c_fluid, _ = get_color_data(s0, args.color)
sc_fluid = ax.scatter(
    s0["x"][fluid_mask] * 1e3,
    s0["y"][fluid_mask] * 1e3,
    s=8, c=c_fluid[fluid_mask], cmap="plasma",
    norm=norm, zorder=3, label="Fluido"
)

# Quiver opcional: vectores de velocidad del fluido (submuestreado para claridad)
if args.quiver:
    # Toma 1 de cada 8 particulas de fluido para no saturar el grafico
    qstep = 8
    xi  = s0["x"][fluid_mask][::qstep] * 1e3
    yi  = s0["y"][fluid_mask][::qstep] * 1e3
    vxi = s0["vx"][fluid_mask][::qstep]
    vyi = s0["vy"][fluid_mask][::qstep]
    qv  = ax.quiver(xi, yi, vxi, vyi,
                    color="white", alpha=0.6,
                    scale=None, scale_units="xy", angles="xy",
                    width=0.003)

# Titulo
title = ax.set_title(f"SPH  paso 0/{len(files)-1}",
                     color="white", fontsize=11)

ax.legend(loc="upper right", fontsize=7,
          facecolor="#333333", labelcolor="white",
          markerscale=2)

# --------------------------------------------------------------------------
# Modo --step: mostrar un solo paso
# --------------------------------------------------------------------------
if args.step is not None:
    step_path = os.path.join(args.outdir, f"state_{args.step:04d}")
    if not os.path.exists(step_path):
        sys.exit(f"[ERROR] No existe '{step_path}'")
    state   = read_state(step_path)
    fmask   = state["id"] < args.nfluid
    bmask   = ~fmask
    cvals,_ = get_color_data(state, args.color)

    ax.cla()
    ax.set_xlim(-0.05*Lx*1e3, 1.05*Lx*1e3)
    ax.set_ylim(-0.05*Ly*1e3, 1.05*Ly*1e3)
    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    ax.set_aspect("equal"); ax.set_facecolor("#111111")
    ax.scatter(state["x"][bmask]*1e3, state["y"][bmask]*1e3,
               s=4, c="#888888", zorder=2)
    ax.scatter(state["x"][fmask]*1e3, state["y"][fmask]*1e3,
               s=8, c=cvals[fmask], cmap="plasma", norm=norm, zorder=3)
    ax.set_title(f"SPH  paso {args.step}", color="white")
    plt.tight_layout()
    plt.show()
    sys.exit(0)

# --------------------------------------------------------------------------
# Modo --frames: guardar PNG individuales
# --------------------------------------------------------------------------
if args.frames:
    frames_dir = "./frames"
    os.makedirs(frames_dir, exist_ok=True)
    print(f"Guardando frames en '{frames_dir}/'...")

    for idx, fpath in enumerate(files):
        state   = read_state(fpath)
        fmask   = state["id"] < args.nfluid
        bmask   = ~fmask
        cvals,_ = get_color_data(state, args.color)

        ax.cla()
        ax.set_xlim(-0.05*Lx*1e3, 1.05*Lx*1e3)
        ax.set_ylim(-0.05*Ly*1e3, 1.05*Ly*1e3)
        ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
        ax.set_aspect("equal"); ax.set_facecolor("#111111")
        ax.scatter(state["x"][bmask]*1e3, state["y"][bmask]*1e3,
                   s=4, c="#888888", zorder=2)
        ax.scatter(state["x"][fmask]*1e3, state["y"][fmask]*1e3,
                   s=8, c=cvals[fmask], cmap="plasma", norm=norm, zorder=3)
        step_num = int(os.path.basename(fpath).replace("state_",""))
        ax.set_title(f"SPH  paso {step_num}", color="white")

        out_png = os.path.join(frames_dir, f"frame_{idx:05d}.png")
        fig.savefig(out_png, dpi=120, facecolor=fig.get_facecolor())

        if idx % 10 == 0:
            print(f"  {idx+1}/{len(files)}  ->  {out_png}")

    print("Frames guardados.")
    print("Para convertir a video con ffmpeg:")
    print(f"  ffmpeg -r 25 -i {frames_dir}/frame_%05d.png "
          f"-vcodec libx264 -pix_fmt yuv420p sph_animation.mp4")
    sys.exit(0)

# --------------------------------------------------------------------------
# Modo default: animacion MP4 con FuncAnimation
# --------------------------------------------------------------------------
def update(frame_idx):
    """Funcion de actualizacion llamada por FuncAnimation en cada frame."""
    state   = read_state(files[frame_idx])
    fmask   = state["id"] < args.nfluid
    bmask   = ~fmask
    cvals,_ = get_color_data(state, args.color)

    # Actualiza posiciones del fluido
    sc_fluid.set_offsets(
        np.column_stack([state["x"][fmask]*1e3, state["y"][fmask]*1e3])
    )
    # Actualiza colores del fluido
    sc_fluid.set_array(cvals[fmask])

    # Actualiza posiciones de la frontera
    sc_bound.set_offsets(
        np.column_stack([state["x"][bmask]*1e3, state["y"][bmask]*1e3])
    )

    # Actualiza vectores de velocidad si se pidio
    if args.quiver:
        xi  = state["x"][fmask][::8]  * 1e3
        yi  = state["y"][fmask][::8]  * 1e3
        vxi = state["vx"][fmask][::8]
        vyi = state["vy"][fmask][::8]
        qv.set_offsets(np.column_stack([xi, yi]))
        qv.set_UVC(vxi, vyi)

    step_num = int(os.path.basename(files[frame_idx]).replace("state_",""))
    title.set_text(f"SPH  paso {step_num}/{len(all_files)-1}")

    return sc_fluid, sc_bound, title

print("Generando animacion MP4...")
plt.tight_layout()

ani = animation.FuncAnimation(
    fig,
    update,
    frames=len(files),    # numero total de frames
    interval=40,          # ms entre frames (25 fps)
    blit=True             # solo redibuja lo que cambia
)

# Escritor de video: usa ffmpeg, que ya esta instalado en el sistema
writer = animation.FFMpegWriter(fps=25, bitrate=1800)
output_video = "sph_animation.mp4"

ani.save(output_video, writer=writer, dpi=120,
         savefig_kwargs={"facecolor": fig.get_facecolor()})

print(f"Video guardado: {output_video}")
