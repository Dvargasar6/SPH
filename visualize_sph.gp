# visualize_sph.gp
# ---------------------------------------------------------------
# Uso:
#   gnuplot visualize_sph.gp                   # abre ventana interactiva
#   gnuplot -e "step=50" visualize_sph.gp      # muestra el paso 50
#   gnuplot -e "make_png=1" visualize_sph.gp   # exporta todos los frames PNG
# ---------------------------------------------------------------

# Paso a mostrar (por defecto el ultimo disponible de la lista)
if (!exists("step")) step = 0
if (!exists("make_png")) make_png = 0

# Numero de particulas de fluido (primeras nfluid lineas de cada archivo)
nfluid = 1600

# Dominio en mm
Lx = 1.0
Ly = 1.0

# ---------------------------------------------------------------
# Funcion auxiliar: construye el nombre del archivo dado el numero
# ---------------------------------------------------------------
statefile(n) = sprintf("./output/state_%04d", n)

# ---------------------------------------------------------------
# Terminal interactiva (wxt es la mas comun en Arch con gnuplot-qt)
# ---------------------------------------------------------------
if (!make_png) {
    set terminal wxt size 700,700 enhanced font "Sans,10" persist
}

# ---------------------------------------------------------------
# Estilo del grafico
# ---------------------------------------------------------------
set size ratio 1
set xlabel "x [mm]"
set ylabel "y [mm]"
set xrange [-0.05 : 1.05]
set yrange [-0.05 : 1.05]
set palette rgbformulae 30,31,32          # paleta "plasma" aproximada
set cbrange [0 : 0.015]                   # rango de velocidad [m/s]
set cblabel "Velocidad [m/s]"
set colorbox

# ---------------------------------------------------------------
# Macro para dibujar un paso: separa fluido y frontera
# Las columnas del archivo son:
#   $1=id  $2=x  $3=y  $4=vx  $5=vy  $6=ax  $7=ay
#   $8=rho $9=mass $10=p $11=c $12=u
# ---------------------------------------------------------------
plot_step(n) = sprintf(\
    "< awk 'NR<=%d' %s", nfluid, statefile(n))

plot_bound(n) = sprintf(\
    "< awk 'NR>%d' %s", nfluid, statefile(n))

# ---------------------------------------------------------------
# Vista de un unico paso
# ---------------------------------------------------------------
if (!make_png) {
    set title sprintf("SPH  paso %d", step)
    plot \
        plot_bound(step) using ($2*1000):($3*1000) \
            with points pt 7 ps 0.3 lc rgb "#888888" title "Frontera", \
        plot_step(step) using ($2*1000):($3*1000):(sqrt($4**2+$5**2)) \
            with points pt 7 ps 0.5 palette title "Fluido (vel)"
}

# ---------------------------------------------------------------
# Exportacion de todos los frames como PNG
# gnuplot -e "make_png=1" visualize_sph.gp
# Luego convierte con:
#   ffmpeg -r 25 -i frames/frame_%05d.png -vcodec libx264 sph_animation.mp4
# ---------------------------------------------------------------
if (make_png) {
    system("mkdir -p frames")
    set terminal pngcairo size 700,700 enhanced font "Sans,10"

    # Detecta cuantos archivos hay usando shell
    n_files = int(system("ls ./output/state_* | wc -l"))
    print sprintf("Exportando %d frames PNG...", n_files)

    do for [i=0:n_files-1] {
        set output sprintf("frames/frame_%05d.png", i)
        set title sprintf("SPH  paso %d", i)
        plot \
            plot_bound(i) using ($2*1000):($3*1000) \
                with points pt 7 ps 0.3 lc rgb "#888888" notitle, \
            plot_step(i) using ($2*1000):($3*1000):(sqrt($4**2+$5**2)) \
                with points pt 7 ps 0.5 palette notitle
    }
    print "Frames guardados en ./frames/"
}
