/*
 * sph.c -- Simulacion SPH 2D del problema de la cavidad con tapa movil
 *          (lid-driven cavity) con busqueda de vecinos por lista enlazada
 *          de celdas (cell linked list).
 *
 * Compilacion:
 *     gcc -O2 -o sph sph.c -lm
 *
 * Uso:
 *     mkdir -p output
 *     ./sph <numero_de_pasos>
 *
 * La fisica es identica a la version de fuerza bruta; lo unico que cambia
 * es como se localizan los vecinos:
 *
 *   - Fuerza bruta: cada particula se compara con las nPart-1 restantes.
 *     Costo O(N^2).
 *   - Celdas: el dominio se divide en una malla de celdas cuadradas de lado
 *     igual al radio de soporte del kernel (kappa*h). Como el kernel se
 *     anula mas alla de ese radio, los vecinos de una particula solo pueden
 *     estar en su propia celda o en las 8 adyacentes. Costo O(N).
 *
 * La malla se representa con dos arreglos (esquema clasico de Hockney &
 * Eastwood): cellHead[c] guarda el indice de la primera particula de la
 * celda c, y cellNext[i] el indice de la siguiente particula de la misma
 * celda que i (-1 marca el final). No se reserva memoria por celda, de modo
 * que no hay un limite artificial de particulas por celda.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define X 0
#define Y 1

#define KAPPA 2.0   /* radio de soporte del kernel en unidades de h */

typedef struct
{
  int id;
  double pos[2];
  double vel[2];
  double accel[2];
  double mass;
  double rho;
  double h;
  double p;
  double c;
  double du;
  double u;
  int *nn;          /* indices globales de los vecinos */
  int nNeighbors;   /* numero de vecinos encontrados */
  int nnCap;        /* capacidad reservada en los arreglos de vecinos */
  double *dx;       /* x_i - x_j por vecino */
  double *dy;       /* y_i - y_j por vecino */
  double *r;        /* distancia por vecino */
  double *W;        /* kernel evaluado por vecino */
  double *dWx;      /* gradiente x del kernel por vecino */
  double *dWy;      /* gradiente y del kernel por vecino */
  int type;         /* 1 = fluido, -1 = frontera */
}Particles;

Particles *part, *auxPart;

int nFluid, nPart;

/* ---- Malla de celdas ---- */
int *cellHead;      /* primera particula de cada celda, -1 si esta vacia */
int *cellNext;      /* siguiente particula de la misma celda, -1 al final */
int nCellsX;        /* numero de celdas en x */
int nCellsY;        /* numero de celdas en y */
int nCells;         /* nCellsX*nCellsY */
double cellSize;    /* lado de la celda = kappa*h */
double domLx;       /* extension del dominio en x */
double domLy;       /* extension del dominio en y */

void ics(int nx, int ny, double dx, double dy, double Lx, double Ly);
double W(double r, double h);
double dW(double r, double dx, double h);
void testKernel(void);
void buildCellGrid(double Lx, double Ly, double hMax);
void buildCellList(void);
void freeCellGrid(void);
void NN(int i);
void test_NN(void);
void density(void);
void eos(void);
void navierStokes(void);
void viscosity(double dx);
void boundaryInteraction(double dx);
void meanVelocity(void);
void acceleration(double dx);
void drift(double dt);
void kick(double dt);
void printState(char *outfile);

int main(int argc, char *argv[])
{
  int i, nx, ny, counter;

  double Lx, Ly, dx, dy;
  double dt = 5e-5;
  double t, tTotal;
  char outfiles[500];

  if( argc < 2 )
    {
      printf("Uso: %s <numero_de_pasos>\n", argv[0]);
      return 1;
    }

  tTotal = atoi(argv[1])*dt;

  printf("voy a correr durante %d pasos, un tiempo total de %lf s\n",atoi(argv[1]),tTotal);

  nx = 40;
  ny = 40;
  Lx = 1e-3;
  Ly = 1e-3;

  dx = Lx/nx;   /* espaciamiento entre particulas en x */
  dy = Ly/ny;   /* espaciamiento entre particulas en y */

  nFluid = nx*ny;

  /* Reserva inicial solo para las particulas de fluido.
     ics() ampliara este bloque con realloc para anadir la frontera. */
  part = (Particles *)malloc((size_t)nFluid*sizeof(Particles));
  if( part==NULL )
    {
      printf("Error alocando part\n");
      exit(1);
    }

  /* Condiciones iniciales: fluido + particulas de frontera */
  ics( nx, ny, dx, dy, Lx, Ly);

  /* Prueba del kernel */
  testKernel();

  /* Malla de celdas. Todas las particulas usan h = dx, de modo que el
     radio de soporte maximo es kappa*dx y ese es el lado de la celda. */
  buildCellGrid( Lx, Ly, dx );

  counter = 0;
  t = 0;

  /* Estado inicial del sistema */
  sprintf(outfiles,"./output/state_%.4d",counter);
  printState(outfiles);

  /* Bucle principal de integracion temporal (leap-frog drift-kick-drift) */
  while( t<=tTotal )
    {
      /* Se reconstruye la lista de celdas: las particulas se movieron en
         el paso anterior, asi que la asignacion celda->particula caduca. */
      buildCellList();

      /* Busqueda de vecinos para cada particula de fluido */
      for( i=0; i<nFluid; i++ )
        NN(i);

      /* Verificacion de la busqueda de vecinos, solo en el primer paso */
      if( counter==0 )
        test_NN();

      /* Densidad por suma de kernels */
      density();

      /* Primer medio paso de deriva (posicion) */
      drift(dt);

      /* Calculo de la aceleracion (presion + viscosidad + frontera + XSPH) */
      acceleration(dx);

      /* Patada: actualizacion de la velocidad con la aceleracion */
      kick(dt);

      /* Segundo medio paso de deriva (posicion) */
      drift(dt);

      t = t + dt;
      counter++;

      /* Volcado del estado del sistema */
      sprintf(outfiles,"./output/state_%.4d",counter);
      printState(outfiles);

      printf("step = %d \n",counter);
    }

  /* Liberacion de la memoria por particula */
  for( i=0; i<nPart; i++ )
    {
      free(part[i].nn);
      free(part[i].dx);
      free(part[i].dy);
      free(part[i].r);
      free(part[i].W);
      free(part[i].dWx);
      free(part[i].dWy);
    }
  freeCellGrid();
  free(part);

  return 0;
}

void ics(int nx, int ny, double dx, double dy, double Lx, double Ly)
{
  int i, j, counter;

  FILE *fFluidIcs, *fbBorder, *frBorder, *ftBorder, *flBorder;
  fFluidIcs = fopen("fluid_ics.output","w");

  /* ---- Condiciones iniciales del fluido ---- */
  counter = 0;
  for( j=0; j<ny; j++)
    {
      for( i=0; i<nx; i++)
	{
	  part[counter].id = counter;
	  part[counter].pos[X] = i*dx+dx/2.0;
	  part[counter].pos[Y] = j*dy+dy/2.0;
	  part[counter].vel[X] = 0.0;
	  part[counter].vel[Y] = 0.0;
	  part[counter].accel[X] = 0.0;
	  part[counter].accel[Y] = 0.0;
	  part[counter].rho = 1000;
	  part[counter].h = dx;
	  part[counter].mass = part[counter].rho*dx*dy;
	  part[counter].p = 0.0;
	  part[counter].c = 0.0;
	  part[counter].du = 0.0;
	  part[counter].u = 357.1;
	  part[counter].nn = NULL;
	  part[counter].nNeighbors = 0;
	  part[counter].nnCap = 0;
	  part[counter].dx = NULL;
	  part[counter].dy = NULL;
	  part[counter].r = NULL;
	  part[counter].W = NULL;
	  part[counter].dWx = NULL;
	  part[counter].dWy = NULL;
	  part[counter].type = 1;
	  counter++;
	}
    }

  /* ---- Condiciones iniciales de la frontera ---- */

  /* velocidad de la tapa movil superior */
  double vBoundary = 1.5e-2;

  int npVirtI = 320;
  int npV = npVirtI/4;

  /* ---- Borde inferior ---- */
  fbBorder = fopen("bottom_border.output","w");

  nPart = nFluid;

  auxPart = NULL;
  auxPart = (Particles *)realloc(part, (size_t)(nPart+npV+1)*sizeof(Particles));
  if(auxPart==NULL)
    {
      printf("error en auxPart\n");
      exit(1);
    }
  else
    {
      part = auxPart;
      auxPart = NULL;
    }

  counter = nPart;

  for( i=0; i<=npV; i++)
    {
      part[counter].id = counter;
      part[counter].pos[X] = i*dx/2.0;
      part[counter].pos[Y] = 0.0;
      part[counter].vel[X] = 0.0;
      part[counter].vel[Y] = 0.0;
      part[counter].accel[X] = 0.0;
      part[counter].accel[Y] = 0.0;
      part[counter].rho = 1000;
      part[counter].h = dx;
      part[counter].mass = part[counter].rho*dx*dy;
      part[counter].p = 0.0;
      part[counter].c = 0.0;
      part[counter].du = 0.0;
      part[counter].u = 357.1;
      part[counter].nn = NULL;
      part[counter].nNeighbors = 0;
      part[counter].nnCap = 0;
      part[counter].dx = NULL;
      part[counter].dy = NULL;
      part[counter].r = NULL;
      part[counter].W = NULL;
      part[counter].dWx = NULL;
      part[counter].dWy = NULL;
      part[counter].type = -1;
      counter++;
    }

  for( i=nPart; i<nPart+npV+1; i++)
    {
      fprintf(fbBorder,"%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
	      part[i].id,
	      part[i].pos[X],part[i].pos[Y],
	      part[i].vel[X],part[i].vel[Y],
	      part[i].accel[X],part[i].accel[Y],
	      part[i].rho,part[i].mass,
	      part[i].p,part[i].c,part[i].u);
    }
  fclose(fbBorder);

  /* ---- Borde derecho ---- */
  frBorder = fopen("right_border.output","w");

  nPart = counter;

  auxPart = NULL;
  auxPart = (Particles *)realloc(part,(size_t)(nPart+npV-1)*sizeof(Particles));
  if(auxPart==NULL)
    {
      printf("error en auxPart\n");
      exit(1);
    }
  else
    {
      part = auxPart;
      auxPart = NULL;
    }

  for( i=0; i<npV-1; i++)
    {
      part[counter].id = counter;
      part[counter].pos[X] = Lx;
      part[counter].pos[Y] = dy/2.0 + i*dy/2.0;
      part[counter].vel[X] = 0.0;
      part[counter].vel[Y] = 0.0;
      part[counter].accel[X] = 0.0;
      part[counter].accel[Y] = 0.0;
      part[counter].rho = 1000;
      part[counter].h = dx;
      part[counter].mass = part[counter].rho*dx*dy;
      part[counter].p = 0.0;
      part[counter].c = 0.0;
      part[counter].du = 0.0;
      part[counter].u = 357.1;
      part[counter].nn = NULL;
      part[counter].nNeighbors = 0;
      part[counter].nnCap = 0;
      part[counter].dx = NULL;
      part[counter].dy = NULL;
      part[counter].r = NULL;
      part[counter].W = NULL;
      part[counter].dWx = NULL;
      part[counter].dWy = NULL;
      part[counter].type = -1;
      counter++;
    }

  for( i=nPart; i<nPart+npV-1; i++)
    {
      fprintf(frBorder,"%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
	      part[i].id,
	      part[i].pos[X],part[i].pos[Y],
	      part[i].vel[X],part[i].vel[Y],
	      part[i].accel[X],part[i].accel[Y],
	      part[i].rho,part[i].mass,
	      part[i].p,part[i].c,part[i].u);
    }
  fclose(frBorder);

  /* ---- Borde superior (tapa movil) ---- */
  ftBorder = fopen("top_border.output","w");

  nPart = counter;

  auxPart = NULL;
  auxPart = (Particles *)realloc(part,(size_t)(nPart+npV+1)*sizeof(Particles));
  if(auxPart==NULL)
    {
      printf("error en auxPart\n");
      exit(1);
    }
  else
    {
      part = auxPart;
      auxPart = NULL;
    }

  for( i=0; i<=npV; i++)
    {
      part[counter].id = counter;
      part[counter].pos[X] = i*dx/2.0;
      part[counter].pos[Y] = Ly;
      part[counter].vel[X] = vBoundary;   /* tapa que se mueve en +X */
      part[counter].vel[Y] = 0.0;
      part[counter].accel[X] = 0.0;
      part[counter].accel[Y] = 0.0;
      part[counter].rho = 1000;
      part[counter].h = dx;
      part[counter].mass = part[counter].rho*dx*dy;
      part[counter].p = 0.0;
      part[counter].c = 0.0;
      part[counter].du = 0.0;
      part[counter].u = 357.1;
      part[counter].nn = NULL;
      part[counter].nNeighbors = 0;
      part[counter].nnCap = 0;
      part[counter].dx = NULL;
      part[counter].dy = NULL;
      part[counter].r = NULL;
      part[counter].W = NULL;
      part[counter].dWx = NULL;
      part[counter].dWy = NULL;
      part[counter].type = -1;
      counter++;
    }

  for( i=nPart; i<nPart+npV+1; i++)
    {
      fprintf(ftBorder,"%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
	      part[i].id,
	      part[i].pos[X],part[i].pos[Y],
	      part[i].vel[X],part[i].vel[Y],
	      part[i].accel[X],part[i].accel[Y],
	      part[i].rho,part[i].mass,
	      part[i].p,part[i].c,part[i].u);
    }
  fclose(ftBorder);

  /* ---- Borde izquierdo ---- */
  flBorder = fopen("left_border.output","w");

  nPart = counter;

  auxPart = NULL;
  auxPart = (Particles *)realloc(part,(size_t)(nPart+npV-1)*sizeof(Particles));
  if(auxPart==NULL)
    {
      printf("error en auxPart\n");
      exit(1);
    }
  else
    {
      part = auxPart;
      auxPart = NULL;
    }

  for( i=0; i<npV-1; i++)
    {
      part[counter].id = counter;
      part[counter].pos[X] = 0.0;
      part[counter].pos[Y] = dy/2.0 + i*dy/2.0;
      part[counter].vel[X] = 0.0;
      part[counter].vel[Y] = 0.0;
      part[counter].accel[X] = 0.0;
      part[counter].accel[Y] = 0.0;
      part[counter].rho = 1000;
      part[counter].h = dx;
      part[counter].mass = part[counter].rho*dx*dy;
      part[counter].p = 0.0;
      part[counter].c = 0.0;
      part[counter].du = 0.0;
      part[counter].u = 357.1;
      part[counter].nn = NULL;
      part[counter].nNeighbors = 0;
      part[counter].nnCap = 0;
      part[counter].dx = NULL;
      part[counter].dy = NULL;
      part[counter].r = NULL;
      part[counter].W = NULL;
      part[counter].dWx = NULL;
      part[counter].dWy = NULL;
      part[counter].type = -1;
      counter++;
    }

  for( i=nPart; i<nPart+npV-1; i++)
    {
      fprintf(flBorder,"%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
	      part[i].id,
	      part[i].pos[X],part[i].pos[Y],
	      part[i].vel[X],part[i].vel[Y],
	      part[i].accel[X],part[i].accel[Y],
	      part[i].rho,part[i].mass,
	      part[i].p,part[i].c,part[i].u);
    }
  fclose(flBorder);

  /* ---- Volcado de todas las particulas ---- */
  nPart = counter;

  for( i=0; i<nPart; i++)
    {
      fprintf(fFluidIcs,"%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
	      part[i].id,
	      part[i].pos[X],part[i].pos[Y],
	      part[i].vel[X],part[i].vel[Y],
	      part[i].accel[X],part[i].accel[Y],
	      part[i].rho,part[i].mass,
	      part[i].p,part[i].c,part[i].u);
    }
  fclose(fFluidIcs);
}

double W(double r, double h)
{
  double R = r/h;
  double alpha = 15.0/(7.0*M_PI*h*h);   /* normalizacion 2D del spline cubico */

  if( (R >= 0.0) && (R < 1.0) )
    return alpha*((2.0/3.0) - R*R + 0.5*R*R*R);

  if( (R >= 1.0) && (R <= 2.0) )
    return alpha*((1.0/6.0)*(2.0-R)*(2.0-R)*(2.0-R));

  return 0.0;
}

double dW(double r, double dx, double h)
{
  double R = r/h;
  double alpha = 15.0/(7.0*M_PI*h*h);

  /* Proteccion contra division por cero cuando dos particulas coinciden */
  if( r < 1e-12 )
    return 0.0;

  if( (R >= 0.0) && (R < 1.0) )
    return alpha*(-2.0 + 1.5*R)*dx/(h*h);

  if( (R >= 1.0) && (R <= 2.0) )
    return alpha*(-0.5*(2.0-R)*(2.0-R))*dx/(h*h*R);

  return 0.0;
}

void testKernel(void)
{
  double r, w, dw;

  FILE *fKernelTest;
  fKernelTest = fopen("kernel_test.output","w");

  for( r=-3.0; r<=3.0; r = r + 0.1)
    {
      w = W( fabs(r), 1.0);
      dw = dW( fabs(r), r/sqrt(3.0), 1.0);

      fprintf(fKernelTest,"%16.10lf %16.10lf %16.10lf\n",r,w,dw);
    }

  fclose(fKernelTest);
}

/* ------------------------------------------------------------------
   Malla de celdas
   ------------------------------------------------------------------ */

/* Reserva la malla. El lado de la celda es el radio de soporte del kernel,
   kappa*hMax, de modo que cualquier vecino de una particula cae forzosamente
   en su celda o en una de las 8 adyacentes. Celdas mas pequenas obligarian a
   revisar mas de 3x3 celdas; mas grandes desperdiciarian comparaciones. */
void buildCellGrid(double Lx, double Ly, double hMax)
{
  cellSize = KAPPA*hMax;

  domLx = Lx;
  domLy = Ly;

  nCellsX = (int)ceil(Lx/cellSize);
  nCellsY = (int)ceil(Ly/cellSize);

  if( nCellsX < 1 ) nCellsX = 1;
  if( nCellsY < 1 ) nCellsY = 1;

  nCells = nCellsX*nCellsY;

  cellHead = (int *)malloc((size_t)nCells*sizeof(int));
  cellNext = (int *)malloc((size_t)nPart*sizeof(int));

  if( cellHead==NULL || cellNext==NULL )
    {
      printf("Error alocando la malla de celdas\n");
      exit(1);
    }

  printf("malla de celdas: %d x %d = %d celdas de lado %.3e m\n",
	 nCellsX, nCellsY, nCells, cellSize);
}

void freeCellGrid(void)
{
  free(cellHead);  cellHead = NULL;
  free(cellNext);  cellNext = NULL;
}

/* Indice de celda en x de una coordenada, recortado al rango valido.
   El recorte es necesario porque una particula de fluido puede salirse del
   dominio; al quedar en la celda del borde sigue encontrando a sus vecinos,
   ya que cualquier particula a menos de cellSize de ella tambien se recorta
   a esa misma franja. El filtro por distancia descarta el resto. */
static int cellIndexX(double x)
{
  int ix = (int)floor(x/cellSize);
  if( ix < 0 ) ix = 0;
  if( ix > nCellsX-1 ) ix = nCellsX-1;
  return ix;
}

static int cellIndexY(double y)
{
  int iy = (int)floor(y/cellSize);
  if( iy < 0 ) iy = 0;
  if( iy > nCellsY-1 ) iy = nCellsY-1;
  return iy;
}

/* Reparte las nPart particulas (fluido y frontera) entre las celdas.
   Se recorre al reves para que cada lista quede en orden creciente de indice,
   lo que hace la salida reproducible. Costo O(N). */
void buildCellList(void)
{
  int i, c;

  for( c=0; c<nCells; c++ )
    cellHead[c] = -1;

  for( i=nPart-1; i>=0; i-- )
    {
      c = cellIndexY(part[i].pos[Y])*nCellsX + cellIndexX(part[i].pos[X]);
      cellNext[i] = cellHead[c];
      cellHead[c] = i;
    }
}

/* ------------------------------------------------------------------
   Busqueda de vecinos
   ------------------------------------------------------------------ */

/* Amplia los arreglos de vecinos de la particula i hasta admitir 'needed'
   entradas. La capacidad se duplica en lugar de crecer de uno en uno, de
   modo que el numero de realloc por paso es logaritmico y no lineal. */
static void growNeighborArrays(int i, int needed)
{
  int newCap;
  void *aux;

  if( needed <= part[i].nnCap )
    return;

  newCap = (part[i].nnCap==0) ? 16 : part[i].nnCap;
  while( newCap < needed )
    newCap = 2*newCap;

  aux = realloc(part[i].nn,  (size_t)newCap*sizeof(int));
  if(aux==NULL) { printf("Error ampliando vecinos de %d\n",i); exit(1); }
  part[i].nn = (int *)aux;

  aux = realloc(part[i].dx,  (size_t)newCap*sizeof(double));
  if(aux==NULL) { printf("Error ampliando vecinos de %d\n",i); exit(1); }
  part[i].dx = (double *)aux;

  aux = realloc(part[i].dy,  (size_t)newCap*sizeof(double));
  if(aux==NULL) { printf("Error ampliando vecinos de %d\n",i); exit(1); }
  part[i].dy = (double *)aux;

  aux = realloc(part[i].r,   (size_t)newCap*sizeof(double));
  if(aux==NULL) { printf("Error ampliando vecinos de %d\n",i); exit(1); }
  part[i].r = (double *)aux;

  aux = realloc(part[i].W,   (size_t)newCap*sizeof(double));
  if(aux==NULL) { printf("Error ampliando vecinos de %d\n",i); exit(1); }
  part[i].W = (double *)aux;

  aux = realloc(part[i].dWx, (size_t)newCap*sizeof(double));
  if(aux==NULL) { printf("Error ampliando vecinos de %d\n",i); exit(1); }
  part[i].dWx = (double *)aux;

  aux = realloc(part[i].dWy, (size_t)newCap*sizeof(double));
  if(aux==NULL) { printf("Error ampliando vecinos de %d\n",i); exit(1); }
  part[i].dWy = (double *)aux;

  part[i].nnCap = newCap;
}

/* Busqueda de los vecinos de la particula i recorriendo unicamente el bloque
   de 3x3 celdas centrado en su celda. Requiere que buildCellList() se haya
   ejecutado con las posiciones actuales. */
void NN(int i)
{
  int j, ix, iy, cx, cy, c, nNeighbors;
  double xij, yij, rij, hij, rcut;

  ix = cellIndexX(part[i].pos[X]);
  iy = cellIndexY(part[i].pos[Y]);

  nNeighbors = 0;

  for( cy=iy-1; cy<=iy+1; cy++ )
    {
      if( cy<0 || cy>nCellsY-1 )
	continue;

      for( cx=ix-1; cx<=ix+1; cx++ )
	{
	  if( cx<0 || cx>nCellsX-1 )
	    continue;

	  c = cy*nCellsX + cx;

	  /* Recorrido de la lista enlazada de la celda c */
	  for( j=cellHead[c]; j!=-1; j=cellNext[j] )
	    {
	      if( i==j )
		continue;

	      hij = 0.5*(part[i].h + part[j].h);   /* h promediado */
	      rcut = KAPPA*hij;

	      xij = part[i].pos[X] - part[j].pos[X];
	      yij = part[i].pos[Y] - part[j].pos[Y];

	      /* Descarte barato antes de la raiz cuadrada */
	      if( fabs(xij) > rcut || fabs(yij) > rcut )
		continue;

	      rij = sqrt( xij*xij + yij*yij );

	      if( rij <= rcut + 1e-12 )            /* dentro del soporte */
		{
		  nNeighbors++;
		  growNeighborArrays(i, nNeighbors);

		  part[i].nn[nNeighbors-1]  = j;
		  part[i].dx[nNeighbors-1]  = xij;
		  part[i].dy[nNeighbors-1]  = yij;
		  part[i].r[nNeighbors-1]   = rij;
		  part[i].W[nNeighbors-1]   = W( rij, hij );
		  part[i].dWx[nNeighbors-1] = dW( rij, xij, hij );
		  part[i].dWy[nNeighbors-1] = dW( rij, yij, hij );
		}
	    }
	}
    }

  part[i].nNeighbors = nNeighbors;
}

void test_NN(void)
{
  int i,j,k;

  FILE *fTestNN;
  fTestNN = fopen("NN_test.output","w");
  srand(time(NULL));

  for( k=0; k<20; k++)
    {
      i = rand() % nFluid;

      printf("testing for particle %d\n",i);
      printf("with %d neighbors\n",part[i].nNeighbors);

      fprintf(fTestNN,"%16d %16.10lf %16.10lf\n",
	      part[i].id,
	      part[i].pos[X],
	      part[i].pos[Y]);

      for( j=0; j<part[i].nNeighbors; j++ )
	fprintf(fTestNN,"%16d %16.10lf %16.10lf\n",
		part[i].nn[j],
		part[part[i].nn[j]].pos[X],
		part[part[i].nn[j]].pos[Y]);
      fprintf(fTestNN,"\n");
    }
  fclose(fTestNN);
}

void density(void)
{
  int i, j;
  double wii, norm;

  for( i=0; i<nFluid; i++ )
    {
      /* autodensidad */
      wii = W( 0.0, part[i].h );

      part[i].rho = part[i].mass*wii;
      for( j=0; j<part[i].nNeighbors; j++ )
	part[i].rho = part[i].rho + part[part[i].nn[j]].mass*part[i].W[j];

      /* normalizacion (Shepard) de la densidad */
      norm = (part[i].mass/part[i].rho)*wii;
      for( j=0; j<part[i].nNeighbors; j++ )
	norm = norm + (part[part[i].nn[j]].mass/part[part[i].nn[j]].rho)*part[i].W[j];

      part[i].rho = part[i].rho/norm;
    }

  printf("density computed\n");
}

void eos(void)
{
  int i;

  for( i=0; i<nPart; i++ )
    {
      part[i].c = 0.01;
      part[i].p = part[i].c*part[i].c*part[i].rho;
    }
}

void navierStokes(void)
{
  int i, j, k;
  double pij, vdw;

  /* velocidad del sonido y presion */
  eos();

  for( i=0; i<nFluid; i++ )
    {
      part[i].accel[X] = part[i].accel[Y] = 0.0;
      part[i].du = 0.0;

      for( k=0; k<part[i].nNeighbors; k++ )
	{
	  j = part[i].nn[k];

	  pij = ( part[i].p/(part[i].rho*part[i].rho) )
	    + ( part[j].p/(part[j].rho*part[j].rho) );
	  part[i].accel[X] = part[i].accel[X] - part[j].mass*pij*part[i].dWx[k];
	  part[i].accel[Y] = part[i].accel[Y] - part[j].mass*pij*part[i].dWy[k];

	  vdw = (part[i].vel[X]-part[j].vel[X])*part[i].dWx[k]
	    + (part[i].vel[Y]-part[j].vel[Y])*part[i].dWy[k];
	  part[i].du = part[i].du + 0.5*part[j].mass*pij*vdw;
	}
    }

  printf("acceleration computed\n");
}

void viscosity(double dx)
{
  int i, j, k;

  double xij, yij, vxij, vyij, vijrij, vdw;
  double hij, cij, phiij, rhoij, Piij;
  double alphapi = 1.0;
  double betapi = 1.0;
  double eps = dx;
  double eps2 = 0.01*eps*eps;

  for( i=0; i<nFluid; i++ )
    {
      for( k=0; k<part[i].nNeighbors ; k++ )
	{
	  j = part[i].nn[k];

	  xij = part[i].pos[X] - part[j].pos[X];
	  yij = part[i].pos[Y] - part[j].pos[Y];
	  vxij = part[i].vel[X] - part[j].vel[X];
	  vyij = part[i].vel[Y] - part[j].vel[Y];
	  vijrij = vxij*xij + vyij*yij;

	  if( vijrij < 0.0 )
	    {
	      hij = 0.5*(part[i].h+part[j].h);
	      phiij = (hij*vijrij)/( xij*xij + yij*yij + eps2);
	      cij = 0.5*(part[i].c+part[j].c);
	      rhoij = 0.5*(part[i].rho+part[j].rho);

	      Piij = ( -alphapi*cij*phiij + betapi*phiij*phiij )/( rhoij );

	      part[i].accel[X] = part[i].accel[X] - part[j].mass*Piij*part[i].dWx[k];
	      part[i].accel[Y] = part[i].accel[Y] - part[j].mass*Piij*part[i].dWy[k];

	      vdw = (part[i].vel[X]-part[j].vel[X])*part[i].dWx[k]
		+ (part[i].vel[Y]-part[j].vel[Y])*part[i].dWy[k];
	      part[i].du = part[i].du + 0.5*part[j].mass*Piij*vdw;
	    }
	}
    }
  printf("viscosity computed\n");
}

void boundaryInteraction(double dx)
{
  int i, j;
  int n1 = 12, n2 = 4;
  double r0 = dx/2.0, D = 0.01;
  double xij, yij, rij, PBxij, PByij;

  for( i=0; i<nFluid; i++ )
    {
      for( j=0; j<part[i].nNeighbors; j++ )
	{
	  if( part[part[i].nn[j]].type==-1 )
	    {
	      xij = part[i].pos[X] - part[part[i].nn[j]].pos[X];
	      yij = part[i].pos[Y] - part[part[i].nn[j]].pos[Y];
	      rij = sqrt( xij*xij + yij*yij );

	      if( rij<r0 )
		{
		  PBxij = D*( pow((r0/rij),n1) - pow((r0/rij),n2) )*(xij/(rij*rij));
		  PByij = D*( pow((r0/rij),n1) - pow((r0/rij),n2) )*(yij/(rij*rij));

		  part[i].accel[X] = part[i].accel[X] + PBxij;
		  part[i].accel[Y] = part[i].accel[Y] + PByij;
		}
	    }
	}
    }
  printf("interaction with boundary computed\n");
}

void meanVelocity(void)
{
  int i, j;
  double epsilon = 0.3;
  double vxMean, vyMean;
  double vxij, vyij, rhoij;

  for( i=0; i<nFluid; i++ )
    {
      vxMean = 0.0;
      vyMean = 0.0;

      for( j=0; j<part[i].nNeighbors; j++ )
	{
	  vxij = part[i].vel[X] - part[part[i].nn[j]].vel[X];
	  vyij = part[i].vel[Y] - part[part[i].nn[j]].vel[Y];
	  rhoij = 0.5*(part[i].rho+part[part[i].nn[j]].rho);
	  vxMean = vxMean + (part[part[i].nn[j]].mass/rhoij)*vxij*part[i].W[j];
	  vyMean = vyMean + (part[part[i].nn[j]].mass/rhoij)*vyij*part[i].W[j];
	}

      part[i].vel[X] = part[i].vel[X] - epsilon*vxMean;
      part[i].vel[Y] = part[i].vel[Y] - epsilon*vyMean;
    }
}

void acceleration(double dx)
{
  /* aceleracion y cambio de energia (presion) */
  navierStokes();

  /* contribucion de la viscosidad artificial */
  viscosity(dx);

  /* interaccion repulsiva con la frontera */
  boundaryInteraction(dx);

  /* correccion XSPH de la velocidad media */
  meanVelocity();

  printf("acceleration computed\n");
}

void drift(double dt)
{
  int i;

  for( i=0; i<nFluid; i++ )
    {
      part[i].pos[X] = part[i].pos[X] + 0.5*dt*part[i].vel[X];
      part[i].pos[Y] = part[i].pos[Y] + 0.5*dt*part[i].vel[Y];
      part[i].u = part[i].u + 0.5*dt*part[i].du;
    }
}

void kick(double dt)
{
  int i;

  for( i=0; i<nFluid; i++ )
    {
      part[i].vel[X] = part[i].vel[X] + dt*part[i].accel[X];
      part[i].vel[Y] = part[i].vel[Y] + dt*part[i].accel[Y];
    }
}

void printState(char *outfile)
{
  int i;

  FILE *fState;
  fState = fopen(outfile,"w");

  /* Si el directorio ./output no existe, fopen devuelve NULL. */
  if( fState==NULL )
    {
      printf("No se pudo abrir %s. Cree el directorio ./output\n", outfile);
      exit(1);
    }

  for( i=0; i<nPart; i++)
    {
      fprintf(fState,"%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
	      part[i].id,
	      part[i].pos[X],part[i].pos[Y],
	      part[i].vel[X],part[i].vel[Y],
	      part[i].accel[X],part[i].accel[Y],
	      part[i].rho,part[i].mass,
	      part[i].p,part[i].c,part[i].u);
    }

  fclose(fState);
}
