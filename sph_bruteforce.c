#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
/* Se elimina <malloc.h>: no es portable y es redundante.
   malloc/free/realloc estan declarados en <stdlib.h>.
   En glibc (Arch) <malloc.h> existe, pero su uso esta desaconsejado. */

#define X 0
#define Y 1

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
  int nNeighbors;   /* numero de vecinos */
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

double calcular_distancia_entre_puntos(double x1, double x2, double y1, double y2);
void ics(int nx, int ny, double dx, double dy, double Lx, double Ly);
double W(double r, double h);
double dW(double r, double dx, double h);
void testKernel(void);
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

  /* Se requiere el numero de pasos como primer argumento */
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
      exit(0);
    }

  /* Condiciones iniciales: fluido + particulas de frontera */
  ics( nx, ny, dx, dy, Lx, Ly);

  /* Prueba del kernel */
  testKernel();

  counter = 0;
  t = 0;

  /* Estado inicial del sistema */
  sprintf(outfiles,"./output/state_%.4d",counter);
  printState(outfiles);

  /* Bucle principal de integracion temporal (leap-frog drift-kick-drift) */
  while( t<=tTotal )
    {
      /* Busqueda de vecinos para cada particula de fluido.
         NN ahora recorre TODAS las particulas (incluida la frontera),
         de modo que la tapa movil puede arrastrar el fluido. */
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
      exit(0);
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
      exit(0);
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
      exit(0);
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
      exit(0);
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

/* Busqueda de vecinos cercanos.

   CORRECCION PRINCIPAL: la version original tenia el cuerpo del bucle de
   busqueda desacoplado del 'for' (una llave mal colocada), por lo que la
   busqueda se ejecutaba una sola vez con un indice invalido. Ademas el
   esquema de celdas usaba coordenadas de 0 a 100 mientras que las particulas
   viven en [0, 1e-3], de modo que TODAS caian en la misma celda.

   Aqui se hace una busqueda directa O(N^2) sobre las nPart particulas
   (incluida la frontera), con el criterio fisico r_ij <= kappa*h_ij.
   Es correcta y robusta; la unica desventaja es el costo cuadratico. */
void NN(int i)
{
  double kappa = 2.0;          /* radio de soporte en unidades de h */
  double xij, yij, rij, hij;
  double *auxDouble;
  int j, *auxInt, nNeighbors;

  /* Se libera la memoria de la iteracion anterior */
  free(part[i].nn);   part[i].nn  = NULL;
  free(part[i].dx);   part[i].dx  = NULL;
  free(part[i].dy);   part[i].dy  = NULL;
  free(part[i].r);    part[i].r   = NULL;
  free(part[i].W);    part[i].W   = NULL;
  free(part[i].dWx);  part[i].dWx = NULL;
  free(part[i].dWy);  part[i].dWy = NULL;

  nNeighbors = 0;

  for( j=0; j<nPart; j++ )
    {
      if( i==j )
	continue;

      hij = 0.5*(part[i].h + part[j].h);             /* h promediado */
      xij = part[i].pos[X] - part[j].pos[X];
      yij = part[i].pos[Y] - part[j].pos[Y];
      rij = sqrt( xij*xij + yij*yij );

      if( rij <= kappa*hij + 1e-12 )                 /* dentro del soporte */
	{
	  nNeighbors++;

	  /* Se amplian los arreglos de vecinos uno a uno (realloc).
	     Se conserva el patron original basado en variables auxiliares. */
	  auxInt = (int *)realloc(part[i].nn,(size_t)(nNeighbors)*sizeof(int));
	  part[i].nn = auxInt;  auxInt = NULL;

	  auxDouble = (double *)realloc(part[i].dx,(size_t)(nNeighbors)*sizeof(double));
	  part[i].dx = auxDouble;  auxDouble = NULL;

	  auxDouble = (double *)realloc(part[i].dy,(size_t)(nNeighbors)*sizeof(double));
	  part[i].dy = auxDouble;  auxDouble = NULL;

	  auxDouble = (double *)realloc(part[i].r,(size_t)(nNeighbors)*sizeof(double));
	  part[i].r = auxDouble;  auxDouble = NULL;

	  auxDouble = (double *)realloc(part[i].W,(size_t)(nNeighbors)*sizeof(double));
	  part[i].W = auxDouble;  auxDouble = NULL;

	  auxDouble = (double *)realloc(part[i].dWx,(size_t)(nNeighbors)*sizeof(double));
	  part[i].dWx = auxDouble;  auxDouble = NULL;

	  auxDouble = (double *)realloc(part[i].dWy,(size_t)(nNeighbors)*sizeof(double));
	  part[i].dWy = auxDouble;  auxDouble = NULL;

	  /* Se almacenan los datos del vecino j */
	  part[i].nn[nNeighbors-1]  = j;
	  part[i].dx[nNeighbors-1]  = xij;
	  part[i].dy[nNeighbors-1]  = yij;
	  part[i].r[nNeighbors-1]   = rij;
	  part[i].W[nNeighbors-1]   = W( rij, hij );
	  part[i].dWx[nNeighbors-1] = dW( rij, xij, hij );
	  part[i].dWy[nNeighbors-1] = dW( rij, yij, hij );
	}
    }

  part[i].nNeighbors = nNeighbors;
}

double calcular_distancia_entre_puntos(double x1, double x2, double y1, double y2)
{
  double distancia;
  distancia = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
  return distancia;
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

  /* Si el directorio ./output no existe, fopen devuelve NULL.
     Se verifica para no escribir sobre un puntero invalido. */
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