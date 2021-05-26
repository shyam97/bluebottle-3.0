#include "bluebottle.h"
// #include <float.h>
// #include <math.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include <stddef.h>
// #include <string.h>
// #include <time.h>
// #include <mpi.h>
// #include <sys/time.h>

#ifndef TRACER_H
#define TRACER_H

void tracer_init(void);
void tracer_execute(double dt);
void tracer_exit(void);

double randgen(void);
double *randomizer(double *randvec, double diffusivity, double tstep);
double distance(double *loc, double x, double y, double z);
int sphere_checker(double *loc, double x, double y, double z, double r);
double *sphere_pusher(double *loc, double x, double y, double z, double r);
int wall_checker(double *loc);
double *wall_reflector(int i, double *loc);
double *reflector(double *loc1, double *loc2, double x,
  double y, double z, double r);
double index_finder(double *loc);
double grid_finder(double index, int axis, double gridloc);
double *interpolator(double *loc, double *uout, double *vout, double *wout);
void exporter(void);

typedef struct tracer {
  int tracerid;
  double pos[3];
  int tracertype;
} tracer;

extern tracer *tracer_array;
extern int tracercheck;
extern double diffusivity;
extern double tstep;
extern int foundindex;
extern double *interpvel;
extern int tracernum;
extern double locx[3];
double dist;

#endif
