#include "bluebottle.h"
#include <sys/stat.h>
// #include <float.h>
// #include <math.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include <stddef.h>
// #include <string.h>
// #include <time.h>
// #include <mpi.h>
// #include <sys/time.h>

#ifndef _TRACER_H
#define _TRACER_H

#define TRACER_DIR "tracer"
#define DEBUG_DIR "debug"
#define PARTICLEINFO_DIR "particleinfo"

void tracer_start(void);
void tracer_execute(double dt);
void tracer_exit(void);

void tracer_init(int i);
double randgen(void);
double *randomizer(double *randvec, double diffusivity, double tstep);
double distance(double *loc, double x, double y, double z);
double *sphere_pusher(double *loc, double x, double y, double z, double r);
int wall_checker(double *loc);
int domain_checker(double *loc);
double *wall_reflector(int i, double *loc);
double *reflector(double *loc1, double *loc2, double x,
  double y, double z, double r);
double *intersector(double *loc1, double *loc2, double x,
    double y, double z, double r);
double index_finder(double *loc);
double grid_finder(double index, int axis, double gridloc);
double *interpolator(double *loc, double *uout, double *vout, double *wout);
void exporter(void);
void tracer_delete(int deleteid);
void tracer_add(void);
void ratekeeper(void);
void maintainer(void);

typedef struct tracer {
  int tracerid;
  double pos[3];
  int tracertype;
} tracer;

extern tracer *tracer_array;
extern real *ucc;
extern real *vcc;
extern real *wcc;
extern long int tracercheck;
extern double diffusivity;
extern double tstep;
extern int foundindex;
extern double interpvel[3];
extern int tracernum;
extern double locx[3];
extern double dist_temp;
extern int idflag;
extern int farwalltreatment;
extern int fixednumber;
extern float fixedrate;
extern float timecheck;
extern int tracertreatment;
extern int tracerdebug;
extern int exportfreq;
extern int exportcount;
extern int readfile;

#endif
