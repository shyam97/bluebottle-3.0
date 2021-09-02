#include "bluebottle.h"
#include "tracer.h"

/*---------------------------------------------------------------------------*/
/* Global variables */
/*---------------------------------------------------------------------------*/

tracer *tracer_array;
real *ucc;
real *vcc;
real *wcc;
double diffusivity;
double tstep;
int foundindex = 0;
double interpvel[3];
int tracernum;
double locx[3];
double dist_temp;
int idflag=0;
int farwalltreatment;
int exportfreq;
int fixednumber;
float fixedrate;
float timecheck = 0;
int tracertreatment;
int exportcount;
int readfile;

/*---------------------------------------------------------------------------*/
/* Initialize tracer simulation and arrays */
/*---------------------------------------------------------------------------*/

void tracer_start(void)
{
  printf("\nInitializing tracer arrays... ");
  readfile = 1;
  tracercheck = 0;
  exportfreq = 100;
  exportcount = 0;
  diffusivity = 1E-4;
  farwalltreatment = 0; //0 for continuing motion, 1 for delete, 2 for stopping at end
  tracertreatment = 0; // 0 for 2type, 1 for fixednumber, 2 for fixedrate

  fixednumber = 500;

  fixedrate = 0.01;
  timecheck = ttime + fixedrate;

  srand((unsigned int)time(NULL));

  tracernum = 2000;
  tracer_array = (tracer*) malloc(tracernum * sizeof(tracer));

  if (readfile == 0){
  for (int i=0; i<tracernum; i++) {
    tracer_init(i);
  }
  }
  else if (readfile ==1){
    char fname[30];
    FILE * fp;
    sprintf(fname, "%s/%s/tracer.csv", ROOT_DIR, INPUT_DIR);   
    fp = fopen (fname, "r");

    for (int i=0;i<tracernum;i++){
        fscanf(fp, "%d, %lf, %lf, %lf, %d", &tracer_array[i].tracerid, 
        &tracer_array[i].pos[0], &tracer_array[i].pos[1], &tracer_array[i].pos[2],
        &tracer_array[i].tracertype);
    }
    fclose(fp);
    printf("\nTracer locations imported successfully.\n");
  }

  struct stat st = {0};
  char buf[CHAR_BUF_SIZE];
  sprintf(buf, "%s/%s/%s", ROOT_DIR, OUTPUT_DIR, TRACER_DIR);
  if (stat(buf, &st) == -1) {
    mkdir(buf, 0700);
  }

  sprintf(buf, "%s/%s/%s", ROOT_DIR, OUTPUT_DIR, PARTICLEINFO_DIR);
  if (stat(buf, &st) == -1) {
    mkdir(buf, 0700);
  }

  // sprintf(buf, "%s/%s/%s", ROOT_DIR, OUTPUT_DIR, DEBUG_DIR);
  // if (stat(buf, &st) == -1) {
  //   mkdir(buf, 0700);
  // }

  FILE *rundata;
  char rundataname[100];
  sprintf(rundataname, "%s/%s/domaindata.csv",ROOT_DIR, OUTPUT_DIR);
  rundata = fopen(rundataname,"w+");
  fprintf(rundata,"%f,%f,%f,%f,%f,%f\n", dom[rank].xs, dom[rank].xe, dom[rank].ys,
          dom[rank].ye, dom[rank].zs, dom[rank].ze);
  fclose(rundata);

  exporter();
  printf("done.\n");

}

/*---------------------------------------------------------------------------*/
/* Remove tracer array from memory */
/*---------------------------------------------------------------------------*/

void tracer_exit(void)
{
  free(tracer_array);
}

/*---------------------------------------------------------------------------*/
/* Initialise each tracer*/
/*---------------------------------------------------------------------------*/

void tracer_init(int i)
{

  tracer_array[i].tracerid = i;

  // tracer_array[i].pos[0] = dom[rank].xs;// + dom[rank].xl*randgen()/2;// + dom[rank].xl*randgen()/6;
  tracer_array[i].pos[0] = dom[rank].xs + (i%2)*(dom[rank].xe - dom[rank].xs);
  //tracer_array[i].pos[0] = dom[rank].xs + 2*dom[rank].dx;
  // tracer_array[i].pos[0] = xlocation[i];

  // tracer_array[i].pos[1] = dom[rank].ys + (dom[rank].yl/2) + (i-9.5)*(dom[rank].yl/40);
  // tracer_array[i].pos[1] = dom[rank].ys;
  tracer_array[i].pos[1] = dom[rank].ys + dom[rank].yl*randgen()/2;// + dom[rank].yl*randgen()/8;
  // tracer_array[i].pos[1] = ylocation[i];

  tracer_array[i].pos[2] = dom[rank].zs + dom[rank].zl*randgen()/2;
  // tracer_array[i].pos[2] = dom[rank].zs;
  // tracer_array[i].pos[2] = zlocation[i];

  tracer_array[i].tracertype = i%2;
}

/*---------------------------------------------------------------------------*/
/* Delete a tracer from current use */
/*---------------------------------------------------------------------------*/

void tracer_delete(int deleteid)
{
  for (int d=deleteid; d<tracernum-1; d++) {
    tracer_array[d].tracerid = tracer_array[d+1].tracerid;
    tracer_array[d].tracertype = tracer_array[d+1].tracertype;

    for (int e=0; e<3; e++) {
      tracer_array[d].pos[e] = tracer_array[d+1].pos[e];
    }
  }

  tracernum -= 1;
}

/*---------------------------------------------------------------------------*/
/* Add new tracer to the simulation */
/*---------------------------------------------------------------------------*/

void tracer_add(void)
{
  tracer_init(tracernum);
  tracernum+=1;
}

/*---------------------------------------------------------------------------*/
/* Export tracer and particle data as csv */
/*---------------------------------------------------------------------------*/

void exporter(void)
{
  FILE *tracercsv;
  char filename[100];
  sprintf(filename, "%s/%s/%s/tracer-%d.csv", ROOT_DIR, OUTPUT_DIR, TRACER_DIR,
          exportcount);
  tracercsv = fopen(filename,"w+");
  for (int i=0; i<tracernum; i++) {
    fprintf(tracercsv, "%d, %f, %f, %f, %d\n", tracer_array[i].tracerid,
    tracer_array[i].pos[0], tracer_array[i].pos[1], tracer_array[i].pos[2],
    tracer_array[i].tracertype);
  }
  fclose(tracercsv);

  FILE *particleinfocsv;
  char particlename[100];
  sprintf(particlename, "%s/%s/%s/particle-%d.csv",ROOT_DIR, OUTPUT_DIR,
          PARTICLEINFO_DIR, exportcount);
  particleinfocsv = fopen(particlename,"w+");
  for (int i=0; i<NPARTS; i++) {
    fprintf(particleinfocsv, "%d, %f, %f, %f, %f\n", parts[i].N, parts[i].x,
            parts[i].y, parts[i].z, parts[i].r);
  }
  fclose(particleinfocsv);

  FILE *runinfo;
  char runname[100];
  sprintf(runname, "%s/%s/tracerinfo.csv", ROOT_DIR, OUTPUT_DIR);
  runinfo = fopen(runname,"a+");
  fprintf(runinfo, "%d, %f, %f, %d\n", exportcount, ttime, dt, nparts);
  fclose(runinfo);

}

/*---------------------------------------------------------------------------*/
/* Maintain number of particles close to inlet */
/*---------------------------------------------------------------------------*/
void maintainer(void)
{
  int tracercount = 0;
  // float safespace = dom[rank].xl/10;
  float safespace = dom[rank].dx;

  for (int a=0; a<tracernum; a++) {
    if (tracer_array[a].pos[0] < dom[rank].xs + safespace) {
      tracercount += 1;
    }
  }

  if (tracercount < fixednumber) {
    tracer_add();
  }
}

/*---------------------------------------------------------------------------*/
/* Keep on adding tracers at fixed rate */
/*---------------------------------------------------------------------------*/
void ratekeeper(void)
{
  if (ttime >= timecheck) {
    tracer_add();
    timecheck += fixedrate;
  }
}

/*---------------------------------------------------------------------------*/
/* Generate random number from -1 to 1 */
/*---------------------------------------------------------------------------*/

double randgen(void)
{
  double randomnum = rand();
  randomnum = fmod(randomnum,1e6);
  randomnum = randomnum/5e5;
  randomnum = randomnum - 1;
  return randomnum;
}

/*---------------------------------------------------------------------------*/
/* Generate the Brownian motion term */
/*---------------------------------------------------------------------------*/

double *randomizer(double *randvec, double diffusivity, double tstep)
{
  double magnitude = pow((6*diffusivity*tstep),0.5);
  double angle1, angle2;
  angle1 = PI * randgen();
  angle2 = PI * randgen();

  double randx = magnitude * cos(angle1) * sin(angle2);
  double randy = magnitude * sin(angle1) * sin(angle2);
  double randz = magnitude * cos(angle2);

  randvec[0] = randx;
  randvec[1] = randy;
  randvec[2] = randz;

  return randvec;
}

/*---------------------------------------------------------------------------*/
/* Calculate distance between two points */
/*---------------------------------------------------------------------------*/

double distance(double *loc, double x, double y, double z)
{
  dist_temp = pow((loc[0] - x),2) + pow((loc[1]-y),2) + pow((loc[2]-z),2);
  dist_temp = pow(dist_temp,0.5);
  return dist_temp;
}

/*---------------------------------------------------------------------------*/
/* Push tracer radially out of a particle before start */
/*---------------------------------------------------------------------------*/

double *sphere_pusher(double *loc, double x, double y, double z, double r)
{

  if (distance(loc, x, y, z) >= r) {
    return loc;
  }
  else {
    // printf("Tracer is inside a particle.\n");
    // printf("Location of particle = (%f, %f, %f).\n",x,y,z);
    // printf("Location of tracer before = (%f, %f, %f).\n",loc[0],loc[1],loc[2]);
    // printf("Distance = %f, r=%f.",distance(loc,x,y,z),r);

    // for (int i=0; i<3; i++) {
    //   printf("%f, ", loc[i]);
    // }

    double locp[3];

    locp[0] = x;
    locp[1] = y;
    locp[2] = z;

    double dn[3];
    double dnmag = distance(loc,x,y,z);

    for (int i=0; i<3; i++) {
      dn[i] = (loc[i] - locp[i])/dnmag;
    }

    double delta = 10;
    double locx_copy[3];

    for (int i=0;i<3;i++) {
      locx[i] = loc[i];
    }

    while (delta<=1E9){

      for (int i=0;i<3;i++) {
        locx_copy[i] = locx[i];
      }

      for (int i=0;i<3;i++) {
        locx[i] += dn[i]/delta;
      }

      if (distance(locx, x, y, z) >= r) {

        for (int i=0;i<3;i++) {
          locx[i] = locx_copy[i];
        }

        delta *= 10;
      }
    }

    for (int i=0;i<3;i++) {
      locx[i] += dn[i]/delta*10;// + dn[i]*(r-dnmag);
      //locx[i] += dn[i]*0.01*r;
    }

    // printf("Location of tracer after = (%f, %f, %f).\n",locx[0],locx[1],locx[2]);
    // printf("Distance = %f, r=%f.",distance(loc,x,y,z),r);
    if (distance(locx, x, y, z) < r) {
      printf("Tracer is still inside a particle.\n");
    }

    return locx;
  }
}

/*---------------------------------------------------------------------------*/
/* Check if tracer is inside a particle */
/*---------------------------------------------------------------------------*/

int particle_checker(double *loc) {
  int flag=0;
  for (int i=0;i<NPARTS;i++) {
    if (distance(loc,parts[i].x,parts[i].y, parts[i].z) < parts[i].r) {
      flag=1;
    }
  }
  return flag;
}

/*---------------------------------------------------------------------------*/
/* Check if tracer is close to a boundary */
/*---------------------------------------------------------------------------*/

int wall_checker(double *loc)
{
  double min_bound, max_bound;
  int count = 1;

  min_bound = dom[rank].xs + dom[rank].xl/(2*dom[rank].xn);
  max_bound = dom[rank].xe - dom[rank].xl/(2*dom[rank].xn);

  if (loc[0] <= min_bound || loc[0] >= max_bound) {
    count *= 2;
  }

  min_bound = dom[rank].ys + dom[rank].yl/(2*dom[rank].yn);
  max_bound = dom[rank].ye - dom[rank].yl/(2*dom[rank].yn);

  if (loc[1] <= min_bound || loc[1] >= max_bound) {
    count *= 3;
  }

  min_bound = dom[rank].zs + dom[rank].zl/(2*dom[rank].zn);
  max_bound = dom[rank].ze - dom[rank].zl/(2*dom[rank].zn);

  if (loc[2] <= min_bound || loc[2] >= max_bound) {
    count *= 5;
  }

  return count;
}

/*---------------------------------------------------------------------------*/
/* Check if tracer is out of a boundary */
/*---------------------------------------------------------------------------*/

int domain_checker(double *loc)
{
  double min_bound, max_bound;
  int count = 1;

  min_bound = dom[rank].xs;
  max_bound = dom[rank].xe;

  if (loc[0] < min_bound || loc[0] > max_bound) {
    count *= 2;
  }

  min_bound = dom[rank].ys;
  max_bound = dom[rank].ye;

  if (loc[1] < min_bound || loc[1] > max_bound) {
    count *= 3;
  }

  min_bound = dom[rank].zs;
  max_bound = dom[rank].ze;

  if (loc[2] < min_bound || loc[2] > max_bound) {
    count *= 5;
  }

  return count;
}

/*---------------------------------------------------------------------------*/
/* Move a tracer if it crosses a wall */
/*---------------------------------------------------------------------------*/

double *wall_reflector(int i, double *loc)
{
  double x = loc[0];
  double y = loc[1];
  double z = loc[2];

  if (x < dom[rank].xs) {
    // x = 2*dom[rank].xs - x;
    x = dom[rank].xe - (dom[rank].xs - x);
    // tracer_array[i].tracertype = 0;
    // printf("Tracer %d reflected xs.\n", i);
  }

  if (x > dom[rank].xe) {

    if (farwalltreatment==0){
    // x = 2*dom[rank].xe - x;
    x = dom[rank].xs + (x - dom[rank].xe);
    // tracer_array[i].tracertype = 1;
    // tracer_array[i].tracertype = 1;
    }
    if (farwalltreatment==1) {
      exporter();
      tracer_delete(i);
    }

    if (farwalltreatment==2) {
      tracer_array[i].tracertype = 1;
    }

    // printf("Tracer %d reflected xe.\n", i);
  }

  if (y < dom[rank].ys) {
    y = 2*dom[rank].ys - y;
    // y = dom[rank].ye - (dom[rank].ys - y);
    tracer_array[i].tracertype = 0;
    // printf("Tracer %d reflected ys.\n", i);
  }

  if (y > dom[rank].ye) {
    y = 2*dom[rank].ye - y;
    // y = dom[rank].ys + (y - dom[rank].ye);
    tracer_array[i].tracertype = 1;
    // printf("Tracer %d reflected ye.\n", i);
  }

  if (z < dom[rank].zs) {
    // z = 2*dom[rank].zs - z;
    z = dom[rank].ze - (dom[rank].zs - z);
    // printf("Tracer %d reflected zs.\n", i);
  }

  if (z > dom[rank].ze) {
    // z = 2*dom[rank].ze - z;
    z = dom[rank].zs + (z - dom[rank].ze);
    // printf("Tracer %d reflected ze.\n", i);
  }

  loc[0] = x;
  loc[1] = y;
  loc[2] = z;

  return loc;
}

/*---------------------------------------------------------------------------*/
/* Reflect a tracer on the surface of a particle */
/*---------------------------------------------------------------------------*/

double *reflector(double *loc1, double *loc2, double x,
  double y, double z, double r)
{
  if (distance(loc2, x, y, z) >= r) {
    return loc2;
  }

  else {

    double deltax = loc2[0] - loc1[0];
    double deltay = loc2[1] - loc1[1];
    double deltaz = loc2[2] - loc1[2];

    double delta = 10;

    double lucy[3];
    double lucy_copy[3];
    lucy[0] = loc1[0];
    lucy[1] = loc1[1];
    lucy[2] = loc1[2];

    while (delta<100000){

      lucy_copy[0] = lucy[0];
      lucy_copy[1] = lucy[1];
      lucy_copy[2] = lucy[2];

      lucy[0] += deltax/delta;
      lucy[1] += deltay/delta;
      lucy[2] += deltaz/delta;

      if (distance(lucy, x, y, z) <= r) {
        lucy[0] = lucy_copy[0];
        lucy[1] = lucy_copy[1];
        lucy[2] = lucy_copy[2];
        delta *= 10;
      }
    }

    double dsmag = distance(loc2, lucy[0],lucy[1],lucy[2]);
    double dnmag = distance(lucy, x, y, z);

    double locp[3];
    locp[0] = x;
    locp[1] = y;
    locp[2] = z;

    double di[3], dn[3], dotvec;
    dotvec = 0;

    for (int i=0;i<3;i++) {
      di[i] = (loc2[i] - lucy[i])/dsmag;
      dn[i] = (lucy[i] - locp[i])/dnmag;
      dotvec += di[i]*dn[i];
    }

    for (int i=0; i<3; i++) {
      locx[i] = (di[i] - 2*dotvec*dn[i])*dsmag + lucy[i];
    }

    return locx;
  }
}

double *intersector(double *loc1, double *loc2, double x,
  double y, double z, double r) {

  if (distance(loc2,x,y,z)>=r) {
    return loc1;
  }

  else {

  double deltax = loc2[0] - loc1[0];
  double deltay = loc2[1] - loc1[1];
  double deltaz = loc2[2] - loc1[2];

  double delta = 10;

  double lucy[3];
  double lucy_copy[3];
  lucy[0] = loc1[0];
  lucy[1] = loc1[1];
  lucy[2] = loc1[2];

  while (delta<100000){

    lucy_copy[0] = lucy[0];
    lucy_copy[1] = lucy[1];
    lucy_copy[2] = lucy[2];

    lucy[0] += deltax/delta;
    lucy[1] += deltay/delta;
    lucy[2] += deltaz/delta;

    if (distance(lucy, x, y, z) <= r) {
      lucy[0] = lucy_copy[0];
      lucy[1] = lucy_copy[1];
      lucy[2] = lucy_copy[2];
      delta *= 10;
    }
  }

  locx[0] = lucy[0];
  locx[1] = lucy[1];
  locx[2] = lucy[2];

  return locx;
  }
}

/*---------------------------------------------------------------------------*/
/* Find the 1D index of the cell-centered grid given tracer location */
/*---------------------------------------------------------------------------*/

double index_finder(double *loc) {
  double xi, yi, zi;
  xi = floor((loc[0] - dom[rank].xs - dom[rank].dx/2)/dom[rank].dx);
  yi = floor((loc[1] - dom[rank].ys - dom[rank].dy/2)/dom[rank].dy);
  zi = floor((loc[2] - dom[rank].zs - dom[rank].dz/2)/dom[rank].dz);
  if (xi<0) { xi =0; }
  if (yi<0) { yi =0; }
  if (zi<0) { zi =0; }
  foundindex = xi + yi * dom[rank].xn
  + zi * dom[rank].xn * dom[rank].yn;
  return foundindex;
}

/*---------------------------------------------------------------------------*/
/* Find the location of the cell-centered grid given the 1D index */
/*---------------------------------------------------------------------------*/

double grid_finder(double index, int axis, double gridloc) {

  double xxi = fmod(index , (double) dom[rank].xn);
  double yyi = fmod(((index - xxi)/dom[rank].xn) , (double) (dom[rank].yn));
  double zzi = ((index - xxi - yyi*dom[rank].xn))/(dom[rank].xn*dom[rank].yn);

  if (axis==0) {
    gridloc = (xxi+0.5)*dom[rank].dx + dom[rank].xs;
    return gridloc;
  }

  else if (axis==1) {
    gridloc = (yyi+0.5)*dom[rank].dy + dom[rank].ys;
    return gridloc;
  }

  else if (axis==2) {
    gridloc = (zzi+0.5)*dom[rank].dz + dom[rank].zs;
    return gridloc;
  }

  else {
    return 0;
  }
}

/*---------------------------------------------------------------------------*/
/* Find the velocity of the fluid at the position of the tracer */
/*---------------------------------------------------------------------------*/

double *interpolator(double *loc, double *ucc, double *vcc, double *wcc) {

  double min_bound, max_bound;
  double gridloc=0;

  switch (wall_checker(loc)) {
    case 1:
    foundindex = index_finder(loc);
    double c000, c001, c010, c011, c100, c101, c110, c111, cx, cy, cz;

    cx = (loc[0] - grid_finder(foundindex,0, gridloc))/dom[rank].dx;
    cy = (loc[1] - grid_finder(foundindex,1, gridloc))/dom[rank].dy;
    cz = (loc[2] - grid_finder(foundindex,2, gridloc))/dom[rank].dz;

    c000 = ucc[foundindex];
    c001 = ucc[foundindex + 1];
    c010 = ucc[foundindex + dom[rank].xn];
    c011 = ucc[foundindex + dom[rank].xn + 1];
    c100 = ucc[foundindex + dom[rank].xn*dom[rank].yn];
    c101 = ucc[foundindex + dom[rank].xn*dom[rank].yn + 1];
    c110 = ucc[foundindex + dom[rank].xn*dom[rank].yn + dom[rank].xn];
    c111 = ucc[foundindex + dom[rank].xn*dom[rank].yn + dom[rank].xn + 1];

    double tempx1 = (1-cx)*c000 + cx*c001;
    double tempx2 = (1-cx)*c010 + cx*c011;
    double tempx3 = (1-cx)*c100 + cx*c101;
    double tempx4 = (1-cx)*c110 + cx*c111;
    double tempy1 = (1-cy)*tempx1 + cy*tempx2;
    double tempy2 = (1-cy)*tempx3 + cy*tempx4;
    interpvel[0]  = (1-cz)*tempy1 + cz*tempy2;

    c000 = vcc[foundindex];
    c001 = vcc[foundindex + 1];
    c010 = vcc[foundindex + dom[rank].xn];
    c011 = vcc[foundindex + dom[rank].xn + 1];
    c100 = vcc[foundindex + dom[rank].xn*dom[rank].yn];
    c101 = vcc[foundindex + dom[rank].xn*dom[rank].yn + 1];
    c110 = vcc[foundindex + dom[rank].xn*dom[rank].yn + dom[rank].xn];
    c111 = vcc[foundindex + dom[rank].xn*dom[rank].yn + dom[rank].xn + 1];

    tempx1 = (1-cx)*c000 + cx*c001;
    tempx2 = (1-cx)*c010 + cx*c011;
    tempx3 = (1-cx)*c100 + cx*c101;
    tempx4 = (1-cx)*c110 + cx*c111;
    tempy1 = (1-cy)*tempx1 + cy*tempx2;
    tempy2 = (1-cy)*tempx3 + cy*tempx4;
    interpvel[1]  = (1-cz)*tempy1 + cz*tempy2;

    c000 = wcc[foundindex];
    c001 = wcc[foundindex + 1];
    c010 = wcc[foundindex + dom[rank].xn];
    c011 = wcc[foundindex + dom[rank].xn + 1];
    c100 = wcc[foundindex + dom[rank].xn*dom[rank].yn];
    c101 = wcc[foundindex + dom[rank].xn*dom[rank].yn + 1];
    c110 = wcc[foundindex + dom[rank].xn*dom[rank].yn + dom[rank].xn];
    c111 = wcc[foundindex + dom[rank].xn*dom[rank].yn + dom[rank].xn + 1];

    tempx1 = (1-cx)*c000 + cx*c001;
    tempx2 = (1-cx)*c010 + cx*c011;
    tempx3 = (1-cx)*c100 + cx*c101;
    tempx4 = (1-cx)*c110 + cx*c111;
    tempy1 = (1-cy)*tempx1 + cy*tempx2;
    tempy2 = (1-cy)*tempx3 + cy*tempx4;
    interpvel[2]  = (1-cz)*tempy1 + cz*tempy2;

    break;

    case 2:

    min_bound = dom[rank].xs + dom[rank].xl/(2*dom[rank].xn);
    max_bound = dom[rank].xe - dom[rank].xl/(2*dom[rank].xn);

    if (loc[0] <= min_bound) {
      foundindex = index_finder(loc);
      interpvel[0] = ucc[foundindex];
      interpvel[1] = vcc[foundindex];
      interpvel[2] = wcc[foundindex];
    }

    if(loc[0] >= max_bound) {
      foundindex = index_finder(loc);
      interpvel[0] = ucc[foundindex];
      interpvel[1] = vcc[foundindex];
      interpvel[2] = wcc[foundindex];
    }
    break;

    case 3:

    min_bound = dom[rank].ys + dom[rank].yl/(2*dom[rank].yn);
    max_bound = dom[rank].ye - dom[rank].yl/(2*dom[rank].yn);

    if (loc[0] <= min_bound) {
      foundindex = index_finder(loc);
      interpvel[0] = ucc[foundindex];
      interpvel[1] = vcc[foundindex];
      interpvel[2] = wcc[foundindex];
    }

    if(loc[0] >= max_bound) {
      foundindex = index_finder(loc);
      interpvel[0] = ucc[foundindex];
      interpvel[1] = vcc[foundindex];
      interpvel[2] = wcc[foundindex];
    }
    break;

    case 5:

    min_bound = dom[rank].zs + dom[rank].zl/(2*dom[rank].zn);
    max_bound = dom[rank].ze - dom[rank].zl/(2*dom[rank].zn);

    if (loc[0] <= min_bound) {
      foundindex = index_finder(loc);
      interpvel[0] = ucc[foundindex];
      interpvel[1] = vcc[foundindex];
      interpvel[2] = wcc[foundindex];
    }

    if(loc[0] >= max_bound) {
      foundindex = index_finder(loc);
      interpvel[0] = ucc[foundindex];
      interpvel[1] = vcc[foundindex];
      interpvel[2] = wcc[foundindex];
    }
    break;

    case 6:
    foundindex = index_finder(loc);
    interpvel[0] = ucc[foundindex];
    interpvel[1] = vcc[foundindex];
    interpvel[2] = wcc[foundindex];
    break;

    case 10:
    foundindex = index_finder(loc);
    interpvel[0] = ucc[foundindex];
    interpvel[1] = vcc[foundindex];
    interpvel[2] = wcc[foundindex];
    break;

    case 15:
    foundindex = index_finder(loc);
    interpvel[0] = ucc[foundindex];
    interpvel[1] = vcc[foundindex];
    interpvel[2] = wcc[foundindex];
    break;

    case 30:
    foundindex = index_finder(loc);
    interpvel[0] = ucc[foundindex];
    interpvel[1] = vcc[foundindex];
    interpvel[2] = wcc[foundindex];
    break;
  }

  return interpvel;
}

/*---------------------------------------------------------------------------*/
/* Main simulation */
/*---------------------------------------------------------------------------*/

void tracer_execute(double dt)
{
  cuda_dom_pull();
  cuda_part_pull();
  tracercheck++;

  double randvec[3];
  double *correction1;
  double *correction2;
  double *correction3;
  double *correction4;
  double *correction5;
  double *correction6;
  double correctionr[3];
  double loc_t[3];

  ucc = (real*) malloc(dom[rank].Gcc.s3 * sizeof(real));
  vcc = (real*) malloc(dom[rank].Gcc.s3 * sizeof(real));
  wcc = (real*) malloc(dom[rank].Gcc.s3 * sizeof(real));
  for (int k = dom[rank].Gcc._ks; k <= dom[rank].Gcc._ke; k++) {
    for (int j = dom[rank].Gcc._js; j <= dom[rank].Gcc._je; j++) {
      for (int i = dom[rank].Gcc._is; i <= dom[rank].Gcc._ie; i++) {
        int C = GCC_LOC(i - DOM_BUF, j - DOM_BUF, k - DOM_BUF,
                          dom[rank].Gcc.s1, dom[rank].Gcc.s2);

        int Cfx_w = GFX_LOC(i - 1, j, k, dom[rank].Gfx.s1b, dom[rank].Gfx.s2b);
        int Cfx = GFX_LOC(i, j, k, dom[rank].Gfx.s1b, dom[rank].Gfx.s2b);
        int Cfx_e = GFX_LOC(i + 1, j, k, dom[rank].Gfx.s1b, dom[rank].Gfx.s2b);
        int Cfx_ee = GFX_LOC(i + 2, j, k, dom[rank].Gfx.s1b, dom[rank].Gfx.s2b);

        int Cfy_s =  GFY_LOC(i, j - 1, k, dom[rank].Gfy.s1b, dom[rank].Gfy.s2b);
        int Cfy = GFY_LOC(i, j, k, dom[rank].Gfy.s1b, dom[rank].Gfy.s2b);
        int Cfy_n =  GFY_LOC(i, j + 1, k, dom[rank].Gfy.s1b, dom[rank].Gfy.s2b);
        int Cfy_nn = GFY_LOC(i, j + 2, k, dom[rank].Gfy.s1b, dom[rank].Gfy.s2b);

        int Cfz_b = GFZ_LOC(i, j, k - 1, dom[rank].Gfz.s1b, dom[rank].Gfz.s2b);
        int Cfz = GFZ_LOC(i, j, k, dom[rank].Gfz.s1b, dom[rank].Gfz.s2b);
        int Cfz_t = GFZ_LOC(i, j, k + 1, dom[rank].Gfz.s1b, dom[rank].Gfz.s2b);
        int Cfz_tt = GFZ_LOC(i, j, k + 2, dom[rank].Gfz.s1b, dom[rank].Gfz.s2b);

        ucc[C] = -0.0625*u[Cfx_w] + 0.5625*u[Cfx] + 0.5625*u[Cfx_e]
                    -0.0625*u[Cfx_ee];
        vcc[C] = -0.0625*v[Cfy_s] + 0.5625*v[Cfy] + 0.5625*v[Cfy_n]
                    -0.0625*v[Cfy_nn];
        wcc[C] = -0.0625*w[Cfz_b] + 0.5625*w[Cfz] + 0.5625*w[Cfz_t]
                    -0.0625*w[Cfz_tt];
      }
    }
  }

  // FILE *debugcsv;
  // char debugname[100];
  // sprintf(debugname, "%s/%s/%s/tracer-%d.csv", ROOT_DIR, OUTPUT_DIR, DEBUG_DIR,
  //         tracercheck);
  // debugcsv = fopen(debugname,"w+");

  for (int i=0; i<tracernum; i++) {

    if (farwalltreatment<2 || tracer_array[i].tracertype==0) {

      // fprintf(debugcsv, "%d, %f, %f, %f\n", i, tracer_array[i].pos[0],
      //      tracer_array[i].pos[1],tracer_array[i].pos[2]);


      // if (particle_checker(tracer_array[i].pos)==1) {
      //   printf("Tracer %d is inside a particle before pushing.\n",i);
      // }

      // 1. Sphere pusher
      while (particle_checker(tracer_array[i].pos)==1) {
        for (int j=0; j<NPARTS; j++) {
          correction1 = sphere_pusher(tracer_array[i].pos,parts[j].x,
            parts[j].y, parts[j].z, parts[j].r);

          for (int k=0; k<3; k++) {
            tracer_array[i].pos[k] = correction1[k];
          }
        }

        correction5 = wall_reflector(i,tracer_array[i].pos);

        for (int j=0; j<3; j++) {
          tracer_array[i].pos[j] = correction5[j];
        }

      }

      // if (domain_checker(tracer_array[i].pos)>1) {
      //   printf("Tracer %d is outside the domain after pushing.\nLocation of tracer = (%f,%f,%f).\n\n",
      //           i,tracer_array[i].pos[0],tracer_array[i].pos[1],tracer_array[i].pos[2]);
      // }

      // if (particle_checker(tracer_array[i].pos)==1) {
      //   printf("Tracer %d is inside a particle after pushing.\n",i);
      // }

      // fprintf(debugcsv, "%f, %f, %f\n", tracer_array[i].pos[0],
      //      tracer_array[i].pos[1],tracer_array[i].pos[2]);

     // 2. Randomizer
      correction2 = randomizer(randvec, diffusivity, dt);

      for (int j=0; j<3; j++) {
        loc_t[j] = tracer_array[i].pos[j] + correction2[j];
      }

      // fprintf(debugcsv, "%f, %f, %f\n", loc_t[0], loc_t[1],loc_t[2]);

      // 3. Interpolator
      correction3 = interpolator(tracer_array[i].pos, ucc, vcc, wcc);

      for (int j=0; j<3; j++) {
        loc_t[j] += correction3[j]*dt;
      }

      // fprintf(debugcsv, "%f, %f, %f\n", loc_t[0], loc_t[1],loc_t[2]);


      // if (particle_checker(loc_t)==1) {
      //   printf("Tracer %d is inside a particle after movement.\n",i);
      // }

      for (int k=0;k<3;k++) {
         correctionr[k] = tracer_array[i].pos[k];
      }

      // 4. Sphere reflector

      {

      int stepcount=0;

      while(particle_checker(loc_t)==1) {

        stepcount++;

        if (stepcount>50) {
          break;
        }

        for (int j=0; j<NPARTS; j++) {

          correction4 = reflector(correctionr, loc_t, parts[j].x, parts[j].y,
            parts[j].z, parts[j].r);

          correction6 = intersector(correctionr,loc_t, parts[j].x, parts[j].y,
            parts[j].z, parts[j].r);

          for (int k=0; k<3; k++) {
            correctionr[k] = correction6[k];
            loc_t[k] = correction4[k];
          }
        }
      }

      for (int j=0;j<3;j++) {
        tracer_array[i].pos[j] = loc_t[j];
      }

      // if (particle_checker(tracer_array[i].pos)==1) {
      //   printf("Tracer %d is inside a particle after sphere reflection.\n",i);
      // }

      }

      // fprintf(debugcsv, "%f, %f, %f\n", tracer_array[i].pos[0],
      //      tracer_array[i].pos[1],tracer_array[i].pos[2]);

      //5. Wall reflector
      correction5 = wall_reflector(i,tracer_array[i].pos);

      for (int j=0; j<3; j++) {
        tracer_array[i].pos[j] = correction5[j];
      }

      // fprintf(debugcsv, "%f, %f, %f\n\n", tracer_array[i].pos[0],
      //      tracer_array[i].pos[1],tracer_array[i].pos[2]);

      // if (particle_checker(tracer_array[i].pos)==1) {
      //   printf("Tracer %d is inside a particle after execution.\n",i);
      // }
    }
  }

  if (tracercheck%exportfreq==0 && ttime<600) {
    exportcount++;
    exporter();
    printf("Time=%.2f, n=%ld.\n",ttime,tracercheck);
  }

  if (tracercheck%10==0 && ttime>=600) {
    exportcount++;
    exporter();
    printf("Time=%.2f, n=%ld.\n",ttime,tracercheck);
  }

  // fclose(debugcsv);
  // exporter();

  free(ucc);
  free(vcc);
  free(wcc);
}
