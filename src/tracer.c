#include "bluebottle.h"
#include "tracer.h"

// double randomnum;
// double rnum;

tracer *tracer_array;
double diffusivity;
double tstep;
int foundindex = 0;
double *interpvel;
int tracernum;
double locx[3];
double dist;

void tracer_init(void)
{
  printf("\nInitializing tracer arrays... ");
  tracercheck = 0;
  diffusivity = 6E-9;

  // we have the number of particles as NPARTS
  // we have the particle sizes and location in parts[i].r,x,y,z

  // printf("Number of parts = %d\n", NPARTS);
  //
  // for (int i=0; i<NPARTS; i++) {
  //   printf("N=%d;\t(%f,%f,%f)\n",i+1, parts[i].x, parts[i].y, parts[i].z);
  // }

  srand((unsigned int)time(NULL));

  tracernum = 1000;
  tracer_array = (tracer*) malloc(tracernum * sizeof(tracer));

  for (int i=0; i<tracernum; i++) {
    tracer_array[i].tracerid = i;
    tracer_array[i].pos[0] = dom[rank].xs;
    tracer_array[i].pos[1] = dom[rank].ys + (dom[rank].yl/2)*(randgen()+1);
    tracer_array[i].pos[2] = dom[rank].zs + (dom[rank].zl/2)*(randgen()+1);
    tracer_array[i].tracertype = 0;
  }

  exporter();
  printf("done.\n");
}

void tracer_exit(void)
{
  free(tracer_array);
}

void exporter(void)
{
  FILE *tracercsv;
  char filename[100];
  sprintf(filename, "%s/%s/tracer-%d.csv", ROOT_DIR, OUTPUT_DIR, tracercheck);
  tracercsv = fopen(filename,"w+");
  for (int i=0; i<tracernum; i++) {
    fprintf(tracercsv, "%d, %f, %f, %f\n", tracer_array[i].tracerid,
    tracer_array[i].pos[0], tracer_array[i].pos[1], tracer_array[i].pos[2]);
  }
  // fprintf(tracercsv,"\n");
  fclose(tracercsv);

  FILE *particleinfocsv;
  char particlename[100];
  sprintf(particlename, "%s/%s/particleinfo-%d.csv",ROOT_DIR, OUTPUT_DIR,
          tracercheck);
  particleinfocsv = fopen(particlename,"w+");
  for (int i=0; i<NPARTS; i++) {
    fprintf(particleinfocsv, "%d, %f, %f, %f, %f\n", parts[i].N, parts[i].x,
            parts[i].y, parts[i].z, parts[i].r);
  }
  fclose(particleinfocsv);
}

double randgen(void)
{
  double randomnum = rand();
  randomnum = fmod(randomnum,1e6);
  randomnum = randomnum/5e5;
  randomnum = randomnum - 1;
  return randomnum;
}

double *randomizer(double *randvec, double diffusivity, double tstep)
{
  double magnitude = pow((6*diffusivity*tstep),0.5);
  // printf("%e",magnitude);
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

double distance(double *loc, double x, double y, double z)
{
  dist = pow((loc[0] - x),2) + pow((loc[1]-y),2) + pow((loc[2]-z),2);
  dist = pow(dist,0.5);
  return dist;
}

int sphere_checker(double *loc, double x, double y, double z, double r)
{
  if (distance(loc,x,y,z) <= r)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

double *sphere_pusher(double *loc, double x, double y, double z, double r)
{

  if (distance(loc, x, y, z) >= r) {
    // printf("Did not push.\n");
    return loc;
  }
  else {
    // printf("Starting to push. Original location is [");

    // for (int i=0; i<3; i++) {
    //   printf("%f, ", loc[i]);
    // }
    // printf("].\n Location of sphere is [%f, %f, %f] with r=%f.",x,y,z,r);
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

    // double locx[3];
    double locx_copy[3];

    for (int i=0;i<3;i++) {
      locx[i] = loc[i];
    }

    while (delta<10000){

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
  // printf("Finished pushing. Final location is [\n");
  // for (int i=0; i<3; i++) {
  //   printf("%f, ", locx[i]);
  // }
  // printf("].\n");
  return locx;
  }
}

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

double *wall_reflector(int i, double *loc)
{
  double x = loc[0];
  double y = loc[1];
  double z = loc[2];

  if (x < dom[rank].xs) {
    x = 2*dom[rank].xs - x;
    tracer_array[i].tracertype = 0;
  }

  if (x > dom[rank].xe) {
    x = 2*dom[rank].xe - x;
    tracer_array[i].tracertype = 1;
  }

  if (y < dom[rank].ys) {
    // y = 2*dom[rank].ys - y;
    y = dom[rank].ye - (dom[rank].ys - y);
  }

  if (y > dom[rank].ye) {
    // y = 2*dom[rank].ye - y;
    y = dom[rank].ys + (y - dom[rank].ye);
  }

  if (z < dom[rank].zs) {
    // z = 2*dom[rank].zs - z;
    z = dom[rank].ze - (dom[rank].zs - z);
  }

  if (z > dom[rank].ze) {
    // z = 2*dom[rank].ze - z;
    z = dom[rank].zs + (z - dom[rank].ze);
  }

  loc[0] = x;
  loc[1] = y;
  loc[2] = z;

  return loc;
}

double *reflector(double *loc1, double *loc2, double x,
  double y, double z, double r)
{
  if (distance(loc2, x, y, z) >= r) {
    // printf("Did not reflect.\n");
    return loc2;
  }
  else {
    // printf("Starting sphere reflection.\n");
    double deltax = loc2[0] - loc1[0];
    double deltay = loc2[1] - loc1[1];
    double deltaz = loc2[2] - loc1[2];

    double delta = 10;

    double locy[3];
    double locx_copy[3];
    locy[0] = loc1[0];
    locy[1] = loc1[1];
    locy[2] = loc1[2];

    // printf("So far so good.\n");

    while (delta<100000){

      locx_copy[0] = locy[0];
      locx_copy[1] = locy[1];
      locx_copy[2] = locy[2];

      locy[0] += deltax/delta;
      locy[1] += deltay/delta;
      locy[2] += deltaz/delta;

      if (distance(locy, x, y, z) <= r) {
        locy[0] = locx_copy[0];
        locy[1] = locx_copy[1];
        locy[2] = locx_copy[2];
        delta *= 10;
      }
    }

    // printf("Hey, it crossed the loop.\n");

    double dsmag = distance(loc2, locx[0],locx[1],locx[2]);
    double dnmag = distance(locx, x, y, z);

    double locp[3];
    locp[0] = x;
    locp[1] = y;
    locp[2] = z;

    double di[3], dn[3],dotvec;
    dotvec = 0;

    for (int i=0;i<3;i++) {
      di[i] = (loc2[i] - locx[i])/dsmag;
      dn[i] = (locx[i] - locp[i])/dnmag;
      dotvec += di[i]*dn[i];
    }

    for (int i=0; i<3; i++) {
      loc2[i] = di[i] - 2*dotvec*dn[i];
    }
    // printf("Finished sphere reflection.\n");
    return loc2;
  }
}

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

double grid_finder(double index, int axis, double gridloc) {

  double xxi = fmod(index , (double) dom[rank].xn);
  double yyi = fmod(((index - xxi)/dom[rank].xn) , (double) (dom[rank].yn));
  double zzi = ((index - xxi - yyi*dom[rank].xn))/(dom[rank].xn*dom[rank].yn);

  if (axis==0) {
    gridloc = (xxi+0.5)*dom[rank].dx;
    return gridloc;
  }

  else if (axis==1) {
    gridloc = (yyi+0.5)*dom[rank].dy;
    return gridloc;
  }

  else if (axis==2) {
    gridloc = (zzi+0.5)*dom[rank].dz;
    return gridloc;
  }

  else {
    return 0;
  }
}

double *interpolator(double *loc, double *uout, double *vout, double *wout) {

  double min_bound, max_bound;
  switch (wall_checker(loc)) {
    case 1:
    foundindex = index_finder(loc);
    double c000, c001, c010, c011, c100, c101, c110, c111, cx, cy, cz;
    double gridloc=0;

    cx = (loc[0] - grid_finder(foundindex,0, gridloc))/dom[rank].dx;
    cy = (loc[1] - grid_finder(foundindex,1, gridloc))/dom[rank].dy;
    cz = (loc[2] - grid_finder(foundindex,2, gridloc))/dom[rank].dz;

    c000 = uout[foundindex];
    c001 = uout[foundindex + 1];
    c010 = uout[foundindex + dom[rank].xn];
    c011 = uout[foundindex + dom[rank].xn + 1];
    c100 = uout[foundindex + dom[rank].xn*dom[rank].yn];
    c101 = uout[foundindex + dom[rank].xn*dom[rank].yn + 1];
    c110 = uout[foundindex + dom[rank].xn*dom[rank].yn + dom[rank].xn];
    c111 = uout[foundindex + dom[rank].xn*dom[rank].yn + dom[rank].xn + 1];

    double tempx1 = (1-cx)*c000 + cx*c001;
    double tempx2 = (1-cx)*c010 + cx*c011;
    double tempx3 = (1-cx)*c100 + cx*c101;
    double tempx4 = (1-cx)*c110 + cx*c111;
    double tempy1 = (1-cy)*tempx1 + cy*tempx2;
    double tempy2 = (1-cy)*tempx3 + cy*tempx4;
    interpvel[0]  = (1-cz)*tempy1 + cz*tempy2;

    c000 = vout[foundindex];
    c001 = vout[foundindex + 1];
    c010 = vout[foundindex + dom[rank].xn];
    c011 = vout[foundindex + dom[rank].xn + 1];
    c100 = vout[foundindex + dom[rank].xn*dom[rank].yn];
    c101 = vout[foundindex + dom[rank].xn*dom[rank].yn + 1];
    c110 = vout[foundindex + dom[rank].xn*dom[rank].yn + dom[rank].xn];
    c111 = vout[foundindex + dom[rank].xn*dom[rank].yn + dom[rank].xn + 1];

    tempx1 = (1-cx)*c000 + cx*c001;
    tempx2 = (1-cx)*c010 + cx*c011;
    tempx3 = (1-cx)*c100 + cx*c101;
    tempx4 = (1-cx)*c110 + cx*c111;
    tempy1 = (1-cy)*tempx1 + cy*tempx2;
    tempy2 = (1-cy)*tempx3 + cy*tempx4;
    interpvel[1]  = (1-cz)*tempy1 + cz*tempy2;

    c000 = wout[foundindex];
    c001 = wout[foundindex + 1];
    c010 = wout[foundindex + dom[rank].xn];
    c011 = wout[foundindex + dom[rank].xn + 1];
    c100 = wout[foundindex + dom[rank].xn*dom[rank].yn];
    c101 = wout[foundindex + dom[rank].xn*dom[rank].yn + 1];
    c110 = wout[foundindex + dom[rank].xn*dom[rank].yn + dom[rank].xn];
    c111 = wout[foundindex + dom[rank].xn*dom[rank].yn + dom[rank].xn + 1];

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
      interpvel[0] = uout[foundindex];
      interpvel[1] = vout[foundindex];
      interpvel[2] = wout[foundindex];
    }

    if(loc[0] >= max_bound) {
      foundindex = index_finder(loc);
      interpvel[0] = uout[foundindex];
      interpvel[1] = vout[foundindex];
      interpvel[2] = wout[foundindex];
    }
    break;

    case 3:

    min_bound = dom[rank].ys + dom[rank].yl/(2*dom[rank].yn);
    max_bound = dom[rank].ye - dom[rank].yl/(2*dom[rank].yn);

    if (loc[0] <= min_bound) {
      foundindex = index_finder(loc);
      interpvel[0] = uout[foundindex];
      interpvel[1] = vout[foundindex];
      interpvel[2] = wout[foundindex];
    }

    if(loc[0] >= max_bound) {
      foundindex = index_finder(loc);
      interpvel[0] = uout[foundindex];
      interpvel[1] = vout[foundindex];
      interpvel[2] = wout[foundindex];
    }
    break;

    case 5:

    min_bound = dom[rank].zs + dom[rank].zl/(2*dom[rank].zn);
    max_bound = dom[rank].ze - dom[rank].zl/(2*dom[rank].zn);

    if (loc[0] <= min_bound) {
      foundindex = index_finder(loc);
      interpvel[0] = uout[foundindex];
      interpvel[1] = vout[foundindex];
      interpvel[2] = wout[foundindex];
    }

    if(loc[0] >= max_bound) {
      foundindex = index_finder(loc);
      interpvel[0] = uout[foundindex];
      interpvel[1] = vout[foundindex];
      interpvel[2] = wout[foundindex];
    }
    break;

    case 6:
    foundindex = index_finder(loc);
    interpvel[0] = uout[foundindex];
    interpvel[1] = vout[foundindex];
    interpvel[2] = wout[foundindex];
    break;

    case 10:
    foundindex = index_finder(loc);
    interpvel[0] = uout[foundindex];
    interpvel[1] = vout[foundindex];
    interpvel[2] = wout[foundindex];
    break;

    case 15:
    foundindex = index_finder(loc);
    interpvel[0] = uout[foundindex];
    interpvel[1] = vout[foundindex];
    interpvel[2] = wout[foundindex];
    break;

    case 30:
    foundindex = index_finder(loc);
    interpvel[0] = uout[foundindex];
    interpvel[1] = vout[foundindex];
    interpvel[2] = wout[foundindex];
    break;
  }

  return interpvel;
}

void tracer_execute(double dt)
{
  cuda_dom_pull();
  tracercheck++;

  double randvec[3];
  double *correction;
  // double pushterm[3];
  // double wallterm[3];
  // double reflectterm[3];
  // double randterm[3];
  double loc_t[3];
  // int flag;
  // double *randterm;
  // randterm = randomizer(randvec, diffusivity, tstep);
  // printf("%e,%e,%e\n", randterm[0],randterm[1],randterm[2]);

  // printf("Hi.\n\n");

  real *ucc = malloc(dom[rank].Gcc.s3 * sizeof(real));
  real *vcc = malloc(dom[rank].Gcc.s3 * sizeof(real));
  real *wcc = malloc(dom[rank].Gcc.s3 * sizeof(real));
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

  // printf("Hi.\n\n");

  for (int i=0; i<tracernum; i++) {
    // flag = 0;

    for (int j=0; j<NPARTS; j++) {
      correction = sphere_pusher(tracer_array[i].pos,parts[j].x,
        parts[j].y, parts[j].z, parts[j].r);

      for (int k=0; k<3; k++) {

        // if (tracer_array[i].pos[k]!=correction[k]) { flag = 1;}
        // pushterm[k] = correction[k];
        tracer_array[i].pos[k] = correction[k];
      }
    }

    // printf("Pusher, \n");

    correction = randomizer(randvec, diffusivity, dt);

    for (int j=0; j<3; j++) {
      // randterm[j] = correction[j];
      loc_t[j] = tracer_array[i].pos[j] + correction[j];
    }

    // printf("randomizer, \n");

    for (int j=0; j<NPARTS; j++) {
      correction = reflector(tracer_array[i].pos, loc_t, parts[j].x, parts[j].y,
        parts[j].z, parts[j].r);

      for (int k=0; k<3; k++) {
        // reflectterm[k] = correction[k];
        tracer_array[i].pos[k] = correction[k];
      }
    }

    // printf("reflector, \n");

    correction = wall_reflector(i,tracer_array[i].pos);

    for (int k=0; k<3; k++) {
      // wallterm[k] = correction[k];
      tracer_array[i].pos[k] = correction[k];
    }

    // if (flag==1) { printf("Tracer %d has been pushed.\n",tracer_array[i].tracerid); }
      // printf(tracer_array[i].pos[0])
    // if (tracer_array[i].pos[0]>-5.5) {
    //   printf("Tracer %d\n", tracer_array[i].tracerid);
      // printf("Pushterm = [%f,%f,%f]\n Randterm=[%f,%f,%f]\n Reflectterm=[%f,%f,%f]\n Wallterm=[%f,%f,%f]\n\n",
      //       pushterm[0],pushterm[1],pushterm[2],
      //       randterm[0], randterm[1], randterm[2],
      //       reflectterm[0], reflectterm[1], reflectterm[2],
      //       wallterm[0], wallterm[1], wallterm[2]);
    // }

    // printf("wall reflector done for tracer %d.\n\n",(i+1));

  }

  exporter();

  free(ucc);
  free(vcc);
  free(wcc);

  // free(randterm);
}
