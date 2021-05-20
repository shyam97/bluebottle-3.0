#include "bluebottle.h"
#include "tracer.h"

// double randomnum;
// double rnum;

tracer *tracer_array;
double diffusivity;
double tstep;
double *extents;
int foundindex = 0;
double *interpvel;

void tracer_init(void)
{
  tracercheck = 7;
  diffusivity = 6E-9;

  // we have the domain in extents
  // we have the number of particles as NPARTS
  // we have the particle sizes and location in parts[i].r,x,y,z

  extents = (double*) malloc(6 * sizeof(double));
  extents[0] = dom[rank].xs;
  extents[1] = dom[rank].xe;
  extents[2] = dom[rank].ys;
  extents[3] = dom[rank].ye;
  extents[4] = dom[rank].zs;
  extents[5] = dom[rank].ze;

  printf("Number of parts = %d\n", NPARTS);

  for (int i=0; i<NPARTS; i++) {
    printf("N=%d;\t(%f,%f,%f)\n",i+1, parts[i].x, parts[i].y, parts[i].z);
  }

  srand((unsigned int)time(NULL));

  tracer_array = (tracer*) malloc(100 * sizeof(tracer));
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
  printf("%e",magnitude);
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
  double dist = pow((loc[0] - x),2) + pow((loc[1]-y),2) + pow((loc[2]-z),2);
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
    return loc;
  }
  else {
    double locp[3];

    locp[0] = x;
    locp[1] = y;
    locp[2] = z;

    double dn[3], dnmag = distance(loc,x,y,z);

    for (int i=0; i<3; i++) {
      dn[i] = (loc[i] - locp[i])/dnmag;
    }

    double delta = 10;

    double locx[3];
    double locx_copy[3];

    for (int i=0;i<3;i++) {
      locx[i] = loc[i];
    }

    while (delta<1e9){

      for (int i=0;i<3;i++) {
        locx_copy[i] = locx[i];
      }

      for (int i=0;i<3;i++) {
        locx[i] += dn[i];
      }

      if (distance(locx, x, y, z) <= r) {

        for (int i=0;i<3;i++) {
          locx[i] = locx_copy[i];
        }

        delta *= 10;
      }
    }
  }
  return loc;
}

int wall_checker(double *loc)
{
  double min_bound, max_bound;
  int count = 0;

  min_bound = dom[rank].xs + dom[rank].xl/dom[rank].xn;
  max_bound = dom[rank].xe - dom[rank].xl/dom[rank].xn;

  if (loc[0] <= min_bound || loc[0] >= max_bound) {
    count += 1;
  }

  min_bound = dom[rank].ys + dom[rank].yl/dom[rank].yn;
  max_bound = dom[rank].ye - dom[rank].yl/dom[rank].yn;

  if (loc[1] <= min_bound || loc[1] >= max_bound) {
    count += 1;
  }

  min_bound = dom[rank].zs + dom[rank].zl/dom[rank].zn;
  max_bound = dom[rank].ze - dom[rank].zl/dom[rank].zn;

  if (loc[2] <= min_bound || loc[2] >= max_bound) {
    count += 1;
  }

  return count;
}

double *wall_reflector(double *loc)
{
  double x = loc[0];
  double y = loc[1];
  double z = loc[2];

  if (x < dom[rank].xs) {
    x = 2*dom[rank].xs - x;
  }

  if (x > dom[rank].xe) {
    x = 2*dom[rank].xe - x;
  }

  if (y < dom[rank].ys) {
    y = 2*dom[rank].ys - y;
  }

  if (y > dom[rank].ye) {
    y = 2*dom[rank].ye - y;
  }

  if (z < dom[rank].zs) {
    z = 2*dom[rank].zs - z;
  }

  if (z > dom[rank].ze) {
    z = 2*dom[rank].ze - z;
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
    return loc2;
  }
  else {
    double deltax = loc2[0] - loc1[0];
    double deltay = loc2[1] - loc1[1];
    double deltaz = loc2[2] - loc1[2];

    double delta = 10;

    double locx[3];
    double locx_copy[3];
    locx[0] = loc1[0];
    locx[1] = loc1[1];
    locx[2] = loc1[2];

    while (delta<1e9){

      locx_copy[0] = locx[0];
      locx_copy[1] = locx[1];
      locx_copy[2] = locx[2];

      locx[0] += deltax/delta;
      locx[1] += deltay/delta;
      locx[2] += deltaz/delta;

      if (distance(locx, x, y, z) <= r) {
        locx[0] = locx_copy[0];
        locx[1] = locx_copy[1];
        locx[2] = locx_copy[2];
        delta *= 10;
      }
    }

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

    return loc2;

  }
}

double index_finder(double *loc) {
  double xi, yi, zi;
  xi = floor((loc[0] - dom[rank].xs - dom[rank].dx/2)/dom[rank].dx);
  yi = floor((loc[1] - dom[rank].ys - dom[rank].dy/2)/dom[rank].dy);
  zi = floor((loc[2] - dom[rank].zs - dom[rank].dz/2)/dom[rank].dz);
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

double *interpolate(double *loc, double *uout, double *vout, double *wout) {

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

  return interpvel;

}

void tracer_execute(double dt)
{
  tstep = dt;
  cuda_dom_pull();
  printf("Tracercheck=%d\n",tracercheck);
  // printf("%d,%d\n",dom[rank].Gcc.s1, dom[rank].Gcc.s2);
  tracercheck++;

  double randvec[3];
  double *randterm;
  randterm = randomizer(randvec, diffusivity, tstep);
  printf("%e,%e,%e\n", randterm[0],randterm[1],randterm[2]);


  real *uout = malloc(dom[rank].Gcc.s3 * sizeof(real));
  real *vout = malloc(dom[rank].Gcc.s3 * sizeof(real));
  real *wout = malloc(dom[rank].Gcc.s3 * sizeof(real));
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

        uout[C] = -0.0625*u[Cfx_w] + 0.5625*u[Cfx] + 0.5625*u[Cfx_e]
                    -0.0625*u[Cfx_ee];
        vout[C] = -0.0625*v[Cfy_s] + 0.5625*v[Cfy] + 0.5625*v[Cfy_n]
                    -0.0625*v[Cfy_nn];
        wout[C] = -0.0625*w[Cfz_b] + 0.5625*w[Cfz] + 0.5625*w[Cfz_t]
                    -0.0625*w[Cfz_tt];
      }
    }
  }




  free(uout);
  free(vout);
  free(wout);


  // free(randterm);
}
