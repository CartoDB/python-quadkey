
#include <Python.h>
#include <math.h>


#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
typedef unsigned long long uint64;
typedef unsigned int uint32;

uint64 xyz2quadint(uint64 x, uint64 y) {
    uint64 B[] = {0x5555555555555555, 0x3333333333333333, 0x0F0F0F0F0F0F0F0F, 0x00FF00FF00FF00FF, 0x0000FFFF0000FFFF};
    uint64 S[] = { 1, 2, 4, 8, 16 };

    x = (x | (x << S[4])) & B[4];
    y = (y | (y << S[4])) & B[4];

    x = (x | (x << S[3])) & B[3];
    y = (y | (y << S[3])) & B[3];

    x = (x | (x << S[2])) & B[2];
    y = (y | (y << S[2])) & B[2];

    x = (x | (x << S[1])) & B[1];
    y = (y | (y << S[1])) & B[1];

    x = (x | (x << S[0])) & B[0];
    y = (y | (y << S[0])) & B[0];

    return x | (y << 1);
}

void lonlat2xyz(double lon, double lat, int zoom, uint32* x, uint32* y) {

  lon = MIN(180.0, MAX(-180.0, lon));
  lat = MIN(85.05112878, MAX(-85.05112878, lat));

  double fx = (lon+180.0)/360.0;
  double sinlat = sin(lat * M_PI/180.0);
  double fy = 0.5 - log((1+sinlat)/(1-sinlat)) / (4*M_PI);

  uint32 mapsize = (1 << zoom);
  *x = (uint32)floor(fx*mapsize);
  *y = (uint32)floor(fy*mapsize);
  *x = MIN(mapsize - 1, MAX(0, *x));
  *y = MIN(mapsize - 1, MAX(0, *y));
}

uint64 lonlat2quadint(double lon, double lat) {
    uint32 x, y;
    lonlat2xyz(lon, lat, 31, &x, &y);
    return xyz2quadint(x, y);
}


static PyObject*
xyz2quadint_py(PyObject* self, PyObject* args)
{
    unsigned int x, y; 

    if (!PyArg_ParseTuple(args, "ii", &x, &y))
        return NULL;

    return Py_BuildValue("K", xyz2quadint(x, y));
}

static PyObject*
lonlat2quadint_py(PyObject* self, PyObject* args)
{
    double lon, lat; 

    if (!PyArg_ParseTuple(args, "dd", &lon, &lat))
        return NULL;

    return Py_BuildValue("K", lonlat2quadint(lon, lat));
}

static PyMethodDef QuadkeyMethods[] =
{
     {"xyz2quadint", xyz2quadint_py, METH_VARARGS, "xyz2quadint_py"},
     {"lonlat2quadint", lonlat2quadint_py, METH_VARARGS, "lonlat2quadint"},
     {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initquadkey(void)
{
     (void) Py_InitModule("quadkey", QuadkeyMethods);
}
