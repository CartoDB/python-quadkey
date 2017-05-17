
#include <Python.h>
#include <math.h>


#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
typedef unsigned long long uint64;
typedef unsigned int uint32;

static uint64 B[] = { 0x5555555555555555, 0x3333333333333333, 0x0F0F0F0F0F0F0F0F, 0x00FF00FF00FF00FF, 0x0000FFFF0000FFFF };
static uint64 S[] = { 1, 2, 4, 8, 16 };


uint64 xy2quadint(uint64 x, uint64 y) {

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

#define MAX_ZOOM 31

#define MAX_LONGITUDE 180.0
#define MAX_LATITUDE  85.05112877980659  /* (2*atan(exp(M_PI))*180.0/M_PI - 90.0) */

#define MIN_LONGITUDE (-MAX_LONGITUDE)
#define MIN_LATITUDE (-MAX_LATITUDE)

#define WEBMERCATOR_R 6378137.0

#define XY_SCALE 2147483648.0 /* (double)((uint32)1 << MAX_ZOOM) */
#define INV_XY_SCALE (1.0/XY_SCALE)
#define WM_RANGE (2.0*M_PI*WEBMERCATOR_R)
#define INV_WM_RANGE (1.0/WM_RANGE)

void lonlat2xy(double lon, double lat, int zoom, uint32* x, uint32* y) {

  lon = MIN(MAX_LONGITUDE, MAX(MIN_LONGITUDE, lon));
  lat = MIN(MAX_LATITUDE, MAX(MIN_LATITUDE, lat));

  double fx = (lon+180.0)/360.0;
  double sinlat = sin(lat * M_PI/180.0);
  double fy = 0.5 - log((1+sinlat)/(1-sinlat)) / (4*M_PI);

  uint32 mapsize = (1 << zoom);
  *x = (uint32)floor(fx*mapsize);
  *y = (uint32)floor(fy*mapsize);
  *x = MIN(mapsize - 1, MAX(0, *x));
  *y = MIN(mapsize - 1, MAX(0, *y));
}

void xy2webmercator(uint32 x, uint32 y, double* wm_x, double* wm_y) {
   *wm_x = (x*INV_XY_SCALE - 0.5)*WM_RANGE;
   *wm_y = (0.5 - y*INV_XY_SCALE)*WM_RANGE;
}

void webmercator2xy(double wm_x, double wm_y, uint32* x, uint32* y) {
   *x = (wm_x*INV_WM_RANGE + 0.5)*XY_SCALE;
   *y = (0.5 - wm_y*INV_WM_RANGE)*XY_SCALE;
}

uint64 lonlat2quadint(double lon, double lat) {
    uint32 x, y;
    lonlat2xy(lon, lat, MAX_ZOOM, &x, &y);
    return xy2quadint(x, y);
}

/*
 Input:  integer x y coordinates based on WebMercator in the range [0,2^31)
 Output: 62-bit quadkey value
*/
static PyObject*
xy2quadint_py(PyObject* self, PyObject* args)
{
    unsigned int x, y;

    if (!PyArg_ParseTuple(args, "ii", &x, &y))
        return NULL;

    return Py_BuildValue("K", xy2quadint(x, y));
}

/*
 Input:  integer x y coordinates based on WebMercator in the range [0,2^31)
 Output: web mercator coordinates (SRID 3857)
*/
static PyObject*
xy2webmercator_py(PyObject* self, PyObject* args)
{
    unsigned int x, y;
    double wm_x, wm_y;

    if (!PyArg_ParseTuple(args, "ii", &x, &y))
        return NULL;

    xy2webmercator(x, y, &wm_x, &wm_y);

    return Py_BuildValue("dd", wm_x, wm_y);
}

/*
 Input: web mercator coordinates (SRID 3857)
 Output:  integer x y coordinates based on WebMercator in the range [0,2^31)
*/
static PyObject*
webmercator2xy_py(PyObject* self, PyObject* args)
{
    unsigned int x, y;
    double wm_x, wm_y;

    if (!PyArg_ParseTuple(args, "dd", &wm_x, &wm_y))
        return NULL;

    webmercator2xy(wm_x, wm_y, &x, &y);

    return Py_BuildValue("ii", x, y);
}

/*
 Input: longitude, latitude in WGS84 (SRID 4326) coordinates in degrees
 Output: 62-bit quadkey value
*/
static PyObject*
lonlat2quadint_py(PyObject* self, PyObject* args)
{
    double lon, lat;

    if (!PyArg_ParseTuple(args, "dd", &lon, &lat))
        return NULL;

    return Py_BuildValue("K", lonlat2quadint(lon, lat));
}

/*
 Input:  web mercator coordinates (SRID 3857)
 Output: 62-bit quadkey value
*/
static PyObject*
webmercator2quadint_py(PyObject* self, PyObject* args)
{
    double wm_x, wm_y;
    unsigned int x, y;

    if (!PyArg_ParseTuple(args, "dd", &wm_x, &wm_y))
        return NULL;

    webmercator2xy(wm_x, wm_y, &x, &y);
    return Py_BuildValue("K", xy2quadint(x, y));
}


static PyMethodDef QuadkeyMethods[] =
{
     {"xy2quadint", xy2quadint_py, METH_VARARGS, "xy2quadint_py"},
     {"lonlat2quadint", lonlat2quadint_py, METH_VARARGS, "lonlat2quadint"},
     {"xy2webmercator", xy2webmercator_py, METH_VARARGS, "xy2webmercator"},
     {"webmercator2xy", webmercator2xy_py, METH_VARARGS, "webmercator2xy"},
     {"webmercator2quadint", webmercator2quadint_py, METH_VARARGS, "webmercator2quadint"},
     {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initquadkey(void)
{
     (void) Py_InitModule("quadkey", QuadkeyMethods);
}
