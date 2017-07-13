
#include <Python.h>
#include <math.h>


#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
typedef unsigned long long uint64;
typedef unsigned int uint32;

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
#define WM_MAX (M_PI*WEBMERCATOR_R)

static inline uint64
xy2quadint(uint64 x, uint64 y)
{

    static uint64 B[] = { 0x5555555555555555, 0x3333333333333333, 0x0F0F0F0F0F0F0F0F, 0x00FF00FF00FF00FF, 0x0000FFFF0000FFFF };
    static uint64 S[] = { 1, 2, 4, 8, 16 };

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

static inline void
quadint2xy(uint64 quadint, uint32* result_x, uint32* result_y)
{
    static const uint64 B[] = {
        0x5555555555555555, 0x3333333333333333, 0x0F0F0F0F0F0F0F0F, 0x00FF00FF00FF00FF, 0x0000FFFF0000FFFF,
        0x00000000FFFFFFFF
    };
    static const unsigned int S[] =
        { 0, 1, 2, 4, 8, 16 };

    uint64 x = quadint;
    uint64 y = quadint >> 1;

    x = (x | (x >> S[0])) & B[0];
    y = (y | (y >> S[0])) & B[0];

    x = (x | (x >> S[1])) & B[1];
    y = (y | (y >> S[1])) & B[1];

    x = (x | (x >> S[2])) & B[2];
    y = (y | (y >> S[2])) & B[2];

    x = (x | (x >> S[3])) & B[3];
    y = (y | (y >> S[3])) & B[3];

    x = (x | (x >> S[4])) & B[4];
    y = (y | (y >> S[4])) & B[4];

    x = (x | (x >> S[5])) & B[5];
    y = (y | (y >> S[5])) & B[5];

    *result_x = x;
    *result_y = y;
}

static inline
uint64 tile_prefix_mask(int zoom)
{
    static const uint64 masks[] = {
        0x0ULL,
        0x3000000000000000ULL,
        0x3c00000000000000ULL,
        0x3f00000000000000ULL,
        0x3fc0000000000000ULL,
        0x3ff0000000000000ULL,
        0x3ffc000000000000ULL,
        0x3fff000000000000ULL,
        0x3fffc00000000000ULL,
        0x3ffff00000000000ULL,
        0x3ffffc0000000000ULL,
        0x3fffff0000000000ULL,
        0x3fffffc000000000ULL,
        0x3ffffff000000000ULL,
        0x3ffffffc00000000ULL,
        0x3fffffff00000000ULL,
        0x3fffffffc0000000ULL,
        0x3ffffffff0000000ULL,
        0x3ffffffffc000000ULL,
        0x3fffffffff000000ULL,
        0x3fffffffffc00000ULL,
        0x3ffffffffff00000ULL,
        0x3ffffffffffc0000ULL,
        0x3fffffffffff0000ULL,
        0x3fffffffffffc000ULL,
        0x3ffffffffffff000ULL,
        0x3ffffffffffffc00ULL,
        0x3fffffffffffff00ULL,
        0x3fffffffffffffc0ULL,
        0x3ffffffffffffff0ULL,
        0x3ffffffffffffffcULL,
        0x3fffffffffffffffULL
    };
    if (zoom < 0)
      zoom = 0;
    else if (zoom > MAX_ZOOM)
      zoom = MAX_ZOOM;
    return masks[zoom];
}

static inline uint64
tile_suffix_mask(int zoom)
{
    static const uint64 masks[] = {
        0x3fffffffffffffffULL,
        0xfffffffffffffffULL,
        0x3ffffffffffffffULL,
        0xffffffffffffffULL,
        0x3fffffffffffffULL,
        0xfffffffffffffULL,
        0x3ffffffffffffULL,
        0xffffffffffffULL,
        0x3fffffffffffULL,
        0xfffffffffffULL,
        0x3ffffffffffULL,
        0xffffffffffULL,
        0x3fffffffffULL,
        0xfffffffffULL,
        0x3ffffffffULL,
        0xffffffffULL,
        0x3fffffffULL,
        0xfffffffULL,
        0x3ffffffULL,
        0xffffffULL,
        0x3fffffULL,
        0xfffffULL,
        0x3ffffULL,
        0xffffULL,
        0x3fffULL,
        0xfffULL,
        0x3ffULL,
        0xffULL,
        0x3fULL,
        0xfULL,
        0x3ULL,
        0x0ULL
    };
    if (zoom < 0)
      zoom = 0;
    else if (zoom > MAX_ZOOM)
      zoom = MAX_ZOOM;
    return masks[zoom];
}

static void
lonlat2xy(double lon, double lat, int zoom, uint32* x, uint32* y)
{

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

static void
xy2webmercator(uint32 x, uint32 y, double* wm_x, double* wm_y)
{
   *wm_x = (x*INV_XY_SCALE - 0.5)*WM_RANGE;
   *wm_y = (0.5 - y*INV_XY_SCALE)*WM_RANGE;
}

static void
webmercator2xy(double wm_x, double wm_y, uint32* x, uint32* y)
{
   *x = (wm_x*INV_WM_RANGE + 0.5)*XY_SCALE;
   *y = (0.5 - wm_y*INV_WM_RANGE)*XY_SCALE;
}

static uint64
lonlat2quadint(double lon, double lat)
{
    uint32 x, y;
    lonlat2xy(lon, lat, MAX_ZOOM, &x, &y);
    return xy2quadint(x, y);
}

static void
tile2bbox_scaled(double scale_x, double scale_y, double offset_x, double offset_y, uint64 quadint, int zoom, double* x_min, double* y_min, double* x_max, double* y_max)
{
  unsigned int x, y;
  int zero_bits = MAX_ZOOM - zoom;
  quadint2xy(quadint, &x, &y);
  x >>= zero_bits;
  y >>= zero_bits;
  *x_min = offset_x + (x * 1.0 / (1ull << zoom)) * scale_x;
  *x_max = offset_x + ((x + 1) * 1.0 / (1ull << zoom)) * scale_x;
  *y_min = offset_y + ((y + 1) * 1.0 / (1ull << zoom)) * scale_y;
  *y_max = offset_y + (y * 1.0 / (1ull << zoom)) * scale_y;
}

static void
tile2bbox_webmercator(uint64 quadint, int zoom, double* x_min, double* y_min, double* x_max, double* y_max)
{
  tile2bbox_scaled(WM_RANGE, -WM_RANGE, -WM_MAX, WM_MAX, quadint, zoom, x_min, y_min, x_max, y_max);
}

static void
tile2bbox(uint64 quadint, int zoom, double* lon_min, double* lat_min, double* lon_max, double* lat_max)
{
  double x_min, y_min, x_max, y_max;
  tile2bbox_scaled(1.0, 1.0, -0.5, -0.5, quadint, zoom, &x_min, &y_min, &x_max, &y_max);

  *lon_min = 360.0 * x_min;
  *lon_max = 360.0 * x_max;

  *lat_min = 90.0 - 360.0*atan(exp(-2 * M_PI * (-y_min))) / M_PI;
  *lat_max = 90.0 - 360.0*atan(exp(-2 * M_PI * (-y_max))) / M_PI;
}

static void
tile_center_scaled(double scale_x, double scale_y, double offset_x, double offset_y, uint64 quadint, int zoom, double* cx, double* cy)
{
  unsigned int x, y;
  int zero_bits = MAX_ZOOM - zoom;
  quadint2xy(quadint, &x, &y);
  x >>= zero_bits;
  y >>= zero_bits;

  *cx = offset_x + ((x + x + 1) * 1.0 / (1ull << (zoom + 1))) * scale_x;
  *cy = offset_y + ((y + y + 1) * 1.0 / (1ull << (zoom + 1))) * scale_y;
}

static void
tile_center_webmercator(uint64 quadint, int zoom, double* cx, double* cy)
{
  tile_center_scaled(WM_RANGE, -WM_RANGE, -WM_MAX, WM_MAX, quadint, zoom, cx, cy);
}

static void
tile_center(uint64 quadint, int zoom, double* lon, double* lat)
{
  double x, y;
  tile_center_scaled(1.0, 1.0, -0.5, -0.5, quadint, zoom, &x, &y);

  *lon = 360.0 * x;
  *lat = 90.0 - 360.0*atan(exp(-2 * M_PI * (-y))) / M_PI;
}

static void
tile2range(uint64 quadint, int zoom, uint64* q_min, uint64* q_max)
{
    *q_min = quadint & tile_prefix_mask(zoom);
    *q_max = quadint | tile_suffix_mask(zoom);
}

static void
tile_children(uint64 quadint, int zoom, uint64* q_sw, uint64* q_nw, uint64* q_se, uint64* q_ne)
{
    int bit = ((MAX_ZOOM - zoom) << 1);
    *q_sw = quadint & tile_prefix_mask(zoom);
    *q_nw = *q_sw | (1ull << (bit - 2));
    *q_se = *q_sw | (1ull << (bit - 1));
    *q_ne = *q_se | (1ull << (bit - 2));
}

static uint64
xyz2quadint(uint32 x, uint32 y, int zoom)
{
    int bits = MAX_ZOOM - zoom;
    return xy2quadint(x << bits, y << bits);
}

static void
tile2xy(uint64 quadint, int zoom, uint32* x, uint32* y)
{
  uint32 qx, qy;
  int bits = MAX_ZOOM - zoom;
  quadint2xy(quadint, &qx, &qy);
  *x = qx >> bits;
  *y = qy >> bits;
}

static double
box_area(double xmin, double ymin, double xmax, double ymax)
{
    double w = xmax - xmin;
    double h = ymax - ymin;
    return (w > 0) && (h > 0) ? w*h : 0.0;
}

static void
parse_box(PyObject* box, double* xmin, double* ymin, double* xmax, double* ymax)
{
    if (PyObject_TypeCheck(box, &PyList_Type)) {
        *xmin = PyFloat_AsDouble(PyList_GetItem(box, 0));
        *ymin = PyFloat_AsDouble(PyList_GetItem(box, 1));
        *xmax = PyFloat_AsDouble(PyList_GetItem(box, 2));
        *ymax = PyFloat_AsDouble(PyList_GetItem(box, 3));
    }
    else {
        PyArg_ParseTuple(box, "dddd", xmin, ymin, xmax, ymax);
    }
}

static inline PyObject*
build_box(double xmin, double ymin, double xmax, double ymax)
{
    return Py_BuildValue("dddd", xmin, ymin, xmax, ymax);
}

static double
disjoint_boxlist_area(PyObject* boxlist)
{
   /* boxlist rectangles must be non-overlapping */
   Py_ssize_t i, nboxes = PyList_Size(boxlist);
   PyObject *box;
   double area = 0.0;
   double xmin, ymin, xmax, ymax;
   for (i = 0; i < nboxes; i++) {
       box = PyList_GetItem(boxlist, i);
       parse_box(box, &xmin, &ymin, &xmax, &ymax);
       area += box_area(xmin, ymin, xmax, ymax);
   }
   return area;
}

static void
split_boxes(PyObject* box1, PyObject* box2, PyObject* output_list1, PyObject* output_list2)
{
    double xmin1, ymin1, xmax1, ymax1, xmin2, ymin2, xmax2, ymax2;
    double clip_xmin, clip_ymin, clip_xmax, clip_ymax;
    parse_box(box1, &xmin1, &ymin1, &xmax1, &ymax1);
    parse_box(box2, &xmin2, &ymin2, &xmax2, &ymax2);
    if (xmin2 <= xmax1 && xmax2 >= xmin1 && ymin2 <= ymax1 && ymax2 >= ymin1) {
        /* Intersecting rectangles */
        if (xmax1 >= xmax2 && xmin1 <= xmin2 && ymax1 >= ymax2 && ymin1 <= ymin2) {
            /* box1 contains box2 */
            PyList_Append(output_list1, box1);
        } else if (xmax2 >= xmax1 && xmin2 <= xmin1 && ymax2 >= ymax1 && ymin2 <= ymin1) {
            /* box2 contains box1 */
            PyList_Append(output_list2, box2);
        } else {
            /* will remve box2 from box1, getting two or four rectangles */
            /* TODO: by sorting coordinates we can simplify this */
            if (xmax1 >= xmax2 && xmin1 <= xmin2) {
                /* substract box2 - box1 */
                clip_xmin = xmin2;
                clip_xmax = xmax2;
                if (ymax2 > ymax1) {
                    clip_ymax = ymax2;
                    clip_ymin = ymax1;
                } else {
                    clip_ymax = ymin1;
                    clip_ymin = ymin2;
                }
                PyList_Append(output_list2, build_box(clip_xmin, clip_ymin, clip_xmax, clip_ymax));
                PyList_Append(output_list1, box1);
            } else if (xmax2 >= xmax1 && xmin2 <= xmin1) {
                /* substract box1 - box2 */
                clip_xmin = xmin1;
                clip_xmax = xmax1;
                if (ymax1 > ymax2) {
                    clip_ymax = ymax1;
                    clip_ymin = ymax2;
                } else {
                    clip_ymax = ymin2;
                    clip_ymin = ymin1;
                }
                PyList_Append(output_list1, build_box(clip_xmin, clip_ymin, clip_xmax, clip_ymax));
                PyList_Append(output_list2, box2);
            } else if (ymax1 >= ymax2 && ymin1 <= ymin2) {
                /* substract box2 - box1 */
                clip_ymin = ymin2;
                clip_ymax = ymax2;
                if (xmax2 > xmax1) {
                    clip_xmax = xmax2;
                    clip_xmin = xmax1;
                } else {
                    clip_xmax = xmin1;
                    clip_xmin = xmin2;
                }
                PyList_Append(output_list2, build_box(clip_xmin, clip_ymin, clip_xmax, clip_ymax));
                PyList_Append(output_list1, box1);
            } else if (ymax2 >= ymax1 && ymin2 <= ymin1) {
                /* substract box1 - box2 */
                clip_ymin = ymin1;
                clip_ymax = ymax1;
                if (xmax1 > xmax2) {
                    clip_xmax = xmax1;
                    clip_xmin = xmax2;
                } else {
                    clip_xmax = xmin2;
                    clip_xmin = xmin1;
                }
                PyList_Append(output_list1, build_box(clip_xmin, clip_ymin, clip_xmax, clip_ymax));
                PyList_Append(output_list2, box2);
            } else {
                if (xmin1 < xmin2 && ymin1 < ymin2) {
                   PyList_Append(output_list1, build_box(xmin1, ymin1, xmin2, ymin2));
                   PyList_Append(output_list1, build_box(xmin2, ymin1, xmax1, ymin2));
                   PyList_Append(output_list1, build_box(xmin1, ymin2, xmin2, ymax1));
                   PyList_Append(output_list2, box2);
                } else if (xmin1 < xmin2 && ymax1 > ymax2) {
                   PyList_Append(output_list1, build_box(xmin1, ymin1, xmin2, ymax2));
                   PyList_Append(output_list1, build_box(xmin1, ymax2, xmin2, ymax1));
                   PyList_Append(output_list1, build_box(xmin2, ymax2, xmax1, ymax1));
                   PyList_Append(output_list2, box2);
                } else if (xmax1 > xmax2 && ymax1 > ymax2) {
                   PyList_Append(output_list1, build_box(xmin1, ymax2, xmax2, ymax1));
                   PyList_Append(output_list1, build_box(xmax2, ymax2, xmax1, ymax1));
                   PyList_Append(output_list1, build_box(xmax2, ymin1, xmax1, ymax2));
                   PyList_Append(output_list2, box2);
                } else { // (xmax1 > xmax2 && ymin1 < ymin2)
                   PyList_Append(output_list1, build_box(xmin1, ymin1, xmax2, ymin2));
                   PyList_Append(output_list1, build_box(xmax2, ymin1, xmax1, ymin2));
                   PyList_Append(output_list1, build_box(xmax2, ymin2, xmax1, ymax1));
                   PyList_Append(output_list2, box2);
                }
            }
        }
    } else {
        /* Non-intersecting rectangles */
        PyList_Append(output_list1, box1);
        PyList_Append(output_list2, box2);
    }
}

static double
box_intersection_area(double xmin1, double ymin1, double xmax1, double ymax1, double xmin2, double ymin2, double xmax2, double ymax2)
{
    return box_area(MAX(xmin1, xmin2), MAX(ymin1, ymin2), MIN(xmax1, xmax2), MIN(ymax1, ymax2));
}

static double
disjoint_boxlist_intersection_area(double xmin, double ymin, double xmax, double ymax, PyObject* boxlist)
{
    // TODO: check if this impacts the performance of adaptive_tiling => convert once to array and use arrays instead of Lists
    Py_ssize_t i, nboxes = PyList_Size(boxlist);
    PyObject *intersections = PyList_New(nboxes);
    PyObject *box;
    double area;
    double box_xmin, box_ymin, box_xmax, box_ymax;
    for (i = 0; i < nboxes; i++) {
        box = PyList_GetItem(boxlist, i);
        parse_box(box, &box_xmin, &box_ymin, &box_xmax, &box_ymax);
        box_xmin = MAX(box_xmin, xmin);
        box_ymin = MAX(box_ymin, ymin);
        box_xmax = MIN(box_xmax, xmax);
        box_ymax = MIN(box_ymax, ymax);
        box = build_box(box_xmin, box_ymin, box_xmax, box_ymax);
        PyList_SetItem(intersections, i, box);
    }
    area = disjoint_boxlist_area(intersections);
    Py_DECREF(intersections);
    return area;
}

static PyObject*
make_disjoint_boxes(PyObject* boxes)
{
    // This is O(n^2)
    // could be O(n*log(n)), see https://stackoverflow.com/questions/244452/what-is-an-efficient-algorithm-to-find-area-of-overlapping-rectangles
    // Also this is buggy now: some boxes go away
    Py_ssize_t i, nboxes = PyList_Size(boxes);
    PyObject* result = PyList_New(0);
    PyObject* alist = PyList_New(0);
    PyObject* blist = PyList_New(0);
    Py_ssize_t j, k;

    PyList_Append(result, PyList_GetItem(boxes, 0));

    for (i = 1; i < nboxes; i++) {
        PyList_Append(blist, PyList_GetItem(boxes, i));
        k = 0;
        for (j = 0; j < PyList_Size(result); j++) {
            if (k >= PyList_Size(blist))
                break;
            split_boxes(PyList_GetItem(result, j), PyList_GetItem(blist, k++),alist,blist);
        }
        // result = alist + blist[k:len(blist)]
        PyObject* tmp = result;
        result = alist;
        PyList_SetSlice(result, PyList_Size(result), PyList_Size(result), PyList_GetSlice(blist, k, PyList_Size(blist)));
        PyList_SetSlice(result, PyList_Size(result), PyList_Size(result), PyList_GetSlice(tmp, j, PyList_Size(tmp)));
        alist = tmp;
        PyList_SetSlice(alist, 0, PyList_Size(alist), NULL);
        PyList_SetSlice(blist, 0, PyList_Size(blist), NULL);
    }
    Py_DECREF(alist);
    Py_DECREF(blist);

    return result;
}

static void
append_tile(PyObject* list, uint64 quadint, int zoom)
{
    PyList_Append(list, Py_BuildValue("Ki", quadint, zoom));
}

#define CIRCULAR_INDEX(base, offset, n) ((base + offset) % n)

#define ADAPTIVE_TILING_BUFFER_SIZE 10000

/*
 * Approximate the area of a bounding box with tiles while limiting the error (as an area ratio).
 * When excess tolerance is given, max_error is the maximum relative error of defect area and
 * excess_tolerance is the maximum error of excess area.
 * For intersecting (tile-covering) use small max_error (so that all area is effectively covered), largish excess_tolerance
 * (e.g. 0.5 would include tiles which have at least 50% intersected by the area; but note that 0.75 would include all tiles with 25% of its area)
 * (1.0 would always return the root tile since the error for this tile would always be )
 * For approximating an area with tiles (tiling) use only a max_error and excess_error = 0.0.
*/
static PyObject*
adaptive_tiling(PyObject* boxes, double max_error, double excess_error)
{
    /* boxes have to be disjoint: TODO: splitting boxes to make them disjoint */
    static uint64 buffer[ADAPTIVE_TILING_BUFFER_SIZE]; /* circular buffer */

    PyObject* tiles = PyList_New(0);
    uint64 tile_quadint;
    double tile_xmin, tile_ymin, tile_xmax, tile_ymax;
    uint64 q_sw, q_nw, q_se, q_ne;

    double area, total_area = 0.0, covered_area = 0.0, int_area, tile_area, err;
    double tile_err;
    int zoom;
    int within_tolerance;

    tile_err =  (excess_error > 0.0) ? excess_error : max_error;

    area = disjoint_boxlist_area(boxes);

    size_t candidate_tiles; /* position of candidates in the buffer */
    size_t next_candidates; /* position of to add next zoom level candidates */
    size_t candidates_end;

    /* Start with the root tile as the only candidate */
    zoom = 0;
    candidate_tiles = 0;
    buffer[candidate_tiles] = 0ull;
    next_candidates = 1;

    for (;;) {
        candidates_end = next_candidates;
        while (candidate_tiles != candidates_end) {
            /* pop candidate from circular buffer */
            tile_quadint = buffer[candidate_tiles];
            candidate_tiles = CIRCULAR_INDEX(candidate_tiles, 1, ADAPTIVE_TILING_BUFFER_SIZE);

            tile2bbox_webmercator(tile_quadint, zoom, &tile_xmin, &tile_ymin, &tile_xmax, &tile_ymax);
            int_area = disjoint_boxlist_intersection_area(tile_xmin, tile_ymin, tile_xmax, tile_ymax, boxes);
            tile_area = box_area(tile_xmin, tile_ymin, tile_xmax, tile_ymax);
            if (fabs(int_area - tile_area) < tile_err*tile_area) {
                /* add the candidate tile to the results */
                total_area += tile_area;
                covered_area += int_area;
                append_tile(tiles, tile_quadint, zoom);
            } else if (int_area > 0 && zoom < MAX_ZOOM) {
                /* schedule the tile children as next-level candidates  */
                tile_children(tile_quadint, zoom, &q_sw, &q_nw, &q_se, &q_ne);
                buffer[next_candidates] = q_sw;
                next_candidates = CIRCULAR_INDEX(next_candidates, 1, ADAPTIVE_TILING_BUFFER_SIZE);
                if (next_candidates == candidate_tiles)
                    break; // buffer full; TODO: report error properly
                buffer[next_candidates] = q_nw;
                next_candidates = CIRCULAR_INDEX(next_candidates, 1, ADAPTIVE_TILING_BUFFER_SIZE);
                if (next_candidates == candidate_tiles) {
                    break; // buffer full; TODO: report error properly
                }
                buffer[next_candidates] = q_se;
                next_candidates = CIRCULAR_INDEX(next_candidates, 1, ADAPTIVE_TILING_BUFFER_SIZE);
                if (next_candidates == candidate_tiles)
                    break; // buffer full; TODO: report error properly
                buffer[next_candidates] = q_ne;
                next_candidates = CIRCULAR_INDEX(next_candidates, 1, ADAPTIVE_TILING_BUFFER_SIZE);
                if (next_candidates == candidate_tiles)
                    break; // buffer full; TODO: report error properly
            }
        }
        /* Check finishing conditions */
        if (excess_error > 0.0) {
            // the uncovered area should be within the tolerance (max_error)
            err = fabs(covered_area - area);
        } else {
            // the are difference between the tiles and area should be withing tolerance
            err = fabs(total_area - area);
        }
        within_tolerance = err < max_error*area;
        if (covered_area >= area || within_tolerance ||  zoom >= MAX_ZOOM || candidate_tiles == next_candidates) {
            break;
        }
        /* Iterate to the nexe zoom level */
        zoom += 1;
    };
    return tiles;
}

static void
tiles_intersecting_webmercator_boxes(PyObject* result, PyObject* boxes, uint64 quadint, int zoom, int max_zoom, int mode)
{
    double tile_xmin, tile_ymin, tile_xmax, tile_ymax;
    double int_area;
    double tile_area;
    double area_tol = 42.0 / (1ull << (max_zoom + 1));
    int    significant_intersection = 0;
    tile2bbox_webmercator(quadint, zoom, &tile_xmin, &tile_ymin, &tile_xmax, &tile_ymax);

    int_area = disjoint_boxlist_intersection_area(tile_xmin, tile_ymin, tile_xmax, tile_ymax, boxes);
    if (int_area > 0) {
        tile_area = box_area(tile_xmin, tile_ymin, tile_xmax, tile_ymax);
        if (int_area - tile_area >= -area_tol) {
            /* the box contains the tile; add the tile to the results */
            append_tile(result, quadint, zoom);
        } else if (zoom == max_zoom) {
            if (mode == 1) {
                /* intended for filtering data out of the bbox */
                /* Nomenclature: tile-covering; i.e. covering the area completely with tiles */
                significant_intersection = int_area >= area_tol;
            } else if (mode == 2)
            {
                /* intended for approximating the bbox with tiles; if max_zoom is too low results may be empty */
                /* Nomenclature: tiling; i.e. approximating the area with tiles; small area fraction could be left uncovered */
                significant_intersection = int_area >= 2*tile_area;
            }
            if (significant_intersection) {
                /* intersection is significant add the tile to the results */
                append_tile(result, quadint, zoom);
            } else {
                /* dismiss intersection */
            }
        } else {
            /* Drill down to next level */
            uint64 child_sw, child_nw, child_se, child_ne;
            tile_children(quadint, zoom, &child_sw, &child_nw, &child_se, &child_ne);
            tiles_intersecting_webmercator_boxes(result, boxes, child_sw, zoom + 1, max_zoom, mode);
            tiles_intersecting_webmercator_boxes(result, boxes, child_nw, zoom + 1, max_zoom, mode);
            tiles_intersecting_webmercator_boxes(result, boxes, child_se, zoom + 1, max_zoom, mode);
            tiles_intersecting_webmercator_boxes(result, boxes, child_ne, zoom + 1, max_zoom, mode);
        }
    } else {
        /* No intersection */
    }
}

static PyObject*
tiles_intersecting_webmercator_box_py(PyObject* self, PyObject* args)
{
    double xmin, ymin, xmax, ymax;
    int max_zoom;

    if (!PyArg_ParseTuple(args, "idddd", &max_zoom, &xmin, &ymin, &xmax, &ymax))
        return NULL;

    PyObject* boxes = PyList_New(1);
    PyList_SetItem(boxes, 0, build_box(xmin, ymin, xmax, ymax));

    PyObject* result = PyList_New(0);
    tiles_intersecting_webmercator_boxes(result, boxes, 0, 0, max_zoom, 1);
    return result;
}

static PyObject*
adaptive_tiling_py(PyObject* self, PyObject* args)
{
    double max_err;
    PyObject* boxes;
    if (!PyArg_ParseTuple(args, "dO", &max_err, &boxes))
        return NULL;
     PyObject* disjoint_boxes = make_disjoint_boxes(boxes);
     // Py_DECREF(boxes);

    return adaptive_tiling(disjoint_boxes, max_err, 0.0);
}

static PyObject*
adaptive_tile_covering_py(PyObject* self, PyObject* args)
{
    double max_err, excess_err;
    PyObject* boxes;

    if (!PyArg_ParseTuple(args, "ddO", &max_err, &excess_err, &boxes))
        return NULL;
     PyObject* disjoint_boxes = make_disjoint_boxes(boxes);
     // we could dismiss small area or small (max widht, height) boxes, given the tol. of tiling
     // Py_DECREF(boxes);

    return adaptive_tiling(disjoint_boxes, max_err, excess_err);
}

static PyObject*
approximate_box_by_tiles_py(PyObject* self, PyObject* args)
{
    double xmin, ymin, xmax, ymax;
    int max_zoom;

    if (!PyArg_ParseTuple(args, "idddd", &max_zoom, &xmin, &ymin, &xmax, &ymax))
        return NULL;

    PyObject* boxes = PyList_New(1);
    PyList_SetItem(boxes, 0, build_box(xmin, ymin, xmax, ymax));

    PyObject* result = PyList_New(0);
    tiles_intersecting_webmercator_boxes(result, boxes, 0, 0, max_zoom, 2);
    return result;
}

static PyObject*
tile_covering_py(PyObject* self, PyObject* args)
{
    PyObject* boxes;
    int max_zoom;

    if (!PyArg_ParseTuple(args, "iO", &max_zoom, &boxes))
        return NULL;

     PyObject* disjoint_boxes = make_disjoint_boxes(boxes);
     // Py_DECREF(boxes);

    PyObject* result = PyList_New(0);
    tiles_intersecting_webmercator_boxes(result, disjoint_boxes, 0, 0, max_zoom, 1);
    return result;
}

static PyObject*
tiling_py(PyObject* self, PyObject* args)
{
    PyObject* boxes;
    int max_zoom;

    if (!PyArg_ParseTuple(args, "iO", &max_zoom, &boxes))
        return NULL;

     PyObject* disjoint_boxes = make_disjoint_boxes(boxes);
     // Py_DECREF(boxes);

    PyObject* result = PyList_New(0);
    tiles_intersecting_webmercator_boxes(result, disjoint_boxes, 0, 0, max_zoom, 2);
    return result;
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
 Input: longitude, latitude in WGS84 (SRID 4326) coordinates in degrees
 Output: 62-bit quadkey value, x, y (in mercator projection)
*/
static PyObject*
lonlat2quadintxy_py(PyObject* self, PyObject* args)
{
    double lon, lat;
    uint32 x, y;

    if (!PyArg_ParseTuple(args, "dd", &lon, &lat))
        return NULL;

    lonlat2xy(lon, lat, MAX_ZOOM, &x, &y);

    return Py_BuildValue("KII", xy2quadint(x, y), x, y);
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

/*
 Input: 62-bit quadkey value
 Output:  web mercator bounding box coordinates (SRID 3857)
*/
static PyObject*
tile2bbox_webmercator_py(PyObject* self, PyObject* args)
{
    uint64 quadint;
    int zoom;
    double x_min, y_min, x_max, y_max;

    if (!PyArg_ParseTuple(args, "Ki", &quadint, &zoom))
        return NULL;

    tile2bbox_webmercator(quadint, zoom, &x_min, &y_min, &x_max, &y_max);

    return build_box(x_min, y_min, x_max, y_max);
}

/*
 Input: 62-bit quadkey value
 Output: WGS84  bounding box coordinates (SRID 4326)
*/
static PyObject*
tile2bbox_py(PyObject* self, PyObject* args)
{
    uint64 quadint;
    int zoom;
    double x_min, y_min, x_max, y_max;

    if (!PyArg_ParseTuple(args, "Ki", &quadint, &zoom))
        return NULL;

    tile2bbox(quadint, zoom, &x_min, &y_min, &x_max, &y_max);

    return build_box(x_min, y_min, x_max, y_max);
}

static PyObject*
tile2range_py(PyObject* self, PyObject* args)
{
    uint64 quadint;
    int zoom;
    uint64 q_min, q_max;

    if (!PyArg_ParseTuple(args, "Ki", &quadint, &zoom))
        return NULL;

    tile2range(quadint, zoom, &q_min, &q_max);

    return Py_BuildValue("KK", q_min, q_max);
}

static PyObject*
tile_mask_py(PyObject* self, PyObject* args)
{
    int zoom;

    if (!PyArg_ParseTuple(args, "i", &zoom))
        return NULL;

    return Py_BuildValue("K", tile_prefix_mask(zoom));
}

static PyObject*
tile_center_webmercator_py(PyObject* self, PyObject* args)
{
    uint64 quadint;
    int zoom;
    double x, y;

    if (!PyArg_ParseTuple(args, "Ki", &quadint, &zoom))
        return NULL;

    tile_center_webmercator(quadint, zoom, &x, &y);

    return Py_BuildValue("dd", x, y);
}

static PyObject*
tile_center_py(PyObject* self, PyObject* args)
{
    uint64 quadint;
    int zoom;
    double x, y;

    if (!PyArg_ParseTuple(args, "Ki", &quadint, &zoom))
        return NULL;

    tile_center(quadint, zoom, &x, &y);

    return Py_BuildValue("dd", x, y);
}

static PyObject*
tile_children_py(PyObject* self, PyObject* args)
{
    uint64 quadint;
    int zoom;
    uint64 q_sw, q_nw, q_se, q_ne;

    if (!PyArg_ParseTuple(args, "Ki", &quadint, &zoom))
        return NULL;

    if (zoom < MAX_ZOOM) {
         tile_children(quadint, zoom, &q_sw, &q_nw, &q_se, &q_ne);
         return Py_BuildValue("KKKK", q_sw, q_nw, q_se, q_ne);
    } else {
         return Py_BuildValue("");
    }
}

static PyObject*
xyz2quadint_py(PyObject* self, PyObject* args)
{
    uint32 x, y;
    int zoom;

    if (!PyArg_ParseTuple(args, "IIi", &x, &y, &zoom))
        return NULL;

    return Py_BuildValue("K", xyz2quadint(x, y, zoom));
}

static PyObject*
tile2xyz_py(PyObject* self, PyObject* args)
{
    uint64 quadint;
    int zoom;
    uint32 x, y;

    if (!PyArg_ParseTuple(args, "Ki", &quadint, &zoom))
        return NULL;

    tile2xy(quadint, zoom, &x, &y);
    return Py_BuildValue("IIi", x, y, zoom);
}

static PyObject*
lonlat2xy_py(PyObject* self, PyObject* args)
{
    double lat, lon;
    uint32 x, y;

    if (!PyArg_ParseTuple(args, "dd", &lon, &lat))
        return NULL;

    lonlat2xy(lon, lat, MAX_ZOOM, &x, &y);
    return Py_BuildValue("II", x, y);
}

static PyMethodDef QuadkeyMethods[] =
{
     {"xy2quadint", xy2quadint_py, METH_VARARGS, "xy2quadint_py"},
     {"lonlat2quadint", lonlat2quadint_py, METH_VARARGS, "lonlat2quadint"},
     {"xy2webmercator", xy2webmercator_py, METH_VARARGS, "xy2webmercator"},
     {"webmercator2xy", webmercator2xy_py, METH_VARARGS, "webmercator2xy"},
     {"webmercator2quadint", webmercator2quadint_py, METH_VARARGS, "webmercator2quadint"},
     {"tile2bbox_webmercator", tile2bbox_webmercator_py, METH_VARARGS, "tile2bbox_webmercator"},
     {"tile2bbox", tile2bbox_py, METH_VARARGS, "tile2bbox"},
     {"tile2range", tile2range_py, METH_VARARGS, "tile2range"},
     {"tile_mask", tile_mask_py, METH_VARARGS, "tile_mask"},
     {"tile_center_webmercator", tile_center_webmercator_py, METH_VARARGS, "tile_center_webmercator"},
     {"tile_center", tile_center_py, METH_VARARGS, "tile_center"},
     {"tile_children", tile_children_py, METH_VARARGS, "tile_children"},
     {"xyz2quadint", xyz2quadint_py, METH_VARARGS, "xyz2quadint"},
     {"tile2xyz", tile2xyz_py, METH_VARARGS, "tile2xyz"},
     {"tiles_intersecting_webmercator_box", tiles_intersecting_webmercator_box_py, METH_VARARGS, "tiles_intersecting_webmercator_box"},
     {"approximate_box_by_tiles", approximate_box_by_tiles_py, METH_VARARGS, "approximate_box_by_tiles"},
     {"adaptive_tiling", adaptive_tiling_py, METH_VARARGS, "adaptive_tiling"},
     {"adaptive_tile_covering", adaptive_tile_covering_py, METH_VARARGS, "adaptive_tile_covering"},
     {"tiling", tiling_py, METH_VARARGS, "tiling"},
     {"tile_covering", tile_covering_py, METH_VARARGS, "tile_covering"},
     {"lonlat2xy", lonlat2xy_py, METH_VARARGS, "lonlat2xy"},
     {"lonlat2quadintxy", lonlat2quadintxy_py, METH_VARARGS, "lonlat2quadintxy"},
     {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initquadkey(void)
{
     (void) Py_InitModule("quadkey", QuadkeyMethods);
}
