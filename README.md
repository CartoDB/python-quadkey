
This is a small python library to work with quadkeys in a fast way

## installing

pip install git+https://github.com/CartoDB/python-quadkey.git

### from source

python -m pip install .

## how to use it

```python
import quadkey
print(quadkey.lonlat2quadint(-73.969558715820312, 40.757678985595703)
1013670064444553679
```

### Quadint vs Quadkey

quadints are in base 10, if you need the _quadkey_, transform the resuting quadint into base 4,
then take the first n digits corresponding to the required zoom level on the tile.

## Acknowledgments

The fast technique for quadkey encodig/decoding is taken from https://github.com/yinqiwen/geohash-int

## Available Methods

Note: WebMercator is a cylindrical projection which does not include latitudes at the poles.

### Conversion functions:

`xy2quadint(x, y)`
* Input:  (uint32) x, y coordinates based on WebMercator in the range [0,2^31)
* Output: (uint64) quadint value

`xy2webmercator(x, y)`
* Input:  (uint32) x, y coordinates based on WebMercator in the range [0,2^31)
* Output: (double) web mercator coordinates (SRID 3857)

`quadint2webmercator(quadint)`
* Input:  (uint64) quadint value
* Output: (double) web mercator coordinates (SRID 3857)

`quadint2xy(quadint)`
* Input:  (uint64) quadint value
* Output: (uint32) x, y coordinates based on WebMercator in the range [0,2^31)

`lonlat2quadint(lon, lat)`
* Input:  (double) longitude, latitude in WGS84 (SRID 4326) coordinates in degrees
* Output: (uint64) quadint value

`lonlat2xy(lon, lat)`
* Input:  (double) longitude, latitude in WGS84 (SRID 4326) coordinates in degrees
* Output: (uint32) x, y coordinates based on WebMercator in the range [0,2^31)

`lonlat2quadintxy(lon, lat)`
* Input:  (double) longitude, latitude in WGS84 (SRID 4326) coordinates in degrees
* Output: (uint64) quadint value
* Output: (uint32) x, y coordinates based on WebMercator in the range [0,2^31)

`webmercator2quadint(x, y)`
* Input:  (double) web mercator coordinates (SRID 3857)
* Output: (uint64) quadint value

`webmercator2xy(x, y)`
* Input:  (double) web mercator coordinates (SRID 3857)
* Output: (uint32) x, y coordinates based on WebMercator in the range [0,2^31)


### Functions to handle tiles (defined by a quadint value and a zoom level)

`tile2range(quadint, zoom)`

* Input: quadint value and zoom level
* Output: minimum and maximum quadint values

`tile2bbox(quadint, zoom)`

* Input: quadint value and zoom level
* Output: bounding box defined by minimum longitude, minimum latitude, maximum longitude and maximum latitude (WGS84)

`tile2bbox_webmercator(quadint, zoom)`

* Input: quadint value and zoom level
* Output: bounding box defined by minimum x, minimum y, maximum x, maximum y (Web Mercator)

`tile_center(quadint, zoom)`

* Input: quadint value and zoom level
* Output: longitude and latitude of the tile center (WGS84)

`tile_center_webmercator(quadint, zoom)`

* Input: quadint value and zoom level
* Output: longitude and latitude of the tile center (Webmercator)

`tile_children(quadint, zoom)`

* Input: quadint value and zoom level
* Output: quadint values of the 4 children of the tile of level zoom+1

`tile_mask(zoom)`

* Input: zoom level
* Output: quadint mask to select tiles of the level (with bitwise AND)

`tile_suffix_mask(zoom)`

* Input: zoom level
* Output: quadint mask to select the significant bits from a zoom level (with bitwise AND)

### Conversion of tile representation between quadint, zoom and x, y zoom:

`xyz2quadint(x, y, zoom)`

* Input: x, y, z coordinates of the tile
* Output: quadint value of the tile

`tile2xyz(quadint, zoom)`

* Input: quadint value and zoom level
* Output: x, y, z of the tile

### Obtaning the list of tiles that intersect a web mercator rectangle,
given a maximum tile zoom level.

`tiles_intersecting_webmercator_box(max_zoom, xmin, ymin, xmax, ymax)`
