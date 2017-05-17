
This is a small python library to work with quadkeys in a fast way

## installing

pip install git+https://github.com/CartoDB/python-quadkey.git

## how to use it

```python
import quadkey
print(quadkey.lonlat2quadint(-73.969558715820312, 40.757678985595703)
1013670064444553679
```
## Aknowledgments

The fast technique for quadkey encodig/decoding is taken from https://github.com/yinqiwen/geohash-int

## Available Methods

`lonlat2quadint(lon, lat)`

* Input: longitude, latitude in WGS84 (SRID 4326) coordinates in degrees
* Output: quadint value

`webmercator2quadint(x, y)`

* Input:  web mercator coordinates (SRID 3857)
* Output: quadint value

### TODO

Functions to handle tiles (defined by a quadint value and a zoom level)

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

`tile_mask(quadint, zoom)`

* Input: quadint value and zoom level
* Output: quadint mask to select the tile

`tile_children(quadint, zoom)`

* Input: quadint value and zoom level
* Output: quadint values of the 4 children of the tile of level zoom+1
