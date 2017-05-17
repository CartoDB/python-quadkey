
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
