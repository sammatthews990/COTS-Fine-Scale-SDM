"""Package a generated site KML without requiring the native zip executable."""
import argparse
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

p = argparse.ArgumentParser()
p.add_argument('kml', type=Path)
p.add_argument('kmz', type=Path)
a = p.parse_args()
assert a.kml.is_file() and a.kmz.suffix.lower() == '.kmz'
with ZipFile(a.kmz, 'w', ZIP_DEFLATED) as z:
    z.write(a.kml, arcname='doc.kml')
with ZipFile(a.kmz) as z:
    assert z.testzip() is None and z.namelist() == ['doc.kml']
