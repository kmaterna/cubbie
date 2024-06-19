#!/usr/bin/env python

import requests

# Search - for whatever reason, this needs to be unauthed
session = requests.session()
x = session.get(
    "https://catalogue.dataspace.copernicus.eu/odata/v1/Products?$filter="
    "Collection/Name eq 'SENTINEL-1' and contains(Name,'SLC') and "
    "ContentDate/Start gt 2021-05-03T00:00:00.000Z and "
    "ContentDate/Start lt 2022-05-21T00:00:00.000Z and "
    # "OData.CSC.Intersects(area=geography'SRID=4326;POINT(-0.5319577002158441 28.65487836189358)') and "
    "Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'orbitDirection' and att/OData.CSC.StringAttribute/Value eq 'DESCENDING')"
).json()

print(x)
print(len(x['value']))
