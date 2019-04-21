import csv
import json
from shapely.geometry import Polygon,shape

data = []
with open("data/poly.csv") as ifile:
    reader = csv.DictReader(ifile)
    for row in reader:
        data.append([float(row["x"]), float(row["y"])])

caseShp = Polygon(data)

data = []
with open("data/example.csv") as ifile:
    reader = csv.DictReader(ifile)
    for row in reader:
        data.append([float(row["x"]), float(row["y"])])

caseShp1 = Polygon(data)

chinaJson = json.loads(open("data/china_sim.geojson").read())
chinaShp = shape(chinaJson["features"][0]["geometry"])

provinces = []
provincesJson = json.loads(open("data/provinces.geojson").read())
for prov in provincesJson["features"]:
    shp = shape(prov["geometry"])
    name = prov["properties"]["CNAME"]
    code = prov["properties"]["GB"]
    provinces.append({"code":code,"name":name, "shape":shp})

cities = []
citiesJson = json.loads(open("data/cities.geojson").read())
for city in citiesJson["features"]:
    shp = shape(city["geometry"])
    name = city["properties"]["CNAME"]
    code = city["properties"]["GB"]
    cities.append({"code":code,"name":name, "shape":shp})