import json
import urllib.request
from urllib.error import HTTPError

from pysisyphus.io.sdf import geom_from_sdf


def cid_for_name(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IUPACName/JSON"
    with urllib.request.urlopen(url) as handle:
        text = handle.read()
    json_ = json.loads(text)
    cid = json_["PropertyTable"]["Properties"][0]["CID"]
    return cid


def sdf_for_cid(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    with urllib.request.urlopen(url) as handle:
        sdf = handle.read().decode("utf-8")
    return sdf


def geom_for_pubchem_name(name, **kwargs):
    try:
        cid = cid_for_name(name)
        sdf = sdf_for_cid(cid)
        geom = geom_from_sdf(sdf, **kwargs)
    except HTTPError:
        geom = None
    return geom
