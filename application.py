import os
import time
import io
import csv
from collections import defaultdict
from functools import wraps
from ExpressionMatrix2 import ExpressionMatrix  # noqa: E402

from flask import Flask, jsonify, send_from_directory, request, make_response, render_template, Response
from flask_restful import reqparse
from flask_restful_swagger_2 import Api, swagger, Resource
from flask_cors import CORS
from flask_compress import Compress
from flask_cache import Cache
import pybrake.flask
import numpy as np
from scipy import stats
import boto3

from schemaparse import parse_schema

# CONSTANTS
REACTIVE_LIMIT = 100000
INVALID_CELL_ID = 4294967295


# Local Errors
class QueryStringError(Exception):
    pass


class Vertex:

    def __init__(self, cellId, x, y):
        self.cellId = cellId
        self.x_coord = x
        self.y_coord = y

    def x(self):
        return float(self.x_coord)

    def y(self):
        return float(self.y_coord)


application = Flask(__name__, static_url_path='/templates')
Compress(application)
CORS(application)
cache = Cache(application, config={'CACHE_TYPE': 'simple'})

# SECRETS
SECRET_KEY = os.environ.get("SECRET_KEY", default=None)
if not SECRET_KEY:
    raise ValueError("No secret key set for Flask application")
APP_USERNAME = os.environ.get("APP_USERNAME", default="")
APP_PASSWORD = os.environ.get("APP_PASSWORD", default="")
CONFIG_FILE = os.environ.get("CONFIG_FILE", default="")
CXG_API_BASE = os.environ.get("CXG_API_BASE", default="")
AIRBRAKE_PROJECT_ID = os.environ.get("AIRBRAKE_PROJECT_ID", default="")
AIRBRAKE_PROJECT_KEY = os.environ.get("AIRBRAKE_PROJECT_KEY", default="")
AIRBRAKE_ENVIRONMENT = os.environ.get("AIRBRAKE_ENVIRONMENT", default="")
application.config['PYBRAKE'] = {
    "project_id": AIRBRAKE_PROJECT_ID,
    "project_key": AIRBRAKE_PROJECT_KEY,
    "environment": AIRBRAKE_ENVIRONMENT
}
application = pybrake.flask.init_app(application)

if not CONFIG_FILE:
    raise ValueError("No config file set for Flask application")
# CONFIG
application.config.from_pyfile(CONFIG_FILE, silent=True)
application.config.update(
    SECRET_KEY=SECRET_KEY,
    USERNAME=APP_USERNAME,
    PASSWORD=APP_PASSWORD,
    CXG_API_BASE=CXG_API_BASE,
)
dir_path = os.path.dirname(os.path.realpath(__file__))
application.config.update(
    SCRATCH_DIR=os.path.join(dir_path, application.config["SCRATCH_DIR"]),
    DATA_DIR=os.path.join(dir_path, application.config["DATA_DIR"]),
)

# SWAGGER
api = Api(application, api_version='0.1', produces=["application/json"], title="cellxgene rest api",
          api_spec_url='/api/swagger',
          description='An API connecting ExpressionMatrix2 clustering algorithm to cellxgene')


def s3_size_with_prefix(resource, bucket, prefix):
    bucket = resource.Bucket(bucket)
    size = 0
    for obj in bucket.objects.filter(Prefix=prefix):
        size += obj.size
    return size


def get_dir_size(path):
    total = 0
    for entry in os.scandir(path):
        if entry.is_file():
            total += entry.stat().st_size
        elif entry.is_dir():
            total += get_dir_size(entry.path)
    return total


def download_data_from_s3():
    resource = boto3.resource('s3')
    prefix = os.path.basename(application.config["DATA_DIR"])
    # Download if doesn't exist or only partially exists
    if (not os.path.exists(application.config["DATA_DIR"])) or (
            get_dir_size(application.config["DATA_DIR"]) != s3_size_with_prefix(resource,
                                                                                application.config["AWS_BUCKET"],
                                                                                prefix)):
        download_s3_bucket(resource, application.config["AWS_BUCKET"],
                           application.config["DATA_DIR"], prefix)


def download_s3_bucket(resource, bucket, dest, prefix):
    print("Downloading")
    bucket = resource.Bucket(bucket)
    for obj in bucket.objects.filter(Prefix=prefix):
        filename = obj.key[len(prefix) + 1:]
        full_destination = os.path.join(dest, filename)
        full_destination_dir = os.path.dirname(full_destination)
        if not os.path.exists(full_destination_dir):
            os.makedirs(full_destination_dir)
        with open(full_destination, 'wb') as data:
            bucket.download_fileobj(obj.key, data)


# Download data if not exists
if not os.path.exists(application.config["DATA_DIR"]):
    download_data_from_s3()

# Param parsers
metadata_parser = reqparse.RequestParser()
metadata_parser.add_argument('celllist', type=str, action="append", required=False, help='List of cells by id')
metadata_parser.add_argument('format', type=str, help='Format: json or csv')

expression_parser = reqparse.RequestParser()
expression_parser.add_argument('celllist', type=str, action="append", required=False, help='List of cells by id')
expression_parser.add_argument('genelist', type=str, action="append", required=False, help='List of genes by name')
expression_parser.add_argument('include_unexpressed_genes', type=bool, required=False,
                               help='Include genes with zero expression across cell set')

differential_parser = reqparse.RequestParser()
differential_parser.add_argument('celllist1', type=str, action="append", required=False,
                                 help='First group of cells by id')
differential_parser.add_argument('celllist2', type=str, action="append", required=False,
                                 help='Second group of cells by id')
differential_parser.add_argument('clusters1', type=str, action="append", required=False,
                                 help='First group of cluters by name')
differential_parser.add_argument('clusters2', type=str, action="append", required=False,
                                 help='Second group of clusters by name')
differential_parser.add_argument('num_genes', type=int, required=False, default=20,
                                 help="Number of diff-expressed genes to return")
differential_parser.add_argument('pval', type=float, required=False, default=0.0001,
                                 help="Pval max of diff-expressed genes")

e = ExpressionMatrix(application.config["DATA_DIR"], True)
schema = parse_schema(os.path.join(application.config["DATA_DIR"], application.config["SCHEMA_FILE"]))


# ---- Helper Functions -------
def check_auth(username, password):
    """This function is called to check if a username /
    password combination is valid.
    """
    return username == application.config["USERNAME"] and password == application.config["PASSWORD"]


def authenticate():
    """Sends a 401 response that enables basic auth"""
    return Response(
        'Could not verify your access level for that URL.\n'
        'You have to login with proper credentials', 401,
        {'WWW-Authenticate': 'Basic realm="Login Required"'})


def requires_auth(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        auth = request.authorization
        if application.config["PASSWORD_PROTECT"] and (not auth or not check_auth(auth.username, auth.password)):
            print("auth!!")
            return authenticate()
        return f(*args, **kwargs)

    return decorated


def filter(cell, qs):
    """
    Check if cell's metadata are within the filtered parameters
    :param cell: cell id
    :param qs: query string
    :return: Boolean True if cell within filter, False if rejected by filter
    """
    global e
    keep = True
    for key, value in qs.items():
        if value["variabletype"] == 'categorical':
            if e.getCellMetaDataValue(cell, key) not in value["query"]:
                keep = False
                break
        elif value["variabletype"] == 'continuous':
            if value["query"]["min"] and value["query"]["max"]:
                keep = value["query"]["min"] <= float(e.getCellMetaDataValue(cell, key)) <= value["query"]["max"]
            elif value["query"]["min"]:
                keep = value["query"]["min"] <= float(e.getCellMetaDataValue(cell, key))
            elif value["query"]["max"]:
                keep = float(e.getCellMetaDataValue(cell, key)) <= value["query"]["max"]

            if not keep:
                break
    return keep


def parse_metadata(cells=False):
    """
    Get dictionary representation of metadata for a list of cells (or all cells in data set if cells not set)
    :param cells: list of cells (False for all cells)
    :return: {metadata_key1: value1, metadata_key2: value1, ...}
    """
    global e
    metadata = []
    if cells is False:
        cells = e.getCellSet('AllCells')
    if INVALID_CELL_ID not in cells:
        mdata = e.getCellsMetaData(cells)
        for i in mdata:
            metadata.append({k: cast_value(k, v) for k, v in i})
    return metadata


def convert_variable(datatype, variable):
    """
    Convert variable to number (float/int)
    Used for dataset metadata and for query string
    :param datatype: type to convert to
    :param variable: value of variable
    :return: converted variable
    :raises: ValueError
    """
    try:
        if variable and datatype == "int":
            variable = int(variable)
        elif variable and datatype == "float":
            variable = float(variable)
        return variable
    except ValueError:
        raise


def cast_value(key, value):
    """
    Cast metadata value to proper type from dataset metadata
    This needs to be flexible since metadata is often dirty
    :param key: string key
    :param value: string metadata value
    :return:
    """
    val_type = schema[key]["type"]
    new_val = value
    try:
        if (key == "EM2Cluster") and value.startswith("Unclustered"):
            new_val = "NoCluster"
        else:
            convert_variable(val_type, value)
    except ValueError:
        pass
    if value == "":
        new_val = ""
    return new_val


def parse_exp_data(cells=(), genes=(), limit=0, unexpressed_genes=False):
    """
    Get expression data for set of cells and genes
    :param cells: list of cell names
    :param genes: list of gene names
    :param limit: optional, limit number of genes returned
    :param unexpressed_genes: boolean, filter out genes with no expression
    :return: json: genes: list of genes, cells: expression matrix:, nonzero_gene_count: number of expressed genes
    """
    if cells:
        cells = sort_celllist_by_id(cells)
    if genes:
        genes = sort_genelist_by_id(genes)
    expression = get_expression(cells, genes)
    if not genes:
        genes = all_genes()
    if not cells:
        cells = all_cells()
    if not unexpressed_genes:
        expression, genes = remove_unexpressed_genes(expression, genes)
    if limit and len(genes) > limit:
        genes = genes[:limit]
        expression = expression[:, :limit]
    cell_data = []
    for idx, cell in enumerate(cells):
        cell_data.append({
            "cellname": cell,
            "e": list(expression[idx]),
        })
    return {
        "genes": genes,
        "cells": cell_data,
        "nonzero_gene_count": int(np.sum(expression.any(axis=1)))
    }


def remove_unexpressed_genes(expression, genes):
    """
    Filter out genes that have no expression data from expression matrix
    :param expression: numpy matrix will expression data
    :param genes: list of genes names
    :return: tuple, filtered expression numpy matrix, list of genes
    """
    genes_expressed = expression.any(axis=0)
    genes = [genes[idx] for idx, val in enumerate(genes_expressed) if val]
    return expression[:, genes_expressed], genes


def get_expression(cells, genes=()):
    """
    Get matrix of expression data
    :param cells: list of cell names, or all cells if empty
    :param genes: list of genes, or all genes if empty
    :return: numpy expression matrix
    """
    global e
    try:
        if not cells:
            cellset = "AllCells"
        else:
            cell_ids = [get_cell_id(name) for name in cells]
            cellset = "ExpressionCellSet"
            e.createCellSet(cellset, cell_ids)
        if not genes:
            expression = e.getDenseExpressionMatrix("AllGenes", cellset)
        else:
            genesetName = "ExpressionGeneSet_{}".format(time.time())
            e.createGeneSetFromGeneNames(genesetName, genes)
            expression = e.getDenseExpressionMatrix(genesetName, cellset)
    except Exception as error:
        raise error
    finally:
        if cells:
            e.removeCellSet(cellset)
        if genes:
            e.removeGeneSet(genesetName)
    return expression


def make_payload(data, errormessage="", errorcode=200):
    """
    Creates JSON respons for requests
    :param data: json data
    :param errormessage: error message
    :param errorcode: http error code
    :return: flask json repsonse
    """
    error = False
    if errormessage:
        error = True
    return make_response(jsonify({
        "data": data,
        "status": {
            "error": error,
            "errormessage": errormessage,
        }
    }), errorcode)


def make_csv(data, filename="export.csv"):
    """
    Creates csv response for requests
    :param data: json data
    :param errormessage: error message
    :param errorcode: http error code
    :return: flask json repsonse
    """
    output = make_response(data)
    output.headers["Content-Disposition"] = "attachment; filename={}".format(filename)
    output.headers["Content-type"] = "text/csv"
    return output


def parse_querystring(qs):
    """
    Parses filter query and transforms to json
    :param qs: dictionary query string
    :return: reformatted query string
    """
    query = {}
    for key in qs:
        if key in ["_nograph", "_includeisolated", "_noexpression"]:
            continue
        value = qs.getlist(key)
        if key not in schema:
            raise QueryStringError("Error: key {} not in metadata schema".format(key))
        query[key] = schema[key]
        if query[key]["variabletype"] == "categorical":
            query[key]["query"] = convert_variable(query[key]["variabletype"], value)
        elif query[key]["variabletype"] == "continuous":
            value = value[0]
            try:
                min, max = value.split(",")
            except ValueError:
                raise QueryStringError("Error: min,max format required for range for key {}, got {}".format(key, value))
            if min == "*":
                min = None
            if max == "*":
                max = None
            try:
                query[key]["query"] = {
                    "min": convert_variable(query[key]["type"], min),
                    "max": convert_variable(query[key]["type"], max)
                }
            except ValueError:
                raise QueryStringError(
                    "Error: expected type {} for key {}, got {}".format(query[key]["type"], key, value))
    return query


def get_metadata_ranges(metadata):
    """
    Parses through all metadata to get available values in current dataset
    :param metadata: dictionary of metadata
    :return: dictionary {categorical_key: { "options": [list of options]}, continous_key:
    {"range": {"min": min_val, "max": max_val}}}
    """
    options = {}
    for s in schema:
        if schema[s]["include"]:
            options[s] = {}
            if schema[s]["variabletype"] == "categorical":
                # Default dict
                options[s]["options"] = defaultdict(int)
            else:
                options[s]["range"] = {
                    "min": None,
                    "max": None
                }

    for cell in metadata:
        for s in schema:
            if schema[s]["include"]:
                try:
                    datum = cell[s]
                    if schema[s]["type"] == "int":
                        datum = int(datum)
                    elif schema[s]["type"] == "float":
                        datum = float(datum)
                    if schema[s]["variabletype"] == "categorical":
                        if (s.startswith("EM2Cluster")) and datum.startswith("Unclustered"):
                            datum = "NoCluster"
                        options[s]["options"][datum] += 1
                    elif schema[s]["variabletype"] == "continuous":
                        if not options[s]["range"]["min"] or datum < options[s]["range"]["min"]:
                            options[s]["range"]["min"] = datum
                        if not options[s]["range"]["max"] or datum > options[s]["range"]["max"]:
                            options[s]["range"]["max"] = datum
                except KeyError:
                    pass
                except ValueError:
                    pass
                except TypeError:
                    pass
    for key, value in dict(options).items():
        if ("options" in value and not value["options"]) or ("range" in value and not value["range"]["max"]):
            del options[key]
    return options


def cells_from_query(qs, output_cellset):
    """
    Get cells matching user's filter
    :param qs: dictionary, query string from user
    :param output_cellset: string, name for cellset
    :return: list of cell ids
    """
    global e
    cellset_from_query(qs, output_cellset)
    cellidlist = e.getCellSet(output_cellset)
    return cellidlist


def cellset_from_query(qs, output_cellset):
    """
    Creates a cell set matching user's filter. **Important** must delete cellset when done
    :param qs: dictionary, query string from user
    :param output_cellset: string, name for cellset
    :return: output cellset name
    """
    global e
    filters = []
    for key, value in qs.items():
        if value["variabletype"] == 'categorical':
            category_filter = []
            if value["query"] == [""]:
                filtername = "{}".format(key)
                e.createCellSetUsingMetaData(filtername, key, "^$", True)
                category_filter.append(filtername)
            elif key.startswith("EM2Cluster") and "NoCluster" in value["query"]:
                filtername = "{}".format(key)
                e.createCellSetUsingMetaData(filtername, key, "Unclustered.*", True)
                category_filter.append(filtername)
            else:
                for idx, item in enumerate(value["query"]):
                    queryval = item
                    filtername = "{}_{}".format(key, idx)
                    e.createCellSetUsingMetaData(filtername, key, queryval, False)
                    category_filter.append(filtername)
            category_output_filter = "{}_out".format(key)
            e.createCellSetUnion(",".join(category_filter), category_output_filter)
            filters.append(category_output_filter)
            for cellset in category_filter:
                e.removeCellSet(cellset)
        elif value["variabletype"] == 'continuous':
            filtername = "{}".format(key)
            filters.append(filtername)
            if value["query"]["min"] and value["query"]["max"]:
                e.createCellSetUsingNumericMetaDataBetween(filtername, key, value["query"]["min"],
                                                           value["query"]["max"])
            elif value["query"]["min"]:
                e.createCellSetUsingNumericMetaDataGreaterThan(filtername, key, value["query"]["min"])
            elif value["query"]["max"]:
                e.createCellSetUsingNumericMetaDataLessThan(filtername, key, value["query"]["max"])
            else:
                filters.pop(filters.index(filtername))
    e.createCellSetIntersection(",".join(filters), output_cellset)
    for cellset in filters:
        e.removeCellSet(cellset)
    return output_cellset


@cache.memoize(86400000)
def create_graph(output_cellset, includeisolated=False, graphname="cellgraph"):
    """
    Computes graphlayout for given cellset
    :param output_cellset: string, name of EM2 cellset
    :param includeisolated: bool, include vertices with no edges
    :param graphname: string, name to save EM2 graph
    :return: list of tuples [(label, x, y), ...], list of kept cellids
    """
    global e
    # Graph using EM2 or use metadata coordinates
    if application.config["GRAPH_EM2"]:
        graph, cellidlist = create_graph_em2(output_cellset, includeisolated=False, graphname="cellgraph")
    else:
        graph, cellidlist = create_graph_tsne(output_cellset)
    return graph, cellidlist


def create_graph_em2(output_cellset, includeisolated=False, graphname="cellgraph"):
    """
    Computes graphlayout for given cellset using Expression Matrix 2
    :param output_cellset: string, name of EM2 cellset
    :param includeisolated: bool, include vertices with no edges
    :param graphname: string, name to save EM2 graph
    :return: list of tuples [(label, x, y), ...], list of kept cellids
    """
    e.createCellGraph(graphname, output_cellset, application.config["HIGH_INFO_NAME"],
                      keepIsolatedVertices=includeisolated)
    e.computeCellGraphLayout(graphname)
    vertices = e.getCellGraphVertices(graphname)
    normalized_verticies = normalize_verticies(vertices)
    graph = []
    # reset cellids if create the graph
    cellidlist = []
    for i in range(len(normalized_verticies["labels"])):
        graph.append((get_cell_name(normalized_verticies["labels"][i]), normalized_verticies["x"][i],
                      normalized_verticies["y"][i]))
        cellidlist.append(normalized_verticies["labels"][i])
    return graph, cellidlist


def create_graph_tsne(output_cellset):
    """
    Layout for given cellset usig precomputed coordinates
    :param output_cellset: string, name of EM2 cellset
    :return: list of tuples [(label, x, y), ...], list of kept cellids
    """
    # get metadata for each cell
    cellidlist = e.getCellSet(output_cellset)
    vertices = [Vertex(get_cell_name(cid), e.getCellMetaDataValue(cid, "tSNE_1"),
                       e.getCellMetaDataValue(cid, "tSNE_2")) for cid in cellidlist]
    normalized_verticies = normalize_verticies(vertices)
    graph = []
    for i in range(len(normalized_verticies["labels"])):
        graph.append((normalized_verticies["labels"][i], normalized_verticies["x"][i],
                      normalized_verticies["y"][i]))
    return graph, cellidlist


def normalize_verticies(verticies):
    """
    normalize x,y verticies for graph
    :param verticies: list of vertex objects
    :return: json: {"labels": list of point labels, "x":
        list of normalized x positions, "y": list of normalized y positions}
    """
    labels = []
    x = []
    y = []
    for v in verticies:
        labels.append(v.cellId)
        x.append(v.x())
        y.append(v.y())
    if x:
        x = normalize(x)
    if y:
        y = normalize(y)
    return {
        "labels": labels,
        "x": x,
        "y": y
    }


def normalize(v):
    """
    run normalization, zero safe
    :param v: vector to normalize
    :return: normalized vector
    """
    vec = np.array(v)
    vec = (vec - min(vec)) / (max(vec) - min(vec))
    return vec


def all_genes():
    """
    Get list of all genenames in dataset
    :return: list of gene names
    """
    global e
    return [e.geneName(i) for i in range(e.geneCount())]


def all_cells():
    """
    Get list of all cell names in dataset
    :return: list of cell names
    """
    global e
    return [get_cell_name(i) for i in e.getCellSet('AllCells')]


def diffexp(expression_1, expression_2, pval=0.001, num_genes=20):
    """
    get top expressed genes from two different cell sets (uses t-test)
    :param expression_1: numpy expression array cell set 1
    :param expression_2: numpy expression array cell set 2
    :param pval: stats limit
    :param num_genes: number of genes to limit results to
    :return: Top genes and mean expression, pvalue, and average difference
     between cell set1 to cell set 2 for each gene in both cell sets
     sorted by t-test value
    """
    global e
    diff_exp = stats.ttest_ind(expression_1, expression_2)
    set1 = np.logical_and(diff_exp.pvalue < pval, diff_exp.statistic > 0)
    set2 = np.logical_and(diff_exp.pvalue < pval, diff_exp.statistic < 0)
    stat1 = diff_exp.statistic[set1]
    stat2 = diff_exp.statistic[set2]
    sort_set1 = np.argsort(stat1)[::-1]
    sort_set2 = np.argsort(stat2)
    pval1 = diff_exp.pvalue[set1][sort_set1]
    pval2 = diff_exp.pvalue[set2][sort_set2]
    mean_ex1_set1 = np.mean(expression_1[:, set1], axis=0)[sort_set1]
    mean_ex2_set1 = np.mean(expression_2[:, set1], axis=0)[sort_set1]
    mean_ex1_set2 = np.mean(expression_1[:, set2], axis=0)[sort_set2]
    mean_ex2_set2 = np.mean(expression_2[:, set2], axis=0)[sort_set2]
    mean_diff1 = mean_ex1_set1 - mean_ex2_set1
    mean_diff2 = mean_ex1_set2 - mean_ex2_set2
    genes_cellset_1 = np.array([e.geneName(idx) for idx, val in enumerate(set1) if val])[sort_set1]
    genes_cellset_2 = np.array([e.geneName(idx) for idx, val in enumerate(set2) if val])[sort_set2]
    return {
        "celllist1": {
            "topgenes": genes_cellset_1.tolist()[:num_genes],
            "mean_expression_cellset1": mean_ex1_set1.tolist()[:num_genes],
            "mean_expression_cellset2": mean_ex2_set1.tolist()[:num_genes],
            "pval": pval1.tolist()[:num_genes],
            "ave_diff": mean_diff1.tolist()[:num_genes]
        },
        "celllist2": {
            "topgenes": genes_cellset_2.tolist()[:num_genes],
            "mean_expression_cellset1": mean_ex1_set2.tolist()[:num_genes],
            "mean_expression_cellset2": mean_ex2_set2.tolist()[:num_genes],
            "pval": pval2.tolist()[:num_genes],
            "ave_diff": mean_diff2.tolist()[:num_genes]
        },
    }


def get_clusters(list_of_clusters):
    """
    Get list of cells that are members of the cluster(s)
    :param list_of_clusters: list of strings cluster names
    :return: list of cell names
    """
    global e
    cluster_filter = "|".join(list_of_clusters)
    cellset_name = "cluster_set_{}".format(cluster_filter)
    cluster_filter = "^({})$".format(cluster_filter)
    e.createCellSetUsingMetaData(cellset_name, application.config["CLUSTER_METADATA_KEY"], cluster_filter, True)
    cellids = e.getCellSet(cellset_name)
    e.removeCellSet(cellset_name)
    return [get_cell_name(i) for i in cellids]


def get_cell_name(cellid):
    """
    Get cell name from cell id
    :param cellid: int cell id
    :return: string name of cell
    """
    global e
    return e.getCellMetaDataValue(cellid, "CellName")


def get_cell_id(cellname):
    """
    Get cell id from cell name
    :param cellname: string cell name
    :return: int id of cell
    """
    global e
    return e.cellIdFromString(cellname)


def get_gene_id(genename):
    """
    Get gene id from gene name
    :param genename: string gene name
    :return: int id of gene
    """
    global e
    return e.geneIdFromName(genename)


def sort_celllist_by_id(celllist):
    """
    Sorts list of cell names by their assigned cellid in EM2 (needed because idx info is lost
    when creating cellset)
    :param celllist:
    :return: cellist
    """
    celllist_ids = np.array([get_cell_id(cell) for cell in celllist])
    sort_idx = np.argsort(celllist_ids)
    celllist = np.array(celllist)[sort_idx]
    return celllist.tolist()


def sort_genelist_by_id(genelist):
    """
    Sorts list of gene names by their assigned geneid in EM2 (needed because idx info is lost
    when creating geneset)
    :param genelist:
    :return: geneist
    """
    genelist_ids = np.array([get_gene_id(gene) for gene in genelist])
    sort_idx = np.argsort(genelist_ids)
    genelist = np.array(genelist)[sort_idx]
    return genelist.tolist()


# ---- Traditional Routes -----
# CellxGene application
@application.route('/')
@requires_auth
def index():
    url_base = application.config["CXG_API_BASE"]
    dataset_title = application.config["DATASET_TITLE"]
    return render_template("index.html", prefix=url_base, datasetTitle=dataset_title)


# renders swagger documentation
@application.route('/swagger')
@requires_auth
def swag():
    return render_template("swagger.html")


@application.route('/favicon.png')
def icon():
    return send_from_directory('templates', 'favicon.png')


# --- Restful Routes ---------
class MetadataAPI(Resource):
    def __init__(self):
        self.parser = metadata_parser

    @swagger.doc({
        'summary': 'Json document containing the metadata for list of cells',
        'tags': ['metadata'],
        'parameters': [
            {
                'name': 'body',
                'in': 'body',
                'schema': {
                    "example": {
                        "celllist": ["1001000173.G8", "1001000173.D4"],
                    }
                }
            }
        ],
        "responses": {
            '200': {
                'description': 'json with list of cell metadata k:v pairs',
                'examples': {
                    'application/json': {
                        "data": {
                            "cell_metadata": [
                                {
                                    "Cell": "1001000173.G8",
                                    "Class": "Neoplastic",
                                    "Cluster_2d": "11",
                                    "Cluster_2d_color": "#8C564B",
                                    "Cluster_CNV": "1",
                                    "Cluster_CNV_color": "#1F77B4",
                                    "ERCC_reads": "152104",
                                    "ERCC_to_non_ERCC": "0.562454470489481",
                                    "Genes_detected": "1962",
                                    "Location": "Tumor",
                                    "Location.color": "#FF7F0E",
                                    "Multimapping_reads_percent": "2.67",
                                    "Neoplastic": "Neoplastic",
                                    "Non_ERCC_reads": "270429",
                                    "Sample.name": "BT_S2",
                                    "Sample.name.color": "#AEC7E8",
                                    "Sample.type": "Glioblastoma",
                                    "Sample.type.color": "#1F77B4",
                                    "Selection": "Unpanned",
                                    "Selection.color": "#98DF8A",
                                    "Splice_sites_AT.AC": "102",
                                    "Splice_sites_Annotated": "122397",
                                    "Splice_sites_GC.AG": "761",
                                    "Splice_sites_GT.AG": "125741",
                                    "Splice_sites_non_canonical": "56",
                                    "Splice_sites_total": "126660",
                                    "Total_reads": "1741039",
                                    "Unique_reads": "1400382",
                                    "Unique_reads_percent": "80.43",
                                    "Unmapped_mismatch": "2.15",
                                    "Unmapped_other": "0.18",
                                    "Unmapped_short": "14.56",
                                    "housekeeping_cluster": "2",
                                    "housekeeping_cluster_color": "#AEC7E8",
                                    "recluster_myeloid": "NA",
                                    "recluster_myeloid_color": "NA"
                                },
                                {
                                    "Cell": "1001000173.D4",
                                    "Class": "Regular",
                                    "Cluster_2d": "2",
                                    "Cluster_2d_color": "#AEC7E8",
                                    "Cluster_CNV": "2",
                                    "Cluster_CNV_color": "#AEC7E8",
                                    "ERCC_reads": "244319",
                                    "ERCC_to_non_ERCC": "1.44247380073801",
                                    "Genes_detected": "515",
                                    "Location": "Tumor",
                                    "Location.color": "#FF7F0E",
                                    "Multimapping_reads_percent": "0.84",
                                    "Neoplastic": "Regular",
                                    "Non_ERCC_reads": "169375",
                                    "Sample.name": "BT_S2",
                                    "Sample.name.color": "#AEC7E8",
                                    "Sample.type": "Glioblastoma",
                                    "Sample.type.color": "#1F77B4",
                                    "Selection": "Unpanned",
                                    "Selection.color": "#98DF8A",
                                    "Splice_sites_AT.AC": "203",
                                    "Splice_sites_Annotated": "103763",
                                    "Splice_sites_GC.AG": "683",
                                    "Splice_sites_GT.AG": "105878",
                                    "Splice_sites_non_canonical": "93",
                                    "Splice_sites_total": "106857",
                                    "Total_reads": "1229919",
                                    "Unique_reads": "1081271",
                                    "Unique_reads_percent": "87.91",
                                    "Unmapped_mismatch": "0.59",
                                    "Unmapped_other": "0",
                                    "Unmapped_short": "10.66",
                                    "housekeeping_cluster": "2",
                                    "housekeeping_cluster_color": "#AEC7E8",
                                    "recluster_myeloid": "NA",
                                    "recluster_myeloid_color": "NA"
                                }
                            ]
                        },
                        "status": {
                            "error": False,
                            "errormessage": "",
                        }
                    }
                }
            },
            '400': {
                'description': 'Cell ids selected not available',
            }
        }
    })
    def post(self):
        args = self.parser.parse_args()
        cell_list = args.celllist
        if cell_list:
            cell_list = [e.cellIdFromString(i) for i in cell_list]
        metadata = parse_metadata(cell_list)
        if cell_list and len(metadata) < len(cell_list):
            return make_payload([], "Some cell ids not available", 400)
        return make_payload({"cell_metadata": metadata})


class ExpressionAPI(Resource):
    def __init__(self):
        self.parser = expression_parser

    @swagger.doc({
        'summary': 'Json with gene list and expression data by cell, limited to first 40 cells',
        'tags': ['expression'],
        'parameters': [
            {
                'name': 'include_unexpressed_genes',
                'description': "Include genes that have 0 expression across all cells in set",
                'in': 'path',
                'type': 'bool',
            }
        ],
        'responses': {
            '200': {
                'description': 'Json for heatmap',
                'examples': {
                    'application/json': {
                        "data": {
                            "cells": [
                                {
                                    "cellname": "1/2-SBSRNA4",
                                    "e": [0, 0, 214, 0, 0]
                                },
                            ],
                            "genes": [
                                "1001000173.G8",
                                "1001000173.D4",
                                "1001000173.B4",
                                "1001000173.A2",
                                "1001000173.E2"
                            ],
                            "nonzero_gene_count": 2857
                        },
                        "status": {
                            "error": False,
                            "errormessage": ""
                        }
                    }
                }
            }
        }
    })
    def get(self):
        # TODO add choice for limit (pagination?)
        args = self.parser.parse_args()
        unexpressed_genes = args.get('include_unexpressed_genes', False)
        data = parse_exp_data(limit=40, unexpressed_genes=unexpressed_genes)
        return make_payload(data)

    @swagger.doc({
        'summary': 'Json with gene list and expression data by cell',
        'tags': ['expression'],
        'parameters': [
            {
                'name': 'body',
                'in': 'body',
                "schema": {
                    "example": {
                        "celllist": ["1001000173.G8", "1001000173.D4"],
                        "genelist": ["1/2-SBSRNA4", "A1BG", "A1BG-AS1", "A1CF", "A2LD1", "A2M", "A2ML1", "A2MP1",
                                     "A4GALT"],
                        "include_unexpressed_genes": True,
                    }

                }
            },
        ],
        'responses': {
            '200': {
                'description': 'Json for expressiondata',
                'examples': {
                    'application/json': {
                        "data": {
                            "cells": [
                                {
                                    "cellname": "1001000173.D4",
                                    "e": [0, 0]
                                },
                                {
                                    "cellname": "1001000173.G8",
                                    "e": [0, 0]
                                }
                            ],
                            "genes": [
                                "ABCD4",
                                "ZWINT"
                            ],
                            "nonzero_gene_count": 2857
                        },
                        "status": {
                            "error": False,
                            "errormessage": ""
                        }

                    }
                }
            },
            '400': {
                'description': 'Required parameter missing/incorrect',
            }
        }
    })
    def post(self):
        args = self.parser.parse_args()
        cell_list = args.get('celllist', [])
        gene_list = args.get('genelist', [])
        unexpressed_genes = args.get('include_unexpressed_genes', False)
        if not (cell_list) and not (gene_list):
            return make_payload([], "must include celllist and/or genelist parameter", 400)
        data = parse_exp_data(cell_list, gene_list, unexpressed_genes=unexpressed_genes)
        if cell_list and len(data['cells']) < len(cell_list):
            return make_payload([], "Some cell ids not available", 400)
        if (gene_list and len(data['genes']) < len(gene_list)) and unexpressed_genes:
            return make_payload([], "Some genes not available", 400)
        return make_payload(data)


# class ClusterAPI(Resource):
#     def get(self):
#         # Retrieve cluster
#         run = randint(0, 10000)
#         graphname = "cellgraph_{}".format(run)
#         e = ExpressionMatrix(application.config['DATA_DIR'], True)
#         args = graph_parser.parse_args()
#         args = graph_parser.parse_args()
#         e.createCellGraph(graphname, args.cellsetname, args.similarpairsname, args.similaritythreshold,
#                           args.connectivity)
#         clustername = "clusters_{}".format(run)
#         params = ClusterGraphCreationParameters()
#         e.createClusterGraph(graphname, params, clustername)
#         clusters = e.getClusterGraphVertices(clustername)
#         clusters = [i for i in clusters]
#         cells = {}
#         for i in clusters:
#             cells[i] = [e.getCellMetaDataValue(cid, "CellName") for cid in e.getClusterCells(clustername, i)]
#         # TODO calculate expression values for these
#         # TODO does the graph really need to be created each time?
#         # TODO parameterize graph with cell sets?
#
#         return make_payload({clustername: clustername}, errorcode=201)


class DifferentialExpressionAPI(Resource):
    def __init__(self):
        self.parser = differential_parser

    @swagger.doc({
        'summary': 'Get the top expressed genes for two cell sets. Calculated using t-test',
        'tags': ['expression'],
        'parameters': [
            {
                'name': 'body',
                'in': 'body',
                'schema': {
                    "example": {
                        "celllist1": ["1001000176.C12", "1001000176.C7", "1001000177.F11"],
                        "celllist2": ["1001000012.D2", "1001000017.F10", "1001000033.C3", "1001000229.D4"],
                        "num_genes": 5,
                        "pval": 0.000001,
                    },
                }
            }
        ],
        "responses": {
            '200': {
                'description': 'top expressed genes for cellset1, cellset2',
                'examples': {
                    'application/json': {
                        "data": {
                            "celllist1": {
                                "ave_diff": [
                                    432.0132935431362,
                                    12470.5623982637,
                                    957.0246880086814
                                ],
                                "mean_expression_cellset1": [
                                    438.6185567010309,
                                    13315.536082474227,
                                    1076.5773195876288
                                ],
                                "mean_expression_cellset2": [
                                    6.605263157894737,
                                    844.9736842105264,
                                    119.55263157894737
                                ],
                                "pval": [
                                    3.8906598089944563e-35,
                                    1.9086226376018916e-25,
                                    7.847480544069826e-21
                                ],
                                "topgenes": [
                                    "TMSB10",
                                    "FTL",
                                    "TMSB4X"
                                ]
                            },
                            "celllist2": {
                                "ave_diff": [
                                    -6860.599158979924,
                                    -519.1314432989691,
                                    -10278.328269126423
                                ],
                                "mean_expression_cellset1": [
                                    2.8350515463917527,
                                    0.6185567010309279,
                                    23.09278350515464
                                ],
                                "mean_expression_cellset2": [
                                    6863.434210526316,
                                    519.75,
                                    10301.421052631578
                                ],
                                "pval": [
                                    4.662891833748732e-44,
                                    3.6278087029927103e-37,
                                    8.396825170618402e-35
                                ],
                                "topgenes": [
                                    "SPARCL1",
                                    "C1orf61",
                                    "CLU"
                                ]
                            }
                        },
                        "status": {
                            "error": False,
                            "errormessage": ""
                        }
                    }
                }
            }
        }
    })
    def post(self):
        args = self.parser.parse_args()
        # TODO allow different calculation methods
        cell_list_1 = args.get('celllist1', [])
        cell_list_2 = args.get('celllist2', [])
        cluster_list_1 = args.get('clusters1', [])
        cluster_list_2 = args.get('clusters2', [])
        num_genes = args.get("num_genes")
        pval = args.get('pval')
        if cluster_list_1 and cluster_list_2:
            cell_list_1 = get_clusters(cluster_list_1)
            cell_list_2 = get_clusters(cluster_list_2)
        if not ((cell_list_1 and cell_list_2)):
            return make_payload([],
                                "must include either (cellist1 and cellist2) or (clusters1 and clusters2) parameters",
                                400)
        cell_list_1 = sort_celllist_by_id(cell_list_1)
        expression_1 = get_expression(cell_list_1, all_genes())
        expression_2 = get_expression(cell_list_2, all_genes())
        data = diffexp(expression_1, expression_2, pval, num_genes)
        return make_payload(data)


class CellsAPI(Resource):
    @swagger.doc({
        'summary': 'filter based on metadata fields to get a subset cells, expression data, and metadata',
        'tags': ['cells'],
        'description': "Cells takes query parameters defined in the schema retrieved from the /initialize enpoint. "
                       "<br>For categorical metadata keys filter based on `key=value` <br>"
                       " For continuous metadata keys filter by `key=min,max`<br> Either value "
                       "can be replaced by a \*. To have only a minimum value `key=min,\*`  To have only a maximum "
                       "value `key=\*,max` <br>Graph data (if retrieved) is normalized"
                       " To only retrieve cells that don't have a value for the key filter by `key`",
        'parameters': [
            {
                'name': '_nograph',
                'description': "Do not calculate and send back graph (graph is sent by default)",
                'in': 'path',
                'type': 'boolean',
            },
            {
                'name': '_includeisolated',
                'description': "Include all cells in graph even if cluster is too small, "
                               "this does nothing if _nograph is true.",
                'in': 'path',
                'type': 'boolean',
            }
        ],

        'responses': {
            '200': {
                'description': 'initialization data for UI',
                'examples': {
                    'application/json': {
                        "data": {
                            "badmetadatacount": 0,
                            "cellcount": 0,
                            "cellids": ["..."],
                            "metadata": [
                                {
                                    "CellName": "1001000173.G8",
                                    "Class": "Neoplastic",
                                    "Cluster_2d": "11",
                                    "Cluster_2d_color": "#8C564B",
                                    "Cluster_CNV": "1",
                                    "Cluster_CNV_color": "#1F77B4",
                                    "ERCC_reads": "152104",
                                    "ERCC_to_non_ERCC": "0.562454470489481",
                                    "Genes_detected": "1962",
                                    "Location": "Tumor",
                                    "Location.color": "#FF7F0E",
                                    "Multimapping_reads_percent": "2.67",
                                    "Neoplastic": "Neoplastic",
                                    "Non_ERCC_reads": "270429",
                                    "Sample.name": "BT_S2",
                                    "Sample.name.color": "#AEC7E8",
                                    "Sample.type": "Glioblastoma",
                                    "Sample.type.color": "#1F77B4",
                                    "Selection": "Unpanned",
                                    "Selection.color": "#98DF8A",
                                    "Splice_sites_AT.AC": "102",
                                    "Splice_sites_Annotated": "122397",
                                    "Splice_sites_GC.AG": "761",
                                    "Splice_sites_GT.AG": "125741",
                                    "Splice_sites_non_canonical": "56",
                                    "Splice_sites_total": "126660",
                                    "Total_reads": "1741039",
                                    "Unique_reads": "1400382",
                                    "Unique_reads_percent": "80.43",
                                    "Unmapped_mismatch": "2.15",
                                    "Unmapped_other": "0.18",
                                    "Unmapped_short": "14.56",
                                    "housekeeping_cluster": "2",
                                    "housekeeping_cluster_color": "#AEC7E8",
                                    "recluster_myeloid": "NA",
                                    "recluster_myeloid_color": "NA"
                                },
                            ],
                            "reactive": True,
                            "graph": [
                                [
                                    "1001000173.G8",
                                    0.93836,
                                    0.28623
                                ],

                                [
                                    "1001000173.D4",
                                    0.1662,
                                    0.79438
                                ]

                            ],
                            "status": {
                                "error": False,
                                "errormessage": ""
                            }

                        },
                    }
                },
            },

            '400': {
                'description': 'bad query params',
            }
        }
    })
    def get(self):
        global e
        # TODO Allow random downsampling
        e = ExpressionMatrix(application.config["DATA_DIR"], True)
        data = {
            "reactive": False,
            "cellids": [],
            "metadata": [],
            "cellcount": 0,
            "badmetadatacount": 0,
            "graph": [],
            "ranges": {},
            "expression": []
        }
        # Parse args
        try:
            args = request.args
            nograph = False
            includeisolated = False
            if "_nograph" in args:
                nograph = bool(args["_nograph"])
            if "_includeisolated" in args:
                includeisolated = bool(args["_includeisolated"])
            qs = parse_querystring(request.args)
        except QueryStringError as e:
            return make_payload({}, str(e), 400)
        cells_metadata = []
        output_cellset = "outcellset_{}".format(request.query_string)
        cellidlist = []
        if len(qs):
            cellidlist = cells_from_query(qs, output_cellset)
        else:
            output_cellset = "AllCells"
            cellidlist = e.getCellSet(output_cellset)
        try:
            graph = None
            if not nograph and len(cellidlist) <= REACTIVE_LIMIT:
                graph, cellidlist = create_graph(output_cellset, includeisolated, "cellgraph")
            cells_metadata = parse_metadata(cellidlist)
            ranges = get_metadata_ranges(cells_metadata)
            data["ranges"] = ranges
            data["cellcount"] = len(cells_metadata)
        finally:
            if output_cellset != "AllCells":
                e.removeCellSet(output_cellset)
            data["reactive"] = True
            data["metadata"] = cells_metadata
            data["cellids"] = [m["CellName"] for m in cells_metadata]
            data["graph"] = graph
        return make_payload(data)


class EgressAPI(Resource):

    @swagger.doc({
        'summary': 'Download expression data for filter based on metadata fields',
        'tags': ['download'],
        'description': "Cells takes query parameters defined in the schema retrieved from the /initialize enpoint. "
                       "<br>For categorical metadata keys filter based on `key=value` <br>"
                       " For continuous metadata keys filter by `key=min,max`<br> Either value "
                       "can be replaced by a \*. To have only a minimum value `key=min,\*`  To have only a maximum "
                       "value `key=\*,max` <br>Graph data (if retrieved) is normalized"
                       " To only retrieve cells that don't have a value for the key filter by `key`",
        'responses': {
            '200': {
                'description': 'csv of cell by gene expression counts',
            },

            '400': {
                'description': 'bad query params or bad download',
            }
        }
    })
    def get(self):
        global e
        e = ExpressionMatrix(application.config["DATA_DIR"], True)
        try:
            qs = parse_querystring(request.args)
        except QueryStringError as error:
            return make_payload({}, str(error), 400)
        output_cellset = "outcellset_{}".format(request.query_string)
        if not len(qs):
            output_cellset = "AllCells"
        try:
            # create cellset
            cells_from_query(qs, output_cellset)
            expression = e.getDenseExpressionMatrix("AllGenes", output_cellset)
            expression = expression.transpose()
            genes = all_genes()
            cellids = e.getCellSet(output_cellset)
            # transform to file
            si = io.StringIO()
            writer = csv.writer(si)
            writer.writerow([""] + [get_cell_name(cid) for cid in cellids])
            for idx, row in enumerate(expression):
                writer.writerow([genes[idx]] + row.tolist())
            return make_csv(si.getvalue(), "expression_matrix.csv")
        except Exception as error:
            print("ERROR creating csv", error)
            return make_payload([], "Error: creating csv download", 400)
        finally:
            if output_cellset != "AllCells":
                e.removeCellSet(output_cellset)


class InitializeAPI(Resource):
    @swagger.doc({
        'summary': 'get metadata schema, ranges for values, and cell count to initialize cellxgene app',
        'tags': ['initialize'],
        'parameters': [],
        'responses': {
            '200': {
                'description': 'initialization data for UI',
                'examples': {
                    'application/json': {
                        "data": {
                            "cellcount": 3589,
                            "options": {
                                "Sample.type": {
                                    "options": {
                                        "Glioblastoma": 3589
                                    }
                                },
                                "Selection": {
                                    "options": {
                                        "Astrocytes(HEPACAM)": 714,
                                        "Endothelial(BSC)": 123,
                                        "Microglia(CD45)": 1108,
                                        "Neurons(Thy1)": 685,
                                        "Oligodendrocytes(GC)": 294,
                                        "Unpanned": 665
                                    }
                                },
                                "Splice_sites_AT.AC": {
                                    "range": {
                                        "max": 1025,
                                        "min": 152
                                    }
                                },
                                "Splice_sites_Annotated": {
                                    "range": {
                                        "max": 1075869,
                                        "min": 26
                                    }
                                }
                            },
                            "schema": {
                                "CellName": {
                                    "displayname": "Name",
                                    "type": "string",
                                    "variabletype": "categorical"
                                },
                                "Class": {
                                    "displayname": "Class",
                                    "type": "string",
                                    "variabletype": "categorical"
                                },
                                "ERCC_reads": {
                                    "displayname": "ERCC Reads",
                                    "type": "int",
                                    "variabletype": "continuous"
                                },
                                "ERCC_to_non_ERCC": {
                                    "displayname": "ERCC:Non-ERCC",
                                    "type": "float",
                                    "variabletype": "continuous"
                                },
                                "Genes_detected": {
                                    "displayname": "Genes Detected",
                                    "type": "int",
                                    "variabletype": "continuous"
                                }
                            },
                            "genes": ["1/2-SBSRNA4", "A1BG", "A1BG-AS1"]

                        },
                        "status": {
                            "error": False,
                            "errormessage": ""
                        }
                    }
                }
            }
        }
    })
    def get(self):
        metadata = parse_metadata()
        options = get_metadata_ranges(metadata)
        genes = all_genes()
        return make_payload({
            "schema": schema,
            "ranges": options,
            "cellcount": len(metadata),
            "reactivelimit": REACTIVE_LIMIT,
            "genes": genes,
        })


api.add_resource(MetadataAPI, "/api/v0.1/metadata")
api.add_resource(ExpressionAPI, "/api/v0.1/expression")
api.add_resource(InitializeAPI, "/api/v0.1/initialize")
api.add_resource(CellsAPI, "/api/v0.1/cells")
api.add_resource(DifferentialExpressionAPI, "/api/v0.1/diffexpression")
api.add_resource(EgressAPI, "/api/v0.1/egress")

if __name__ == "__main__":
    application.run(host='0.0.0.0', debug=True)
