import json

from flask import (Flask, jsonify, send_from_directory,
                   request, make_response, render_template, Response)
from flask_restful import reqparse
from flask_restful_swagger_2 import Api, swagger, Resource
from flask_cors import CORS
from flask_compress import Compress
from flask_cache import Cache
from werkzeug.datastructures import Headers
import numpy as np
import pandas as pd
from scipy import stats

# SCANPY
import scanpy.api as sc
import pandas.api.types


# CONSTANTS
REACTIVE_LIMIT = 100000
INVALID_CELL_ID = 4294967295

class Float32JSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.float32):
            return float(obj)
        elif isinstance(obj, np.integer):
            return int(obj)
        return json.JSONEncoder.default(self, obj)

class QueryStringError(Exception):
    pass



application = Flask(__name__, static_url_path='/templates')
Compress(application)
CORS(application)
cache = Cache(application, config={'CACHE_TYPE': 'simple'})

# SWAGGER
api = Api(application, api_version='0.1', produces=["application/json"], title="cellxgene rest api",
          api_spec_url='/api/swagger',
          description='An API connecting ExpressionMatrix2 clustering algorithm to cellxgene')


ADATA = sc.read_h5ad("/home/ubuntu/cellxgene/pmbc3k.h5ad")

def _metadata_fields(_adata):
    return _adata.obs.columns

def _is_metadata_categorical(_adata, field):
    metadata_type_dict = _adata.obs.dtypes.to_dict()
    return pandas.api.types.is_categorical_dtype(metadata_type_dict[field])

def _is_metadata_continuous(_adata, field):
    return not _is_metadata_categorical(_adata, field)


def parse_metadata(_adata):
    """
    Get dictionary representation of metadata for a list of cells (or all cells in
    data set if cells not set)

    :param cells: list of cells (False for all cells)
    :return: {metadata_key1: value1, metadata_key2: value1, ...}
    """

    metadata = _adata.obs.to_dict(orient="records")
    for i, idx in enumerate(_adata.obs.index):
        metadata[i]["CellName"] = idx
    return metadata


def get_metadata_ranges(_adata):
    """
    Parses through all metadata to get available values in current dataset
    :param metadata: dictionary of metadata
    :return: dictionary {categorical_key: { "options": [list of options]}, continous_key:
    {"range": {"min": min_val, "max": max_val}}}
    """
    metadata_ranges = {}

    for field in _metadata_fields(_adata):
        if _is_metadata_categorical(_adata, field):
            metadata_ranges[field] = {"options": _adata.obs.groupby(field).size().to_dict()}
        else:
            metadata_ranges[field] = {
                "range": {"min": _adata.obs[field].min(),
                          "max": _adata.obs[field].max()}}

    metadata_ranges["CellName"] = {"options": {i: 1 for i in _adata.obs.index}}

    return metadata_ranges


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
    # Questionable
    data = json.loads(json.dumps(data, cls=Float32JSONEncoder))
    return make_response(jsonify({
        "data": data,
        "status": {
            "error": error,
            "errormessage": errormessage,
        }
    }), errorcode)


def all_genes(_adata):
    """
    Get list of all genenames in dataset
    :return: list of gene names
    """
    return _adata.var.index.tolist()


def all_cells(_adata):
    """
    Get list of all cell names in dataset
    :return: list of cell names
    """
    return _adata.obs.index.tolist()

def query_adata(_adata, query_string):

    # Start with all true
    cell_idx = np.ones((len(all_cells(_adata)),), dtype=bool)

    for key in query_string:
        if key in ["_nograph", "_includeisolated", "_noexpression"]:
            continue
        if key not in _metadata_fields(_adata):
            raise QueryStringError("Error: key {} not in metadata schema".format(key))

        value = query_string.getlist(key)

        if _is_metadata_categorical(_adata, key):
            key_idx = np.in1d(getattr(_adata.obs, key), value)
            cell_idx = np.logical_and(cell_idx, key_idx)
        else:
            value = value[0]
            try:
                min_, max_ = value.split(",")
            except ValueError:
                raise QueryStringError("Error: min,max format required for range for key {}, got {}".format(key, value))
            if min_ != "*":
                key_idx = np.array((getattr(_adata.obs, key) >= min_).data)
                cell_idx = np.logical_and(cell_idx, key_idx)
            if max_ != "*":
                key_idx = np.array((getattr(_adata.obs, key) <= min_).data)
                cell_idx = np.logical_and(cell_idx, key_idx)

    return _adata[cell_idx, :]

def create_graph(_adata, graph_method="umap"):

    # Run the graph method
    getattr(sc.tl, graph_method)(_adata)

    # scanpy stores the coordinates in one of the metadata tables
    df_raw = pd.DataFrame(data=_adata.obsm["X_" + graph_method], index=_adata.obs_names)

    # cellxgene wants values between 0 and 1 I think
    df_norm = (df_raw - df_raw.min()) / (df_raw.max() - df_raw.min())

    return np.hstack((df_norm.index.values.reshape(df_norm.index.shape[0], 1), df_norm.values)).tolist()


def get_expression(_adata, cells, genes=()):
    """
    Get matrix of expression data
    :param cells: list of cell names, or all cells if empty
    :param genes: list of genes, or all genes if empty
    :return: numpy expression matrix
    """
    print(cells)
    print(genes)

    if cells:
        cells_idx = np.in1d(_adata.obs_names, cells)
    else:
        cells_idx = np.ones((len(all_cells(_adata)),), dtype=bool)
    
    if genes:
        genes_idx = np.in1d(_adata.var_names, genes)
    else:
        genes_idx = np.ones((len(all_genes(_adata)),), dtype=bool)
    
    subsetted_adata = _adata[cells_idx,:][:,genes_idx]

    return subsetted_adata

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
    diff_exp = stats.ttest_ind(expression_1.X, expression_2.X)
    set1 = np.logical_and(diff_exp.pvalue < pval, diff_exp.statistic > 0)
    set2 = np.logical_and(diff_exp.pvalue < pval, diff_exp.statistic < 0)
    stat1 = diff_exp.statistic[set1]
    stat2 = diff_exp.statistic[set2]
    sort_set1 = np.argsort(stat1)[::-1]
    sort_set2 = np.argsort(stat2)
    pval1 = diff_exp.pvalue[set1][sort_set1]
    pval2 = diff_exp.pvalue[set2][sort_set2]
    mean_ex1_set1 = np.mean(expression_1.X[:, set1], axis=0)[sort_set1]
    mean_ex2_set1 = np.mean(expression_2.X[:, set1], axis=0)[sort_set1]
    mean_ex1_set2 = np.mean(expression_1.X[:, set2], axis=0)[sort_set2]
    mean_ex2_set2 = np.mean(expression_2.X[:, set2], axis=0)[sort_set2]
    mean_diff1 = mean_ex1_set1 - mean_ex2_set1
    mean_diff2 = mean_ex1_set2 - mean_ex2_set2
    genes_cellset_1 = np.array([expression_1.var_names[idx] for idx, val in enumerate(set1) if val])[sort_set1]
    genes_cellset_2 = np.array([expression_2.var_names[idx] for idx, val in enumerate(set2) if val])[sort_set2]
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

def parse_exp_data(_adata, cells=(), genes=(), limit=0, unexpressed_genes=False):
    """
    Get expression data for set of cells and genes
    :param cells: list of cell names
    :param genes: list of gene names
    :param limit: optional, limit number of genes returned
    :param unexpressed_genes: boolean, filter out genes with no expression
    :return: json: genes: list of genes, cells: expression matrix:, nonzero_gene_count: number of expressed genes
    """

    expression = get_expression(_adata, cells, genes).X

    if len(expression.shape) == 1:
        expression = expression.reshape(expression.shape[0], 1)

    if not genes:
        genes = all_genes(_adata)
    if not cells:
        cells = all_cells(_adata)

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

# ---- Traditional Routes -----
# CellxGene application
@application.route('/')
def index():
    url_base = application.config["CXG_API_BASE"]
    dataset_title = application.config["DATASET_TITLE"]
    return render_template("index.html", prefix=url_base, datasetTitle=dataset_title)


# renders swagger documentation
@application.route('/swagger')
def swag():
    return render_template("swagger.html")


@application.route('/favicon.png')
def icon():
    return send_from_directory('templates', 'favicon.png')

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
        num_genes = args.get("num_genes")
        pval = args.get('pval')

        if not ((cell_list_1 and cell_list_2)):
            return make_payload([],
                                "must include either (cellist1 and cellist2) or (clusters1 and clusters2) parameters",
                                400)
        expression_1 = get_expression(ADATA, cell_list_1)
        expression_2 = get_expression(ADATA, cell_list_2)
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
            queried_adata = query_adata(ADATA, request.args)
        except QueryStringError as e:
            return make_payload({}, str(e), 400)

        graph = None

        if not nograph and len(all_cells(queried_adata)) <= REACTIVE_LIMIT:
            graph = create_graph(queried_adata)

        cells_metadata = parse_metadata(queried_adata)
        ranges = get_metadata_ranges(queried_adata)

        data["ranges"] = ranges
        data["cellcount"] = len(cells_metadata)

        data["reactive"] = True
        data["metadata"] = cells_metadata
        data["cellids"] = [m["CellName"] for m in cells_metadata]
        data["graph"] = graph

        return make_payload(data)


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
        options = get_metadata_ranges(ADATA)
        genes = all_genes(ADATA)
        return make_payload({
            "schema": None,
            "ranges": options,
            "cellcount": len(all_cells(ADATA)),
            "reactivelimit": REACTIVE_LIMIT,
            "genes": genes,
        })

expression_parser = reqparse.RequestParser()
expression_parser.add_argument('celllist', type=str, action="append", required=False, help='List of cells by id')
expression_parser.add_argument('genelist', type=str, action="append", required=False, help='List of genes by name')
expression_parser.add_argument('include_unexpressed_genes', type=bool, required=False,
                               help='Include genes with zero expression across cell set')

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
        data = parse_exp_data(ADATA, limit=40, unexpressed_genes=unexpressed_genes)
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
        data = parse_exp_data(ADATA, cell_list, gene_list, unexpressed_genes=unexpressed_genes)
        if cell_list and len(data['cells']) < len(cell_list):
            return make_payload([], "Some cell ids not available", 400)
        if (gene_list and len(data['genes']) < len(gene_list)) and unexpressed_genes:
            return make_payload([], "Some genes not available", 400)
        return make_payload(data)

api.add_resource(InitializeAPI, "/api/v0.1/initialize")
api.add_resource(CellsAPI, "/api/v0.1/cells")
api.add_resource(DifferentialExpressionAPI, "/api/v0.1/diffexpression")
api.add_resource(ExpressionAPI, "/api/v0.1/expression")

if __name__ == "__main__":
    application.run(host='0.0.0.0', debug=True, port=4002)
