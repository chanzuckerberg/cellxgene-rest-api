import csv

from functools import wraps
from flask import Flask, jsonify, redirect, url_for, send_file, request
from flask_restful import reqparse
from flask_restful_swagger_2 import Api, swagger, Resource
from flask_dance.contrib.google import make_google_blueprint, google
from flask_cors import CORS

application = Flask(__name__)
CORS(application)

application.config.from_pyfile('app.cfg', silent=True)
api = Api(application, api_version='0.1', api_spec_url='/api/swagger')
blueprint = make_google_blueprint(
    client_id=application.config["GOOGLE_CLIENT_ID"],
    client_secret=application.config["GOOGLE_CLIENT_SECRET"],
    scope=["profile", "email"]
)
application.register_blueprint(blueprint, url_prefix="/login")

metadata_parser = reqparse.RequestParser()
metadata_parser.add_argument('celllist', type=str, action="append", required=False, help='List of cells by id')
metadata_parser.add_argument('format', type=str, help='Format: json or csv')


# ---- Helper Functions -------
def parse_metadata(cell_ids=False):
    with open(application.config["GBM_DIR"] + "GBM_metadata.csv") as fi:
        reader = csv.reader(fi)
        metadata = []
        header = next(reader)
        for row in reader:
            if cell_ids and row[0] not in cell_ids:
                continue
            metadata.append({k: v for k, v in zip(header, row)})
    return {"cell_metadata": metadata}


def parse_exp_data(limit=40):
    with open(application.config["GBM_DIR"] + "GBM_data-noERCC.csv") as fi:
        reader = csv.reader(fi)
        data = []
        header = next(reader)
        for idx, row in enumerate(reader):
            if limit and idx == limit:
                break
            data.append({
                "cellname": row[0],
                "e": [int(i) for i in row[1:]]
            })
    return {
        "genes": header[1:],
        "cells": data
    }


def make_payload(data, errormessage=""):
    error = False
    if errormessage:
        error = True
    return jsonify({
        "data": data,
        "status": {
            "error": error,
            "errormessage": errormessage,
        }
    })


# ---- Decorators -------------
def login_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        if not google.authorized:
            return redirect(url_for("google.login"))
        return f(*args, **kwargs)

    return decorated_function


# ---- Traditional Routes -----
@application.route('/')
def index():
    return "Clutering API"


# @application.route('/glogin')
# def glogin():
#     if not google.authorized:
#         return redirect(url_for("google.login"))
#     resp = google.get("/oauth2/v2/userinfo")
#     assert resp.ok, resp.text
#     return "You are {email} on Google".format(email=resp.json()["email"])

# @application.route('/')
# def index():
#     if not google.authorized:
#         return redirect(url_for("google.login"))
#     resp = google.get("/oauth2/v2/userinfo")
#     assert resp.ok, resp.text
#     return "You are {email} on Google".format(email=resp.json()["email"])

# @application.route("/test")
# @login_required
# def test():
#     print(google)
#     if not google.authorized:
#         return "Not Authorized"
#     else:
#         return "You are OK"

# --- Restful Routes ---------
class MetadataAPI(Resource):
    # TODO Swagger documentation

    def __init__(self):
        self.parser = metadata_parser

    @swagger.doc({
        'description': 'Returns a either a csv file or a json document containing the metadata for all cells',
        'parameters': [
            {
                'name': 'format',
                'description': 'format of file',
                'in': 'query',
                'type': 'string'
            }
        ],
        'responses': {
            '200': {
                'description': 'A csv file or json file',
            }
        }
    })
    def get(self):
        args = self.parser.parse_args()
        if args.format and args.format.lower() == 'json':
            metadata = parse_metadata()
            return make_payload(metadata)
        return send_file(application.config["GBM_DIR"] + "GBM_metadata.csv", mimetype="text/csv")

    @swagger.doc({
        'description': 'json document containing the metadata for list of cells',
        'parameters': [
            {
                'name': 'celllist',
                'description': 'list of cellid strings',
                'in': 'body',
                'type': 'list of strings',
                "schema": {
                    "cellids": []
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
            }
        }
    })
    def post(self):
        args = self.parser.parse_args()
        cell_list = args.celllist
        metadata = parse_metadata(cell_list)
        if len(metadata['cell_metadata']) < len(cell_list):
            return make_payload([], "Some cell ids not available")
        return make_payload(metadata)


class HeatmapAPI(Resource):
    # TODO Swagger documentation

    @swagger.doc({
        'description': 'Json with gene list and expression data by cell, limited to first 40 cells',
        'parameters': [],
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
                            ]
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
        data = parse_exp_data()
        return make_payload(data)


api.add_resource(MetadataAPI, "/api/v0.1/metadata")
api.add_resource(HeatmapAPI, "/api/v0.1/heatmap")

if __name__ == "__main__":
    application.run(debug=True)
