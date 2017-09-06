import csv, sys, os

from functools import wraps
from flask import Flask, jsonify, redirect, url_for, send_file, request
from flask_restful import reqparse
from flask_restful_swagger_2 import Api, swagger, Resource
from flask_dance.contrib.google import make_google_blueprint, google
from flask_cors import CORS

sys.path.insert(0, "/home/ubuntu/ExpressionMatrix2Test")
from ExpressionMatrix2 import *

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

graph_parser = reqparse.RequestParser()
graph_parser.add_argument('cellsetname', type=str, default='AllCells', required=False, help="Named cell set")
graph_parser.add_argument('similarpairsname', type=str, default='ExactHighInformationGenes', required=False, help="Named setof pairs")
graph_parser.add_argument('similaritythreshold', type=float, required=False, default=0.3, help="Threshold between 0-1")
graph_parser.add_argument('connectivity', type=int, required=False, default=20, help="Maximum connectivity, default is 20")

metadata_parser = reqparse.RequestParser()
metadata_parser.add_argument('celllist', type=str, action="append", required=False, help='List of cells by id')
metadata_parser.add_argument('format', type=str, help='Format: json or csv')

expression_parser = reqparse.RequestParser()
expression_parser.add_argument('celllist', type=str, action="append", required=False, help='List of cells by id')
expression_parser.add_argument('genelist', type=str, action="append", required=False, help='List of genes by name')


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


def parse_exp_data(cells=[], genes=[], limit=0):
	with open(application.config["GBM_DIR"] + "GBM_data-noERCC.csv") as fi:
		reader = csv.reader(fi)
		cell_list = next(reader)
		cell_idxs = []
		cell_data = {}
		gene_list = []
		for idx, el in enumerate(cell_list[1:]):
			if not cells or \
					(cells and el in cells):
				cell_idxs.append(idx + 1)

		if limit:
			if cell_idxs and len(cell_idxs) > limit:
				cell_idxs = cell_idxs[:limit + 1]
			else:
				cell_idxs = [i for i in range(1, limit)]
		for idx in cell_idxs:
			cell_data[cell_list[idx]] = []
		for idx, row in enumerate(reader):
			gene = row[0]
			if not genes or (gene in genes):
				gene_list.append(gene)
			else:
				continue

			for idx in cell_idxs:
				cell_data[cell_list[idx]].append(int(row[idx]))
		cell_data = [{"cellname": k, "e": v} for k, v in cell_data.items()]
	return {
		"genes": gene_list,
		"cells": cell_data
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
#	 if not google.authorized:
#		 return redirect(url_for("google.login"))
#	 resp = google.get("/oauth2/v2/userinfo")
#	 assert resp.ok, resp.text
#	 return "You are {email} on Google".format(email=resp.json()["email"])

# @application.route('/')
# def index():
#	 if not google.authorized:
#		 return redirect(url_for("google.login"))
#	 resp = google.get("/oauth2/v2/userinfo")
#	 assert resp.ok, resp.text
#	 return "You are {email} on Google".format(email=resp.json()["email"])

# @application.route("/test")
# @login_required
# def test():
#	 print(google)
#	 if not google.authorized:
#		 return "Not Authorized"
#	 else:
#		 return "You are OK"

# --- Restful Routes ---------
class MetadataAPI(Resource):
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
		if cell_list and len(metadata['cell_metadata']) < len(cell_list):
			return make_payload([], "Some cell ids not available")
		return make_payload(metadata)


class ExpressionAPI(Resource):
	def __init__(self):
		self.parser = expression_parser

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
		data = parse_exp_data(limit=40)
		return make_payload(data)

	@swagger.doc({
		'description': 'Json with gene list and expression data by cell',
		'parameters': [
			{
				'name': 'celllist',
				'description': 'list of cellid strings',
				'in': 'body',
				'type': 'list of strings',
				"schema": {
					"cellids": []
				}
			},
			{
				'name': 'genelist',
				'description': 'list of gene name strings',
				'in': 'body',
				'type': 'list of strings',
				"schema": {
					"genelist": []
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
	def post(self):
		args = self.parser.parse_args()
		cell_list = args.get('celllist', [])
		gene_list = args.get('genelist', [])
		if not(cell_list) and not(gene_list):
			return make_payload([], "must include either celllist and/or genelist parameter")
		data = parse_exp_data(cell_list, gene_list)
		if cell_list and len(data['cells']) < len(cell_list):
			return make_payload([], "Some cell ids not available")
		if gene_list and len(data['genes']) < len(gene_list):
			return make_payload([], "Some genes not available")
		return make_payload(data)

class GraphAPI(Resource):
	def __init__(self):
		self.parser = graph_parser

	@swagger.doc({
		'description': 'computes the graph for a named cell set',
		'parameters': [
			{
				'name': 'cellsetname',
				'description': 'Named cell set',
				'in': 'query',
				'type': 'str',
			},
			{
				'name': 'similarpairsname',
				'description': 'Named setof pairs',
				'in': 'query',
				'type': 'str',
			},
			{
				'name': 'similaritythreshold',
				'description': 'Threshold between 0-1',
				'in': 'query',
				'type': 'float',
			},
			{
				'name': 'connectivity',
				'description': 'Maximum connectivity',
				'in': 'query',
				'type': 'int',
			}
		],
		'responses': {
			'200': {
				'description': 'xy points for cell graph',
				'examples': {
					'application/json': {
						"data": [
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

					}
				}
			}
		}
	})
	def get(self):
		args = self.parser.parse_args()
		os.chdir(application.config['SCRATCH_DIR'])
		e = ExpressionMatrix(application.config['DATA_DIR'])
		e.createCellGraph('AllCells', args.cellsetname, args.similarpairsname, args.similaritythreshold, args.connectivity)
		e.computeCellGraphLayout('AllCells')
		vertices = e.getCellGraphVertices('AllCells')
		data = [[e.getCellMetaData(v.cellId, 'CellName'), v.x(), v.y()] for v in vertices]
		return make_payload(data)


api.add_resource(MetadataAPI, "/api/v0.1/metadata")
api.add_resource(ExpressionAPI, "/api/v0.1/expression")
api.add_resource(GraphAPI, "/api/v0.1/graph")

if __name__ == "__main__":
	application.run(debug=True)
