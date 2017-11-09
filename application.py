import sys, os
from collections import defaultdict
from random import randint
from re import escape

from functools import wraps
from flask import Flask, jsonify, send_from_directory, request, make_response, render_template
from flask_restful import reqparse
from flask_restful_swagger_2 import Api, swagger, Resource
from flask_cors import CORS
import numpy as np

from schemaparse import parse_schema

application = Flask(__name__, static_url_path='/templates')
CORS(application)

REACTIVE_LIMIT = 5000
# SECRETS
SECRET_KEY = os.environ.get("SECRET_KEY", default=None)
if not SECRET_KEY:
	raise ValueError("No secret key set for Flask application")
CXG_API_BASE = os.environ.get("CXG_API_BASE", default=None)
if not CXG_API_BASE:
	raise ValueError("No api base")

application.config.from_pyfile('app.cfg', silent=True)
application.config.update(
	CXG_API_BASE=CXG_API_BASE,
	SECRET_KEY=SECRET_KEY
)
dir_path = os.path.dirname(os.path.realpath(__file__))
application.config.update(
	GBM_DIR=os.path.join(dir_path, application.config["GBM_DIR"]),
	SCRATCH_DIR=os.path.join(dir_path, application.config["SCRATCH_DIR"]),
	DATA_DIR=os.path.join(dir_path, application.config["DATA_DIR"]),
	SCHEMAS_DIR=os.path.join(dir_path, application.config["SCHEMAS_DIR"]),
	EM2_DIR=os.path.join(dir_path, application.config["EM2_DIR"]),
)

api = Api(application, api_version='0.1', produces=["application/json"], title="cellxgene rest api",
		  api_spec_url='/api/swagger',
		  description='A API connecting ExpressionMatrix2 clustering algorithm to cellxgene')

sys.path.insert(0, application.config["EM2_DIR"])
from ExpressionMatrix2 import ExpressionMatrix, CellIdList, invalidCellId, ClusterGraphCreationParameters

graph_parser = reqparse.RequestParser()
graph_parser.add_argument('cellsetname', type=str, default='AllCells', required=False, help="Named cell set")
graph_parser.add_argument('similarpairsname', type=str, default='ExactHighInformationGenes', required=False,
						  help="Named setof pairs")
graph_parser.add_argument('similaritythreshold', type=float, required=False, default=0.3, help="Threshold between 0-1")
graph_parser.add_argument('connectivity', type=int, required=False, default=20,
						  help="Maximum connectivity, default is 20")

metadata_parser = reqparse.RequestParser()
metadata_parser.add_argument('celllist', type=str, action="append", required=False, help='List of cells by id')
metadata_parser.add_argument('format', type=str, help='Format: json or csv')

expression_parser = reqparse.RequestParser()
expression_parser.add_argument('celllist', type=str, action="append", required=False, help='List of cells by id')
expression_parser.add_argument('genelist', type=str, action="append", required=False, help='List of genes by name')
expression_parser.add_argument('include_unexpressed_genes', type=bool, required=False,
							   help='Include genes with zero expression across cell set')

cluster_parser = reqparse.RequestParser()
cluster_parser.add_argument('clustername', type=str, required=True, help="Name of cell graph")


# ---- Helper Functions -------
def stringstringpairs2dict(ssp):
	return {i.first: i.second for i in ssp}


def parse_metadata(cells):
	e = ExpressionMatrix(application.config["DATA_DIR"], True)

	mdict = []
	if invalidCellId not in cells:
		metadata = e.getCellsMetaData(cells)
		mdict = [stringstringpairs2dict(m) for m in metadata]
	return {"cell_metadata": mdict}


def parse_exp_data(cells, genes=(), limit=0, unexpressed_genes=False):
	e = ExpressionMatrix(application.config["DATA_DIR"], True)

	if not genes:
		genes = [e.geneName(gid) for gid in range(e.geneCount())]
	expression = np.zeros([len(genes), len(list(cells))])
	for idx, gene in enumerate(genes):
		expression[idx] = list(e.getCellsExpressionCountFromGeneName(cells, gene))
	if not unexpressed_genes:
		genes_expressed = expression.any(axis=1)
		genes = [genes[idx] for idx, val in enumerate(genes_expressed) if val]
		expression = expression[genes_expressed]

	if limit and len(genes) > limit:
		genes = genes[:limit]
		expression = expression[:limit]
	cell_data = []
	for idx, cid in enumerate(cells):
		cell_data.append({
			"cellname": e.getCellMetaDataValue(cid, "CellName"),
			"e": list(expression[:, idx]),
		})
	return {
		"genes": genes,
		"cells": cell_data,
		"nonzero_gene_count": int(np.sum(expression.any(axis=1)))
	}


def toCellIDsCpp(cells=()):
	e = ExpressionMatrix(application.config["DATA_DIR"], True)
	cell_number_ids = CellIdList()
	if cells:
		for cell_id in cells:
			cell_number_ids.append(e.cellIdFromString(cell_id))
	else:
		cell_number_ids = e.getCellSet('AllCells')
	return cell_number_ids


def make_payload(data, errormessage="", errorcode=200):
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


def parse_querystring(qs):
	schema = parse_schema(os.path.join(application.config["SCHEMAS_DIR"], "test_data_schema.json"))
	query = {}
	for key in qs:
		if key == "_nograph":
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


def get_metadata_ranges(schema, metadata):
	options = {}
	for s in schema:
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
			try:
				datum = cell[s]
				if schema[s]["type"] == "int":
					datum = int(datum)
				elif schema[s]["type"] == "float":
					datum = float(datum)
				if schema[s]["variabletype"] == "categorical":
					options[s]["options"][datum] += 1
				elif schema[s]["variabletype"] == "continuous":
					if not options[s]["range"]["min"] or datum < options[s]["range"]["min"]:
						options[s]["range"]["min"] = datum
					if not options[s]["range"]["max"] or datum > options[s]["range"]["max"]:
						options[s]["range"]["max"] = datum
			except KeyError:
				pass
			except TypeError:
				pass

	return options


def normalize_verticies(verticies):
	labels = []
	x = []
	y = []
	for v in verticies:
		labels.append(v.cellId)
		x.append(v.x())
		y.append(v.y())
	x = normalize(x)
	y = normalize(y)
	return {
		"labels": labels,
		"x": x,
		"y": y
	}


def normalize(v):
	vec = np.array(v)
	norm = np.linalg.norm(vec, ord=np.inf)
	if norm == 0:
		return vec
	return vec / norm


class QueryStringError(Exception):
	pass


def convert_variable(datatype, variable):
	try:
		if variable and datatype == "int":
			variable = int(variable)
		elif variable and datatype == "float":
			variable = float(variable)
		return variable
	except ValueError:
		print("Bad conversion")
		raise


# ---- Traditional Routes -----
@application.route('/')
def index():
	url_base = application.config["CXG_API_BASE"]
	split_base = url_base.split("/")
	prefix = '/'.join(split_base[:-2]) + "/"
	version = "{}/".format(split_base[-2])
	return render_template("index.html", prefix=prefix, version=version)


@application.route('/swagger')
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
		metadata = parse_metadata(toCellIDsCpp(cell_list))
		if cell_list and len(metadata['cell_metadata']) < len(cell_list):
			return make_payload([], "Some cell ids not available", 400)
		return make_payload(metadata)


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

		data = parse_exp_data(cells=toCellIDsCpp(), limit=40, unexpressed_genes=unexpressed_genes)
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
		data = parse_exp_data(toCellIDsCpp(cell_list), gene_list, unexpressed_genes=unexpressed_genes)
		if cell_list and len(data['cells']) < len(cell_list):
			return make_payload([], "Some cell ids not available", 400)
		if gene_list and len(data['genes']) < len(gene_list):
			return make_payload([], "Some genes not available", 400)
		return make_payload(data)


class GraphAPI(Resource):
	def __init__(self):
		self.parser = graph_parser

	@swagger.doc({
		'summary': 'computes the graph for a named cell set',
		'tags': ['graph'],
		'parameters': [
			{
				'name': 'cellsetname',
				'description': "Named cell set ex. 'AllCells'",
				'in': 'path',
				'type': 'string',
			},
			{
				'name': 'similarpairsname',
				'description': "Named set of pairs ex. 'ExtractHighInformationGenes'",
				'in': 'path',
				'type': 'string',
			},
			{
				'name': 'similaritythreshold',
				'description': 'Threshold between 0-1 ex. 0.3',
				'in': 'path',
				'type': 'number',
			},
			{
				'name': 'connectivity',
				'description': 'Maximum connectivity ex. 20',
				'in': 'path',
				'type': 'number',
			}
		],
		'responses': {
			'200': {
				'description': 'xy points for cell graph',
				'examples': {
					'application/json': {
						"data": {
							"1001000173.G8":
								[
									0.93836,
									0.28623
								],

							"1001000173.D4":
								[
									0.1662,
									0.79438
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
		args = self.parser.parse_args()
		os.chdir(application.config['DATA_DIR'])
		e = ExpressionMatrix(application.config['DATA_DIR'], True)
		run = randint(0, 10000)
		graphname = "cellgraph_{}".format(run)
		e.createCellGraph(graphname, args.cellsetname, args.similarpairsname, args.similaritythreshold,
						  args.connectivity)
		e.computeCellGraphLayout(graphname)
		vertices = e.getCellGraphVertices(graphname)
		data = {e.getCellMetaDataValue(v.cellId, 'CellName'): [v.x(), v.y()] for v in vertices}
		return make_payload(data)


class ClusterAPI(Resource):
	def get(self):
		# Retrieve cluster
		run = randint(0, 10000)
		graphname = "cellgraph_{}".format(run)
		e = ExpressionMatrix(application.config['DATA_DIR'], True)
		args = graph_parser.parse_args()
		e.createCellGraph(graphname, args.cellsetname, args.similarpairsname, args.similaritythreshold,
						  args.connectivity)
		clustername = "clusters_{}".format(run)
		params = ClusterGraphCreationParameters()
		e.createClusterGraph(graphname, params, clustername)
		clusters = e.getClusterGraphVertices(clustername)
		clusters = [i for i in clusters]
		cells = {}
		for i in clusters:
			cells[i] = [e.getCellMetaDataValue(cid, "CellName") for cid in e.getClusterCells(clustername, i)]
		# TODO calculate expression values for these
		# TODO does the graph really need to be created each time?
		# TODO parameterize graph with cell sets?

		return make_payload({clustername: clustername}, errorcode=201)


class CellsAPI(Resource):
	@swagger.doc({
		'summary': 'filter based on metadata fields to get a subset cells, expression data, and metadata',
		'tags': ['cells'],
		'description': "Cells takes query parameters definied in the schema retrieved from the /initialize enpoint. <br>For categorical metadata keys filter based on `key=value` <br>"
					   " For continous  metadata keys filter by `key=min,max`<br> Either value can be replaced by a \*. To have only a minimum value `key=min,\*`  To have only a maximum value `key=\*,max` <br>Graph data (if retrieved) is normalized",
		'parameters': [{
			'name': '_nograph',
			'description': "Do not calculate and send back graph (graph is sent by default)",
			'in': 'path',
			'type': 'boolean',
		}],

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
		# TODO Allow random downsampling

		e = ExpressionMatrix(application.config['DATA_DIR'], True)
		data = {
			"reactive": False,
			"cellids": [],
			"metadata": [],
			"cellcount": 0,
			"badmetadatacount": 0,
			"ranges": {}
		}
		try:
			args = request.args
			nograph = False
			if "_nograph" in args:
				nograph = bool(args["_nograph"])
			qs = parse_querystring(request.args)
		except QueryStringError as e:
			return make_payload({}, str(e), 400)
		schema = parse_schema(os.path.join(application.config["SCHEMAS_DIR"], "test_data_schema.json"))
		metadata = parse_metadata(toCellIDsCpp())['cell_metadata']
		bad_metadata_count = 0
		keptcells = []
		run = randint(0, 10000)
		output_cellset = "out_{}".format(run)
		if len(qs):
			# TODO catch errors
			error = False
			filters = []
			for key, value in qs.items():
				if value["variabletype"] == 'categorical':
					category_filter = []
					for idx, item in enumerate(value["query"]):
						queryval = escape(item)
						filtername = "{}_{}_{}".format(key, run, idx)
						e.createCellSetUsingMetaData(filtername, key, queryval, False)
						category_filter.append(filtername)
					category_output_filter = "{}_out_{}".format(key, run)
					e.createCellSetUnion(",".join(category_filter), category_output_filter)
					filters.append(category_output_filter)
					for cellset in category_filter:
						e.removeCellSet(cellset)
				elif value["variabletype"] == 'continuous':
					filtername = "{}_{}".format(key, run)
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
			keptcells = parse_metadata(e.getCellSet(output_cellset))["cell_metadata"]
			for cellset in filters:
				e.removeCellSet(cellset)
		else:
			keptcells = metadata
			output_cellset = "AllCells"
		ranges = get_metadata_ranges(schema, keptcells)
		data["cellcount"] = len(keptcells)
		graph = None
		if not nograph:
			default_graph_params = {
				"similarpairsname": "ExactHighInformationGenes",
				"similaritythreshold": 0.3,
				"connectivity": 20,
			}
			os.chdir(application.config['DATA_DIR'])
			graphname = "cellgraph_{}".format(run)
			e.createCellGraph(graphname, output_cellset, default_graph_params["similarpairsname"],
							  default_graph_params["similaritythreshold"],
							  default_graph_params["connectivity"])
			e.computeCellGraphLayout(graphname)
			vertices = e.getCellGraphVertices(graphname)
			normalized_verticies = normalize_verticies(vertices)
			graph = [
				(e.getCellMetaDataValue(normalized_verticies["labels"][i], 'CellName'), normalized_verticies["x"][i],
				 normalized_verticies["y"][i])
				for i in range(len(normalized_verticies["labels"]))]
		if output_cellset != "AllCells":
			e.removeCellSet(output_cellset)

		if data["cellcount"] <= REACTIVE_LIMIT:
			data["reactive"] = True
			data["metadata"] = keptcells
			data["cellids"] = [m["CellName"] for m in keptcells]
			data["ranges"] = ranges
			# data["expression"] = parse_exp_data2(data["cellids"])
			data["badmetadatacount"] = bad_metadata_count
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
	def get(self):
		schema = parse_schema(os.path.join(application.config["SCHEMAS_DIR"], "test_data_schema.json"))
		metadata = parse_metadata(toCellIDsCpp())["cell_metadata"]

		options = get_metadata_ranges(schema, metadata)
		return make_payload({"schema": schema, "options": options, "cellcount": len(metadata)})


api.add_resource(MetadataAPI, "/api/v0.1/metadata")
api.add_resource(ExpressionAPI, "/api/v0.1/expression")
api.add_resource(GraphAPI, "/api/v0.1/graph")
api.add_resource(InitializeAPI, "/api/v0.1/initialize")
api.add_resource(CellsAPI, "/api/v0.1/cells")
api.add_resource(ClusterAPI, "/api/v0.1/cluster")

if __name__ == "__main__":
	application.run(host='0.0.0.0', debug=True)
