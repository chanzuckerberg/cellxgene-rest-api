import sys, os, time
from collections import defaultdict
# import atexit

from flask import Flask, jsonify, send_from_directory, request, make_response, render_template
from flask_restful import reqparse
from flask_restful_swagger_2 import Api, swagger, Resource
from flask_cors import CORS
import numpy as np
from scipy import stats

from schemaparse import parse_schema

# CONSTANTS
REACTIVE_LIMIT = 100000
INVALID_CELL_ID = 4294967295

# Local Errors
class QueryStringError(Exception):
	pass

application = Flask(__name__, static_url_path='/templates')
CORS(application)

# SECRETS
SECRET_KEY = os.environ.get("SECRET_KEY", default=None)
if not SECRET_KEY:
	raise ValueError("No secret key set for Flask application")
CXG_API_BASE = os.environ.get("CXG_API_BASE", default=None)
if not CXG_API_BASE:
	raise ValueError("No api base")

# CONFIG
application.config.from_pyfile('app-human-mtg.cfg', silent=True)
application.config.update(
	CXG_API_BASE=CXG_API_BASE,
	SECRET_KEY=SECRET_KEY
)
dir_path = os.path.dirname(os.path.realpath(__file__))
application.config.update(
	SCRATCH_DIR=os.path.join(dir_path, application.config["SCRATCH_DIR"]),
	DATA_DIR=os.path.join(dir_path, application.config["DATA_DIR"]),
	EM2_DIR=os.path.join(dir_path, application.config["EM2_DIR"]),
)

# SWAGGER
api = Api(application, api_version='0.1', produces=["application/json"], title="cellxgene rest api",
		  api_spec_url='/api/swagger',
		  description='An API connecting ExpressionMatrix2 clustering algorithm to cellxgene')

# Load expression matrix libray only after location path is configured
sys.path.insert(0, application.config["EM2_DIR"])
from ExpressionMatrix2 import ExpressionMatrix

# Param parsers
graph_parser = reqparse.RequestParser()
graph_parser.add_argument('cellsetname', type=str, default='AllCells', required=False, help="Named cell set")
graph_parser.add_argument('similarpairsname', type=str, default='LshHighInfo', required=False,
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

differential_parser = reqparse.RequestParser()
differential_parser.add_argument('celllist1', type=str, action="append", required=False, help='First group of cells by id')
differential_parser.add_argument('celllist2', type=str, action="append", required=False, help='Second group of cells by id')
differential_parser.add_argument('clusters1', type=str, action="append", required=False, help='First group of cluters by name')
differential_parser.add_argument('clusters2', type=str, action="append", required=False, help='Second group of clusters by name')
differential_parser.add_argument('num_genes', type=int, required=False, default=20,
						  help="Number of diff-expressed genes to return")

e = ExpressionMatrix(application.config["DATA_DIR"], True)
schema = parse_schema(os.path.join(application.config["DATA_DIR"], application.config["SCHEMA_FILE"]))

# @atexit.register
# def goodbye():
# 	print("You are now leaving the Python sector.")

# ---- Helper Functions -------
def filter(cell, qs):
	"""
	Check if cell's metadata are within the filtered parameters
	:param cell: cell id
	:param qs: query string
	:return: Boolean True if cell within filter, False if rejected by filter
	"""
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
	metadata = []
	if cells == False:
		cells = e.getCellSet('AllCells')
	if INVALID_CELL_ID not in cells:
		mdata = e.getCellsMetaData(cells)
		for i in mdata:
			metadata.append({k:cast_value(k,v) for k,v in i})
	return metadata


def cast_value(key, value):
	"""
	Cast metadata value to proper type
	:param key: string key
	:param value: string metadata value
	:return:
	"""
	val_type = schema[key]["type"]
	new_val = value
	try:
		if val_type == "int":
			new_val = int(value)
		elif val_type == "float":
			new_val = float(value)
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
	expression = get_expression(cells, genes)
	if not genes:
		genes = all_genes()
	if not cells:
		cells = all_cells()
	if not unexpressed_genes:
		expression, genes = remove_unexpressed_genes(expression, genes)
	if limit and len(genes) > limit:
		genes = genes[:limit]
		expression = expression[:,:limit]
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


def get_expression(cells, genes = ()):
	"""
	Get matrix of expression data
	:param cells: list of cell names, or all cells if empty
	:param genes: list of genes, or all genes if empty
	:return: numpy expression matrix
	"""
	# TODO speed this up
	if not cells:
		cellset = "AllCells"
	else:
		cell_ids = [get_cell_id(name) for name in cells]
		cellset = "ExpressionCellSet"
		e.createCellSet(cellset, cell_ids)
	if not genes:
		e.getDenseExpressionMatrix("AllGenes", cellset)
	else:
		genesetName = "ExpressionGeneSet_{}".format(time.time())
		e.createGeneSetFromGeneNames(genesetName, genes)
		expression = e.getDenseExpressionMatrix(genesetName, cellset)
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
	:return: dictionary {categorical_key: { "options": [list of options]}, continous_key: {"range": {"min": min_val, "max": max_val}}}
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
						if (s == "EM2Cluster") and datum.startswith("Unclustered"):
							continue
						options[s]["options"][datum] += 1
					elif schema[s]["variabletype"] == "continuous":
						if not options[s]["range"]["min"] or datum < options[s]["range"]["min"]:
							options[s]["range"]["min"] = datum
						if not options[s]["range"]["max"] or datum > options[s]["range"]["max"]:
							options[s]["range"]["max"] = datum
				except KeyError as e:
					pass
				except ValueError:
					pass
				except TypeError:
					pass
	for key, value in dict(options).items():
		if ("options" in value and not value["options"]) or ("range" in value and not value["range"]["max"]):
			del options[key]

	return options


def normalize_verticies(verticies):
	"""
	normalize x,y verticies for graph
	:param verticies: list of vertex objects
	:return: json: {"labels": list of point labels, "x": list of normalized x positions, "y": list of normalized y positions}
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
	norm = np.linalg.norm(vec, ord=np.inf)
	if norm == 0:
		return vec
	return vec / norm


def all_genes():
	"""
	Get list of all genenames in dataset
	:return: list of gene names
	"""
	return [e.geneName(i) for i in range(e.geneCount())]


def all_cells():
	"""
	Get list of all cell names in dataset
	:return: list of cell names
	"""
	return [get_cell_name(i) for i in e.getCellSet('AllCells')]


def diffexp(expression_1, expression_2, num_genes=20):
	"""
	get top expressed genes from two different cell sets (uses t-test)
	:param expression_1: numpy expression array cell set 1
	:param expression_2: numpy expression array cell set 2
	:param num_genes: number of cells to return
	:return: Top genes and mean expression for each gene in both cell sets
	"""
	# TODO make sure genes aren't nans
	diff_exp = stats.ttest_ind(expression_1, expression_2)
	set1 = np.argsort(diff_exp.statistic)[::-1]
	set1 = np.roll(set1, -np.count_nonzero(np.isnan(diff_exp.statistic)))[:num_genes]
	set2 = np.argsort(diff_exp.statistic)[:num_genes]
	mean_ex1_set1 = np.mean(expression_1[:, set1], axis=0)
	mean_ex2_set1 = np.mean(expression_2[:, set1], axis=0)
	mean_ex1_set2 = np.mean(expression_1[:, set2], axis=0)
	mean_ex2_set2 = np.mean(expression_2[:, set2], axis=0)


	genes_cellset_1 = [e.geneName(gid) for gid in set1]
	genes_cellset_2 = [e.geneName(gid) for gid in set2]
	return {
		"celllist1": {
			"topgenes": genes_cellset_1,
			"mean_expression_cellset1": mean_ex1_set1.tolist(),
			"mean_expression_cellset2": mean_ex2_set1.tolist(),
		},
		"celllist2": {
			"topgenes": genes_cellset_2,
			"mean_expression_cellset1": mean_ex1_set2.tolist(),
			"mean_expression_cellset2":mean_ex2_set2.tolist(),
		},
	}


# TODO do i need this one?
def convert_variable(datatype, variable):
	"""
	Convert variable to number (float/int)
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
		print("Bad conversion")
		raise


def get_clusters(list_of_clusters):
	"""
	Get list of cells that are members of the cluster(s)
	:param list_of_clusters: list of strings cluster names
	:return: list of cell names
	"""
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
	return e.getCellMetaDataValue(cellid, "CellName")

def get_cell_id(cellname):
	"""
	Get cell id from cell name
	:param cellname: string cell name
	:return: int id of cell
	"""
	return e.cellIdFromString(cellname)


# ---- Traditional Routes -----
# CellxGene application
@application.route('/')
def index():
	url_base = application.config["CXG_API_BASE"]
	split_base = url_base.split("/")
	prefix = '/'.join(split_base[:-2]) + "/"
	version = "{}/".format(split_base[-2])
	return render_template("index.html", prefix=prefix, version=version)


# renders swagger documentation
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
		#print("get expression args", args)
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
		#print("post expression args", args)
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
# 	def get(self):
# 		# Retrieve cluster
# 		run = randint(0, 10000)
# 		graphname = "cellgraph_{}".format(run)
# 		e = ExpressionMatrix(application.config['DATA_DIR'], True)
# 		args = graph_parser.parse_args()
# 		args = graph_parser.parse_args()
# 		e.createCellGraph(graphname, args.cellsetname, args.similarpairsname, args.similaritythreshold,
# 						  args.connectivity)
# 		clustername = "clusters_{}".format(run)
# 		params = ClusterGraphCreationParameters()
# 		e.createClusterGraph(graphname, params, clustername)
# 		clusters = e.getClusterGraphVertices(clustername)
# 		clusters = [i for i in clusters]
# 		cells = {}
# 		for i in clusters:
# 			cells[i] = [e.getCellMetaDataValue(cid, "CellName") for cid in e.getClusterCells(clustername, i)]
# 		# TODO calculate expression values for these
# 		# TODO does the graph really need to be created each time?
# 		# TODO parameterize graph with cell sets?
#
# 		return make_payload({clustername: clustername}, errorcode=201)


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
					"example": {"clusters1": ["6"], "clusters2": ["3"]}
				}
			}
		],
		"responses": {
			'200': {
				'description': 'top expressed genes for cellset1, cellset2',
				'examples': {
					'application/json': {
						"data":  {
							"celllist1": {
								"mean_expression_cellset1": [
									3339.285714285714,
									4806.642857142857,
									241.07142857142858,
									619.0714285714286,
									203.42857142857142
								],
								"mean_expression_cellset2": [
									93.55555555555556,
									1667.3333333333333,
									68.22222222222223,
									2.2222222222222223,
									0
								],
								"topgenes": [
									"SAT1",
									"FTL",
									"RPS16",
									"SRGN",
									"CTSS"
								]
							},
							"celllist2": {
								"mean_expression_cellset1": [
									95.78571428571429,
									47.642857142857146,
									54.357142857142854,
									47.142857142857146,
									500.35714285714283
								],
								"mean_expression_cellset2": [
									1523.5555555555557,
									203.44444444444446,
									752.8888888888889,
									240.22222222222223,
									3110.222222222222
								],
								"topgenes": [
									"CLU",
									"COX7C",
									"CKB",
									"LOC550643",
									"TUBA1A"
								]
							}
						},
						"status": {
							"error": False,
							"errormessage": "",
						}
					}
				}
			},
		}
	})
	def post(self):
		args = self.parser.parse_args()
		# TODO allow different calculation methods
		# get 2 cell sets
		cell_list_1 = args.get('celllist1', [])
		cell_list_2 = args.get('celllist2', [])
		cluster_list_1 = args.get('clusters1', [])
		cluster_list_2 = args.get('clusters2', [])
		num_genes = args.get("num_genes")
		if cluster_list_1 and cluster_list_2:
			cell_list_1 = get_clusters(cluster_list_1)
			cell_list_2 = get_clusters(cluster_list_2)
		if not ((cell_list_1 and cell_list_2)):
			return make_payload([], "must include either (cellist1 and cellist2) or (clusters1 and clusters2) parameters", 400)
		# TODO else get reverse
		expression_1 = get_expression(cell_list_1, all_genes())
		expression_2 = get_expression(cell_list_2, all_genes())
		# expression_1_average = np.mean(expression_1, axis=1)
		# expression_2_average = np.mean(expression_2, axis=1)
		# t-test between cell sets
		# TODO make sure genes aren't nans
		data = diffexp(expression_1, expression_2, num_genes)
		# TODO also return statistic value
		# return top expressed genes
		return make_payload(data)


class CellsAPI(Resource):
	@swagger.doc({
		'summary': 'filter based on metadata fields to get a subset cells, expression data, and metadata',
		'tags': ['cells'],
		'description': "Cells takes query parameters definied in the schema retrieved from the /initialize enpoint. <br>For categorical metadata keys filter based on `key=value` <br>"
					   " For continous  metadata keys filter by `key=min,max`<br> Either value can be replaced by a \*. To have only a minimum value `key=min,\*`  To have only a maximum value `key=\*,max` <br>Graph data (if retrieved) is normalized"
					   " To only retrieve cells that don't have a value for the key filter by `key`",
		'parameters': [{
			'name': '_nograph',
			'description': "Do not calculate and send back graph (graph is sent by default)",
			'in': 'path',
			'type': 'boolean',
		},
			{
				'name': '_includeisolated',
				'description': "Include all cells in graph even if cluster is too small, this does nothing if _nograph is true.",
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
		try:
			args = request.args
			#print("get cell args", args)
			nograph = False
			includeisolated = False
			if "_nograph" in args:
				nograph = bool(args["_nograph"])

			if "_includeisolated" in args:
				includeisolated = bool(args["_includeisolated"])
			qs = parse_querystring(request.args)
		except QueryStringError as e:
			return make_payload({}, str(e), 400)
		keptcells = []
		output_cellset = "outcellset2"
		if len(qs):
			# TODO catch errors
			error = False
			filters = []
			for key, value in qs.items():
				if value["variabletype"] == 'categorical':
					category_filter = []
					if value["query"] == [""]:
						filtername = "{}".format(key)
						e.createCellSetUsingMetaData(filtername, key, "^$", True)
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
			cellidlist = e.getCellSet(output_cellset)
			for cellset in filters:
				e.removeCellSet(cellset)
		else:
			output_cellset = "AllCells"
			cellidlist = e.getCellSet(output_cellset)
		try:
			graph = None
			if not nograph and len(cellidlist) <= REACTIVE_LIMIT:
				graphname = "cellgraph"
				e.createCellGraph(graphname, output_cellset, 'LshHighInfo', keepIsolatedVertices=includeisolated)
				e.computeCellGraphLayout(graphname)
				vertices = e.getCellGraphVertices(graphname)
				normalized_verticies = normalize_verticies(vertices)
				graph = []
				# reset cellids if create the graph
				cellidlist = []
				for i in range(len(normalized_verticies["labels"])):
					graph.append((get_cell_name(i), normalized_verticies["x"][i],
								  normalized_verticies["y"][i]))
					cellidlist.append(normalized_verticies["labels"][i])

			keptcells = parse_metadata(cellidlist)
			ranges = get_metadata_ranges(keptcells)
			data["ranges"] = ranges
			data["cellcount"] = len(keptcells)

		finally:
			if output_cellset != "AllCells":
				e.removeCellSet(output_cellset)
			data["reactive"] = True
			data["metadata"] = keptcells
			data["cellids"] = [m["CellName"] for m in keptcells]
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
		metadata = parse_metadata()
		options = get_metadata_ranges(metadata)
		return make_payload({"schema": schema, "options": options, "cellcount": len(metadata)})


api.add_resource(MetadataAPI, "/api/v0.1/metadata")
api.add_resource(ExpressionAPI, "/api/v0.1/expression")
api.add_resource(InitializeAPI, "/api/v0.1/initialize")
api.add_resource(CellsAPI, "/api/v0.1/cells")
api.add_resource(DifferentialExpressionAPI, "/api/v0.1/diffexpression")

if __name__ == "__main__":
	application.run(host='0.0.0.0', debug=True)
