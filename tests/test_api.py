import unittest
import requests


class EndPoints(unittest.TestCase):
	"""Test Case for endpoints"""

	def setUp(self):
		# Local
		self.url_base = "http://localhost:5000/api/v0.1/"
		# Dev
		# self.url_base = "http://api-dev.clustering.czi.technology/api/v0.1/"
		# Prod
		# self.url_base = "http://api.clustering.czi.technology/api/v0.1/"
		self.session = requests.Session()

	def test_cell_list(self):
		url = "{base}{endpoint}".format(base=self.url_base, endpoint="metadata")
		result = self.session.post(url, data={"celllist": ["1001000012.C3"]})
		assert result.status_code == 200
		result_json = result.json()
		assert 'data' in result_json
		assert len(result_json['data']['cell_metadata'])
		assert result_json['data']['cell_metadata'][0]["CellName"] == "1001000012.C3"

	def test_cell_list_error(self):
		url = "{base}{endpoint}".format(base=self.url_base, endpoint="metadata")
		result = self.session.post(url, data={"celllist": ["A"]})
		assert result.status_code == 400
		result_json = result.json()
		assert result_json['status']['error']

	def test_heatmap(self):
		url = "{base}{endpoint}".format(base=self.url_base, endpoint="expression")
		result = self.session.get(url)
		assert result.status_code == 200
		result_json = result.json()
		assert not result_json['status']['error']
		assert len(result_json['data']['genes'])
		assert len(result_json["data"]['cells'])
		assert len(result_json["data"]['cells'][0]['e'])

	def test_heatmap_post_celllist(self):
		url = "{base}{endpoint}".format(base=self.url_base, endpoint="expression")
		result = self.session.post(url, data={"celllist": ["1001000012.C3"]})
		assert result.status_code == 200
		result_json = result.json()
		assert 'data' in result_json
		assert len(result_json['data']['genes'])
		assert len(result_json["data"]['cells']) == 1
		assert len(result_json["data"]['cells'][0]['e'])
		assert result_json["data"]['nonzero_gene_count']
		result = self.session.post(url, data={"celllist": ["1001000012.C3"], "include_unexpressed_genes": True})
		assert result.status_code == 200
		result_json_unexpressed = result.json()
		assert len(result_json["data"]['genes']) < len(result_json_unexpressed["data"]['genes'])


	def test_heatmap_post_genelist(self):
		url = "{base}{endpoint}".format(base=self.url_base, endpoint="expression")
		result = self.session.post(url, data={"genelist": ["ABCD4"]})
		assert result.status_code == 200
		result_json = result.json()
		assert 'data' in result_json
		assert len(result_json['data']['genes']) == 1
		assert len(result_json["data"]['cells'])
		assert len(result_json["data"]['cells'][0]['e']) == 1

	def test_heatmap_post_genelist_celllist(self):
		url = "{base}{endpoint}".format(base=self.url_base, endpoint="expression")
		result = self.session.post(url, data={"celllist": ["1001000012.C3"], "genelist": ["AASS"]})
		print(result.status_code)
		assert result.status_code == 200
		result_json = result.json()
		assert 'data' in result_json
		assert len(result_json['data']['genes']) == 1
		assert len(result_json["data"]['cells']) == 1
		assert len(result_json["data"]['cells'][0]['e']) == 1

	def test_graph_endpoint(self):
		url = "{base}{endpoint}?{params}".format(base=self.url_base, endpoint="graph", params="&".join(
			["cellsetname=AllCells", "similarpairsname=ExactHighInformationGenes", "similaritythreshold=0.3",
			 "connectivity=20"]))
		result = self.session.get(url)
		assert result.status_code == 200

	def test_initialize(self):
		url = "{base}{endpoint}".format(base=self.url_base, endpoint="initialize")
		result = self.session.get(url)
		assert result.status_code == 200
		result_json = result.json()

		assert result_json["data"]["cellcount"] > 0

	def test_cells(self):
		url = "{base}{endpoint}?{params}".format(base=self.url_base, endpoint="cells", params="&".join(
			["Class[]=Neoplastic", "ERCC_reads=150000,160000"]))
		result = self.session.get(url)
		assert result.status_code == 200
		result_json = result.json()
		assert result_json["data"]["cellcount"] > 0
	
	def test_cells_failure(self):
		url = "{base}{endpoint}?{params}".format(base=self.url_base, endpoint="cells", params="&".join(
			["Class[]=Neoplastic", "ERCC_reads=a,160000"]))
		result = self.session.get(url)
		assert result.status_code == 400
		url = "{base}{endpoint}?{params}".format(base=self.url_base, endpoint="cells", params="&".join(
			["sadfkljds[]=Neoplastic", "ERCC_reads=150000,160000"]))
		result = self.session.get(url)
		assert result.status_code == 400
		url = "{base}{endpoint}?{params}".format(base=self.url_base, endpoint="cells", params="&".join(
			["Class[]=Neoplastic", "ERCC_reads=150000"]))
		result = self.session.get(url)
		assert result.status_code == 400
