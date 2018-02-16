import unittest
import requests
import os


class EndPoints(unittest.TestCase):
    """Test Case for endpoints"""

    def setUp(self):
        # Local
        self.url_base = os.environ["CXG_API_BASE"]
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
        result = self.session.post(url, data={"genelist": ["ABCD4", "AASS"]})
        assert result.status_code == 200
        result_json = result.json()
        assert 'data' in result_json
        assert len(result_json['data']['genes']) == 2
        assert len(result_json["data"]['cells'])
        assert len(result_json["data"]['cells'][0]['e']) == 2

    def test_heatmap_post_genelist_celllist(self):
        url = "{base}{endpoint}".format(base=self.url_base, endpoint="expression")
        result = self.session.post(url, data={"celllist": ["1001000012.C3"], "genelist": ["AASS"]})
        assert result.status_code == 200
        result_json = result.json()
        assert 'data' in result_json
        assert len(result_json['data']['genes']) == 1
        assert len(result_json["data"]['cells']) == 1
        assert len(result_json["data"]['cells'][0]['e']) == 1

    def test_initialize(self):
        url = "{base}{endpoint}".format(base=self.url_base, endpoint="initialize")
        result = self.session.get(url)
        assert result.status_code == 200
        result_json = result.json()

        assert result_json["data"]["cellcount"] > 0

    def test_cells(self):
        url = "{base}{endpoint}".format(base=self.url_base, endpoint="cells")
        result = self.session.get(url)
        assert result.status_code == 200
        url = "{base}{endpoint}?{params}".format(base=self.url_base, endpoint="cells", params="&".join(
            ["Selection=Unpanned", "ERCC_reads=150000,170000"]))
        result = self.session.get(url)
        assert result.status_code == 200
        result_json = result.json()
        assert result_json["data"]["cellcount"] > 0
        assert len(result_json["data"]["graph"])
        assert len(result_json["data"]["graph"]) == result_json["data"]["cellcount"]

    def test_cells_nograph(self):
        url = "{base}{endpoint}?{params}".format(base=self.url_base, endpoint="cells", params="&".join(
            ["Selection=Unpanned", "ERCC_reads=150000,170000", "_nograph=True"]))
        result = self.session.get(url)
        result_json = result.json()
        assert result.status_code == 200
        assert result_json["data"]["graph"] is None

    def test_cells_failure(self):
        url = "{base}{endpoint}?{params}".format(base=self.url_base, endpoint="cells", params="&".join(
            ["Selection=Unpanned", "ERCC_reads=a,160000"]))
        result = self.session.get(url)
        assert result.status_code == 400
        url = "{base}{endpoint}?{params}".format(base=self.url_base, endpoint="cells", params="&".join(
            ["sadfkljds=Unpanned", "ERCC_reads=150000,160000"]))
        result = self.session.get(url)
        assert result.status_code == 400
        url = "{base}{endpoint}?{params}".format(base=self.url_base, endpoint="cells", params="&".join(
            ["Selection=Unpanned", "ERCC_reads=150000"]))
        result = self.session.get(url)
        assert result.status_code == 400

    def test_diff_expression(self):
        url = "{base}{endpoint}".format(base=self.url_base, endpoint="diffexpression")
        result = self.session.post(url, data={
            "celllist1": [
                "1001000010.C5",
                "1001000010.C8",
                "1001000010.D5",
                "1001000010.D9",
                "1001000010.E8",
                "1001000010.F1",
                "1001000010.F10",
                "1001000010.F7",
                "1001000010.F8",
                "1001000010.G2",
                "1001000010.G3",
                "1001000010.H3",
                "1001000012.A1",
                "1001000012.A1",
                "1001000012.A10"
            ],
            "celllist2": [
                "1001000012.H4",
                "1001000012.H5",
                "1001000012.H6",
                "1001000012.H7",
                "1001000012.H8",
                "1001000012.H9",
                "1001000014.B2",
                "1001000014.C1",
                "1001000014.C3"
            ],
            "pval": 0.01,
            "num_genes": 10
        })
        assert result.status_code == 200
        result_json = result.json()
        assert 'data' in result_json
        assert len(result_json['data']["celllist1"]['topgenes']) > 0
        result = self.session.post(url, data={"clusters1": ["6"], "clusters2": ["3"]})
        assert result.status_code == 200
        result_json = result.json()
        assert 'data' in result_json
        assert len(result_json['data']["celllist1"]['topgenes']) > 0
