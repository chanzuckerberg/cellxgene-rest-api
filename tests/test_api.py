import unittest
import requests
import json


class MetaData(unittest.TestCase):
    """Test Case for metadata endpoint"""

    def setUp(self):
        self.url_base = "http://localhost:5000/api/v0.1/"
        self.session = requests.Session()

    def test_get(self):
        url = "{base}{endpoint}".format(base=self.url_base, endpoint="metadata")
        result = self.session.get(url)
        assert result.status_code == 200
        assert len(result.text)

    def test_get_json(self):
        url = "{base}{endpoint}?{params}".format(base=self.url_base, endpoint="metadata", params="format=json")
        result = self.session.get(url)
        assert result.status_code == 200
        result_json = result.json()
        assert 'data' in result_json
        assert len(result_json['data']['cell_metadata'])

    def test_cell_list(self):
        url = "{base}{endpoint}".format(base=self.url_base, endpoint="metadata")
        result = self.session.post(url, data={"celllist": ["1001000012.C3"]})
        assert result.status_code == 200
        result_json = result.json()
        assert 'data' in result_json
        assert len(result_json['data']['cell_metadata'])
        assert result_json['data']['cell_metadata'][0]["Cell"] == "1001000012.C3"

    def test_cell_list_error(self):
        url = "{base}{endpoint}".format(base=self.url_base, endpoint="metadata")
        result = self.session.post(url, data={"celllist": ["A"]})
        assert result.status_code == 200
        result_json = result.json()
        assert result_json['status']['error']
