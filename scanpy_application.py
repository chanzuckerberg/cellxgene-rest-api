import json

from flask import (Flask, jsonify, send_from_directory,
                   request, make_response, render_template, Response, stream_with_context)
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

import scanpy_rest_api

# CONSTANTS
REACTIVE_LIMIT = 100000
INVALID_CELL_ID = 4294967295


scanpy_rest_api.initialize("/home/ubuntu/result.h5ad")

app = Flask(__name__, static_url_path='/templates')
Compress(app)
CORS(app)
cache = Cache(app, config={'CACHE_TYPE': 'simple'})

# SWAGGER
api = Api(app, api_version='0.1', produces=["application/json"], title="cellxgene rest api",
          api_spec_url='/api/swagger',
          description='An API connecting ExpressionMatrix2 clustering algorithm to cellxgene')

# ---- Traditional Routes -----
# CellxGene application
@app.route('/')
def index():
    url_base = app.config["CXG_API_BASE"]
    dataset_title = app.config["DATASET_TITLE"]
    return render_template("index.html", prefix=url_base, datasetTitle=dataset_title)


# renders swagger documentation
@app.route('/swagger')
def swag():
    return render_template("swagger.html")


@app.route('/favicon.png')
def icon():
    return send_from_directory('templates', 'favicon.png')

class CellsAPI(Resource):
    def get(self):
        return Response(scanpy_rest_api.generate_cells_get(request.args),
                        content_type="application/json")
                        

class InitializeAPI(Resource):
    def get(self):
        return Response(scanpy_rest_api.generate_initialize_get(), mimetype="application/json")


api.add_resource(InitializeAPI, "/api/v0.1/initialize")
api.add_resource(CellsAPI, "/api/v0.1/cells")

if __name__ == "__main__":
    app.run(host='0.0.0.0', port=4002)
