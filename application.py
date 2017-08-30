import csv

from functools import wraps
from flask import Flask, jsonify, redirect, url_for, send_file, request
from flask_restful import Api, Resource, reqparse
from flask_restful_swagger import swagger
from flask_dance.contrib.google import make_google_blueprint, google
from flask_cors import CORS


application = Flask(__name__)
CORS(application)

application.config.from_pyfile('app.cfg', silent=True)
api = swagger.docs(Api(application), apiVersion="1.0")
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
            metadata.append({k:v for k,v in zip(header,row)})
    return {"cell_metadata": metadata}


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

    def get(self):
        args = self.parser.parse_args()
        if args.format and args.format.lower() == 'json':
            metadata = parse_metadata()
            return make_payload(metadata)
        return send_file(application.config["GBM_DIR"] + "GBM_metadata.csv", mimetype="text/csv")

    def post(self):
        args = self.parser.parse_args()
        cell_list = args.celllist
        metadata = parse_metadata(cell_list)
        if len(metadata['cell_metadata']) < len(cell_list):
            return make_payload([], "Some cell ids not available")
        return make_payload(metadata)


api.add_resource(MetadataAPI, "/api/v0.1/metadata")

if __name__ == "__main__":
    application.run(debug=True)
