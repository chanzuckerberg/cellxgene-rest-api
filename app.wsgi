import sys, os
sys.path.insert(0, '/var/www/html/clustering_flask_api')

def execfile(filename):
    globals = dict( __file__ = filename )
    exec( open(filename).read(), globals )

activate_this = os.path.join( '/home/ubuntu/clustering_flask_api', 'venv/bin', 'activate_this.py' )
execfile( activate_this )

def application(environ, start_response):
    for key in ['SECRET_KEY']:
        os.environ[key] = environ.get(key, '')
    from application import application as _application

    return _application(environ, start_response)
