import sys, os
sys.path.insert(0, '/var/www/html/clustering_api_layer')

def execfile(filename):
    globals = dict( __file__ = filename )
    exec( open(filename).read(), globals )

activate_this = os.path.join( '/home/ubuntu/clustering_api_layer', 'venv/bin', 'activate_this.py' )
execfile( activate_this )

from application import application as application
