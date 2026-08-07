# server.py
import http.server
import socketserver
from urllib.parse import urlparse
import os

PORT = 8000
BASE_DIR = os.path.expanduser("~/Desktop/reax_tools")

class RewriteHandler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        self.directory = BASE_DIR
        super().__init__(*args, **kwargs)
    
    def translate_path(self, path):
        if path.startswith('/reax_tools/'):
            path = path[12:]
        elif path == '/reax_tools':
            path = '/'
        return super().translate_path(path)

with socketserver.TCPServer(("", PORT), RewriteHandler) as httpd:
    print(f"Serving {BASE_DIR} at http://localhost:{PORT}")
    print(f"http://localhost:{PORT}/")
    httpd.serve_forever()
