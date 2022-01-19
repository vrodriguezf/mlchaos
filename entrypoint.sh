#!/bin/bash --login
set -e
pip install -e ~/work
exec "$@"
