# Dockerfile for Debian GCC8 env April 2022
FROM debian:buster 

RUN \
apt-get update && \
apt-get install -y build-essential git xsltproc \
python2 libpython3.7 libpython3.7-dev python3-pip python3-venv python3-mysqldb python3-tk python3-mysqldb-dbg \
man tcllib environment-modules libfabric1 libfabric-dev gcc-8 libblas-dev liblapack-dev libopenmpi3

# General python packages required by test suite, excitingtools and Jupyter tutorials
RUN \
pip3 install --upgrade --force pillow==7.0.0 wheel>=0.35.0 numpy>=1.14.5 matplotlib>=2.2.0   

# Test suite dependencies
RUN \
pip3 install termcolor lxml pytest pyyaml
