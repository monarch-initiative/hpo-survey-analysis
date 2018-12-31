#!/bin/sh
set -e

# Create a directory to store the c binary and python source files
mkdir redland-bindings

# Get the gnupg key
wget https://www.dajobe.org/gnupg.asc
gpg --batch --import gnupg.asc

# Download raptor
wget http://download.librdf.org/source/raptor2-2.0.15.tar.gz
wget http://download.librdf.org/source/raptor2-2.0.15.tar.gz.asc
gpg --no-default-keyring --verify raptor2-2.0.15.tar.gz.asc raptor2-2.0.15.tar.gz
tar -zxf raptor2-2.0.15.tar.gz && cd raptor2-2.0.15
./configure --prefix=/usr/local/
make install && cd /build

# Download rasqal
wget http://download.librdf.org/source/rasqal-0.9.33.tar.gz
wget http://download.librdf.org/source/rasqal-0.9.33.tar.gz.asc
gpg --no-default-keyring --verify rasqal-0.9.33.tar.gz.asc rasqal-0.9.33.tar.gz
tar -zxf rasqal-0.9.33.tar.gz && cd rasqal-0.9.33
./configure --prefix=/usr/local/
make install && cd /build

# Download redland
wget http://download.librdf.org/source/redland-1.0.17.tar.gz
wget http://download.librdf.org/source/redland-1.0.17.tar.gz.asc
gpg --no-default-keyring --verify redland-1.0.17.tar.gz.asc redland-1.0.17.tar.gz
tar -zxf redland-1.0.17.tar.gz && cd redland-1.0.17
./configure --prefix=/usr/local/
make install &&  cd /build

# Download bindings
wget http://download.librdf.org/source/redland-bindings-1.0.17.1.tar.gz
wget http://download.librdf.org/source/redland-bindings-1.0.17.1.tar.gz.asc
gpg --no-default-keyring --verify redland-bindings-1.0.17.1.tar.gz.asc redland-bindings-1.0.17.1.tar.gz
tar -zxf redland-bindings-1.0.17.1.tar.gz && cd redland-bindings-1.0.17.1
./configure --with-python=python3 && cd python
make install DESTDIR=/build && cd /build

cp /build/usr/local/lib/python3.7/site-packages/RDF.py redland-bindings/
cp /build/usr/local/lib/python3.7/site-packages/*Redland* redland-bindings/

tar -zcf redland-bindings.tar.gz redland-bindings
