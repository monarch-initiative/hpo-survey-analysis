#!/bin/sh
set -e

# Get the version of glibc and gcc and replace in dockerfile
GLIBC_VERSION=$(apt-cache policy libc6 | grep Installed | sed 's/  Installed: //')
GCC_VERSION=$(apt-cache policy gcc | grep Installed | sed 's/  Installed: //')
UBUNTU_RELEASE=$(lsb_release -rs)

sed s/{LIBC_VERSION}/$GLIBC_VERSION/ Dockerfile-base > Dockerfile
sed -i s/{GCC_VERSION}/$GCC_VERSION/ Dockerfile
sed -i s/{UBUNTU_RELEASE}/$UBUNTU_RELEASE/ Dockerfile

