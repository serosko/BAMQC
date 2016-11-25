#!/bin/bash
while [[ $# > 1 ]]
do
key="$1"

case $key in
    -v|--version)
    GCC_VERSION="$2"
    shift
    ;;
    -r|--repository)
	REPOSITORY="$2"
	shift
    ;;
    *)
    # unknown option
    ;;
esac
shift
done

if [ -z "${GCC_VERSION}" ]; then
	GCC_VERSION="4.9"
fi

if [ -z "${REPOSITORY}" ]; then
	REPOSITORY="ubuntu-toolchain-r/test"
fi

set -e

sudo add-apt-repository -y ppa:${REPOSITORY}
sudo apt-get update -qq

sudo apt-get install -qq g++-${GCC_VERSION}
sudo apt-get install -qq gcc-${GCC_VERSION}
sudo ln -s -v -f $(which g++-${GCC_VERSION}) /usr/bin/g++
sudo ln -s -v -f $(which gcc-${GCC_VERSION}) /usr/bin/gcc

set +e
