#! /bin/bash

echo "cleaning PFNET ..."
find . -name \*~ -delete
find . -name aclocal.m4 -delete
find . -name autom4te.cache -type d -exec rm -rf {} +
find . -name compile -delete
find . -name config.guess -delete
find . -name config.sub -delete
find . -name configure -delete
find . -name depcomp -delete
find . -name install.sh -delete
find . -name ltmain.sh -delete
find . -name m4 -type d -exec rm -rf {} +
find . -name Makefile.in -delete
find . -name missing -delete
find . -name test-driver -delete
rm -r build
