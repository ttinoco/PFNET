echo "cleaning..."
find ./pfnet -name \*.so -delete
find ./pfnet -name \*.pyc -delete
rm -f ./pfnet/cpfnet.c
rm -rf build
rm -rf lib/pfnet
