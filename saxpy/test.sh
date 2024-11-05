#!/bin/bash

make

echo "Most SAXPY benchmarks from https://github.com/bennylp/saxpy-benchmark"

echo ""
echo "Python loop"
python saxpy_loop.py

echo ""
echo "CPU"
./saxpy_cpu

echo ""
echo "Python w/ Numpy"
python saxpy_numpy.py

#echo ""
#echo "Python w/ mxnet"
#python saxpy_mxnet.py

echo ""
echo "Python w/ Pandas"
python saxpy_pandas.py

echo ""
echo "Python w/ PyCUDA"
python saxpy_pycuda.py

echo ""
echo "Python w/ PyOpenCL"
python saxpy_pyocl.py

echo ""
echo "OpenMP"
./saxpy_omp

echo ""
echo "OpenCL"
./saxpy_ocl

echo ""
echo ""
echo "OpenACC"
./saxpy_acc

echo ""
echo "CUDA"
./saxpy_cuda

echo ""
echo "CUBLAS"
./saxpy_cublas

