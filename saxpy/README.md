# SAXPY GPGPU Benchmark

SAXPY (Single Precision A * X Plus Y) is basically:
```python
  for i=0 to N:
     Y[i] = A * X[i] + Y[i]
```

This is a subset of the benchmarks from https://github.com/bennylp/saxpy-benchmark
with some modifications and an added OpenACC benchmark.

In order to compile most C++ and CUDA examples, first execute
```bash
module load NVHPC GCC/12.3.0 Python
```
then execute the following to install necessary Pyhton libraries
```bash
python -m pip install --user numpy pandas mxnet pycuda pyopencl
```

The `mxnet` and `pyopencl` versions of the benchmark currently do not work. Feel free to propose fixes.