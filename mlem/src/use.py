import numpy as np
from cffi import FFI

ffi = FFI()
ffi.cdef("float sum_floats(float*, size_t);")
mlem = ffi.dlopen("../../target/debug/libmlem.so")

def sumit(data):
    array = np.array(data, dtype=np.float32)
    pointer = ffi.cast('float *', array.ctypes.data)
    return mlem.sum_floats(pointer, len(array))


data = [1,2,3,4]
total = sumit(data)
print(f'Sum of {data} = {total}')
assert total == 10.0
