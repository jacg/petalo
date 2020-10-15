from cffi import FFI
ffi = FFI()
ffi.cdef("int square(int);")
C = ffi.dlopen("../../target/debug/libcffi.so")
for i in range(10):
    print(C.square(i))
