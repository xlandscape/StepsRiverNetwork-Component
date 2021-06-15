
from numba import vectorize, float64
from numba import cfunc
import numpy as np
import time

i = 100000000
x = np.arange(i).reshape(10, int(i/10))
x = x.astype(np.float64)
@vectorize([float64(float64, float64)])
def f(x, y):
    return x + y


@cfunc("float64(float64, float64)")
def add(x, y):
    return x + y


# DO NOT REPORT THIS... COMPILATION TIME IS INCLUDED IN THE EXECUTION TIME!
start = time.time()
f(x,x)
end = time.time()
print("Elapsed (with compilation) = %s" % (end - start))
# NOW THE FUNCTION IS COMPILED, RE-TIME IT EXECUTING FROM CACHE
start = time.time()
f(x,x)
end = time.time()
print("Elapsed (after compilation) = %s" % (end - start))

#start = time.time()
#add.ctypes(x,x)
#end = time.time()
#print("Elapsed (after compilation) = %s" % (end - start))
#start = time.time()
#add.ctypes(x,x)
#end = time.time()
#print("Elapsed (after compilation) = %s" % (end - start))

start = time.time()
y = x+x
end = time.time()
print("Elapsed (after compilation) = %s" % (end - start))