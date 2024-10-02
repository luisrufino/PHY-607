from sort_alg import bubble
from sort_alg import heap_sort
import numpy as np
import time
#import pandas as pd
import matplotlib.pyplot as plt

# sizes = []
# times = []
# for i in range (1,11):
#     size =  i*i*50
#     randomlist = np.random.randint(1, 1000, size=size)
#     start = time.time()
#     sorted = bubble(randomlist)
#     elapsedtime = (time.time() - start)
#     sizes.append(size)
#     times.append(elapsedtime)
#     print("The time it took to do bubble-sort of randomly sorted",size, "integers is",elapsedtime )
#
# print("Times array",times)
# print("Sizes array",sizes)
#
# values = {'n': sizes, 'Time results': times}
# df = pd.DataFrame(values)
# print(df)
#
# plot = df.plot(x='n', y='Time results', kind='line')
# plot.set_xlabel("size of array")
# plot.set_ylabel("time elapsed")
# plot.set_title('size vs time elapsed')
# plt.show()

t = []
del_t = []
size = np.arange(1,1000,10)
for sizes in size:
    print(f"Size: {sizes}")
    tmp_t = [] ## Contains the average time of bouble and heap for each size
    for run in range(10):
        r_list = np.random.randint(1,100, size = sizes)
        start = time.time()
        sort_list = bubble(r_list)
        end = time.time()
        tot_b = end - start

        start = time.time()
        sort_list = heap_sort(r_list)
        end = time.time()
        tot_heap = end - start
        tmp_time = [tot_b, tot_heap]
        tmp_t.append(tmp_time)
    tmp_t = np.array(tmp_t)
    avg_b = np.mean(tmp_t[:,0])
    abs_h = np.mean(tmp_t[:,1])
    avg_t = [avg_b, abs_h]
    std_b = np.std(tmp_t[:,0])
    std_h = np.std(tmp_t[:,1])
    std_t = [std_b, std_h]
    del_t.append(std_t)
    t.append(avg_t)
t = np.array(t)
del_t = np.array(del_t)
plt.errorbar(size, t[:,0], yerr=del_t[:,0], label ='bubble sort', fmt = '-', ecolor="red")
plt.errorbar(size, t[:,1], yerr=del_t[:,1], label ='heap sort', fmt = '-', ecolor="red")
plt.xlabel('Size of array')
plt.ylabel('Time of execuation')
plt.legend()
plt.title('Comnparison between bouble and heap sort')
plt.savefig("comparison.png")


