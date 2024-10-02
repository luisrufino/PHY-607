from sort_alg import bubble
import numpy as np
import time
import pandas as pd
import matplotlib.pyplot as plt

sizes = []
times = []
for i in range (1,11):
    size =  i*i*50
    randomlist = np.random.randint(1, 1000, size=size)
    start = time.time()
    sorted = bubble(randomlist)
    elapsedtime = (time.time() - start)
    sizes.append(size)
    times.append(elapsedtime)
    print("The time it took to do bubble-sort of randomly sorted",size, "integers is",elapsedtime )

print("Times array",times)
print("Sizes array",sizes)

values = {'n': sizes, 'Time results': times}
df = pd.DataFrame(values)
print(df)

plot = df.plot(x='n', y='Time results', kind='line')
plot.set_xlabel("size of array")
plot.set_ylabel("time elapsed")
plot.set_title('size vs time elapsed')
plt.show()
