import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("TkAgg")
import numpy as np

total_times_100 = []
pcat_times_100 = []
total_times_500 = []
pcat_times_500 = []
test_nums = [1,2,3,4] # number of threads

plt.scatter(test_nums,pcat_times_100,c='black')
plt.scatter(test_nums,total_times_100,c='black')
plt.scatter(test_nums,pcat_times_500,c='black')
plt.scatter(test_nums,total_times_500,c='black')

plt.plot(np.unique(test_nums), np.poly1d(np.polyfit(test_nums, total_times_100, 1))(np.unique(test_nums)),label='100 Samps - Total') # 100 samps , total
plt.plot(np.unique(test_nums), np.poly1d(np.polyfit(test_nums, pcat_times_100, 1))(np.unique(test_nums)),label='100 Samps - PCAT') # 100 samps, pcat
plt.plot(np.unique(test_nums), np.poly1d(np.polyfit(test_nums, total_times_500, 1))(np.unique(test_nums)),label='500 Samps - Total') # 500 samps, total
plt.plot(np.unique(test_nums), np.poly1d(np.polyfit(test_nums, pcat_times_500, 1))(np.unique(test_nums)),label='500 Samps - PCAT') # 500 samps, pcat

plt.title('Multi-Threading Tests on PCAT-SPIRE')
plt.xlabel('Num of Threads')
plt.ylabel('Run Time [sec]')
plt.legend()
plt.show()
