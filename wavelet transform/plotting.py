import numpy as np
import matplotlib.pyplot as plt


file = open("cmake-build-debug\out1.txt")

coeffs = np.zeros(4096)
line = file.readline()

nums = line.split()
for i in range(4096):
    coeffs[i] = float(nums[i])

# x = np.linspace(0, 2048, 2048)
# print(len(x))
#
# plt.plot(x, coeffs[0:2048])
#
# plt.show()
#
# plt.plot(x, coeffs[2048:4096])
#
# plt.show()

x = np.linspace(0, 4096, 4096)

plt.plot(x, coeffs)

plt.show()