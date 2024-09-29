import numpy as np

# # filter coeff for a bank of filter
# b0 = np.array([1, 2, 3])
# b1 = 0.4*np.ones(3)
# b2 = np.array([3, 5, 6])

# # Concatenate b0, b1, b2 into a matrix B where each row is [b0, b1, b2]
# B = np.column_stack([b0, b1, b2])

# # print('----------------')
# # print(B)
# # print('----------------')

# sos = np.concatenate([B, B], axis=1)

# print(sos)
# # print(sos.shape[0])

# for i in range(sos.shape[0]):  # normalization of each raw by its a0 coefficient
#     sos[i,:] /= sos[i,3]

# print(sos)


A = 10**(-25/40)
Q = 10
wc = 2*np.pi*50

mu = (2*A**4-2*(A*Q)**2-1)/(A*Q)**2

print(mu)

delta_w_theo = wc/np.sqrt(2)*(np.sqrt(-mu+np.sqrt((mu-2)*(mu+2)))-np.sqrt((-mu-np.sqrt((mu-2)*(mu+2)))))

print(delta_w_theo)

print(10**(-3/40))
print(0.5**0.25)