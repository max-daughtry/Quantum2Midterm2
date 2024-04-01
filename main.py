import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

eps = 0.0000001

data = np.loadtxt('Li1Pi_u-potential.dat')

x = data[:,0]
V = data[:,1]

m = 6.941 / (5.4857990*10**-4)

min_i = int(np.where(V == min(V))[0])
min_x = x[min_i]
min_V = V[min_i]

E = []
index_bounds = []
for n in range(5):
    V_test_max = max(V[min_i:])
    V_test_min = min_V
    V_test = V_test_max

    res = np.pi
    target = np.pi*(n+0.5)

    iteration = 0
    while abs(res - target) > eps:
        V_test = (V_test_max + V_test_min) / 2

        i = min_i
        best_left_i = i
        delta_min = abs(V[0] - min_V)
        while i >= 0:
            delta = abs(V_test - V[i])
            if delta < delta_min:
                delta_min = delta
                best_left_i = i
            i -= 1

        i = min_i
        best_right_i = i
        delta_min = abs(V[0] - min_V)
        while i <= np.where(V == max(V[min_i:]))[0]:
            delta = abs(V_test - V[i])
            if delta < delta_min:
                delta_min = delta
                best_right_i = i
            i += 1
        

        f = np.sqrt(2 * m * (V_test - V[best_left_i+1:best_right_i]))

        res = spi.trapz(f, x[best_left_i+1:best_right_i])

        if (res > target):
            V_test_max = V_test
        elif (res < target):
            V_test_min = V_test
        else:
            break

        if iteration > 100:
            break

        iteration += 1
        
    print()
    print("n: ", n)
    print("Energy: ", V_test)
    print("error: ", abs(res-target))
    E.append(V_test)
    index_bounds.append((best_left_i, best_right_i))

print()
print()
for i in range(len(E)-1):
    print(abs(E[i+1]-E[i]))

l=20
u=200

plt.figure(figsize=(16/1.5,9/1.5))
plt.title("Energy levels", fontsize=18)
plt.ylabel("Energy (a.u.)", fontsize=15)
plt.xlabel("Internuclear distance (a.u.)", fontsize=15)
plt.scatter(x[l:u],V[l:u], label='V(x)',color='blue')
for (n, (e, i)) in enumerate(zip(E,index_bounds)):
    plt.plot(x[i[0]:i[1]],[e] * (i[1]-i[0]), color='red')
plt.legend()
plt.savefig('levels.png')
# plt.show()