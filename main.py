# Alright so this is basically just a dumb binary search
# over the energy levels

# You use sqrt(2m(V_g - V(x))) and integrate that where V_g is your guess and V(x) is the potential equation (left hand side of wkb)
# You set the right hand side h\pi(n+1/2) to whatever n value you want
# Compare the left side to the right side, depending on whether the lef side is bigger or smaller you adjust your upper
# and lower bounds for your energy guess until you converge on some energy and then that should be your energy

# Be prepared to see some interesting techniques but I didn't really know how else to go about it

import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt

# Convergence parameter
eps = 0.0000001

# Load in the file
data = np.loadtxt('Li1Pi_u-potential.dat')

# Select the data for the x and V
x = data[:,0]
V = data[:,1]

# Mass of Na in a.u. (the value in a.m.u converted to Hartree Fock value...)
m = 6.941 / (5.4857990*10**-4)

# Index of minimum energy value
min_i = int(np.where(V == min(V))[0])
# X value corresponding to minimum energy value
min_x = x[min_i]
# The actual minimum energy value
min_V = V[min_i]

# List of the energies for each energy level (the answers)
E = []

# used for debugging pretty much
index_bounds = []

# Loops over how many energy levels you want to do
for n in range(5):

    # maximum bound of binary search values, it's the value of the maximum level of the crest (the maximum on the right side of the minimum potential value)
    V_test_max = max(V[min_i:])
    # minimum bound, it's just the minimum potential
    V_test_min = min_V
    # the value we're testing (the value between max and min)
    V_test = (V_test_max+V_test_min)/2

    # The result of the integral
    res = None
    # the target value/the right hand side of wkb
    target = np.pi*(n+0.5)

    # Counter to break the loop if it doesn't converge
    iteration = 0

    # Loop until it converges (or until a certain max number of iterations)
    while abs(res - target) > eps:

        V_test = (V_test_max + V_test_min) / 2

        # Ok so this is where the weirdness really starts.
        # I needed a way to get the corresponding x-values for each V_test
        # I did this by starting at the index of the minimum energy value,
        # then moving to the left until a potential value is found that is closest to V_test.
        # Then I take that best possible index (best_left_i) and use the corresponding x-value as the lower bound of the integration
        i = min_i
        best_left_i = i
        delta_min = abs(V[0] - min_V)
        while i >= 0:
            delta = abs(V_test - V[i])
            if delta < delta_min:
                delta_min = delta
                best_left_i = i
            i -= 1

        # everything I said above except now it's to the right and it's the upper bound of the integration
        i = min_i
        best_right_i = i
        delta_min = abs(V[0] - min_V)
        while i <= np.where(V == max(V[min_i:]))[0]:
            delta = abs(V_test - V[i])
            if delta < delta_min:
                delta_min = delta
                best_right_i = i
            i += 1
        
        # integrand
        f = np.sqrt(2 * m * (V_test - V[best_left_i+1:best_right_i]))
        # integral
        res = spi.trapz(f, x[best_left_i+1:best_right_i])

        # The binary search part
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