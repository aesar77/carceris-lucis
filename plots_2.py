import numpy as np
import matplotlib.pyplot as plt


data1 = np.loadtxt("output2A.csv", delimiter=",")
data2 = np.loadtxt("output2B.csv", delimiter=",")
data3 = np.loadtxt("output2C.csv", delimiter=",")
data4 =  np.loadtxt("output2D.csv", delimiter=",")
data5 = np.loadtxt("output2E.csv", delimiter=",")

orbits = ["Orbit A", "Orbit B", "Orbit C", "Orbit D", "Orbit E"]

r0 = [1+np.sqrt(2), 1+np.sqrt(3), 3, 1+2*np.sqrt(2), 2]

dat = [data1, data2, data3, data4, data5]

errors = []
times = []

fig = plt.figure(figsize=plt.figaspect(0.2))

for d in range(len(dat)):
    data = dat[d]
    t = data[:, 0]
    r = data[:, 1]
    th = data[:, 2]
    phi = data[:, 3]
    xx = []
    yy = []
    zz = []
    a = 1
    
    errors.append(abs(r - r0[d])/r0[d])
    times.append(t)

    for i in range(len(r)):
        x = np.sqrt(r[i]**2 + a**2) * np.sin(th[i])*np.cos(phi[i])
        y = np.sqrt(r[i]**2 + a**2) * np.sin(th[i])*np.sin(phi[i])
        z = r[i] * np.cos(th[i])
        xx.append(x)
        yy.append(y)
        zz.append(z)



    u = np.linspace(0, np.pi, 30)
    v = np.linspace(0, 2 * np.pi, 30)

    xxx = np.outer(np.sin(u), np.sin(v))*0.5
    yyy = np.outer(np.sin(u), np.cos(v))*0.5
    zzz = np.outer(np.cos(u), np.ones_like(v))*0.5

    ax = fig.add_subplot(1, 5, d+1, projection='3d')
    ax.plot_surface(xxx, yyy, zzz, color="black")
    ax.plot(xx, yy, zz, color="orange")
    plt.legend(["BH", orbits[d]], loc="upper center")

    

plt.show()

for e in range(len(errors)):

    
    plt.xlabel("t/M")
    plt.ylabel(r"|r-$r_0$|/$r_0$")
    plt.plot(times[e], errors[e])
    plt.legend(orbits)

plt.ylim(0, 0.1)
plt.show()
