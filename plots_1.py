import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt("output1-a.csv", delimiter=",")
data2 = np.loadtxt("output1-b.csv", delimiter=",")

d = [data1, data2]
title = [r"$r_{obs}$ = 10 M, $\xi$ = 0.48857874", r"$r_{obs}$ = 20 M, $\xi$ = 0.24904964"]
title 

for j in range(2):
    data = d[j]
    t = data[:, 0]
    r = data[:, 1]
    phi = data[:, 3]
    xx = []
    yy = []

    for i in range(len(r)):
        x = r[i]*np.cos(phi[i])
        y = r[i]*np.sin(phi[i])
        xx.append(x)
        yy.append(y)

    plt.plot(xx,yy, color = "blue")
    c=plt.Circle((0, 0), 2, color='black')
    plt.gca().add_artist(c)
    plt.title(title[j], loc="left")
    intr = np.interp(xx[-1], xx, yy)
    plt.title(r"Interp. at $r_{obs}$ = " + "{:.4f}".format(intr), loc="right" )
    plt.xlabel("x/M")
    plt.ylabel("y/M")
    plt.show()

    print(intr)
