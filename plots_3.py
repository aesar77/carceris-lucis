import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
matplotlib.rcParams.update({'font.size': 22})

dat1= np.loadtxt("output_raytracing1-a.csv", delimiter=",")
dat2= np.loadtxt("output_raytracing1-b.csv", delimiter=",")
dat3= np.loadtxt("output_raytracing2-a.csv", delimiter=",")
dat4= np.loadtxt("output_raytracing2-b.csv", delimiter=",")
dat5= np.loadtxt("output_raytracing3-a.csv", delimiter=",")
dat6= np.loadtxt("output_raytracing3-b.csv", delimiter=",")

fig1 = plt.figure(figsize=plt.figaspect(0.2))
ax = fig1.add_subplot(1, 2, 1)
ax.imshow(dat1, cmap="hot", extent=(-25, 25, -25, 25))
plt.xlabel(r"$x_{sc}/M$")
plt.ylabel(r"$y_{sc}/M$")
ax = fig1.add_subplot(1, 2, 2)
ax.imshow(dat2, cmap="hot", extent=(-25, 25, -25, 25))
plt.xlabel(r"$x_{sc}/M$")
plt.ylabel(r"$y_{sc}/M$")
plt.show()

fig2 = plt.figure(figsize=plt.figaspect(0.2))
ax2 = fig2.add_subplot(1, 2, 1)
ax2.imshow(dat3, cmap="hot", extent=(-25, 25, -25, 25))
plt.xlabel(r"$x_{sc}/M$")
plt.ylabel(r"$y_{sc}/M$")
ax2 = fig2.add_subplot(1, 2, 2)
ax2.imshow(dat4, cmap="hot", extent=(-25, 25, -25, 25))
plt.xlabel(r"$x_{sc}/M$")
plt.ylabel(r"$y_{sc}/M$")
plt.show()


fig3 = plt.figure(figsize=plt.figaspect(0.2))
ax2 = fig3.add_subplot(1, 2, 1)
ax2.imshow(dat5, cmap="hot", extent=(-25, 25, -25, 25))
plt.xlabel(r"$x_{sc}/M$")
plt.ylabel(r"$y_{sc}/M$")
ax2 = fig3.add_subplot(1, 2, 2)
ax2.imshow(dat6, cmap="hot", extent=(-25, 25, -25, 25))
plt.xlabel(r"$x_{sc}/M$")
plt.ylabel(r"$y_{sc}/M$")
plt.show()
# for i in range(len(temp)):
#     if (temp[i] != 0):
#         plt.scatter(x[i], y[i])

# plt.show()

# plt.plot(xx,yy,z)
# plt.show()
# t = data[:, 0]
# r = data[:, 1]
# th = data[:, 2]
# phi = data[:, 3]
# xx = []
# yy = []
# zz = []

# for i in range(len(r)):
#     x = np.sqrt(r[i]**2 + a**2) * np.sin(th[i])*np.cos(phi[i])
#     y = np.sqrt(r[i]**2 + a**2) * np.sin(th[i])*np.sin(phi[i])
#     z = r[i] * np.cos(th[i])
#     xx.append(x)
#     yy.append(y)
#     zz.append(z)



# plt.plot(xx,yy, color = "blue")
# plt.title(r"$r_{obs}$ = 10M, $\xi$ = 0.48857874")
# plt.xlabel("x/M")
# plt.ylabel("y/M")
# plt.show()

# ax = plt.figure().add_subplot(projection='3d')
# ax.plot(xx, yy, z, label='parametric curve')
# ax.legend()

# plt.show()

