import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

points = []

for q in [2,3,4,5,6,8,9,10]:
    points.append((0.4,0.9,q))

for c in [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.95,0.99]:
    points.append((0.4,c,7))

for c in [0.37,0.34,0.31,0.28,0.26,0.23,0.21,0.19,0.17,0.14,0.12,0.10,0.08,0.06,0.04,0.02,0.44,0.49,0.54,0.55,0.57,0.59,0.61,0.63,0.65,0.67,0.70]:
    points.append((c,0.9,7))

#for c in [0.1,0.4,0.1,0.4,0.1,0.4,0.1,0.4,0.1,0.4,0.1,0.4,0.1,0.4,0.1,0.4,0.1,0.4]:
#    points.append((c,0.9,7))

x = [p[0] for p in points]
y = [p[1] for p in points]
z = [p[2] for p in points]

ax.set_xlabel(r'$\chi_{NS}$')
ax.set_ylabel(r'$\chi_{BH}$')
ax.set_zlabel('q')

ax.set_xlim(0.0,5.0)
ax.set_ylim(0,1.0)
ax.set_zlim(1,10)

ax.scatter(x,y,z,c="blue",marker="o",s=20)
ax.scatter(0.4,0.9,7,c="red",marker='o',s=100)
plt.show()
plt.savefig("3dparam.png")
