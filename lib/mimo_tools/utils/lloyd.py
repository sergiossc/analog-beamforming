import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import uuid

n = 10
k = 5

# class who defines a point
class Point:
    def __init__(self, pos_x, pos_y):
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.id = uuid.uuid4()

    def get_coordinates(self):
        return np.array([self.pos_x, self.pos_y])

# setting point
points = []
for i in range(n):
    x = np.random.rand()
    y = np.random.rand()
    p = Point(x,y)
    points.append(p)


# getting some clusters(random) from points
clusters = []
for i in range(k):
    cluster = points[np.random.randint(n)]
    clusters.append(cluster)


# getting the closest cluster from each point
for p in points:
    print('>>>point: ', p.id)
    d_min = min([ norm(p.get_coordinates() - c.get_coordinates()) for c in clusters ])
    print('---d_min: ', d_min)
        


#plt.plot(x, y, '*') 
##plt.plot(clusters_x, clusters_y, '*') 
#plt.show()
