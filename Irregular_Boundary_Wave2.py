import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm

# Lets 1st define our grid
h = 1/(2**3) # only use (1/(2**k)), k = 1,2,3,4,5.... due to inconsistency in numpy arange (this is the space step)
T = 300
x = np.arange(0, 10+h, h)
y = np.arange(0, 10+h, h)
X, Y = np.meshgrid(x,y)


# Now lets define our domain function
def domain(x, y):
    listdomx = []
    listdomy = []
    for i in range(0, len(X)):
        for j in range(0, len(Y)):
            if ((X[i][j] - 5)**2 + (Y[i][j] - 5)**2) < (3.5**2):
                listdomx.append(X[i][j])
                listdomy.append(Y[i][j])
    
    return(listdomx, listdomy)
    
xx, yy = domain(X, Y)

# We need to now get rid of bad boundary points i.e points with only 1 neighbouring point
def bad_points(x, y, h):
    bad = []
    g = list(zip(x, y))
    for i in range(0, len(g)):
        points = [(g[i][0]+h, g[i][1]), (g[i][0]-h, g[i][1]), (g[i][0], g[i][1]+h),(g[i][0], g[i][1]-h)]
        count = 0
        for j in range(0, len(points)):
            if points[j] in g:
                count += 1
        
        if count == 1:
            bad.append(g[i])
        return(bad)

deleted = bad_points(xx, yy, h)

old_coors = list(zip(xx, yy))
new_coors = [i for i in old_coors if i not in deleted]
xx = []
yy = []
for i in range(0, len(new_coors)):
    xx.append(new_coors[i][0])
    yy.append(new_coors[i][1])

# Now we append the boundary points
def boundary(x, y, h):
    bdry = []
    g = list(zip(x, y))
    for i in range(0, len(g)):
        points = [(g[i][0]+h, g[i][1]), (g[i][0]-h, g[i][1]), (g[i][0], g[i][1]+h),(g[i][0], g[i][1]-h)]
        count = 0
        for j in range(0, len(points)):
            if points[j] in g:
                count += 1
        
        if count != 4:
            bdry.append(g[i])
    return(bdry)

bdrys = boundary(xx,yy,h) # coordinates of the boundary points
xxx = []
yyy = []
for i in range(0, len(bdrys)):
    xxx.append(bdrys[i][0])
    yyy.append(bdrys[i][1])

# Lets color the rest of the grid black (function), this also returns all points not within the boundary
def black(x, y):
    listdomx = []
    listdomy = []
    for i in range(0, len(X)):
        for j in range(0, len(Y)):
            if ((X[i][j] - 5)**2 + (Y[i][j] - 5)**2) > (3.5**2):
                listdomx.append(X[i][j])
                listdomy.append(Y[i][j])
    
    return(listdomx, listdomy)

# To plot the grid using filled contour and color the rest of the grid black
# =============================================================================
# blackx, blacky = black(X, Y)
# z = np.random.randn(len(xx))
# plt.scatter(blackx, blacky, c='black')
# plt.tricontourf(xx,yy,z)
# plt.scatter(blackx, blacky, c='black')
# =============================================================================


# =============================================================================
# #plt.scatter(xx, yy)
# plt.scatter(xxx,yyy, c ='r')
# plt.axis("Equal")
# =============================================================================

# Now we solve the 2D wave equation within the irregualr boundary
blackx, blacky = black(X, Y)
non_eva = list(zip(blackx,blacky)) # Points that are not evaluated 
eva = list(zip(xx,yy)) # coordinates of all points to solve for (evaluated)
# bdrys contains the boundary points
u = np.zeros((len(y),len(x))) # initialize the displacement values for each coordinate point
r = 0.5 #r = cdt/dx
c = 1
dt = (r*h)/c

prev_u = u.copy()
next_u = u.copy()

# Now we iterate through time
t = 0
while t < T:
    
    t += dt
    prev_u = u.copy()
    u = next_u.copy()
    
    # Now initaite the source at some coordinate point(s)
    sauce = (5,7) # Coordinate point of source
    sauce2 = (5,3)
    # Now locate the index of this point
    for i in range(0, len(y)):
        for j in range(0, len(x)):
            if (X[i][j],Y[i][j]) == sauce:
                u[i][j] = (dt**2)*20*np.sin(30*np.pi*t/20) # Sine wave source
            if (X[i][j],Y[i][j]) == sauce2:
                u[i][j] = (dt**2)*20*np.sin(30*np.pi*t/20) # Sine wave source
                
    # 1st lets set all points outside the boundary to be zero
    for i in range(0, len(y)):
        for j in range(0, len(x)):
            if (X[i][j],Y[i][j]) in non_eva:
                u[i][j] = 0
                
    # Next we set the boundary points to be zero (this is dirichlet)
    for i in range(0, len(y)):
        for j in range(0, len(x)):
            if (X[i][j],Y[i][j]) in bdrys:
                next_u[i][j] = 0
    
    # Now we solve for the rest of the points
    for i in range(0, len(y)):
        for j in range(0, len(x)):
            if (X[i][j],Y[i][j]) in eva and (X[i][j],Y[i][j]) not in bdrys:
                next_u[i][j] = 2*u[i][j] - prev_u[i][j] + (r**2)*(u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 4*u[i][j])
    
    # Now we plot
    plt.scatter(blackx, blacky, c='black')
    plt.contourf(X,Y,u, vmin = -0.04, vmax = 0.04)
    plt.scatter(blackx, blacky, c='black')
    plt.axis('Equal')
    plt.xlim([0,9])
    plt.ylim([0,9])
    plt.show()
# =============================================================================
#     fig = plt.figure(figsize = [12,8])
#     ax = fig.gca(projection = '3d')
#     ax.plot_surface(X,Y,u, cmap = cm.coolwarm, vmin = -0.04, vmax = 0.04)
#     ax.set_zlim([-0.04,0.04])
#     plt.show()
# =============================================================================
