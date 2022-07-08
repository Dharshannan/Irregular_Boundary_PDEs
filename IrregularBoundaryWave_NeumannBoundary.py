import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm

# Lets 1st define our grid
h = 1/(2**4) # only use (1/(2**k)), k = 1,2,3,4,5.... due to inconsistency in numpy arange (this is the space step)
T = 300
x = np.arange(-2, 2+h, h)
y = np.arange(-2, 2+h, h)
X, Y = np.meshgrid(x,y)


# Now lets define our domain function
def domain(x, y):
    listdomx = []
    listdomy = []
    for i in range(0, len(X)):
        for j in range(0, len(Y)):
            if (2.8*x[i][j]**2*(x[i][j]**2*(2.5*x[i][j]**2+y[i][j]**2-2)+1.2*y[i][j]**2*(y[i][j]*(3*y[i][j]-0.75)-6.0311)+3.09)+0.98*y[i][j]**2*((y[i][j]**2-3.01)*y[i][j]**2+3)-1.005) < 0:
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
            if (2.8*x[i][j]**2*(x[i][j]**2*(2.5*x[i][j]**2+y[i][j]**2-2)+1.2*y[i][j]**2*(y[i][j]*(3*y[i][j]-0.75)-6.0311)+3.09)+0.98*y[i][j]**2*((y[i][j]**2-3.01)*y[i][j]**2+3)-1.005) > 0:
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
    sauce = (0,0) # Coordinate point of source
    sauce2 = (0,1)
    # Now locate the index of this point
    for i in range(0, len(y)):
        for j in range(0, len(x)):
            if (X[i][j],Y[i][j]) == sauce:
                u[i][j] = (dt**2)*20*np.sin(100*np.pi*t/20) # Sine wave source
            if (X[i][j],Y[i][j]) == sauce2:
                u[i][j] = (dt**2)*20*np.sin(100*np.pi*t/20)
                
    # 1st lets set all points outside the boundary to be zero
    for i in range(0, len(y)):
        for j in range(0, len(x)):
            if (X[i][j],Y[i][j]) in non_eva:
                u[i][j] = 0
    
    # Now we solve for the rest of the points
    for i in range(0, len(y)):
        for j in range(0, len(x)):
            
            # Solve for nuemann (du/dn = 0) boundary points
            if (X[i][j],Y[i][j]) in bdrys:
                pts = [(i+1,j),(i-1,j),(i,j+1),(i,j-1)]
                if (X[pts[0]],Y[pts[0]]) not in eva:
                    pts[0] = pts[1]
                if (X[pts[1]],Y[pts[1]]) not in eva:
                    pts[1] = pts[0]
                if (X[pts[2]],Y[pts[2]]) not in eva:
                    pts[2] = pts[3]
                if (X[pts[3]],Y[pts[3]]) not in eva:
                    pts[3] = pts[2]
                    
                next_u[i,j] = 2*u[i,j] - prev_u[i,j] + (r**2)*(u[pts[0]] + u[pts[1]] + u[pts[2]] + u[pts[3]] - 4*u[i,j])
                
            # Solve for rest of the points
            if (X[i][j],Y[i][j]) in eva and (X[i][j],Y[i][j]) not in bdrys:
                next_u[i][j] = 2*u[i][j] - prev_u[i][j] + (r**2)*(u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 4*u[i][j])
    
    # Now we plot
    plt.scatter(blackx, blacky, c='black')
    plt.contourf(X,Y,u, vmin = -0.0175, vmax = 0.0175)
    plt.scatter(blackx, blacky, c='black')
    plt.axis('Equal')
    plt.xlim([-1.5,1.5])
    plt.ylim([-1.5,1.5])
    plt.show()
