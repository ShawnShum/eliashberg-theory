import numpy as np
from numpy.lib.function_base import delete
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cbook
from matplotlib import cm
from matplotlib.colors import LightSource
from scipy import integrate
from scipy.optimize import fsolve
from scipy.interpolate import make_interp_spline
import seaborn as sns
import mpl_toolkits.mplot3d
temp = 0.005
freq_cut = 8
om_E = 1
g = 4
angle = np.pi/16*np.arange(32)

FS_x_quarter = []
FS_y_quarter = []
FS_x = []
FS_y = []
for i in range(5):

    def f(x):
        return -2*(np.cos(x)+np.cos(np.tan(i*np.pi/16)*x))+0.24

    sol1_fsolve = fsolve(f, [1])
    y_sol1_fsolve = np.tan(i*np.pi/16)*(sol1_fsolve)
    FS_x.append(sol1_fsolve[0])
    FS_y.append(y_sol1_fsolve[0])

FS_x_quarter = FS_x + sorted(FS_y, reverse=True)
FS_y_quarter = FS_y + sorted(FS_x, reverse=False)

del FS_x_quarter[4]
del FS_y_quarter[4]

FS_x_half = np.vstack((np.array(FS_x_quarter), (-1 * np.array(FS_y_quarter))))
FS_y_half = np.vstack((np.array(FS_y_quarter), (np.array(FS_x_quarter))))
FS_x_half = np.reshape(FS_x_half, (18, 1))
FS_y_half = np.reshape(FS_y_half, (18, 1))

FS_x_half = np.delete(FS_x_half, 8)
FS_y_half = np.delete(FS_y_half, 8)

FS_x = np.vstack((np.array(FS_x_half), (-1 * np.array(FS_x_half))))
FS_y = np.vstack((np.array(FS_y_half), (-1 * np.array(FS_y_half))))
FS_x = np.reshape(FS_x, (34, 1))
FS_y = np.reshape(FS_y, (34, 1))

FS_x = np.delete(FS_x, 17)
FS_y = np.delete(FS_y, 17)
FS_x = np.delete(FS_x, 0)
FS_y = np.delete(FS_y, 0)

# fermion velocity
FS_v = np.sqrt(4*np.sin(FS_x)**2+4*np.sin(FS_y)**2)
''' x_new = np.linspace(0, 32, 200)
y_smooth = make_interp_spline(range(32), 1/FS_v)(x_new)
plt.plot(x_new, y_smooth, color='b', marker='o', linestyle='-')
 '''

# fermion surface plot
for N in range(100):
    E = N/100

    def f(x, y):
        return 1

    def h(x):
        return np.arccos(E-np.cos(x))
    v, err = integrate.dblquad(f, 0, np.arccos(E-1), lambda x: 0, h)
    Ratio = 4*v/(4*np.pi**2)
    if (np.abs(Ratio-7/16) < 0.0015):
        print(Ratio)
        print(-2*E)
        step = 0.01
'''         a = np.arange(-np.pi, np.pi, step)
        b = np.arange(-np.pi, np.pi, step)
        X, Y = np.meshgrid(a, b)
        Z = -2*(np.cos(X)+np.cos(Y))
        plt.axis("equal")
        plt.title('The Brillouin zone patching scheme')
        plt.xlabel("$k_x$")
        plt.ylabel("$k_y$")
        plt.ylim(-3.5, 3.5)
        plt.xlim(-3.5, 3.5)
        contour1 = plt.contour(
            X, Y, Z, [-2*E, 0.0], label='n', colors='b')
        plt.scatter(FS_x, FS_y, marker='x', color='r')
        plt.show()
        plt.savefig('fermion_surface_patch.jpg', dpi=500, bbox_inches='tight') '''
           


# Eliashberg kernel

def matsubara(idx, temp):
    return (2*idx+1)*np.pi*temp


mu = -0.24
b = -1*mu/2


# f, ax= plt.subplots(figsize = (14, 10))
# sns.heatmap(kernel_g, linewidths = 0.05, ax = ax)
# ax.set_title('convolution kernel of g')
# ax.invert_yaxis()
# f.savefig('kernel_g.jpg', dpi=100, bbox_inches='tight')
# self_energy

# Eliashberg


def matsubara(idx, temp):
    return (2*idx+1)*np.pi*temp


def num_freq(freq_cut, temp):
    return int(2*np.floor(freq_cut/(2*np.pi*temp)))


def phonon_pro(om_E, mf_1, mf_2):
    return 1.0*om_E**2/(om_E**2 + np.power((mf_1-mf_2), 2))


def mat_freq_list(num_freq, temp):
    return matsubara(np.arange(num_freq)-num_freq/2, temp)


phonon = np.zeros((num_freq(freq_cut, temp), num_freq(freq_cut, temp)))
mf_list1 = mat_freq_list(num_freq(freq_cut, temp), temp)
for i in range(num_freq(freq_cut, temp)):
    for j in range(num_freq(freq_cut, temp)):
        phonon[i][j] = 1.0*om_E**2 / \
            (om_E**2 + np.power((mf_list1[i]-mf_list1[j]), 2))

''' f, ax= plt.subplots(figsize = (14, 10))
sns.heatmap(phonon, linewidths = 0.05, ax = ax)
ax.set_title('phonon propagator')
ax.invert_yaxis()
f.savefig('phonon propagator.jpg', dpi=100, bbox_inches='tight') 
 '''

kernel_p_out = np.ones(32)


for j in range(32):
    kernel_p_out[j] = 0
    for i in range(32):
        if (i == 31):
            dl = np.sqrt((FS_x[i]-FS_x[0])**2+(FS_y[i]-FS_y
                                               [0])**2)
        else:
            dl = np.sqrt((FS_x[i]-FS_x[i+1])**2+(FS_y[i]-FS_y
                                                 [i+1])**2)
        kernel_p_out[j] = 1/np.pi**2*dl*g**2*(np.sin(FS_x[j])+np.sin(
            FS_y[j])-np.sin(FS_x[i])-np.sin(FS_y[i]))**2/FS_v[i] + kernel_p_out[j]
x_new = np.linspace(0, 32, 200)
y_smooth = make_interp_spline(range(32), kernel_p_out)(x_new)
plt.plot(x_new, y_smooth, color='b', marker='o', linestyle='-')


kernel = np.zeros((32, num_freq(freq_cut, temp), num_freq(freq_cut, temp)))
for i in range(32):
    for j in range(num_freq(freq_cut, temp)):
        for k in range(num_freq(freq_cut, temp)):
            kernel[i][j][k] = kernel_p_out[i]*phonon[j][k]

''' f, ax= plt.subplots(figsize = (14, 10))
sns.heatmap(kernel[:,:,int(num_freq(freq_cut, temp)/2)], linewidths = 0.05, ax = ax)
ax.set_title('convolution kernel of g')
ax.invert_yaxis()
f.savefig('kernel.jpg', dpi=100, bbox_inches='tight')  
 '''
""" a, b = np.meshgrid(range(32), range(num_freq(freq_cut, temp)))
c = kernel_p_out[a]*phonon[b]
ax = plt.subplot(111, projection='3d')
ax.plot_surface(a,b,c, cstride=1,  cmap='rainbow',
                linewidth=0, antialiased=False, shade=False)
f.savefig('kernel_3d.jpg', dpi=100, bbox_inches='tight')
 """

# Interation of eliashberg function

max_iter = 300  
tol = 0.00001  
Gap_function_out = np.ones((32, num_freq(freq_cut, temp)))
Gap_function_in = np.ones((32, num_freq(freq_cut, temp)))
''' for i in range(32):
    Gap_function_in[i, :] = (np.cos(FS_x[i])-np.cos(FS_y[i]))*np.ones(num_freq(freq_cut, temp)) '''
''' for i in range(32):
    Gap_function_in[i, :] = (np.sin(FS_x[i])*np.sin(FS_y[i]))*np.ones(num_freq(freq_cut, temp)) '''
# Gap_function_in = 0.1*np.ones((32, num_freq(freq_cut, temp)))  
convoluted_function = np.ones((32, num_freq(freq_cut, temp)))

for k in range(32):  

    for i in range(max_iter):  
        for j in range(num_freq(freq_cut, temp)):  
            convoluted_function[k, :] = (Gap_function_in[k, :]-mf_list1/mf_list1[j]
                                         * Gap_function_in[k][j])/np.sqrt(mf_list1**2+Gap_function_in[k, :]**2)
            Gap_function_out[k][j] = np.pi*temp * \
                np.dot(kernel[k, j, :], convoluted_function[k, :])
        diff = np.linalg.norm(
            Gap_function_in[k, :] - Gap_function_out[k, :]) 
        Gap_function_in[k, :] = Gap_function_out[k, :]
        if (diff < tol):
            print('\nat i = ', i, 'p = ', k,
            ' the convergence is achieved')
            break
x_new = np.linspace(0, 32, 200)
y_smooth = make_interp_spline(range(32), Gap_function_out[:, int(
    num_freq(freq_cut, temp)/2)])(x_new)
plt.plot(x_new, y_smooth, color='b', marker='o', linestyle='-')
plt.title('Gap function at T = %.3f' %temp )
plt.xlabel("$p$")
plt.ylabel("$\Delta$")
plt.show()


''' f, ax = plt.subplots(figsize=(14, 10))
sns.heatmap(Gap_function_out, linewidths=0.05, ax=ax)
ax.set_title('Gap_function')
ax.invert_yaxis() 
 '''



''' 
def myplot(matrix1):
    xsize, ysize = matrix1.shape
    x = np.arange(0, ysize, 1)
    y = np.arange(0, xsize, 1)
    xs, ys = np.meshgrid(x, y)
    z1 = matrix1
    fig=plt.figure(figsize = (14, 10))
    ax1 = fig.add_subplot(2, 1, 1, projection='3d')
    ax1.plot_surface(xs, ys, z1,color="blue",alpha=0.5,cmap='gnuplot',rstride=1, cstride=10)
    ax1.contourf(xs, ys, z1, zdir='z', offset=-0.5, cmap='gnuplot')
    ax1.contourf(xs, ys, z1, zdir='x', offset=-10, cmap='gnuplot')
    ax1.set_xlim(-10, 328)
    ax1.set_ylim(0, 40)
    ax1.set_zlim(-0.5, 1.2835)
    plt.title('Gap function at T = %.3f' %temp )
    plt.ylabel("$p$")
    plt.xlabel("$\omega_n$")
    plt.tight_layout
    fig.savefig('gap_function_3D.jpg', dpi=500, bbox_inches='tight')

myplot(Gap_function_out)
 '''

 