import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.animation as animation


ani = 'y'
static = ''


save_fig = ''
num_saves = 100

#animation parameters
loop = 'y'
if save_fig:
    loop = 'y'


dx = 1e-10
dt = 1e-26
T  = 6e-15

iterations = int( np.ceil(T/dt) )

L = 110e-9
boundaries = [47e-9, 53e-9, 57e-9,63e-9]

DL = 1e-9
beta_peaks = [0, 0.5, 0.5]
beta_prime_peaks = [0, 0.7, 0.7]
D_peaks = [0.005, 0.003, 0.003]
lambda_sf_peaks = [20e-9, 5e-9, 5e-9]
lambda_J_peaks = [1e-20, 4e-9, 4e-9]
lambda_phi_peaks = [1e-20, 4e-9, 4e-9]
Mx_peaks = [0, 0, 0]
My_peaks = [0, 0, 0]
Mz_peaks = [0, 1, 0]
j_e_peak = 1e15
minf = 4e7
teq = 1e-50


iterations = int(iterations)

animation_speed = 1000

def ax_filler(ax, intensity=0.15):
    color_intensities = 1 - (beta/np.max(beta)* intensity)
    for i in range(len(Uz)):
        if beta[i] > 0:
            ax.add_patch(Rectangle((i*dx*1e9,-1e20),dx*1e9,2e20, facecolor=str(color_intensities[i]), edgecolor="none"))
    return ax

def p(y):
    plt.close()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x,y)
    ax_filler(ax)
    ax.set_xlabel('X [nm]')
    fig.show()

def diffuse_calc(x, non_mag, mag1, mag2):
    mag1 = mag1 - non_mag
    mag2 = mag2 - non_mag

    mid_points = [boundaries[i+1]-boundaries[i] for i in range(len(boundaries)-1)]
    y1 = mag1/2*( np.tanh( ( x - boundaries[0] ) / (DL/2) ) + 1 )
    y2 = -mag1/2*( np.tanh( ( x - boundaries[0] - mid_points[0] ) / (DL/2) ) + 1 ) + mag1
    y3 = mag2/2*( np.tanh( ( x - boundaries[1] - mid_points[1] ) / (DL/2) ) + 1 )
    y4 = -mag2/2*( np.tanh( ( x - boundaries[2] - mid_points[2] ) / (DL/2) ) + 1 ) + mag2
    return (y1+y2+y3+y4)-(y1[0]+y2[0]+y3[0]+y4[0])+non_mag

def spin_current(m, M):
    jm = np.zeros(len(m))
    diff_m = np.gradient(m, dx)
    for i in range(len(m)):
        first_term = beta[i] * M[i] * lin_increase_je(i*dt, teq, j_e_peak)
        second_term = diff_m[i]
        third_term = beta[i] * beta_prime[i] * M[i] * M[i] * diff_m[i]
        jm[i] = first_term - 2 * D[i] * (second_term - third_term)
    return jm


def lin_increase_je(t, teq, je):
    if t < teq:
        return t*je/teq
    else:
        return je


x = np.arange(0,L,dx)
beta = diffuse_calc(x, *beta_peaks)
beta_prime = diffuse_calc(x, *beta_prime_peaks)
D = diffuse_calc(x, *D_peaks)
lambda_sf = diffuse_calc(x, *lambda_sf_peaks)
lambda_J = diffuse_calc(x, *lambda_J_peaks)
lambda_phi = diffuse_calc(x, *lambda_phi_peaks)
Mx = diffuse_calc(x, Mx_peaks[0], Mx_peaks[1], Mx_peaks[2])
My = diffuse_calc(x, My_peaks[0], My_peaks[1], My_peaks[2])
Mz = diffuse_calc(x, Mz_peaks[0], Mz_peaks[1], Mz_peaks[2])

fact1 = dt/(lambda_J * lambda_J)
fact2 = dt/(lambda_phi * lambda_phi)
fact3 = dt/(lambda_sf*lambda_sf)


jmx = np.zeros(len(x))
jmy = np.zeros(len(x))
jmz = np.zeros(len(x))

Ux = np.zeros(len(x))
Uy = np.zeros(len(x))
Uz = np.zeros(len(x))


def solver(Ux, Uy, Uz, jmx, jmy, jmz):
    Ux_new = np.zeros(len(x))
    Uy_new = np.zeros(len(x))
    Uz_new = np.zeros(len(x))

    grad_jmx = np.gradient(jmx,dx)
    grad_jmy = np.gradient(jmy,dx)
    grad_jmz = np.gradient(jmz,dx)

    # Ux terms
    for i in range(len(x)):
        first_term = Ux[i]
        second_term = dt*grad_jmx[i]
        third_term = fact1[i] * ( -My[i]*Uz[i] + Mz[i]*Uy[i] )
        fourth_term = fact2[i] * (-Mx[i]*My[i]*Uy[i] - Mx[i]*Mz[i]*Uz[i] + Ux[i]*(My[i]*My[i] + Mz[i]*Mz[i]) )
        fifth_term = fact3[i] * (Ux[i] - minf*Mx[i])
        Ux_new[i] = first_term - (second_term + third_term + fourth_term + fifth_term)

    # Uy terms
    for i in range(len(x)):
        first_term = Uy[i]
        second_term = dt*grad_jmy[i]
        third_term = fact1[i] * ( Mx[i]*Uz[i] - Mz[i]*Ux[i] )
        fourth_term = fact2[i] * (Mx[i]*Mx[i]*Uy[i] - Mx[i]*My[i]*Ux[i] - My[i]*Mz[i]*Uz[i] + Mz[i]*Mz[i]*Uy[i])
        fifth_term = fact3[i] * (Uy[i] - minf*My[i])
        Uy_new[i] = first_term - (second_term + third_term + fourth_term + fifth_term)

    # Uz terms
    for i in range(len(x)):
        first_term = Uz[i]
        second_term = dt*grad_jmz[i]
        third_term = fact1[i] * ( -Mx[i]*Uy[i] + My[i]*Ux[i] )
        fourth_term = fact2[i] * (Mx[i]*Mx[i]*Uz[i] - Mx[i]*Mz[i]*Ux[i] + My[i]*My[i]*Uz[i] - My[i]*Mz[i]*Uy[i] )
        fifth_term = fact3[i] * (Uz[i] - minf*Mz[i])
        Uz_new[i] = first_term - (second_term + third_term + fourth_term + fifth_term)

    jmx = spin_current(Ux_new, Mx)
    jmy = spin_current(Uy_new, My)
    jmz = spin_current(Uz_new, Mz)

    return Ux_new, Uy_new, Uz_new, jmx, jmy, jmz


def plot_setup():
    fig = plt.figure(facecolor='white', figsize=(16, 8))
    ax = plt.subplot2grid((6,2), (0,0), rowspan=5)
    ax_jm = plt.subplot2grid((6,2), (0,1), rowspan=5)
    ax_params = plt.subplot2grid((6,2), (5,0), colspan=2)
    params = [beta_peaks, beta_prime_peaks, D_peaks, lambda_sf_peaks, lambda_phi_peaks, lambda_J_peaks, Mx_peaks, My_peaks, Mz_peaks]
    ax_params.axis('off')
    param_data = [['Non-Magnet']+["%.2g"%i[0] for i in params], ['F1']+["%.2g"%i[1] for i in params], ['F2']+["%.2g"%i[2] for i in params] ]
    headers = ['',r'$\beta$', r"$\beta$'",r"$D$", r"$\lambda_{sf}$", r"$\lambda_{\phi}$", r"$\lambda_{J}$", r"$M_x$", r"$M_y$", r"$M_z$"]
    table = ax_params.table(cellText = param_data, loc='center', colLabels=headers, cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(15)
    table.scale(1, 2.5)
    
    
    ax_filler(ax)
    ax.set_title("Spin Accumuation -Z", fontsize=14)
    ax.set_xlabel(r"X [nm]")
    ax_jm.set_xlabel(r"X [nm]")
    #ax.set_ylabel(r"Uz [$C m ^{-3}$]")
    #ax_jm.set_ylabel(r"$j_{m_{z}}$ [$A m^{-2}$]")
    ax.set_xlim([0,L*1e9])
    
    ax_filler(ax_jm)
    ax_jm.set_title("Spin Current -Z", fontsize=14)
    
    ax.set_xlim([25, 74])
    ax_jm.set_xlim([25, 75])
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig, ax, ax_jm, ax_params

if ani == 'y': 
    fig, ax, ax_jm, ax_params = plot_setup()
    counter = 0
    counter2 = 1
    
    line1, = ax.plot([0], [0], 'b-')
    line2, = ax_jm.plot(x*1e9, jmz)
    
    ylims1 = [-1, 1]
    ylims2 = [-1, 1]

    def ani_plot(i):
        global counter, counter2
        global line1, line2, ylims1, ylims2
        global jmz, jmx, jmy
        global Uz, Ux, Uy
        
        if save_fig == 'y':
            save_steps = int(iterations/num_saves)
            if save_steps == 0:
                save_steps = 1
            if (counter2 == num_saves+1):
                global save_fig 
                save_fig = ''
            elif (i%(save_steps) == 0):
                fig.savefig("./img/%i.png"%counter2, format='png')
                print("Saved figure %i/%i"%(counter2, num_saves))
                counter2 += 1
        
        fig.suptitle("Diffusion Length = %.2gnm, Current time = %.2g s"%(DL*1e9, (i+1)*dt), fontsize=16)
        counter += 1
        if i == 0 and loop == 'y':
            jmx = np.zeros(len(x))
            jmy = np.zeros(len(x))
            jmz = np.zeros(len(x))
            
            Ux = np.zeros(len(x))
            Uy = np.zeros(len(x))
            Uz = np.zeros(len(x))
         
        Ux, Uy, Uz, jmx, jmy, jmz = solver(Ux, Uy, Uz, jmx, jmy, jmz)
        
        max_Uz = np.max(Uz)
        min_Uz = np.min(Uz)
        if max_Uz > ylims1[1] or min_Uz < ylims1[0]:
            ylims1[0] = 2*min_Uz
            ylims1[1] = 2*max_Uz
            ax.set_ylim(ylims1)
            
        max_jmz = np.max(jmz)
        min_jmz = np.min(jmz)
        if max_jmz > ylims2[1] or min_jmz < ylims2[0]:
            ylims2[0] = 2*min_jmz
            ylims2[1] = 2*max_jmz
            ax_jm.set_ylim(ylims2)
            
        line1.set_data(x*1e9, Uz )
        line2.set_data(x*1e9, jmz )
        plt.draw()
        return line1, line2
    
    line_ani = animation.FuncAnimation(fig, ani_plot, iterations, interval=animation_speed, blit=False)
    
    plt.show()
    plt.show()
if static == 'y':
    
    fit = np.polyfit([1000,100000],[10,200] , 1)
    num_displayed_iterations = int(np.polyval(fit, iterations))
    for i in range(int(iterations)):
        if (i%(int(iterations/num_displayed_iterations)+1) == 0):
            print("Iteration = "+str(i)+"/"+str(iterations))
        Ux, Uy, Uz, jmx, jmy, jmz = solver(Ux, Uy, Uz, jmx, jmy, jmz)
        
    fig, ax_static, ax_static1, ax_params = plot_setup()
    
    ax_static.plot(x*1e9,Ux, 'r-')
    ax_static.plot(x*1e9,Uy, 'g-')
    ax_static.plot(x*1e9,np.gradient(Uz, dx), 'b-')
    
    ax_static1.plot(x*1e9,jmz, 'b-')
    
    fig.suptitle("Diffusion Length = %.2gnm, Current time = %.2g s"%(DL*1e9, T), fontsize=16)
    ax_static.set_ylim([np.min(np.gradient(Uz, dx))*1.1, np.max(np.gradient(Uz, dx))*1.1])
    ax_static1.set_ylim([np.min(jmz)*1.1, np.max(jmz)*1.1])
    plt.show()
