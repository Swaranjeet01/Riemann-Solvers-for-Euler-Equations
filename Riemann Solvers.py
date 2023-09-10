from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
# import sys

# orig_stdout = sys.stdout
# fp = open('Assignment2_Output.txt', 'w')
# sys.stdout = fp




# Heat capacity ratio
gamma=1.4

# Input data for pressure, veloctiy and density
W_L=[[460.894,19.5975,5.99924]]
W_R=[[46.0950,-6.19633,5.99242]]

x=int(input("Enter position to find final pressure, velocity and density:"))
t=int(input("Enter time to find final pressure velocity and density:"))

# Solving for all five cases
for i in range(len(W_L)):
    
    print("\n"+"*"*70,"\n")
    print("Question",i+1,":\n")
    
    # Assigning pressure, velocity and density
    [P_L,u_L,rho_L]=W_L[i]
    [P_R,u_R,rho_R]=W_R[i]
    
    # Sonic speed
    a_L=(gamma*P_L/ rho_L)**0.5
    a_R=(gamma*P_R/ rho_R)**0.5
    
    # Initial guess for P*
    P_s=0.5*(P_L+P_R)
    if P_L==P_R: P_s=1e-6
    
    iter=0
    
    # N-R iteration to find P*
    while(iter<1):
        if(P_s>P_L):
            f_L=(P_s-P_L)*sqrt((2.0/((gamma+1.0)*rho_L))/(P_s+((gamma-1.0)/(gamma+1.0))*P_L))
            d_fL=sqrt((2.0/((gamma+1.0)*rho_L))/(P_s+((gamma-1.0)/(gamma+1.0))*P_L)) * \
                (1.0-(P_s-P_L)/(2.0*(P_s+P_L*(gamma-1.0)/(gamma+1.0))))
    
        else:
            f_L=((2.0*a_L)/(gamma-1.0))*((P_s/P_L)**((gamma-1.0)/(2.0*gamma)) - 1.0)
            d_fL=(1/(rho_L*a_L)) * (P_s/P_L)**(-(gamma+1.0)/(2.0*gamma))
        if(P_s>P_R):
            f_R=(P_s-P_R)*sqrt((2.0/((gamma+1.0)*rho_R))/(P_s+((gamma-1.0)/(gamma+1.0))*P_R))
            d_fR=sqrt((2.0/((gamma+1.0)*rho_R))/(P_s+((gamma-1.0)/(gamma+1.0))*P_R)) * \
                (1.0-(P_s-P_R)/(2.0*(P_s+P_R*(gamma-1.0)/(gamma+1.0))))
        else:
            f_R=((2.0*a_R)/(gamma-1.0))*((P_s/P_R)**((gamma-1.0)/(2.0*gamma)) - 1.0)
            d_fR=(1/(rho_R*a_R)) * (P_s/P_R)**(-(gamma+1.0)/(2.0*gamma))
        dP=-(f_L+f_R+u_R-u_L)/(d_fL+d_fR)
        iter+=1
        if(abs(dP)<1e-6):
            break
        
        # Updating P*
        P_s=P_s+dP
    
    # Calculating U*
    u_s=0.5*(u_L+u_R)+0.5*(f_R-f_L)
    
    print("Pressure in star region, P* =",round(P_s,6))
    print("Contact discontinuity velocity, U* =",round(u_s,6),"\n")
    
    # Density rho*L and rho*R and shock and rarefaction speed
    if (P_s>P_R):
        rho_s_R=rho_R*((P_s/P_R + (gamma-1)/(gamma+1))/((gamma-1)/(gamma+1)*P_s/P_R + 1))
        s_R=u_R+a_R*((gamma+1)/(2*gamma)*P_s/P_R + (gamma-1)/(2*gamma))**0.5
        print("Right side shock: \nShock speed =",round(s_R,6))
        print("rho*_R =",round(rho_s_R,6),"\n")
    else:
        rho_s_R=rho_R*(P_s/P_R)**(1/gamma)
        a_s_R=(gamma*P_s/rho_s_R)**0.5    
        r_RH=u_R+a_R
        r_RT=u_s+a_s_R
        print("Right side rarefaction: \nHead velocity =",round(r_RH,6),"\nTail velocity =",round(r_RT,6))
        print("rho*_R =",round(rho_s_R,6),"\n")
    
    if (P_s>P_L):
        rho_s_L=rho_L*((P_s/P_L + (gamma-1)/(gamma+1))/((gamma-1)/(gamma+1)*P_s/P_L + 1))
        s_L=u_L-a_L*((gamma+1)/(2*gamma)*P_s/P_L + (gamma-1)/(2*gamma))**0.5
        print("Left side shock: \nShock speed =",round(s_L,6))
        print("rho*_L =",round(rho_s_L,6))
    else:
        rho_s_L=rho_L*(P_s/P_L)**(1/gamma)
        a_s_L=(gamma*P_s/rho_s_L)**0.5
        r_LH=u_L-a_L
        r_LT=u_s-a_s_L
        print("Left side rarefaction: \nHead velocity =",round(r_LH,6),"\nTail velocity =",round(r_LT,6))
        print("rho*_L =",round(rho_s_L,6))
    
    
    
    
    
    # Finding W depending on region (x,t)
    s=x/t
    W=[0]*3
    if (s<u_s):
        if(P_s>P_L):
            if(s<s_L):
                W=[P_L,u_L,rho_L]
            else:
                W=[P_s,u_s,rho_s_L]
        else:
            if(s<r_LH):
                W=[P_L,u_L,rho_L]
            else:
                if(s>r_LT):
                    W=[P_s,u_s,rho_s_L]
                else:
                    rho_L_fan=rho_L*(2/(gamma+1) + (gamma-1)/((gamma+1)*a_L)*(u_L-s))**(2/(gamma-1))
                    u_L_fan=2/(gamma+1)*(a_L + (gamma-1)/2*u_L + s)
                    P_L_fan=P_L*(2/(gamma+1) + (gamma-1)/((gamma+1)*a_L)*(u_L-s))**(2*gamma/(gamma-1))
                    W=[P_L_fan,u_L_fan,rho_L_fan]
    else:
        if(P_s>P_R):
            if(s>s_R):
                W=[P_R,u_R,rho_R]
            else:
                W=[P_s,u_s,rho_s_R]
        else:
            if(s>r_RH):
                W=[P_R,u_R,rho_R]
            else:
                if(s<r_RT):
                    W=[P_s,u_s,rho_s_R]
                else:
                    rho_R_fan=rho_R*(2/(gamma+1) - (gamma-1)/((gamma+1)*a_R)*(u_R-s))**(2/(gamma-1))
                    u_R_fan=2/(gamma+1)*(-a_R + (gamma-1)/2*u_R + s)
                    P_R_fan=P_R*(2/(gamma+1) - (gamma-1)/((gamma+1)*a_R)*(u_R-s))**(2*gamma/(gamma-1))
                    W=[P_R_fan,u_R_fan,rho_R_fan]
    
    print("\nAt (x,t)=(%.2f,%.2f):"%(x,t),"\nPressure, velocity, density=(%.6f,%.6f,%.6f)"%(W[0],W[1],W[2]))
    
    
    
    
    
    # Plotting contact discontinuity, shock and/or rarefaction waves
    fig,ax=plt.subplots(dpi=300)
    t_plot=np.linspace(0,5,100)
    x_star_plot=[u_s*i for i in t_plot]
    ax.plot(x_star_plot,t_plot,'k-.',label='Contact discontinuity')
    if(P_s>P_R):
        x_rs_plot=[s_R*i for i in t_plot]
        ax.plot(x_rs_plot,t_plot,'b',label='Right side shock')
    else:
        x_r_rh_plot=[r_RH*i for i in t_plot]
        x_r_rt_plot=[r_RT*i for i in t_plot]
        ax.plot(x_r_rh_plot,t_plot,'--',color='blue',label='Right side rarefaction head')
        ax.plot(x_r_rt_plot,t_plot,'--',color='dodgerblue',label='Right side rarefaction tail')
    
    if(P_s>P_L):
        x_ls_plot=[s_L*i for i in t_plot]
        ax.plot(x_ls_plot,t_plot,'r',label='Left side shock')
    else:
        x_r_lh_plot=[r_LH*i for i in t_plot]
        x_r_lt_plot=[r_LT*i for i in t_plot]
        ax.plot(x_r_lh_plot,t_plot,'--',color='red',label='Left side rarefaction head')
        ax.plot(x_r_lt_plot,t_plot,'--',color='orange',label='Left side rarefaction tail')
    ax.set_ylabel("Time")
    ax.set_xlabel("Distance")
    ax.set_title("Question "+str(i+1))
    ax.grid(linewidth=0.5)
    ax.legend(loc='lower center', bbox_to_anchor=(0.5,-0.35),
          ncol=2,fancybox=True)
# sys.stdout = orig_stdout
# fp.close()