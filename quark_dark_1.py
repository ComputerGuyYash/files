from pylab import *
from cmath import *
from scipy.optimize import fsolve
e0q = 1.01
e0d = 393961.4974
B_eff = 57.55602 * 8.96260385*(10**(-7))/ e0q
Ro = 1.476
dr = 0.0001
R = 0.000000001
Ar = -1.70419484e-05
Anr = 1.69251932e+00
Asr = 2.57248611e+00
def e_bar_q(p):
    e_bar = (3*p) + (4*B_eff)
    return e_bar

def e_bar_d(x):
    e_den = Ar*(x**(3/4)) + Anr*(x**(3/5)) + Asr*x
    return e_den

    
def dp_dr_q(p_q, p_d, m_q, m_d, r, r_ec):
    dp_dr_gr_dm = (-1*Ro*(m_d+m_q)*e_bar_q(p_q)/(r*r))*(1 + (p_q/e_bar_q(p_q)))*(1+(4*pi*r*r*r*((p_d*e0d) + (e0q*p_q))/(m_q+m_d)))*((1-(2*Ro*(m_d+m_q)/r))**(-1))
   
    return dp_dr_gr_dm

def dp_dr_d(p_q, p_d, m_q, m_d, r, r_ec):
    dp_dr_gr_dm = r_ec*(-1*Ro*(m_q+m_d)*e_bar_d(p_d)/(r*r))*(1 + (p_d/e_bar_d(p_d)))*(1+(4*pi*r*r*r*((p_d*e0d) + (e0q*p_q))/(m_q+m_d)))*((1-(2*Ro*(m_q+m_d)/r))**(-1))
    
    return dp_dr_gr_dm


    
    
def dm_dr_q(p_q, p_d, m_q, m_d, r, r_ec):
    dm_dr = 4*pi*r*r*(e_bar_q(p_q))*e0q
   
    return dm_dr
def dm_dr_d(p_q, p_d, m_q, m_d, r, r_ec):
    dm_dr = 4*pi*r*r*(e_bar_d(p_d))*e0d*r_ec
    
    return dm_dr

def rk4_gr(p_q, p_d, m_q, m_d, r, r_ec):
    
    q_p_1 = dp_dr_q(p_q, p_d, m_q, m_d, r, r_ec)
    q_m_1 = dm_dr_q(p_q, p_d, m_q, m_d, r, r_ec)
    d_p_1 = dp_dr_d(p_q, p_d, m_q, m_d, r, r_ec)
    d_m_1 = dm_dr_d(p_q, p_d, m_q, m_d, r, r_ec)
    
    q_p_2 = dp_dr_q(p_q + (dr*q_p_1/2), p_d + (dr*d_p_1/2), m_q + (dr*q_m_1/2), m_d + (dr*d_m_1/2), r + (dr/2), r_ec)
    q_m_2 = dm_dr_q(p_q + (dr*q_p_1/2), p_d + (dr*d_p_1/2), m_q + (dr*q_m_1/2), m_d + (dr*d_m_1/2), r + (dr/2), r_ec)
    d_p_2 = dp_dr_d(p_q + (dr*q_p_1/2), p_d + (dr*d_p_1/2), m_q + (dr*q_m_1/2), m_d + (dr*d_m_1/2), r + (dr/2), r_ec)
    d_m_2 = dm_dr_d(p_q + (dr*q_p_1/2), p_d + (dr*d_p_1/2), m_q + (dr*q_m_1/2), m_d + (dr*d_m_1/2), r + (dr/2), r_ec)
   
    q_p_3 = dp_dr_q(p_q + (dr*q_p_2/2), p_d + (dr*d_p_2/2), m_q + (dr*q_m_2/2), m_d + (dr*d_m_2/2), r + (dr/2), r_ec)
    q_m_3 = dm_dr_q(p_q + (dr*q_p_2/2), p_d + (dr*d_p_2/2), m_q + (dr*q_m_2/2), m_d + (dr*d_m_2/2), r + (dr/2), r_ec)
    d_p_3 = dp_dr_d(p_q + (dr*q_p_2/2), p_d + (dr*d_p_2/2), m_q + (dr*q_m_2/2), m_d + (dr*d_m_2/2), r + (dr/2), r_ec)
    d_m_3 = dm_dr_d(p_q + (dr*q_p_2/2), p_d + (dr*d_p_2/2), m_q + (dr*q_m_2/2), m_d + (dr*d_m_2/2), r + (dr/2), r_ec)
    
    q_p_4 = dp_dr_q(p_q + (dr*q_p_3), p_d + (dr*d_p_3), m_q + (dr*q_m_3), m_d + (dr*d_m_3), r + (dr), r_ec)
    q_m_4 = dm_dr_q(p_q + (dr*q_p_3), p_d + (dr*d_p_3), m_q + (dr*q_m_3), m_d + (dr*d_m_3), r + (dr), r_ec)
    d_p_4 = dp_dr_d(p_q + (dr*q_p_3), p_d + (dr*d_p_3), m_q + (dr*q_m_3), m_d + (dr*d_m_3), r + (dr), r_ec)
    d_m_4 = dm_dr_d(p_q + (dr*q_p_3), p_d + (dr*d_p_3), m_q + (dr*q_m_3), m_d + (dr*d_m_3), r + (dr), r_ec)
    
    
    
    q_p = (1/6)*(q_p_1 + 2*q_p_2 + 2*q_p_3 + q_p_4)
    d_p = (1/6)*(d_p_1 + 2*d_p_2 + 2*d_p_3 + d_p_4)
    q_m = (1/6)*(q_m_1 + 2*q_m_2 + 2*q_m_3 + q_m_4)
    d_m = (1/6)*(d_m_1 + 2*d_m_2 + 2*d_m_3 + d_m_4)
    
    
    return (p_q + q_p*dr, p_d + d_p*dr, m_q + q_m*dr, m_d + d_m*dr, r + dr, r_ec)


P_C =  ( 8.96260385*(10**(-4)))/e0q

def star_solver(p_qt, p_dt):
    p_q = p_qt*(P_C)
    p_d = p_dt
    m_q = (4*pi*R*R*R/3)*(e_bar_q(p_q))*(e0q)
    m_d = (4*pi*R*R*R/3)*(e_bar_d(p_d))*(e0d)
    r_ = R
    r_ec = 1
    state = [(p_q, p_d, m_q, m_d, r_, r_ec)]
    count = 1
    while (p_q.imag == 0 and isnan(p_q) == 0):
        state.append(rk4_gr(p_q, p_d, m_q, m_d, r_, r_ec))
        p_q = state[len(state)-1][0]
        p_d = state[len(state)-1][1]
        m_q = state[len(state)-1][2]
        m_d = state[len(state)-1][3]
        r_ = state[len(state)-1][4]
    p_d_Img = state[len(state)-4][1]
    p_q = p_qt*(P_C)
    p_d = p_dt
    m_q = (4*pi*R*R*R/3)*(e_bar_q(p_q))*(e0q)
    m_d = (4*pi*R*R*R/3)*(e_bar_d(p_d))*(e0d)
    r_ = R
    r_ec = 1
    state_final = [(p_q, p_d, m_q, m_d, r_, r_ec)]
    count = 1
    rcount = 1
    ko = 1/dr
    ko = int(ko)
    while (p_q > 0):

        if rcount % ko == 0:
            pass
        rcount = rcount + 1
        if p_d < p_d_Img:
            count = count +1
            r_ec = 10**(-9*(count))
        state_final.append(rk4_gr(p_q, p_d, m_q, m_d, r_, r_ec))
        p_q = state_final[len(state_final)-1][0]
        p_d = state_final[len(state_final)-1][1]
        m_q = state_final[len(state_final)-1][2]
        m_d = state_final[len(state_final)-1][3]
        r_ = state_final[len(state_final)-1][4]
    return (m_q + m_d, r_, m_d/m_q)

m_t_1_final = []
m_d_1_final = []
dm_1_frac = []
radius_1_final = []
fin_1_count = 1
m_d_final_change = []
r_d_final_change = []
p_d_final_change = []
dmfrac_d_final_change = []
p_q_final_change = []
p_inp_ar = linspace(0.001, 10, 300)
out_count = 1
for p_dm_ in linspace(4.72776964e-12, 4.72776964e-8, 300):
    for p_inp in p_inp_ar:
        print(fin_1_count, out_count)
        fin_1_count = fin_1_count + 1
        x = star_solver(p_inp, p_dm_)
        m_t_1_final.append(x[0])
        dm_1_frac.append(x[2])
        radius_1_final.append(x[1])
    m_d_final_change.append(max(m_t_1_final))
    r_d_final_change.append(radius_1_final[m_t_1_final.index(max(m_t_1_final))])
    p_d_final_change.append(p_dm_)
    dmfrac_d_final_change.append(dm_1_frac[m_t_1_final.index(max(m_t_1_final))])
    p_q_final_change.append(p_inp_ar[m_t_1_final.index(max(m_t_1_final))])
    m_t_1_final = []
    m_d_1_final = []
    dm_1_frac = []
    radius_1_final = []
    fin_1_count = 1
   
    out_count = out_count + 1
    
import pickle
with open('dm_q_md_1.pkl', 'wb') as f:
    pickle.dump([m_d_final_change, r_d_final_change, dmfrac_d_final_change, p_q_final_change], f)



    