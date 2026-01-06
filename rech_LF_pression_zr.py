import numpy as np

gamma = 1.4


"""def pression_droite(rho_droite, Zl, gamma=1.4):
    #Pression (gaz parfait)
    rho_l, u_l, P_l = Zl
    return P_l*(2*rho_droite + (1-gamma)*(rho_droite + rho_l)) / (2*rho_l + (gamma-1)*(rho_droite + rho_l))
"""
def pression_droite(rho_r, Zl, gamma=1.4): #forme canonique,  merci chatGPT
    rho_l, u_l, P_l = Zl
    return P_l * (
        (gamma + 1)*rho_r - (gamma - 1)*rho_l
    ) / (
        (gamma + 1)*rho_l - (gamma - 1)*rho_r
    )

def calcul_J(P_r, P_l, rho_r, rho_l):
    J2 = (P_r - P_l) / (1.0 / rho_l - 1.0 / rho_r) #J²
    if J2 < 0:
        raise ValueError("J² négatif : choc non admissible")
    print ("J=", np.sqrt(J2))
    return np.sqrt(J2)

def Zr_fct_de_Zl(Zl, rho_r, gamma):
    rho_l, u_l, P_l = Zl
    P_r= pression_droite(rho_r, Zl, gamma)
    J=calcul_J(P_r, P_l, rho_r, rho_l)
    u_r = u_l - J * (1.0 / rho_l - 1.0 / rho_r) #on calcule Ur
    return rho_r, u_r, P_r


def conservatives(Z, gamma=1.4):
    rho, u, P = Z
    e = P / ((gamma - 1) * rho)
    E = e + 0.5 * u**2
    return np.array([rho, rho*u, rho*E])

# Domaine 
"""
Nx = 500
x = np.linspace(-2.5, 2.5, Nx)
dx = x[1] - x[0]

# Etat gauche
Zl = (1.0, -500.0, 1e5)

# Calcul RH de l'état droit
rho_r = 2.0
Zr = Zr_fct_de_Zl(Zl, rho_r, gamma)

# Variables conservatives
Wl = conservatives(Zl, gamma)
Wr = conservatives(Zr, gamma)

# Initialisation W(x,0)
W = np.zeros((3, Nx))
W[:, x < 0.5]  = Wl[:, None]
W[:, x >= 0.5] = Wr[:, None]
"""

def pressure(W, gamma=1.4):
    rho = W[0]
    u   = W[1] / rho
    E   = W[2] / rho
    return (gamma - 1) * rho * (E - 0.5 * u**2)

def flux(W, gamma=1.4):
    rho = W[0]
    u   = W[1] / rho
    P   = pressure(W, gamma)
    return np.array([
        rho * u,
        rho * u**2 + P,
        u * (W[2] + P)
    ])

def lax_friedrichs(W, dx, dt, gamma=1.4):
    Wnew = W.copy()
    for i in range(1, W.shape[1] - 1):
        Wnew[:, i] = 0.5*(W[:, i+1] + W[:, i-1]) \
                     - dt/(2*dx)*(flux(W[:, i+1], gamma)
                                 - flux(W[:, i-1], gamma))
    return Wnew


def calcul_erreur_L1(x, rho_num, Zl, Zr, t):
    rho_l, u_l, P_l = Zl
    rho_r, u_r, P_r = Zr

    dx = x[1] - x[0]
    
    # u - J
    sigma = (rho_r*u_r - rho_l*u_l) / (rho_r - rho_l)
    #solution exacte
    rho_exact = np.where(x < sigma*t, rho_l, rho_r)
    
    erreur = np.sum(np.abs(rho_num - rho_exact)) * dx
    print("sigma =", sigma)
    return erreur

Nx=500
x = np.linspace(-10, 5, Nx)
dx = x[1] - x[0]
  
# Etat gauche
Zl = (1.0, 850.0, 1e5)
    
# calcul de l'état droit
rho_r = 3
Zr = Zr_fct_de_Zl(Zl, rho_r, gamma)
  
# Variables conservatives
Wl = conservatives(Zl, gamma)
Wr = conservatives(Zr, gamma)
    
    # Initialisation W(x,0)
W = np.zeros((3, Nx))
W[:, x < 0.5]  = Wl[:, None]
W[:, x >= 0.5] = Wr[:, None]
    
CFL = 0.45
Tfinal = 0.01
t = 0.0
    
while t < Tfinal:
    
        rho = W[0]
        u   = W[1] / rho
        P   = pressure(W, gamma)
        c   = np.sqrt(gamma * P / rho)
    
        dt = CFL * dx / np.max(np.abs(u) + c)
        if t + dt > Tfinal:
            dt = Tfinal - t
    
        W = lax_friedrichs(W, dx, dt, gamma)
    
        # Conditions aux limites simples
        W[:, 0]  = W[:, 1]
        W[:, -1] = W[:, -2]
    
        t += dt
        
        
rho_num = W[0]
err = calcul_erreur_L1(x, rho_num, Zl, Zr, t)
print("Erreur L1 =", err)


import matplotlib.pyplot as plt

# États analytiques RH
rho_l, u_l, P_l = Zl
rho_r, u_r, P_r = Zr

"""
#Densité
plt.figure()
plt.plot(x, rho) #R-H
plt.axhline(rho_l, linestyle="--") #L-F
plt.axhline(rho_r, linestyle="--")
plt.xlabel("x")
plt.ylabel("densité")
plt.title("Densité R-H (pointillés), L-F")
plt.show()


#Vitesse
plt.figure()
plt.plot(x, u)
plt.axhline(u_l, linestyle="--")
plt.axhline(u_r, linestyle="--")
plt.xlabel("x")
plt.ylabel("vitesse")
plt.title("Vitesse")
plt.show()
"""

#Pression
plt.figure()
plt.plot(x, P)
plt.axhline(P_l, linestyle="--")
plt.axhline(P_r, linestyle="--")
plt.xlabel("x")
plt.ylabel("pression")
plt.title("Pression")
plt.show()

print(Zl, "\n", Zr)