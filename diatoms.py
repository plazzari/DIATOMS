import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.integrate import solve_ivp

# ============================================================
# PARAMETRI
# ============================================================
id = "small_diatoms"  # "large_diatoms" o "small_diatoms"
if id == "large_diatoms":
    S_MAX = 120.0      # micron
    S_MIN = 50.0        # taglia minima: default 20 flusso uscente al bordo sinistro
    SST = 70.0         # sexualization: default 50 size threshold
elif id == "small_diatoms":
    S_MAX = 50.0       # micron
    S_MIN = 20.0        # taglia minima: default 20 flusso uscente al bordo sinistro
    SST = 30.0          # sexualization: default 50 size threshold
else:
    raise ValueError(f"ID '{id}' non riconosciuto. Usa 'large_diatoms' o 'small_diatoms'.")

T_MAX = 360.0      # tempo finale

Ns = 600           # punti in size
Nt = 10000          # punti temporali

s = np.linspace(S_MIN, S_MAX, Ns)
ds = s[1] - s[0]

t = np.linspace(0.0, T_MAX, Nt)
dt = t[1] - t[0]

# ============================================================
# PARAMETRI BIOLOGICI - DA CALIBRARE
# ============================================================


omega = 0.5e-5  # scaling factor per la velocità di crescita/splitting: da calibrare per ottenere tempi realistici (giorni) e dimensioni realistiche (micron)

C = 3.8 * (6 * np.pi) ** (-0.17)  # costante Tang 1995 & Irwin 2006: mu = C * s^(-0.51)
print (f"Costante C per mu(s): {C:.4e}")
K_rho = 0.693 * omega * C           # costante della velocità di splitting: rho = -K_rho * s^(2.5-0.51)

S_crit_delta = 30.0
delta_0 = 0.1
c_delta = 2
b_delta = 0.03
a_delta = - b_delta/(2*S_crit_delta)


beta_0 = 0.01

RUN_NUMERIC = True  # se False, salta il time loop e usa solo la soluzione analitica

T_growth = 5.0    # periodo di alternanza crescita (stessa unità di T_MAX)
T_on     = 5.0   # durata della fase attiva dentro ogni periodo

K_base = 1.0      # baseline carrying capacity
K_lin_max = 0.05   # valore del trend lineare a S_MAX
K_lin_min = 0.01   # valore del trend lineare a S_MIN
K_A1   = 2.5     # ampiezza primo picco
K_s1   = 70.0     # posizione primo picco (micron)
K_sig1 = 3.0      # larghezza primo picco
K_A2   = 25.0     # ampiezza secondo picco
K_s2   = 120.0     # posizione secondo picco (micron)
K_sig2 = 3.0      # larghezza secondo picco

# ============================================================
# FUNZIONI DEL MODELLO
# ============================================================

def mu(s):
    s0 = np.maximum(s, 0.005)  # evita zero per il log
    volume = 4/3 * np.pi * (s0/2)**3
    p_sum = 3.8 * np.power(volume,(-0.17))  # Tang 1995 & Irwin 2006
    return p_sum # mu_max * (c_mu + b_mu * s + a_mu * s**2)

def delta(s):
    s0 = np.maximum(s, 0.005)  # evita zero per il log
    volume = 4/3 * np.pi * (s0/2)**3
    p_srs = 0.063 - 0.008 * np.log10(volume)  # Shimoda 2016 
    return p_srs  # mortalità totale: SRS + mortalità base (non negativa)

def rho(s):
    # splitting is growth * ln(2) / doubling_time ln(2) ~ 0.693
    return - K_rho * s**1.99  # velocità di crescita/splitting: dipende da s^2.5 (volume) ma con esponente ridotto perch´=e tiene conto della splitting rate che dipende da s^(-0.51) (Tang 1995 & Irwin 2006)

def K(s):
    """Carrying capacity: trend lineare decrescente + due picchi gaussiani"""
    S_MIN_K= 20.0  # taglia minima per il carrying capacity (al di sotto di questa taglia, K è costante)
    S_MAX_K = 120.0 # taglia massima per il carrying capacity (al di sopra di questa taglia, K è costante)
    linear = K_lin_min + (K_lin_max - K_lin_min) * (s - S_MIN) / (S_MAX - S_MIN)
    g1 = K_A1 * np.exp(-((s - K_s1) / K_sig1)**2)
    g2 = K_A2 * np.exp(-((s - K_s2) / K_sig2)**2)
    return linear + g1 + g2

def beta(s):
    # sigmoide al posto dello step: transizione morbida attorno a SST
    return beta_0 / (1.0 + np.exp((s - SST) / 2.0))

def t_char(s, s0):
    """Inversa di s_char: dato s(t) e s0, restituisce il tempo t.
    Da ds/dt = -K_rho * s^1.99  =>  t = (s^{-0.99} - s0^{-0.99}) / (0.99 * K_rho)
    Nota: s < s0 (la cellula decresce), quindi t >= 0.
    """
    return (s**(-0.99) - s0**(-0.99)) / (0.99 * K_rho)

# ============================================================
# CONTROLLO CFL
# ============================================================

v = rho(s)
cfl = np.max(np.abs(v)) * dt / ds

if cfl > 1:
    print(f"ATTENZIONE: CFL = {cfl:.2f} > 1. Lo schema esplicito può essere instabile.")
else:
    print(f"CFL = {cfl:.3f}")

# ============================================================
# CONDIZIONE INIZIALE
# ============================================================

def initial_condition(s):
    mean = 100.0
    sigma = 10.0
    return np.exp(-0.5 * ((s - mean) / sigma)**2)

n = np.zeros((Nt, Ns))
n[0, :] = initial_condition(s)

# ============================================================
# TIME LOOP
# ============================================================

base_loss = -delta(s) - beta(s)   # mortalità sempre attiva
mu_s      = mu(s)                 # crescita massima (dipende solo da s)

if RUN_NUMERIC:
    for k in range(Nt - 1):

        # Crescita attiva se siamo nella fase ON del ciclo
        growth_active = (t[k] % T_growth) < T_on
        growth = mu_s + base_loss if growth_active else base_loss

        nk = n[k, :].copy()
        nnew = nk.copy()

        vel = rho(s)

        # Upwind per velocità negativa: usa derivata forward
        for i in range(Ns - 1):
            adv = -vel[i] * (nk[i + 1] - nk[i]) / ds
            reac = growth[i] * nk[i] - nk[i]* nk[i] / K(s[i])  
            nnew[i] = nk[i] + dt * (adv + reac)

        # Bordo sinistro: flusso uscente (outflow) a S_MIN
        # vel < 0 => le cellule escono dal dominio, nessun reingresso
       # nnew[0] = max(nk[0] + dt * (growth[0] * nk[0] - nk[0] * nk[0] / K(s[0])
       #                             - abs(vel[0]) * (nk[0]) / ds), 0.0)

        # Bordo destro: n(t, S_MAX) = integral beta(s) n(t,s) ds
        birth_flux = np.trapz(growth*beta(s) * nk, s)
        nnew[-1] = birth_flux
       # nnew[-1] = birth_flux / abs(rho(S_MAX))

        # Evita densità negative numeriche
        nnew[nnew < 0] = 0.0

        n[k + 1, :] = nnew
else:
    print("Time loop saltato (RUN_NUMERIC=False).")
    n[:,:]  = 0.01  # riempi con NaN per evidenziare che non sono stati calcolati
# ============================================================
# PLOT
# ============================================================

print("Generazione dei plot...", flush=True)

# ── Panel unico: mu/delta/beta  |  soluzione sulla caratteristica ──
s0_char = 100.0
# Soluzione analitica: integrazione stiff (Radau) da S_MAX verso S_MIN
# ODE lungo la caratteristica:
#   dn/ds = [(\mu(s) - \delta(s) - \beta(s)) \, n - n^2/K] / \rho(s)
def f_ode(s_val, y):
    n_val = y[0]
    return [((mu(s_val) - delta(s_val) - beta(s_val)) * n_val - n_val**2 / K(s_val)) / rho(s_val)]

sol = solve_ivp(
    f_ode,
    t_span=(s[-1], s[0]),          # da S_MAX verso S_MIN
    y0=[n[-1, -1]],                # condizione iniziale a S_MAX
    method='Radau',                # solver implicito per ODE stiff
    t_eval=s[::-1],                # valuta su tutti i punti griglia (invertiti)
    rtol=1e-6, atol=1e-10
)
n_analytical = np.maximum(sol.y[0, ::-1], 0.0)  # riporta nell'ordine originale

# Caratteristica che parte da s0 al tempo t=0: ds/dt = rho(s) = -K_rho*s^1.99
# Soluzione: s_char(t) = (s0^{-0.99} + 0.99*K_rho*t)^{-1/0.99}
s_char = (s0_char**(-0.99) + 0.99 * K_rho * t)**(-1.0/0.99)

# Estrai n lungo la caratteristica (solo punti dentro il dominio)
mask = (s_char >= S_MIN) & (s_char <= S_MAX)
s_char_valid = s_char[mask]
n_char_valid = np.array([
    n[k, np.argmin(np.abs(s - s_char[k]))]
    for k in np.where(mask)[0]
])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Panel 1: mu, delta, beta (sinistra) + K(s) (destra)
ax1.plot(s, mu(s),    label=r"$\mu(s)$")
ax1.plot(s, delta(s), label=r"$\delta(s)$")
ax1.plot(s, beta(s),  label=r"$\beta(s)$")
ax1.set_xlabel("Size s [$\mu$m]")
ax1.set_ylabel("Rate")
#ax1.set_title(r"$\mu(s)$, $\delta(s)$, $\beta(s)$, $K(s)$")
ax1.text(-0.02, 1.02, "A", transform=ax1.transAxes, fontsize=14, fontweight='bold')

ax1b = ax1.twinx()
ax1b.plot(s, K(s), color='gray', linestyle='--', label=r"$K(s)$")
ax1b.set_ylabel(r"$K(s)$", color='gray')
ax1b.tick_params(axis='y', labelcolor='gray')

# Unisci le leggende dei due assi
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax1b.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2)

# Panel 2: soluzione numerica, analitica, punti sulla caratteristica
n_eq = np.maximum((mu(s) - delta(s) - beta(s)) * K(s), 0.0)   # equilibrio: n* = (mu-delta-beta)*K
if RUN_NUMERIC:
    ax2.plot(s, n[-1, :], label=r"Numeric $t=t_{fin}$")
ax2.plot(s, n_analytical, '--', label=r"Characteristics")
ax2.plot(s, n_eq, ':', label=r"Stable Equilibrium $n^*=(\mu-\delta-\beta)\,K(s)$")
#ax2.scatter(s_char_valid, n_char_valid, s=10, zorder=5,
#            label=rf"Caratteristica $s_0={s0_char:.0f}\,\mu$m")
ax2.set_xlabel("Size s [$\mu$m]")
ax2.set_ylabel(r"$n(t,s)$")
ax2.set_yscale('log')
#ax2.set_title("Comparison: Numeric vs Characteristics vs Equilibrium")
ax2.text(-0.02, 1.02, "B", transform=ax2.transAxes, fontsize=14, fontweight='bold')
ax2.axvline(SST, color='k', linestyle=':', linewidth=1.2, label=r"$s_{SST}$")
ax2.legend()

# Asse x secondario (in cima): età della diatomea = t_char(s, S_MAX)
def _s_to_age(s_val):
    s_val = np.asarray(s_val, dtype=float)
    return (s_val**(-0.99) - S_MAX**(-0.99)) / (0.99 * K_rho)

def _age_to_s(t_val):
    t_val = np.asarray(t_val, dtype=float)
    return (S_MAX**(-0.99) + 0.99 * K_rho * t_val)**(-1.0 / 0.99)

for _ax in (ax1, ax2):
    ax_top = _ax.secondary_xaxis('top', functions=(_s_to_age, _age_to_s))
    ax_top.set_xlabel("Diatom cohort age [days]")

plt.tight_layout()
plt.savefig(f"rates_and_characteristic_{id}.png", dpi=150)
plt.show()

if RUN_NUMERIC:
    plt.figure(figsize=(8, 5))
    n_plot = np.clip(n.T, a_min=n[n > 0].min() if (n > 0).any() else 1e-15, a_max=None)
    plt.imshow(
        n_plot,
        aspect='auto',
        origin='lower',
        extent=[t[0], t[-1], s[0], s[-1]],
        #norm=LogNorm()
    )
    plt.colorbar(label=r"$n(t,s)$")
    plt.xlabel("Time")
    plt.ylabel("Size s [$\mu$m]")
    plt.title("Diatom size-structured model")
    plt.tight_layout()
    plt.savefig(f"heatmap_{id}.png", dpi=150)
    plt.show()

    plt.figure(figsize=(8, 5))
    for idx in [0, Ns//4, Ns//2, -1]:
        plt.plot(t, n[:, idx], label=f"s = {s[idx]:.1f} $\mu$m")
    plt.xlabel("Time")
    plt.ylabel(r"$n(t,s)$")
    plt.yscale('log')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"time_series_{id}.png", dpi=150)
    plt.show()