import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.integrate import odeint
from scipy.linalg import eig
from scipy.integrate import solve_ivp
from scipy.special import hermite, factorial

# Shooting method for nonlinear equation
def shoot2_nonlin(psi, x, gamma, epsilon):
    return [psi[1], (gamma*np.abs(psi[0])**2 + x**2 - epsilon) * psi[0]]

# Shooting method for linear equation
def shoot2_lin(psi, x, k, epsilon):
    return [psi[1], (k*(x**2) - epsilon) * psi[0]]

def solve_ivp_shoot(x, psi, k, epsilon):
    return [psi[1], (k * (x ** 2) - epsilon) * psi[0]]

def rhs(x, y, E):
    return [y[1], (x**2 - E) * y[0]]

# Solve the ODE using the specified method and tolerance
def solve_ode(method, tol, E, y0, x_span):
    options = {'rtol': tol, 'atol': tol}
    sol = solve_ivp(rhs, x_span, y0, method=method, args=(E,), **options)
    return sol.t  # Return the time steps taken by the solver

# Compute the average step size from the solver's time steps
def compute_average_step_size(sol_t):
    step_sizes = np.diff(sol_t)
    avg_step_size = np.mean(step_sizes)
    return avg_step_size

# Define the Gaussian envelope
def gaussian_envelope(x):
    return np.exp(-x**2 / 2)

# Normalization constant for the n-th wavefunction
def normalization_constant(n):
    return 1.0 / np.sqrt(2**n * factorial(n) * np.sqrt(np.pi))

# Generate the n-th Hermite polynomial wavefunction with normalization
def quantum_harmonic_oscillator_wavefunction(n, x):
    H_n = hermite(n)  # Generate the n-th Hermite polynomial
    N_n = normalization_constant(n)  # Normalization constant
    psi_n = N_n * H_n(x) * gaussian_envelope(x)  # Normalized wavefunction
    return psi_n

'''
Problem (a)
'''
L = 4
xspan = np.arange(-L, L+0.1, 0.1)
tol = 1e-4
k = 1
col = ['r', 'b', 'g', 'c', 'm', 'k']
shoot_efuncs = []
shoot_evals = []
epsilon_start = 0.1
for modes in range(1, 6):
    epsilon = epsilon_start
    depsilon = 0.1
    for _ in range(1000):
        x0 = [1, np.sqrt(L**2-epsilon)]
        y = odeint(shoot2_lin, x0, xspan, args=(k, epsilon))

        if abs(y[-1, 1] + np.sqrt(L**2-epsilon)*y[-1,0]) < tol:
            shoot_evals.append(epsilon)
            break

        if (-1)**(modes + 1) * (y[-1, 1] + np.sqrt(L**2-epsilon)*y[-1,0]) > 0:
            epsilon += depsilon
        else:
            epsilon -= depsilon / 2
            depsilon /= 2
    
    epsilon_start = epsilon + 0.1
    norm = np.trapz(y[:, 0] * y[:, 0], xspan)  # calculate the normalization
    plt.plot(xspan, y[:, 0] / np.sqrt(norm), col[modes - 1])  # plot modes
    shoot_efuncs.append(abs(y[:, 0] / np.sqrt(norm)))
plt.title('Eigenfunctions for Shooting Method')
plt.grid(True)
plt.show()
    
'''
Problem (b)
'''
N = len(xspan)
dx = xspan[1] - xspan[0]
A = np.zeros((N-2,N-2))
A[0, 0] = 2/3 + xspan[1]**2 * dx**2
A[0, 1] = -2/3
for i in range(1, N-3):
    A[i, i] = 2 + (dx**2) * (xspan[i+1]**2)
    A[i, i-1] = -1
    A[i, i+1] = -1
A[78, 77] = -2/3
A[78, 78] = 2/3 + xspan[N-1]**2 * dx**2

D,V = eig(A) 

# Filter real eigenvalues and their corresponding eigenvectors
real_indices = np.where(np.isclose(D.imag, 0))[0]  # Indices of real eigenvalues
D = D[real_indices].real  # Get the real part of the real eigenvalues
V = V[:, real_indices]  # Corresponding eigenvectors

# Sort real eigenvalues and get first five
sorted_indices = np.argsort(D)  # Sort indices by eigenvalue magnitude
first_five_indices = sorted_indices[:5]  # Get the indices of the first five

# Get first five real eigenvalues and corresponding eigenfunctions
direct_evals = D[first_five_indices]/dx**2
# print(A4)
temp_efunc = V[:, first_five_indices]
direct_efuncs = []

for j in range(0, 5):
    temp = temp_efunc[:, j]
    temp_len = len(temp)
    vec = np.append((4 * temp[0] - temp[1]) / (3 + 2*dx*np.sqrt(16-direct_evals[j])), temp)
    vec = np.append(vec, (4 * temp[temp_len-1] - temp[temp_len-2]) / (3 + 2*dx*np.sqrt(16-direct_evals[j])))
    norm = np.trapz(vec**2, xspan)  # calculate the normalization
    plt.plot(xspan, vec / np.sqrt(norm), col[j - 1])  # plot modes
    direct_efuncs.append(abs(vec / np.sqrt(norm)))
plt.grid(True)
plt.title('Eigenfunctions for Direct Method')
plt.show()

# Collect the absolute values of eigenfunctions from both methods
shoot_efunc_abs = np.array(shoot_efuncs)  # Eigenfunctions from shooting method
direct_efunc_abs = np.array(direct_efuncs)  # Eigenfunctions from direct method

# Plot the absolute values of the eigenfunctions from both methods
plt.figure(figsize=(10, 6))

# Plot eigenfunctions from shooting method
for i in range(shoot_efunc_abs.shape[0]):
    plt.plot(xspan, shoot_efunc_abs[i], col[i], label=f'Shooting Mode {i+1}' if i < 5 else '')

# Plot eigenfunctions from direct method
for j in range(direct_efunc_abs.shape[0]):
    plt.plot(xspan, direct_efunc_abs[j], col[j], linestyle='--', label=f'Direct Mode {j+1}' if j < 5 else '')

# Add title and labels
plt.title('Comparison of Eigenfunctions (Absolute Values) from Both Methods')
plt.xlabel('x')
plt.ylabel('Eigenfunction Value (Absolute)')

# Add grid and legend
plt.grid(True)
plt.legend(loc='upper right')

# Show the plot
plt.show()

'''
Problem (c)
'''
# Gamma = 0.05
L = 2
xspan = np.arange(-L, L+0.1, 0.1)
tol = 1e-4
gamma = 0.05
pnonlin_efuncs = []
pnonlin_evals = []
epsilon_start = 0.001
A = 0.01
for modes in range(1, 3):
    epsilon = epsilon_start
    depsilon = 0.1
    
    for _ in range(1000):
        x0 = [A, A*np.sqrt(L**2-epsilon)]
        y = odeint(shoot2_nonlin, x0, xspan, args=(gamma, epsilon))
        bound = y[-1, 1] + np.sqrt(L**2-epsilon)*y[-1,0]
        norm = np.trapz(y[:, 0]**2, xspan)

        if abs(norm - 1) < tol:
            pnonlin_evals.append(epsilon)
            break
        else:
            A = A/np.sqrt(norm)

        if abs(bound) < tol:
            pnonlin_evals.append(epsilon)
            break

        if (-1)**(modes + 1) * (bound) < 0:
            epsilon -= depsilon / 2
            depsilon /= 2
        else:
            epsilon += depsilon
    
    epsilon_start = epsilon + 0.1 # calculate the normalization
    plt.plot(xspan, y[:, 0] / np.sqrt(norm), col[modes - 1])  # plot modes
    pnonlin_efuncs.append(abs(y[:, 0] / np.sqrt(norm)))
pnonlin_efuncs = np.transpose(pnonlin_efuncs)
plt.grid(True)
plt.title('Eigenfuctions for Nonlinear Equation: Gamma = 0.05')
plt.show()

# Gamma = -0.05
gamma = -0.05
A = 0.01
nnonlin_efuncs = []
nnonlin_evals = []
epsilon_start = 0.01
for modes in range(1, 3):
    epsilon = epsilon_start
    depsilon = 0.1
    
    for _ in range(1000):
        x0 = [A, A*np.sqrt(L**2-epsilon)]
        y = odeint(shoot2_nonlin, x0, xspan, args=(gamma, epsilon))
        bound = y[-1, 1] + np.sqrt(L**2-epsilon)*y[-1,0]
        norm = np.trapz(y[:, 0]**2, xspan)
        
        if abs(norm - 1) < tol:
            nnonlin_evals.append(epsilon)
            break
        else:
            A = A/np.sqrt(norm)

        if abs(bound) < tol:
            nnonlin_evals.append(epsilon)
            break

        if (-1)**(modes + 1) * (bound) > 0:
            epsilon += depsilon
        else:
            epsilon -= depsilon / 2
            depsilon /= 2
    
    epsilon_start = epsilon + 0.1
    norm = np.trapz(y[:, 0]**2, xspan)  # calculate the normalization
    plt.plot(xspan, y[:, 0] / np.sqrt(norm), col[modes - 1])  # plot modes
    nnonlin_efuncs.append(abs(y[:, 0] / np.sqrt(norm)))
nnonlin_efuncs = np.transpose(nnonlin_efuncs)
plt.grid(True)
plt.title('Eigenfuctions for Nonlinear Equation: Gamma = -0.05')
plt.show()

'''
Problem (d)
'''
K = 1  # Constant
L = 2  # Length of domain
y0 = [1, np.sqrt(K * L**2 - 1)]  # Initial conditions
x_span = [-L, L]  # Domain for x

# Tolerances and methods
TOL = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
methods = ['RK45', 'RK23', 'Radau', 'BDF']

avg_step_sizes = {method: [] for method in methods}

# Solve the ODE for each method and tolerance
for method in methods:
    for tol in TOL:
        sol_t = solve_ode(method, tol, K, y0, x_span)
        avg_step_size = compute_average_step_size(sol_t)
        avg_step_sizes[method].append(avg_step_size)

# Plot log-log graph of average step size vs. tolerance
slopes = []
for method in methods:
    log_tol = np.log10(TOL)
    log_avg_step_size = np.log10(avg_step_sizes[method])

    # Plotting log-log data
    plt.plot(log_avg_step_size, log_tol, label=method)

    # Compute the slope using polyfit
    slope = np.polyfit(log_avg_step_size, log_tol, 1)[0]
    slopes.append(slope)

# Display the graph
plt.xlabel('log(Average Step Size)')
plt.ylabel('log(Tolerance)')
plt.legend()
plt.title('Step Size vs Tolerance for ODE Solvers')
plt.grid(True)
plt.show()

# Save the slopes as a 4x1 vector
errors = np.transpose(slopes)
# print("Slopes (RK45, RK23, Radau, BDF):", errors)

'''
Problem (e)
'''
# Gets the exact values using solve_ivp
L = 4
xspan = np.arange(-L, L+0.1, 0.1)
xspan2 = (-L,L)
tol = 1e-4
k = 1
eigenvalues_actual = []
epsilon_start = 0.1
A = 0.01
x0 = [0, 1]
for modes in range(1, 6):
    epsilon = epsilon_start
    depsilon = 0.1
    for _ in range(1000):
        sol = solve_ivp(solve_ivp_shoot, xspan2, x0, method='RK45', args=(k, epsilon))
        
        # Check if the solution meets the tolerance boundary condition
        if abs(sol.y[1, -1] + np.sqrt(L ** 2 - epsilon) * sol.y[0, -1]) < tol:
            eigenvalues_actual.append(epsilon)
            break

        # Update epsilon based on the mode and boundary condition
        if (-1) ** (modes + 1) * (sol.y[1, -1] + np.sqrt(L ** 2 - epsilon) * sol.y[0, -1]) > 0:
            epsilon += depsilon
        else:
            epsilon -= depsilon / 2
            depsilon /= 2
    
    epsilon_start = epsilon + 0.1

# Error calculation for shooting and direct method
direct_efuncs = np.transpose(direct_efuncs)
shoot_efuncs = np.transpose(shoot_efuncs)
shoot_error_func = []
shoot_error_val = []
direct_error_func = []
direct_error_val = []
for j in range(5):
    exact_wave = quantum_harmonic_oscillator_wavefunction(j, xspan)
    exact_wave = np.transpose(exact_wave)   
    integrand1 = (abs(shoot_efuncs[:,j]) - abs(exact_wave))**2
    integrand2 = (abs(direct_efuncs[:,j]) - abs(exact_wave))**2
    error_integral1 = integrate.simpson(integrand1, x=xspan)
    shoot_error_func.append(error_integral1)
    error_integral2 = integrate.simpson(integrand2, x=xspan)
    direct_error_func.append(error_integral2)
    a_error = 100 * (abs(shoot_evals[j] - eigenvalues_actual[j]) / eigenvalues_actual[j])
    b_error = 100 * (abs(direct_evals[j] - eigenvalues_actual[j]) / eigenvalues_actual[j])
    shoot_error_val.append(a_error)
    direct_error_val.append(b_error)