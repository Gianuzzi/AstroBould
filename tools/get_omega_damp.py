# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt

# EQUATION
# Omega = Omega0 + dOmegadt * t
# 0 = Omega0 + dOmegadt * tau
# tau = - Omega0 / dOmegadt

# Set actual parameters
P_actual = 7.004 / 24.  # [days]
Omega_actual = 2 * np.pi / P_actual  # [Rad/day]

# Example parameters
Prot0_example = P_actual
Omega0_example = 2 * np.pi / Prot0_example  # [Rad/day]
dOmegadt_example = - Omega0_example * 0.3 / 1e7  # [Rad/day^2]
Omega_end_example = Omega0_example + dOmegadt_example * 1e7  # [Rad/day]
Prot_end_example = 2 * np.pi / Omega_end_example  # [days]

# Calculate the damping time
tau_example = - Omega0_example / dOmegadt_example

# Set time
time_end = 365.25 * 5e4  # [yrs]

# My parameters
as_example = True

if as_example:
    Omega0 = Omega0_example
    Omega_end = Omega_end_example
    dOmegadt = dOmegadt_example

else:

    # Set my parameters
    Omega0 = 0  # [Rad/day]
    Omega_end = Omega_actual  # [Rad/day]

    # Calculate my initial Omega or end Omega
    if Omega0 <= 0:
        Omega0 = Omega_end - dOmegadt_example * time_end
    else:
        Omega_end = Omega0 + dOmegadt_example * time_end

    # Calculate my dOmegadt
    dOmegadt = (Omega_end - Omega0) / time_end


# Calculate my damping time
tau = - Omega0 / dOmegadt

# Calculate change rate
change_rate = Omega0/Omega_actual

# Print results
print(f"dOmegadt = {dOmegadt:.4e} rad/day^2")
print(f"Change rate = {change_rate:.7f} == {change_rate*100:.7f}%")
print(f"Time_end = {time_end:.5e} days")
print("Actual parameters")
print(f" Omega = {Omega_actual:.7f}")
print(f" Prot = {P_actual:.7f} days")
print("My parameters")
print(f" Omega0 = {Omega0:.7f}")
print(f" Prot0 = {2*np.pi/Omega0:.7f} days")
print(f" Omega_end = {Omega_end:.7f}")
print(f" Prot_end = {2*np.pi/Omega_end:.7f} days")
print(f" tau_omega = {tau:.4e} days")
print("Example parameters")
print(f" Omega0 = {Omega0_example:.7f}")
print(f" Prot0 = {2*np.pi/Omega0_example:.7f} days")
print(f" Omega_end = {Omega_end_example:.7f}")
print(f" Prot_end = {Prot_end_example:.7f} days")
print(f" tau_omega_example = {tau_example:.4e} days")

# Plot
time = np.sort(np.concatenate((
    np.logspace(-2, np.log10(time_end), 250),
    np.linspace(0, time_end, 251))
))

Omega = Omega0 + dOmegadt * time
Omega0_example = Omega0_example + dOmegadt_example * time
dOmegadt_exponente = np.floor(np.log10(np.abs(dOmegadt)))
dOmegadt_base = dOmegadt / (10**dOmegadt_exponente)
tau_exponente = np.floor(np.log10(tau))
tau_base = tau / (10**tau_exponente)

plt.figure(dpi=120, figsize=(8, 6))

plt.gca().axhline(Omega0, color='r', ls='--',
                  label=f'Initial: {Omega0:.2f} rad/day')
plt.gca().axhline(Omega_actual, color='k', ls='--',
                  label=f'Actual: {Omega_actual:.2f} rad/day')

plt.plot(time/365.25, Omega0_example, label='Example')
plt.plot(time/365.25, Omega, label='Mine')

# plt.xscale('log')

plt.title('$d\\Omega/dt = %.2f\\times10^{%d}$ rad/day$^2$' % (
    dOmegadt_base, dOmegadt_exponente))

plt.xlabel('Time [years]')
plt.ylabel('Omega [rad/day]')

plt.legend(title="$\\tau_\\Omega = %.2f\\times10^{%d}$ days" % (
    tau_base, tau_exponente))

plt.show()
