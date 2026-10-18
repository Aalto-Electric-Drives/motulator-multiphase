"""
4.4-kW dual-three-phase IPMSM, FVC
==================================

This example simulates sensorless flux-vector control (FVC) of an asymmetrical
dual-three-phase interior permanent-magnet synchronous machine (IPMSM) drive. Each
three-phase set corresponds to the 2.2-kW IPMSM of the public examples, and the sets
are displaced by 30°. The common mode (CM) of the sets is controlled using the public
flux-vector controller, parametrized with the parallel-connected equivalent machine.
The differential mode (DM) is current-controlled, and the current sharing between the
sets is changed during the simulation [#Rub2021]_.

"""

# %%
from math import pi

import motulator.drive.control.sm as sm_control
import numpy as np

from motulator_multiphase import control, model, utils

# %%
# Compute base values of the parallel-connected equivalent machine (just for figures).

n_sets = 2
nom = utils.NominalValues(U=370, I=2 * 4.3, f=75, P=2 * 2.2e3, tau=2 * 14)
base = utils.BaseValues.from_nominal(nom, n_p=3)

# %%
# Configure the system model. The parameters describe one three-phase set, when both
# sets carry equal currents. The leakage inductance `L_ls` of a set is not coupled to
# the other set, and it determines the DM dynamics.

par = model.SynchronousMachinePars(n_p=3, R_s=3.6, L_d=0.036, L_q=0.051, psi_f=0.545)
L_ls = 0.005
machine = model.MultiThreePhaseMachine(par, L_ls=L_ls, n_sets=n_sets)
mechanics = model.MechanicalSystem(J=2 * 0.015)
converter = model.MultiThreePhaseConverter(u_dc=540, n_sets=n_sets)
mdl = model.MultiThreePhaseDrive(machine, mechanics, converter)

# %%
# Configure the control system. The public flux-vector controller is used for the CM.

est_par = utils.to_parallel(par, n_sets)
cfg = sm_control.FluxVectorControllerCfg(i_s_max=1.5 * base.i, sensorless=True)
vector_ctrl = sm_control.FluxVectorController(est_par, cfg)
dm_ctrl = control.DifferentialModeController(L_ls, n_sets=n_sets, alpha_c=2 * pi * 200)
speed_ctrl = sm_control.SpeedController(J=2 * 0.015, alpha_s=2 * pi * 4)
ctrl = control.MultiThreePhaseControlSystem(vector_ctrl, dm_ctrl, speed_ctrl)

# %%
# Set the speed reference, the load torque, and the current-sharing coefficients. In
# the end, the first set carries 75% of the total current.

ctrl.set_speed_ref(lambda t: (t > 0.2) * base.w_M)
mdl.mechanics.set_external_load_torque(lambda t: (t > 0.8) * nom.tau)
ctrl.set_current_sharing(lambda t: [1.5, 0.5] if t > 1.2 else [1.0, 1.0])

# %%
# Simulate and plot the results. The first figure shows the CM quantities in per-unit
# values, and the second one the currents of the sets and the DM.

sim = model.Simulation(mdl, ctrl)
res = sim.simulate(t_stop=1.6)
utils.plot(res, base)
utils.plot_sets(res)
print("Torque at the end (Nm):", np.round(res.mdl.machine.tau_M[-1], 2))

# %%
# .. rubric:: References
#
# .. [#Rub2021] Rubino, Dordevic, Armando, Bojoi, Levi, "A novel matrix transformation
#    for decoupled control of modular multiphase PMSM drives," IEEE Trans. Power
#    Electron., 2021, https://doi.org/10.1109/TPEL.2020.3043083
