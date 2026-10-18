"""
1.1-kW triple-three-phase SPMSM, CVC
====================================

This example simulates current-vector control (CVC) of a nine-phase surface-mounted
permanent-magnet synchronous machine (SPMSM) with three three-phase sets. The sets are
displaced by 15°, which differs from the standard 20° of asymmetrical nine-phase
machines. The machine parameters are from [#Rub2021]_. The common mode (CM) of the sets
is controlled using the public current-vector controller, parametrized with the
parallel-connected equivalent machine, and the differential modes (DMs) are
current-controlled. The current sharing between the sets is changed during the
simulation.

"""

# %%
from math import pi

import motulator.drive.control.sm as sm_control
import numpy as np

from motulator_multiphase import control, model, utils

# %%
# Configure the system model. In the multi-stator model of [#Rub2021]_, the flux
# linkage of the set k is `psi_k = L_ls*i_k + M*sum(i_z) + psi_f`. When all sets carry
# equal currents, the inductance of a set is `L_ls + n_sets*M`, which is used here.
# The model assumes identical sets, so the stator resistance and the leakage
# inductance of the sets 1 and 3 are used (the set 2 has 7.9 Ω and 10.3 mH).

n_sets = 3
R_s, L_ls, M, psi_f = 8.2, 18.5e-3, 10.5e-3, 0.265
L_s = L_ls + n_sets * M
par = model.SynchronousMachinePars(n_p=3, R_s=R_s, L_d=L_s, L_q=L_s, psi_f=psi_f)
theta_sets = np.radians([0, 15, 30])
machine = model.MultiThreePhaseMachine(
    par, L_ls=L_ls, n_sets=n_sets, theta_sets=theta_sets
)
mechanics = model.MechanicalSystem(J=0.0133)
converter = model.MultiThreePhaseConverter(u_dc=450, n_sets=n_sets)
mdl = model.MultiThreePhaseDrive(machine, mechanics, converter)

# %%
# Configure the control system. The current amplitude of each set is limited to 3.5 A,
# as in [#Rub2021]_.

est_par = utils.to_parallel(par, n_sets)
cfg = sm_control.CurrentVectorControllerCfg(i_s_max=n_sets * 3.5, sensorless=False)
vector_ctrl = sm_control.CurrentVectorController(est_par, cfg)
dm_ctrl = control.DifferentialModeController(L_ls, n_sets=n_sets)
speed_ctrl = sm_control.SpeedController(J=0.0133, alpha_s=2 * pi * 4)
ctrl = control.MultiThreePhaseControlSystem(
    vector_ctrl, dm_ctrl, speed_ctrl, theta_sets=theta_sets
)

# %%
# Set the speed reference, the rated load torque of 7.1 Nm, and the current sharing.
# The rated current is 1.4 A (rms) per set. In the end, the sets carry 50%, 30%, and
# 20% of the total current.

ctrl.set_speed_ref(lambda t: (t > 0.1) * 2 * pi * 1000 / 60)
mdl.mechanics.set_external_load_torque(lambda t: (t > 0.5) * 7.1)
ctrl.set_current_sharing(lambda t: [1.5, 0.9, 0.6] if t > 0.8 else n_sets * [1.0])

# %%
# Simulate and plot the results. The first figure shows the CM quantities, where the
# current is the sum of the set currents, and the second one the currents of the sets
# and the DMs.

sim = model.Simulation(mdl, ctrl)
res = sim.simulate(t_stop=1.2)
utils.plot(res)
utils.plot_sets(res)

# %%
# .. rubric:: References
#
# .. [#Rub2021] Rubino, Dordevic, Armando, Bojoi, Levi, "A novel matrix transformation
#    for decoupled control of modular multiphase PMSM drives," IEEE Trans. Power
#    Electron., 2021, https://doi.org/10.1109/TPEL.2020.3043083
