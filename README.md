# motulator-multiphase

**Multi-three-phase machine drives for [motulator](https://github.com/Aalto-Electric-Drives/motulator)**

Documentation: https://aalto-electric-drives.github.io/motulator-multiphase/

This package extends *motulator* with models and control of machines having *n* three-phase winding sets (e.g., dual-three-phase machines), each set supplied by its own converter. The models and the control are based on the decoupled multi-stator (DMS) transformation [[1]](#references), which maps the sets to the common mode (CM) and *n* − 1 differential modes (DMs). The CM is equivalent to a three-phase machine whose sets are connected in parallel, so the machine models and vector controllers of *motulator* (synchronous and induction machines, flux-vector and current-vector control, sensored and sensorless) are reused for it. The DMs are resistive-inductive circuits, controlled by complex-vector PI current controllers, and the current sharing between the sets can be changed during operation.

## Installation

The package is not yet on PyPI. Install it from GitHub in a [virtual environment](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/) (*motulator* is installed as a dependency). On Linux and macOS:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install git+https://github.com/Aalto-Electric-Drives/motulator-multiphase
```

On Windows, the environment is created and activated with `py -m venv .venv` and `.venv\Scripts\activate`. The installation from GitHub requires Git.

For development, clone the repository and install the package with `pip install -e ".[dev]"`.

## Usage

A sensorless flux-vector-controlled dual-three-phase IPMSM drive, whose sets are displaced by 30°:

```python
from math import pi

import motulator.drive.control.sm as sm_control
from motulator_multiphase import control, model, utils

# System model: the parameters describe one three-phase set, when both sets carry
# equal currents, and L_ls is the leakage inductance not coupled to the other set
par = model.SynchronousMachinePars(n_p=3, R_s=3.6, L_d=0.036, L_q=0.051, psi_f=0.545)
machine = model.MultiThreePhaseMachine(par, L_ls=0.005, n_sets=2)
converter = model.MultiThreePhaseConverter(u_dc=540, n_sets=2)
mdl = model.MultiThreePhaseDrive(machine, model.MechanicalSystem(J=0.03), converter)

# Control system: the flux-vector controller of motulator for the CM, parametrized
# with the parallel-connected equivalent machine, and a current controller for the DM
est_par = utils.to_parallel(par, 2)
cfg = sm_control.FluxVectorControllerCfg(i_s_max=18, sensorless=True)
vector_ctrl = sm_control.FluxVectorController(est_par, cfg)
dm_ctrl = control.DifferentialModeController(L_ls=0.005, n_sets=2)
speed_ctrl = sm_control.SpeedController(J=0.03, alpha_s=2 * pi * 4)
ctrl = control.MultiThreePhaseControlSystem(vector_ctrl, dm_ctrl, speed_ctrl)

# Speed reference, load torque, and current sharing (75% in the first set at the end)
ctrl.set_speed_ref(lambda t: (t > 0.2) * 150)
mdl.mechanics.set_external_load_torque(lambda t: (t > 0.8) * 28)
ctrl.set_current_sharing(lambda t: [1.5, 0.5] if t > 1.2 else [1.0, 1.0])

res = model.Simulation(mdl, ctrl).simulate(t_stop=1.6)
utils.plot(res)  # CM quantities, the current being the sum of the set currents
utils.plot_sets(res)  # Set currents and DM currents
```

The sets are displaced by 60°/*n* by default (30° for a dual-three-phase machine); other arrangements are given with `theta_sets`. Complete examples are in the [example gallery](https://aalto-electric-drives.github.io/motulator-multiphase/examples/index.html).

## Scope

- Identical sets, any number *n* ≥ 1; a drive with one set equals the three-phase drive of *motulator*.
- The saturated machine models and GradNet flux maps of *motulator* can be used for the CM, see `utils.to_parallel`.
- Each DM is a resistive-inductive circuit of the stator resistance and the leakage inductance of a set; it is not saturated and produces no torque.
- Gate blocking (standby control) is not supported.

## References

1. Rubino, Dordevic, Armando, Bojoi, Levi, "A novel matrix transformation for decoupled control of modular multiphase PMSM drives," IEEE Trans. Power Electron., 2021, https://doi.org/10.1109/TPEL.2020.3043083

## Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->

<!-- prettier-ignore-start -->

<!-- markdownlint-disable -->

<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/SandroRubino"><img src="https://avatars.githubusercontent.com/u/346736963?v=4?s=50" width="50px;" alt="Sandro Rubino"/><br /><sub><b>Sandro Rubino</b></sub></a><br /><a href="#ideas-SandroRubino" title="Ideas, Planning, & Feedback">🤔</a> <a href="#research-SandroRubino" title="Research">🔬</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://research.aalto.fi/en/persons/marko-hinkkanen"><img src="https://avatars.githubusercontent.com/u/76600872?v=4?s=50" width="50px;" alt="Marko Hinkkanen"/><br /><sub><b>Marko Hinkkanen</b></sub></a><br /><a href="https://github.com/Aalto-Electric-Drives/motulator-multiphase/commits?author=mhinkkan" title="Code">💻</a> <a href="#example-mhinkkan" title="Examples">💡</a></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->

<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
