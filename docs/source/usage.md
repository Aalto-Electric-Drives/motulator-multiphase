# Usage

After {doc}`installation`, the package is used in the same manner as *motulator*: a continuous-time system model, a discrete-time control system, and a simulation object are created, see the [usage of *motulator*](https://aalto-electric-drives.github.io/motulator/usage.html). The package is organized like {mod}`motulator.drive`:

```python
from motulator_multiphase import control, model, utils
```

## Quick Start

```{include} ../../README.md
---
start-after: '## Usage'
end-before: '## Scope'
---
```

## Machine Parameters

The machine parameters, such as {class}`motulator.drive.model.SynchronousMachinePars`, describe one three-phase set, when all sets carry equal currents. Hence, the parameters of a three-phase machine can be used for each set. The leakage inductance `L_ls` of a set, which is not coupled to the other sets, is needed in addition, see {doc}`model`. The controllers of *motulator* are parametrized with the equivalent machine, whose sets are connected in parallel, given by {func}`motulator_multiphase.utils.to_parallel`. This function also scales the magnetic models of saturated machines, e.g., GradNet flux maps.

## Simulation Results

The simulation results are organized as in *motulator*, see the [usage of *motulator*](https://aalto-electric-drives.github.io/motulator/usage.html). In addition, the following signals are available:

- `res.mdl.machine.i_c_ab` and `res.mdl.converter.u_c_ab` are the currents and voltages of the sets in their own coordinates, with the shape `(n_sets, N)`.
- `res.mdl.machine.i_dm_ab` and `res.mdl.machine.psi_dm_ab` are the DM currents and flux linkages in stationary coordinates, with the shape `(n_sets - 1, N)`.
- `res.mdl.machine.i_s_ab`, `res.mdl.machine.psi_s_ab`, and `res.mdl.machine.u_s_ab` are the current, flux linkage, and voltage of the parallel-connected equivalent machine. Hence, the current is the sum of the set currents. These signals allow using the plotting functions of three-phase drives, such as {func}`motulator.drive.utils.plot`.
- `res.ctrl.dm` contains the DM signals, see {class}`motulator_multiphase.control.DifferentialModeSignals`, with the shape `(N, n_sets - 1)`. They are in the same coordinates as the other signals of the vector controller, e.g., in the estimated rotor flux coordinates for induction machines.

The function {func}`motulator_multiphase.utils.plot_sets` plots the current magnitudes of the sets and the DM currents with their references.

## Scope

```{include} ../../README.md
---
start-after: '## Scope'
end-before: '## Reference'
---
```
