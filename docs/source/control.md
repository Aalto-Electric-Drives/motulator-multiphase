# Control

This document describes the control system of the {mod}`motulator_multiphase.control` module, implemented in the {class}`motulator_multiphase.control.MultiThreePhaseControlSystem` class. The control is based on the DMS transformation, see {doc}`model`. The CM is controlled by a vector controller of *motulator*, e.g., {class}`motulator.drive.control.sm.FluxVectorController` or {class}`motulator.drive.control.sm.CurrentVectorController`, configured with the parameters of the parallel-connected equivalent machine given by {func}`motulator_multiphase.utils.to_parallel`. The DMs are current-controlled in the coordinates of the vector controller. Each set has its own pulse-width modulator (PWM). The speed-control and torque-control modes work in the same manner as in the vector control systems of *motulator*, see {external+motulator:doc}`control/common/control_system`.

## Feedback Signals

The measured phase currents of the sets are transformed to the CM and the DMs, see {ref}`dms-transformation`. The vector controller of the CM obtains the current $\iss = n\icm^\mathrm{s}$ of the equivalent machine, i.e., the sum of the set currents, and the CM voltage $\uss = \ucm^\mathrm{s}$. The CM voltage is obtained from the realized voltages of the PWMs of the sets in the same manner as in *motulator*: the observers use the voltage of the ongoing sampling period, while the control laws use the average of the voltages of the previous and ongoing sampling periods. The observer of the vector controller gives the angle $\thetac$ and the angular speed $\omegac$ of the controller coordinates. The DM currents are transformed to these coordinates as

```{math}
---
label: dm_controller_coordinates
---
    \idm{m} = \e^{-\jj\thetac}\idm{m}^\mathrm{s}
```

and the DM voltages in the same manner.

## Current Sharing

The current sharing between the sets is defined by the coefficients $\kappa_k$, which sum up to the number of sets, $\sum_{k=1}^{n} \kappa_k = n$. The reference of the set current is $\kappa_k\icm$, where $\icm = \is/n$ is the CM current. Applying the DMS transformation to these references, the DM current references become

```{math}
---
label: dm_current_ref
---
    \idm{m}^\mathrm{ref} = c_m\icm
```

where $c_m = \sum_{k=1}^{n} T_{mk}\kappa_k$. Equal sharing, $\kappa_k = 1$, results in zero DM current references. For example, the coefficients $\kappa_1 = 1.5$ and $\kappa_2 = 0.5$ of a dual-three-phase machine make the first set carry 75% of the total current. The coefficients are set by {meth}`motulator_multiphase.control.MultiThreePhaseControlSystem.set_current_sharing`, either as constants or as a function of time.

## Differential-Mode Current Control

The DM current controller is implemented in the {class}`motulator_multiphase.control.DifferentialModeController` class. In the controller coordinates rotating at $\omegac$, the DM model, see {ref}`differential-modes`, becomes

```{math}
---
label: dm_model_controller
---
    \Lls\frac{\D\idm{m}}{\D t} = \udm{m} - (\Rs + \jj\omegac\Lls)\idm{m}
```

Each DM is controlled by a flux-linkage-based complex-vector 2DOF PI controller {cite}`Bri1999,Awa2019a`, in the same manner as the current controllers of *motulator*, see {external+motulator:ref}`complex-vector-2dof-pi-controller`. The controller is

```{math}
---
label: dm_controller
---
    \frac{\D\mathbf{u}_{\mathrm{i},m}}{\D t} &= (\ki + \jj\omegac\kt)\left(\psidm{m}^\mathrm{ref} - \psidm{m}\right) \\
    \udm{m}^\mathrm{ref} &= \kt\psidm{m}^\mathrm{ref} - \kp\psidm{m} + \mathbf{u}_{\mathrm{i},m}
```

where $\psidm{m}^\mathrm{ref} = \Lls\idm{m}^\mathrm{ref}$ and $\psidm{m} = \Lls\idm{m}$, and $\mathbf{u}_{\mathrm{i},m}$ is the integral state. The gains are

```{math}
---
label: dm_gains
---
    \kt = \alphac \qquad \kp = \alphac + \alphai \qquad \ki = \alphac\alphai
```

where $\alphac$ is the current-control bandwidth and $\alphai$ is the integral-action bandwidth. If the resistance is neglected, the closed-loop system from $\psidm{m}^\mathrm{ref}$ to $\psidm{m}$ becomes the first-order system $\alphac/(s + \alphac)$. As in *motulator*, the anti-windup of the integrator uses the realized voltage.

If the sampling period $\Ts$ is given, the direct discrete-time design of {class}`motulator.common.control.DiscreteComplexController` is used instead. The vector controller of the CM must then be a direct discrete-time design as well, since the angle advance of the voltage references, explained below, is common to all modes.

## Voltage References

The voltage reference $\usref$ of the CM is given by the vector controller and the DM voltage references by {eq}`dm_controller`, all in the controller coordinates. The voltage references of the sets in their own coordinates are obtained using the inverse DMS transformation as

```{math}
---
label: set_voltage_ref
---
    \mathbf{u}_k^\mathrm{ref} = \e^{\jj(\thetac + k_\mathrm{comp}\Ts\omegac - \vartheta_k)}\left(\usref + n\sum_{m=1}^{n-1} T_{mk}\udm{m}^\mathrm{ref}\right)
```

where the angle advance $k_\mathrm{comp}\Ts\omegac$ compensates for the computational and PWM delays, as in *motulator*. The factor $k_\mathrm{comp}$ is given by the vector controller. The PWM of each set computes the duty ratios from its voltage reference and the measured DC-bus voltage. The limited voltages of the PWMs are transformed back to the CM and the DMs in the controller coordinates. The limited CM voltage is used by the vector controller, and the limited DM voltages are used for the anti-windup in the direct discrete-time design.
