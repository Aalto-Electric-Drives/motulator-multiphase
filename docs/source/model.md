# System Model

This document describes the continuous-time models of the {mod}`motulator_multiphase.model` module. Peak-valued complex space vectors are used, as in *motulator*, see {external+motulator:doc}`model/common/space_vectors`.

## Space Vectors of the Sets

The machine has $n$ three-phase winding sets, $k = 1, \ldots, n$. The phase a of the set $k$ is located at the angle $\vartheta_k$ with respect to the stationary $\alpha$ axis. The space vector of the set current in its own coordinates is

```{math}
---
label: set_vector
---
    \mathbf{i}_k = \frac{2}{3}\left(i_{\mathrm{a}k} + i_{\mathrm{b}k}\e^{\jj 2\pi/3} + i_{\mathrm{c}k}\e^{\jj 4\pi/3}\right)
```

where $i_{\mathrm{a}k}$, $i_{\mathrm{b}k}$, and $i_{\mathrm{c}k}$ are the phase currents of the set. In the common stationary coordinates, the space vector of the set is $\mathbf{i}_k^\mathrm{s} = \e^{\jj\vartheta_k}\mathbf{i}_k$. The voltages and flux linkages are treated in the same manner. By default, the sets are displaced by $\pi/(3n)$, i.e., $\vartheta_k = (k - 1)\pi/(3n)$, see {func}`motulator_multiphase.utils.default_set_angles`. For example, the sets of an asymmetrical dual-three-phase machine are displaced by 30°. Other arrangements are given with the parameter `theta_sets`.

(dms-transformation)=

## Decoupled Multi-Stator Transformation

The decoupled multi-stator (DMS) transformation maps the space vectors of the sets to the common mode (CM) and $n - 1$ differential modes (DMs) {cite}`Rub2021`. The CM is the mean of the sets,

```{math}
---
label: dms_cm
---
    \icm^\mathrm{s} = \frac{1}{n}\sum_{k=1}^{n} \mathbf{i}_k^\mathrm{s}
```

and the DMs are

```{math}
---
label: dms_dm
---
    \idm{m}^\mathrm{s} = \sum_{k=1}^{n} T_{mk}\mathbf{i}_k^\mathrm{s}
```

where $m = 1, \ldots, n - 1$. The coefficients are

```{math}
---
label: dms_coefficients
---
    T_{mk} = \frac{1}{\sqrt{n(n - m)(n - m + 1)}}
    \begin{cases}
        n - m & \text{if } k = m \\
        -1 & \text{if } k > m \\
        0 & \text{if } k < m
    \end{cases}
```

The coefficients of each DM sum up to zero, i.e., the DMs are orthogonal to the CM. The inverse transformation is

```{math}
---
label: dms_inverse
---
    \mathbf{i}_k^\mathrm{s} = \icm^\mathrm{s} + n\sum_{m=1}^{n-1} T_{mk}\idm{m}^\mathrm{s}
```

For example, a dual-three-phase machine has $\icm^\mathrm{s} = (\mathbf{i}_1^\mathrm{s} + \mathbf{i}_2^\mathrm{s})/2$ and $\idm{1}^\mathrm{s} = (\mathbf{i}_1^\mathrm{s} - \mathbf{i}_2^\mathrm{s})/2$, and the set currents are $\mathbf{i}_1^\mathrm{s} = \icm^\mathrm{s} + \idm{1}^\mathrm{s}$ and $\mathbf{i}_2^\mathrm{s} = \icm^\mathrm{s} - \idm{1}^\mathrm{s}$. The transformation matrix, whose first row gives the CM and the other rows the DMs, is given by {func}`motulator_multiphase.utils.dms_matrix`.

## Machine Model

The machine model is implemented in the {class}`motulator_multiphase.model.MultiThreePhaseMachine` class. Its parameters describe one three-phase set, when all sets carry equal currents, i.e., when only the CM is present. Furthermore, the leakage inductance $\Lls$ of a set, which is not coupled to the other sets, is needed for the DMs.

### Common Mode

The CM is equivalent to a three-phase machine whose sets are connected in parallel. The voltage and the flux linkage of this equivalent machine are those of the CM, $\uss = \ucm^\mathrm{s}$ and $\psiss = \psicm^\mathrm{s}$, while its current is the sum of the set currents,

```{math}
---
label: parallel_current
---
    \iss = n\icm^\mathrm{s} = \sum_{k=1}^{n} \mathbf{i}_k^\mathrm{s}
```

Consequently, the resistances and inductances of the equivalent machine are those of a set divided by $n$, and its core-loss conductance is multiplied by $n$. For example, the stator resistance of the equivalent machine is $\Rs/n$, where $\Rs$ is the stator resistance of a set. The magnetic models of saturated machines are scaled accordingly. These parameters are given by {func}`motulator_multiphase.utils.to_parallel`. The equivalent machine is modeled with the three-phase machine models of *motulator*, i.e., {class}`motulator.drive.model.SynchronousMachine` or {class}`motulator.drive.model.InductionMachine`, see {external+motulator:doc}`model/drive/synchronous_machine` and {external+motulator:doc}`model/drive/induction_machine`. The equivalent machine produces the electromagnetic torque $\tauM$ of the whole machine.

(differential-modes)=

### Differential Modes

Each DM is a resistive-inductive circuit, consisting of the stator resistance $\Rs$ and the leakage inductance $\Lls$ of a set. In stationary coordinates, the DM is modeled as

```{math}
---
label: dm_model
---
    \frac{\D\psidm{m}^\mathrm{s}}{\D t} = \udm{m}^\mathrm{s} - \Rs\idm{m}^\mathrm{s}
```

where the flux linkage is $\psidm{m}^\mathrm{s} = \Lls\idm{m}^\mathrm{s}$. The DMs are not saturated, and they produce no torque. The state variables of the machine model are those of the equivalent machine, followed by the DM flux linkages. The set currents in their own coordinates are obtained from {eq}`dms_inverse` as

```{math}
---
label: set_currents
---
    \mathbf{i}_k = \e^{-\jj\vartheta_k}\left(\frac{1}{n}\iss + n\sum_{m=1}^{n-1} T_{mk}\idm{m}^\mathrm{s}\right)
```

## Converter Model

The {class}`motulator_multiphase.model.MultiThreePhaseConverter` class models $n$ three-phase converters sharing a DC bus of constant voltage $\udc$. The voltage of the set $k$ in its own coordinates is

```{math}
---
label: set_voltage
---
    \mathbf{u}_k = \mathbf{q}_k\udc
```

where $\mathbf{q}_k$ is the switching-state space vector of the set. The converter takes a modulator of *motulator*, i.e., {class}`motulator.common.model.ZOH` (default) or {class}`motulator.common.model.CarrierComparison`, which converts the duty ratios of the legs into their gate states, including the computational delay and the optional dead time. As in the three-phase converter model of *motulator*, the switches are ideal except for the dead time. Gate blocking is not supported.

## Drive Model

The {class}`motulator_multiphase.model.MultiThreePhaseDrive` class interconnects the converter, the machine, and a mechanical system model of *motulator*, i.e., {class}`motulator.drive.model.MechanicalSystem`, {class}`motulator.drive.model.TwoMassMechanicalSystem`, or {class}`motulator.drive.model.ExternalRotorSpeed`.
