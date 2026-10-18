# *motulator-multiphase:* Multi-Three-Phase Machine Drives for *motulator*

This package extends [*motulator*](https://aalto-electric-drives.github.io/motulator/) with models and control of machines having $n$ three-phase winding sets (e.g., dual-three-phase machines), each set supplied by its own converter. The models and the control are based on the decoupled multi-stator (DMS) transformation {cite}`Rub2021`, which maps the sets to the common mode (CM) and $n - 1$ differential modes (DMs). The CM is equivalent to a three-phase machine whose sets are connected in parallel. Hence, the machine models, magnetic models, and vector controllers of *motulator* are reused for the CM, while the DMs are modeled and controlled separately. The current sharing between the sets can be changed during operation.

This documentation covers only the multi-three-phase extension. The machine models, control methods, and conventions common with three-phase drives are described in the [documentation of *motulator*](https://aalto-electric-drives.github.io/motulator/).

```{note}
The package is in early development, and its interfaces may change before version 1.0. It currently supports *motulator* 0.14, and the support is extended after each new release of *motulator* has been tested.
```

```{toctree}
---
titlesonly:
caption: Getting Started
name: getting_started
maxdepth: 1
---
installation
usage
```

```{toctree}
---
titlesonly:
caption: Examples
name: examples
maxdepth: 1
---
examples/index
```

```{toctree}
---
titlesonly:
caption: Theory
name: theory
maxdepth: 1
---
model
control
```

```{toctree}
---
titlesonly:
caption: Reference
name: reference
maxdepth: 1
---
API <autoapi/motulator_multiphase/index>
bibliography
motulator <https://aalto-electric-drives.github.io/motulator/>
```
