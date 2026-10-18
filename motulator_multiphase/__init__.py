"""
Multi-three-phase machine drives.

This package provides models and controls of machines with n three-phase winding
sets (e.g., dual-three-phase machines), based on the decoupled multi-stator (DMS)
transformation [#Rub2021]_. The common mode of the winding sets is equivalent to a
three-phase machine whose sets are connected in parallel. Therefore, the three-phase
machine models, magnetic models, and vector controllers of motulator are reused for
the common mode, while the differential modes are modeled and controlled separately.
The package is organized like `motulator.drive`::

    from motulator_multiphase import control, model, utils

References
----------
.. [#Rub2021] Rubino, Dordevic, Armando, Bojoi, Levi, "A novel matrix transformation
   for decoupled control of modular multiphase PMSM drives," IEEE Trans. Power
   Electron., 2021, https://doi.org/10.1109/TPEL.2020.3043083

"""
