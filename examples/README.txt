Examples
========

These examples demonstrate control of multi-three-phase machine drives, in which each three-phase winding set of the machine is supplied by its own converter. The machine models and the control are based on the decoupled multi-stator (DMS) transformation [#Rub2021]_. The common mode of the sets is equivalent to a three-phase machine whose sets are connected in parallel. Hence, the three-phase machine models and vector controllers of motulator are used for the common mode, while the differential modes are current-controlled. The current sharing between the sets can be changed.

.. rubric:: References

.. [#Rub2021] Rubino, Dordevic, Armando, Bojoi, Levi, "A novel matrix transformation for decoupled control of modular multiphase PMSM drives," IEEE Trans. Power Electron., 2021, https://doi.org/10.1109/TPEL.2020.3043083
