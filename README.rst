Data and analysis for Huang et al, 2021
---------------------------------------

This reproduces the EC-MS figures of:

Junheng Huang et al, **Online Electrochemistry−Mass Spectrometry Evaluation of the
Acidic Oxygen Evolution Reaction at Supported Catalysts** ACS Catal., 11, 12745-12753, **2021**.
https://doi.org/10.1021/acscatal.1c03430

Instructions
............

1. Clone this repository.

2. Install ``ixdat`` version >= 0.1.6. The command is:

  ``pip install --upgrade ixdat``

3. Download the data onto your computer from this Dropbox link:
   https://www.dropbox.com/sh/7u5ffi94upkvphc/AADlAWi2QtNJsjxV4YZ9b1yka?dl=0
   **We will make the data available in another way soon**

4. Change the path in the ``data_dir`` variable at the top of each script to the appropriate location on your computer.

5. Run the scripts! `Cite the article <https://doi.org/10.1021/acscatal.1c03430>`_ if you find this useful!

The scripts
...........

- Figure 2: `Figure 2/fig2.py <https://github.com/ScottSoren/Huang2021/blob/main/Figure%202/fig2.py>`_

  This script plots calibrated signals vs potential for the first cycles of the two catalyst types.
  It makes use of equivalent left and right y-axes to show both molar fluxes and partial current
  densities normalized to ruthenium loading.

- Figure 3: `Figure 3/fig3.py <https://github.com/ScottSoren/Huang2021/blob/main/Figure%203/fig3.py>`_

  This script plots calculates and plots the cumulative partial current densities and Faradaic efficiencies for several
  cycles in cyclic voltammatry of the two catalyst types. It demonstrates the use of ixdat's
  CyclicVoltammatry objects.

- Figure 3: `Figure 4/fig4.py <https://github.com/ScottSoren/Huang2021/blob/main/Figure%204/fig4.py>`_

  This script plots calculates and plots the Faradaic efficiencies and cumulative partial current
  densities for constant-current water electrolysis of two catalyst types. The failure
  of the reduced catalyst is quite dramatic!