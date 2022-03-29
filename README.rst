Data and analysis for Huang et al, 2021
---------------------------------------

This reproduces the EC-MS figures of:

Junheng Huang et al, **Online Electrochemistry−Mass Spectrometry Evaluation of the
Acidic Oxygen Evolution Reaction at Supported Catalysts** ACS Catal., 11, 12745-12753, **2021**.
https://doi.org/10.1021/acscatal.1c03430

Instructions
............

1. Clone this repository.

2. Install ``ixdat`` version >= 0.2.0. The command is:

    ``pip install --upgrade ixdat``

   For a version of the repository compatible with ixdat v0.1.x, use the
   `ixdat_v0p1 branch <https://github.com/ScottSoren/Huang2021/tree/ixdat_v0p1>`_.

   (ixdat version >= 0.2.1 is needed for data/data_conversion.py)

3. There are two ways to run it. If you wish to run from our raw data, as it was exported from the potentiostat and mass spectrometr, do the following:

     Switch to the `raw data branch <https://github.com/ScottSoren/Huang2021/tree/raw_data>`_ and follow the instructions.
   
   If you are content to run on a reduced dataset that has been assempled and exported with ``ixdat`` in the script ``data/data_conversion.py``, 
   then the scripts will just run once you've cloed the repository!

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
