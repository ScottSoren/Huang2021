"""This module just saves the calibrations in ixdat's formats.

The sensitivity factors were determined by experiments like that in Figure S6.
"""

from ixdat.techniques.ms import MSCalResult
from ixdat.techniques.ec_ms import ECMSCalibration


A_el = 0.196  # reference potential on RHE scale in [V]
RE_vs_RHE = 0.260  # electrode area in [cm^2]

# for Ru/G-450Red: Fig2. ab, Fig3 abc
calibration_1 = ECMSCalibration(
    ms_cal_results=[
        MSCalResult(mol="O2", mass="M32", cal_type="internal", F=0.02321959502828152),
        MSCalResult(mol="CO2", mass="M44", cal_type="internal", F=0.017006001217645464),
    ],
    RE_vs_RHE=RE_vs_RHE,
    A_el=A_el,
)
calibration_1.export("./calibration_1.ix")

# for RuO2/G-450Ox: Fig2 cd, Fig3 def
calibration_2 = ECMSCalibration(
    ms_cal_results=[
        MSCalResult(mol="O2", mass="M32", cal_type="internal", F=0.0124),
        MSCalResult(mol="CO2", mass="M44", cal_type="internal", F=0.0097),
    ],
    RE_vs_RHE=RE_vs_RHE,
    A_el=A_el,
)
calibration_2.export("./calibration_2.ix")

# Fig4. for both Ru/G-450Red and RuO2/G-450Ox

calibration_3 = ECMSCalibration(
    ms_cal_results=[
        MSCalResult(mol="O2", mass="M32", cal_type="internal", F=0.0495593221803036),
        MSCalResult(mol="CO2", mass="M44", cal_type="internal", F=0.046182030605598154),
    ],
    RE_vs_RHE=RE_vs_RHE,
    A_el=A_el,
)
calibration_3.export("./calibration_3.ix")
