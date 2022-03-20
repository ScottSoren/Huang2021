from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid

from matplotlib import pyplot as plt

from ixdat import Measurement
from ixdat.techniques.ec_ms import ECMSCalibration
from ixdat.constants import FARADAY_CONSTANT

# ---------- import and calibrate the data ------------- #
data_dir = (
    Path(r"~/Dropbox/DATA/Huang2021").expanduser()
    / r"Figure 4\the raw data of MS and EC"
)

ms = Measurement.read(data_dir / "timescan-sem190-LOW.txt", reader="rgasoft")
# ms.plot()  # plot MS data to check if import worked


ec_1 = Measurement.read(data_dir / "01-Ru rgo 450H-50Ag.txt", reader="chi")
# ec_1.plot()  # plot EC data to check if import worked
ecms_1 = ec_1 + ms
# ecms_1.plot()  # overview plot

ec_2 = Measurement.read(data_dir / "02-RuO2 rgo 450oxo-50Ag.txt", reader="chi")
# ec_2.plot()
ecms_2 = ec_2 + ms
# ecms_2.plot()

# define and apply the MS calibration:
ecms_1.add_calibration(ECMSCalibration.read("../calibration/calibration_3.ix"))
# ecms_1.plot(mol_list=["O2", "CO2"])  # overview plot with calibrated data
ecms_2.add_calibration(ECMSCalibration.read("../calibration/calibration_3.ix"))
# ecms_2.plot(mol_list=["O2", "CO2"])  # overview plot with calibrated data

# ------- define the timespans (using the overview plots) --------- #
# first, get the start of the data to t=0
ecms_1.tstamp += ecms_1.t[0]
ecms_2.tstamp += ecms_2.t[0]

tspan = [0, 4200]

ecms_1.set_bg(tspan_bg=[-20, 0])
ecms_2.set_bg(tspan_bg=[-20, 0])

ecms_1.plot(mol_list=["O2", "CO2"], tspan=tspan)
ecms_2.plot(mol_list=["O2", "CO2"], tspan=tspan)

fig_a, axes_a = plt.subplots(nrows=3, sharex=True)
fig_b, ax_b = plt.subplots()
fig_c, ax_c = plt.subplots()
axes_a[2].set_xlabel("time / (s)")
axes_a[0].set_ylabel("FE O2 / (%)")
axes_a[1].set_ylabel("FE CO2 / (%)")
axes_a[2].set_ylabel("U$_{RHE}$ / (V)")
ax_b.set_xlabel("time / (s)")
ax_b.set_ylabel("charge / (mC cm$^{-2}$)")
ax_c.set_xlabel("time / (s)")
ax_c.set_ylabel("charge / (mC cm$^{-2}$)")

for (meas, color_a, ax_bc) in [(ecms_1, "r", ax_b), (ecms_2, "b", ax_c)]:
    t, v = meas.grab("potential", tspan=tspan)
    j = meas.grab_for_t("current", t=t)
    n_dot_O2 = meas.grab_for_t("n_dot_O2", t=t)
    n_dot_CO2 = meas.grab_for_t("n_dot_CO2", t=t)

    j_O2 = n_dot_O2 * 4 * FARADAY_CONSTANT * 1e3 / meas.A_el  # [mA/cm^2]
    FE_O2 = j_O2 / j

    j_CO2 = n_dot_CO2 * 4 * FARADAY_CONSTANT * 1e3 / meas.A_el  # [mA/cm^2]
    FE_CO2 = j_CO2 / j

    axes_a[0].plot(FE_O2 * 100, color=color_a)
    axes_a[1].plot(FE_CO2 * 100, color=color_a)
    axes_a[2].plot(v, color=color_a)

    q = cumulative_trapezoid(j, t, initial=0)  # [mC/cm^2]
    q_O2 = cumulative_trapezoid(j_O2, t, initial=0)  # [mC/cm^2]
    q_CO2 = cumulative_trapezoid(j_CO2, t, initial=0)  # [mC/cm^2]

    ax_bc.plot(t, q, "r")
    ax_bc.plot(t, q_O2, "k")
    ax_bc.plot(t, q_CO2, "brown")
    ax_bc.plot(t, q_O2 + q_CO2, "b--")


fig_a.savefig("fig4a.png")
fig_b.savefig("fig4b.png")
fig_c.savefig("fig4c.png")
