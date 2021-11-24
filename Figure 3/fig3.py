from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid

from matplotlib import pyplot as plt
from matplotlib import gridspec

from ixdat import Measurement
from ixdat.techniques.ec_ms import ECMSCalibration
from ixdat.constants import FARADAY_CONSTANT

# ---------- import and calibrate the data ------------- #
data_dir = (
    Path(r"C:\Users\scott\Dropbox\WORKSPACES\China\Junheng\published data")
    / r"Figure 2\the raw data of MS and EC"
)

ec_1 = Measurement.read(data_dir / "RuG-450Red CV-46cyc.txt", reader="chi")
# ec_1.plot()  # plot EC data to check if import worked
ms_1 = Measurement.read(data_dir / "timescan-RuG-450Red CV-46cyc.txt", reader="rgasoft")
# ms_1.plot()  # plot MS data to check if import worked
ecms_1 = ec_1 + ms_1
# ecms_1.plot()  # overview plot

ec_2 = Measurement.read(data_dir / "RuO2 G-450OX-10cyc.txt", reader="chi")
# ec_2.plot()
ms_2 = Measurement.read(
    data_dir / "timescan-RuO2G-450Ox-CV-10cyc.txt", reader="rgasoft"
)
# ms_2.plot()
ecms_2 = ec_2 + ms_2
# ecms_2.plot()

# define and apply the MS calibration:
ecms_1.calibration = ECMSCalibration.read("../calibration/calibration_1.ix")
# ecms_1.plot(mol_list=["O2", "CO2"])  # overview plot with calibrated data
ecms_2.calibration = ECMSCalibration.read("../calibration/calibration_2.ix")
# ecms_2.plot(mol_list=["O2", "CO2"])  # overview plot with calibrated data

# ------- define the timespans (using the overview plots) --------- #
# first, get the start of the data to t=0
ecms_1.tstamp += ecms_1.t[0]
ecms_2.tstamp += ecms_2.t[0]

tspan_1 = [0, 2100]
tspan_2 = [0, 1600]

ecms_1.set_bg(tspan_bg=[0, 20])
ecms_2.set_bg(tspan_bg=[0, 20])

# ecms_1.plot(mol_list=["O2", "CO2"], tspan=tspan_1)
# ecms_2.plot(mol_list=["O2", "CO2"], tspan=tspan_2)

# ----------- figures a,b, and c (Ru/G-450Red) --------------- #

# prepare the axes using a gridspec:
fig_abc = plt.figure()

gs_abc = gridspec.GridSpec(8, 1, fig_abc)
# gs.update(hspace=0.025)
ax_a_total = plt.subplot(gs_abc[0:1, 0])  # total current axis
ax_a_O2 = plt.subplot(gs_abc[1:2, 0])  # OER current axis
ax_a_CO2 = plt.subplot(gs_abc[2:3, 0])  # carbon oxidation current axis
ax_b = plt.subplot(gs_abc[3:6, 0])  # cumulative charge axis
ax_b_right = ax_b.twinx()  # potential axis
ax_c = plt.subplot(gs_abc[7:8, 0])  # FE axis

fig_abc.set_figheight(fig_abc.get_figwidth() * 1.5)

ax_a_total.set_xlabel("time / s")
ax_a_total.xaxis.set_label_position("top")
ax_c.set_xlabel("time / s")

ax_a_total.set_ylabel("J$_{total}$ / (mA cm$^{-2}$)")
ax_a_O2.set_ylabel("J$_{O2}$ / (mA cm$^{-2}$)")
ax_a_CO2.set_ylabel("J$_{CO2}$ / (mA cm$^{-2}$)")
ax_b.set_ylabel("FE / (%)")
ax_b_right.set_ylabel("charge / (mC cm$^{-2}$)")
ax_c.set_ylabel("U$_{RHE}$ / (V)")
ax_a_total.tick_params(
    axis="x", top=True, bottom=True, labeltop=True, labelbottom=False
)
ax_a_O2.tick_params(
    axis="x", top=True, bottom=True, labeltop=False, labelbottom=False
)
ax_a_CO2.tick_params(
    axis="x", top=True, bottom=True, labeltop=False, labelbottom=False
)


# plot the current and partial currents

t, j = ecms_1.grab("current", tspan=tspan_1)  # in [mA/cm^2]
n_dot_O2 = ecms_1.grab_for_t("n_dot_O2", t)  # in [mol/s]
n_dot_CO2 = ecms_1.grab_for_t("n_dot_CO2", t)  # in [mol/s]
j_O2 = n_dot_O2 * 4 * FARADAY_CONSTANT * 1e3 / ecms_1.A_el  # in [mA/cm^2]
j_CO2 = n_dot_CO2 * 4 * FARADAY_CONSTANT * 1e3 / ecms_1.A_el  # in [mA/cm^2]

ax_a_total.plot(t, j, "r")
ax_a_O2.plot(t, j_O2, "k")
ax_a_CO2.plot(t, j_CO2, "brown")


# plot the faradaic efficiencies:
cv_1 = ecms_1.cut(tspan_1).as_cv()
# iterate cycle each anodic scan through 0.8 V_RHE:
cv_1.redefine_cycle(start_potential=0.8, redox=True)
# cv_1.plot_measurement(J_str="cycle")  # to see where cycle 1 starts
cycle_numbers = list(range(1, 8))
centers = []
for c in cycle_numbers:
    cycle = cv_1[c]
    t_c, j_c = cycle.grab("current")  # in [mA/cm^2]
    Q = np.trapz(j_c, t_c)  # in mC/cm^2
    n_dot_c_O2 = cycle.grab_for_t("n_dot_O2", t_c)  # in [mol/s]
    n_dot_c_CO2 = cycle.grab_for_t("n_dot_CO2", t_c)  # in [mol/s]
    j_c_O2 = n_dot_c_O2 * 4 * FARADAY_CONSTANT * 1e3 / ecms_1.A_el  # in [mA/cm^2]
    j_c_CO2 = n_dot_c_CO2 * 4 * FARADAY_CONSTANT * 1e3 / ecms_1.A_el  # in [mA/cm^2]
    FE_c_O2 = np.trapz(j_c_O2, t_c) / Q   # ratio of O2 current to total current
    FE_c_CO2 = np.trapz(j_c_CO2, t_c) / Q   # ratio of CO2 current to total current
    center = np.mean(t_c)
    centers.append(center)
    width = (t_c[-1] - t_c[0]) / 2
    ax_b.bar(center, FE_c_O2 * 100, color="k", alpha=0.5, width=width)
    ax_b.bar(
        center, FE_c_CO2 * 100, bottom=FE_c_O2 * 100, color="brown", alpha=0.5, width=width
    )

ax_b.set_xticks(centers)
ax_b.set_xticklabels(cycle_numbers)
ax_b.set_xlabel("cycle number")
ax_b.set_ylim([0, 100])

q_total = cumulative_trapezoid(j, t, initial=0)  # [mC/cm^2]
q_O2 = cumulative_trapezoid(j_O2, t, initial=0)  # [mC/cm^2]
q_CO2 = cumulative_trapezoid(j_CO2, t, initial=0)  # [mC/cm^2]

ax_b_right.plot(t, q_total, "r")
ax_b_right.plot(t, q_O2, "k")
ax_b_right.plot(t, q_CO2, "brown")
ax_b_right.tick_params(
    axis="x", top=True, bottom=False, labeltop=False, labelbottom=False
)

v = ecms_1.grab_for_t("potential", t)
ax_c.plot(t, v, "b")

# ----------- figures d,e, and f (RuO2/G-450Ox) --------------- #

# prepare the axes using a gridspec:
fig_def = plt.figure()

gs_def = gridspec.GridSpec(8, 1, fig_def)
# gs.update(hspace=0.025)
ax_d_total = plt.subplot(gs_abc[0:1, 0])  # total current axis
ax_d_O2 = plt.subplot(gs_abc[1:2, 0])  # OER current axis
ax_d_CO2 = plt.subplot(gs_abc[2:3, 0])  # carbon oxidation current axis
ax_e = plt.subplot(gs_abc[3:6, 0])  # cumulative charge axis
ax_e_right = ax_e.twinx()  # potential axis
ax_f = plt.subplot(gs_abc[7:8, 0])  # FE axis

fig_def.set_figheight(fig_def.get_figwidth() * 1.5)

ax_d_total.set_xlabel("time / s")
ax_d_total.xaxis.set_label_position("top")
ax_f.set_xlabel("time / s")

ax_d_total.set_ylabel("J$_{total}$ / (mA cm$^{-2}$)")
ax_d_O2.set_ylabel("J$_{O2}$ / (mA cm$^{-2}$)")
ax_d_CO2.set_ylabel("J$_{CO2}$ / (mA cm$^{-2}$)")
ax_e.set_ylabel("FE / (%)")
ax_e_right.set_ylabel("charge / (mC cm$^{-2}$)")
ax_f.set_ylabel("U$_{RHE}$ / (V)")
ax_d_total.tick_params(
    axis="x", top=True, bottom=True, labeltop=True, labelbottom=False
)
ax_d_O2.tick_params(
    axis="x", top=True, bottom=True, labeltop=False, labelbottom=False
)
ax_d_CO2.tick_params(
    axis="x", top=True, bottom=True, labeltop=False, labelbottom=False
)


# plot the current and partial currents

t, j = ecms_2.grab("current", tspan=tspan_2)  # in [mA/cm^2]
n_dot_O2 = ecms_2.grab_for_t("n_dot_O2", t)  # in [mol/s]
n_dot_CO2 = ecms_2.grab_for_t("n_dot_CO2", t)  # in [mol/s]
j_O2 = n_dot_O2 * 4 * FARADAY_CONSTANT * 1e3 / ecms_2.A_el  # in [mA/cm^2]
j_CO2 = n_dot_CO2 * 4 * FARADAY_CONSTANT * 1e3 / ecms_2.A_el  # in [mA/cm^2]

ax_d_total.plot(t, j, "r")
ax_d_O2.plot(t, j_O2, "k")
ax_d_CO2.plot(t, j_CO2, "brown")


# plot the faradaic efficiencies:
cv_1 = ecms_2.cut(tspan_2).as_cv()
# iterate cycle each anodic scan through 0.8 V_RHE:
cv_1.redefine_cycle(start_potential=0.8, redox=True)
# cv_1.plot_measurement(J_str="cycle")  # to see where cycle 1 starts
cycle_numbers = list(range(1, 6))
centers = []
for c in cycle_numbers:
    cycle = cv_1[c]
    t_c, j_c = cycle.grab("current")  # in [mA/cm^2]
    Q = np.trapz(j_c, t_c)  # in mC/cm^2
    n_dot_c_O2 = cycle.grab_for_t("n_dot_O2", t_c)  # in [mol/s]
    n_dot_c_CO2 = cycle.grab_for_t("n_dot_CO2", t_c)  # in [mol/s]
    j_c_O2 = n_dot_c_O2 * 4 * FARADAY_CONSTANT * 1e3 / ecms_2.A_el  # in [mA/cm^2]
    j_c_CO2 = n_dot_c_CO2 * 4 * FARADAY_CONSTANT * 1e3 / ecms_2.A_el  # in [mA/cm^2]
    FE_c_O2 = np.trapz(j_c_O2, t_c) / Q   # ratio of O2 current to total current
    FE_c_CO2 = np.trapz(j_c_CO2, t_c) / Q   # ratio of CO2 current to total current
    center = np.mean(t_c)
    centers.append(center)
    width = (t_c[-1] - t_c[0]) / 2
    ax_e.bar(center, FE_c_O2 * 100, color="k", alpha=0.5, width=width)
    ax_e.bar(
        center, FE_c_CO2 * 100, bottom=FE_c_O2 * 100, color="brown", alpha=0.5, width=width
    )

ax_e.set_xticks(centers)
ax_e.set_xticklabels(cycle_numbers)
ax_e.set_xlabel("cycle number")
ax_e.set_ylim([0, 100])

q_total = cumulative_trapezoid(j, t, initial=0)  # [mC/cm^2]
q_O2 = cumulative_trapezoid(j_O2, t, initial=0)  # [mC/cm^2]
q_CO2 = cumulative_trapezoid(j_CO2, t, initial=0)  # [mC/cm^2]

ax_e_right.plot(t, q_total, "r")
ax_e_right.plot(t, q_O2, "k")
ax_e_right.plot(t, q_CO2, "brown")
ax_e_right.tick_params(
    axis="x", top=True, bottom=False, labeltop=False, labelbottom=False
)

v = ecms_2.grab_for_t("potential", t)
ax_f.plot(t, v, "b")


# ------ Get the corresponding axes y limts to be equal and save the figs! --------- #

for ax_abc, ax_def in [
    (ax_a_total, ax_d_total), (ax_a_O2, ax_d_O2), (ax_a_CO2, ax_d_CO2),
    (ax_b, ax_e), (ax_b_right, ax_e_right), (ax_c, ax_f)
]:

    ylim_abc = ax_abc.get_ylim()
    ylim_def = ax_def.get_ylim()
    ylim = (min((ylim_abc[0], ylim_def[0])), max((ylim_abc[-1], ylim_def[-1])))
    ax_abc.set_ylim(ylim)
    ax_def.set_ylim(ylim)


fig_abc.savefig("fig3_abc.png")
fig_def.savefig("fig3_def.png")