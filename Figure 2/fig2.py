from pathlib import Path
from matplotlib import pyplot as plt

from ixdat import Measurement
from ixdat.techniques.ec_ms import MSCalResult, ECMSCalibration
from ixdat.constants import FARADAY_CONSTANT

# ---------- import and calibrate the data ------------- #
data_dir = (
    Path(r"C:\Users\scott\Dropbox\WORKSPACES\China\Junheng\published data")
    / r"Figure 2 and 3\the raw data of MS and EC"
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

tspan_1 = [15, 280]
tspan_2 = [15, 300]
# ecms_1.plot(mol_list=["O2", "CO2"], tspan=tspan_1)
# ecms_2.plot(mol_list=["O2", "CO2"], tspan=tspan_2)

LOADING = 0.04e-3 * ecms_1.A_el  # loading in [g], based on 0.04 mg / cm^2


# ------ export the calibrated data ----------- #

# ecms_1.export("reduced_sample.csv", mol_list=["O2", "CO2"])
# ecms_2.export("oxidized_sample.csv", mol_list=["O2", "CO2"])

# ---------- fig 2a and b (Ru/G-450Red) ---------- #

fig_ab, axes_ab = plt.subplots(nrows=3, sharex=True)
axes_ab[0].set_ylabel("O2 cal. signal / (pmol s$^{-1}$)")
axes_ab[1].set_ylabel("CO2 cal. signal / (pmol s$^{-1}$)")
axes_ab[2].set_ylabel("total current / (mA cm${-1}$)")
axes_ab[2].set_xlabel("Potential (V vs RHE)")

t, v = ecms_1.grab("potential", tspan=tspan_1)  # time in [s], potential vs RHE in [V]
t_MS = t + 0  # A delay in the MS data can be applied here (in [s])
#  ^ for example, the average delay for O2 at L=100um is 3.3s [Trimarco2018]
#  ^ needs to be stated clearly and transparently if you choose to do so!
n_dot_O2 = ecms_1.grab_for_t("n_dot_O2", t=t_MS, tspan_bg=[0, 20])  # O2 flux in [mol/s]
n_dot_CO2 = ecms_1.grab_for_t(
    "n_dot_CO2", t=t_MS, tspan_bg=[0, 20]
)  # CO2 flux in [mol/s]
j_total = ecms_1.grab_for_t("raw_current", t=t)  # current in [mA]

axes_ab[0].plot(v, n_dot_O2 * 1e12, "k")
axes_ab[1].plot(v, n_dot_CO2 * 1e12, "brown")
axes_ab[2].plot(v, j_total, "r")

# now we calculate the mass activity.
j_O2 = n_dot_O2 * (4 * FARADAY_CONSTANT)  # OER current in [A]
j_CO2 = n_dot_O2 * (4 * FARADAY_CONSTANT)  # graphene oxidation current in [A]

J_mass_O2 = j_O2 / LOADING  # mass-normalized OER current in [A/g]
J_mass_CO2 = j_CO2 / LOADING  # mass-normalized graphene oxidation current in [A/g]

# We notice that it's proportional to the flux by the factor:
FACTOR = (4 * FARADAY_CONSTANT) / LOADING  # [(A/g) / (mol/s)]

# but we plotted the flux in pmol/s rather than mol/s, so we have to factor that in:
FACTOR_MS = FACTOR * 1e-12  # [(A/g) / (pmol/s)]

# Using that, we add the right y-axes, just scaling them.
right_axes_ab = []
for ax in axes_ab:
    right_ax = ax.twinx()
    right_ax.set_ylim([y * FACTOR_MS for y in ax.get_ylim()])
    right_axes_ab.append(right_ax)

# it's actually a different factor for the current:
FACTOR_J = ecms_1.A_el * 1e-3 / LOADING  # [(A/g) / (mA/cm^2)]
right_axes_ab[2].set_ylim([y * FACTOR_J for y in axes_ab[2].get_ylim()])

right_axes_ab[0].set_ylabel("j$_{O2}$ / (A g$^{-1}$)")
right_axes_ab[1].set_ylabel("j$_{CO2}$ / (A g$^{-1}$)")
right_axes_ab[2].set_ylabel("j$_{total}$ / (A g$^{-1}$)")


fig_ab.set_figheight(fig_ab.get_figwidth() * 2)


# ---------- fig 2c and d (RuO2/G-450Ox)---------- #

fig_cd, axes_cd = plt.subplots(nrows=3, sharex=True)
axes_cd[0].set_ylabel("O2 cal. signal / (pmol s$^{-1}$)")
axes_cd[1].set_ylabel("CO2 cal. signal / (pmol s$^{-1}$)")
axes_cd[2].set_ylabel("total current / (mA cm${-1}$)")
axes_cd[2].set_xlabel("Potential (V vs RHE)")

t, v = ecms_2.grab("potential", tspan=tspan_2)  # time in [s], potential vs RHE in [V]
t_MS = t + 0  # A delay in the MS data can be applied here (in [s])
#  ^ for example, the average delay for O2 at L=100um is 3.3s [Trimarco2018]
#  ^ needs to be stated clearly and transparently if you choose to do so!
n_dot_O2 = ecms_2.grab_for_t("n_dot_O2", t=t_MS, tspan_bg=[0, 20])  # O2 flux in [mol/s]
n_dot_CO2 = ecms_2.grab_for_t(
    "n_dot_CO2", t=t_MS, tspan_bg=[0, 20]
)  # CO2 flux in [mol/s]
j_total = ecms_2.grab_for_t("raw_current", t=t)  # current in [mA]

axes_cd[0].plot(v, n_dot_O2 * 1e12, "k")
axes_cd[1].plot(v, n_dot_CO2 * 1e12, "brown")
axes_cd[2].plot(v, j_total, "r")

# now we calculate the mass activity.
j_O2 = n_dot_O2 * (4 * FARADAY_CONSTANT)  # OER current in [A]
j_CO2 = n_dot_O2 * (4 * FARADAY_CONSTANT)  # graphene oxidation current in [A]

J_mass_O2 = j_O2 / LOADING  # mass-normalized OER current in [A/g]
J_mass_CO2 = j_CO2 / LOADING  # mass-normalized graphene oxidation current in [A/g]

# We notice that it's proportional to the flux by the factor:

FACTOR = (4 * FARADAY_CONSTANT) / LOADING  # [(A/g) / (mol/s)]

# but we plotted the flux in pmol/s rather than mol/s, so we have to factor that in:
FACTOR_MS = FACTOR * 1e-12  # [(A/g) / (pmol/s)]

# Using that, we add the right y-axes, just scaling them.
right_axes_cd = []
for ax in axes_cd:
    right_ax = ax.twinx()
    right_ax.set_ylim([y * FACTOR_MS for y in ax.get_ylim()])
    right_axes_cd.append(right_ax)

# it's actually a different factor for the current:
FACTOR_J = ecms_2.A_el * 1e-3 / LOADING  # [(A/g) / (mA/cm^2)]
right_axes_cd[2].set_ylim([y * FACTOR_J for y in axes_cd[2].get_ylim()])

right_axes_cd[0].set_ylabel("j$_{O2}$ / (A g$^{-1}$)")
right_axes_cd[1].set_ylabel("j$_{CO2}$ / (A g$^{-1}$)")
right_axes_cd[2].set_ylabel("j$_{total}$ / (A g$^{-1}$)")

fig_cd.set_figheight(fig_cd.get_figwidth() * 2)

# ------------------- equalize axes and save figures ------------------------ #
# to make the figures comparable, we set the same y-axis limits on fig_ab and fig_bc:
for i in range(3):
    ylim_ab = axes_ab[i].get_ylim()
    ylim_cd = axes_cd[i].get_ylim()
    ylim = (min((ylim_ab[0], ylim_cd[0])), max((ylim_ab[-1], ylim_cd[-1])))
    axes_ab[i].set_ylim(ylim)
    axes_cd[i].set_ylim(ylim)

    # We also have to rescale the right axes for the MS data...
    right_axes_ab[i].set_ylim([y * FACTOR_MS for y in ylim])
    right_axes_cd[i].set_ylim([y * FACTOR_MS for y in ylim])

# And for the EC data...
right_axes_ab[2].set_ylim([y * FACTOR_J for y in axes_ab[2].get_ylim()])
right_axes_cd[2].set_ylim([y * FACTOR_J for y in axes_cd[2].get_ylim()])

# Now we're ready to save the figures:
fig_ab.savefig("fig2ab.png")
fig_cd.savefig("fig2cd.png")
