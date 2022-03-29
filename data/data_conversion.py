"""This module reads in raw data and exports it as .csv's for use by other scripts."""

from pathlib import Path
from ixdat import Measurement

# Change this path to the location of the raw data on your own computer.
DATA_DIRECTORY = Path(r"~/Dropbox/DATA/Huang2021").expanduser()

f23 = DATA_DIRECTORY / "Figure 2 and 3/the raw data of MS and EC"
f4 = DATA_DIRECTORY / "Figure 4/the raw data of MS and EC"

for (name, folder, raw_EC_data, raw_MS_data, tspan) in [
    ("fig3ab", f23, "RuG-450Red CV-46cyc.txt", "timescan-RuG-450Red CV-46cyc.txt", [0, 2500]),
    ("fig3cd", f23, "RuO2 G-450OX-10cyc.txt", "timescan-RuO2G-450Ox-CV-10cyc.txt", [0, 2500]),
    ("fig4b", f4, "01-Ru rgo 450H-50Ag.txt", "timescan-sem190-LOW.txt", [-50, 4200]),
    ("fig4c", f4, "02-RuO2 rgo 450oxo-50Ag.txt", "timescan-sem190-LOW.txt", [-50, 4200])
]:
    ec = Measurement.read(folder / raw_EC_data, reader="chi")
    ms = Measurement.read(folder / raw_MS_data, reader="rgasoft")
    ecms = ec + ms
    ecms.tstamp += ecms.t[0]  # necessary as CHI records the end and not start timestamp
    ecms.plot(tspan=tspan)  # a sanity check here.
    ecms.export(name + ".csv", tspan=tspan)
