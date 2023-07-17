import sys
import os
import pandas as pd

from pathlib import Path

from ..lib import config
from .future_land_use_simulator_dialog import FutureLandUseSimulatorDialog


def dist_populate(gui):
    reg_name = gui.dlg.cb_region.currentText()
    print(reg_name)

    if reg_name == "N/A":
        gui.dlg.cb_district.addItems(config.test)
    elif reg_name == "Ahafo":
        gui.dlg.cb_district.addItems(config.ahafo)
    elif reg_name == "Ashanti":
        gui.dlg.cb_district.clear()
        gui.dlg.cb_district.addItems(config.ashanti)
    elif reg_name == "Bono":
        gui.dlg.cb_district.clear()
        gui.dlg.cb_district.addItems(config.bono)
    elif reg_name == "Bono East":
        gui.dlg.cb_district.clear()
        gui.dlg.cb_district.addItems(config.bono_east)
    elif reg_name == "Central":
        gui.dlg.cb_district.clear()
        gui.dlg.cb_district.addItems(config.central)
    elif reg_name == "Eastern":
        gui.dlg.cb_district.clear()
        gui.dlg.cb_district.addItems(config.eastern)
    elif reg_name == "Greater Accra":
        gui.dlg.cb_district.clear()
        gui.dlg.cb_district.addItems(config.accra)
    elif reg_name == "North East":
        gui.dlg.cb_district.clear()
        gui.dlg.cb_district.addItems(config.ne)
    elif reg_name == "Northern":
        gui.dlg.cb_district.clear()
        gui.dlg.cb_district.addItems(config.north)
    elif reg_name == "Oti":
        gui.dlg.cb_district.clear()
        gui.dlg.cb_district.addItems(config.oti)
    elif reg_name == "Savannah":
        gui.dlg.cb_district.clear()
        gui.dlg.cb_district.addItems(config.sav)
    elif reg_name == "Upper East":
        gui.dlg.cb_district.clear()
        gui.dlg.cb_district.addItems(config.ue)
    elif reg_name == "Upper West":
        gui.dlg.cb_district.clear()
        gui.dlg.cb_district.addItems(config.uw)
    elif reg_name == "Volta":
        gui.dlg.cb_district.clear()
        gui.dlg.cb_district.addItems(config.volta)
    elif reg_name == "Western":
        gui.dlg.cb_district.clear()
        gui.dlg.cb_district.addItems(config.western)
    elif reg_name == "Western North":
        gui.dlg.cb_district.clear()
        gui.dlg.cb_district.addItems(config.west_north)
    else:
        gui.dlg.cb_district.clear()


def csv_calc(gui):
    root_dir = (Path(__file__).parents[1]).joinpath('data/population.csv')
    yr = gui.dlg.cb_year.currentText()

    proj_pop_field = gui.dlg.le_proj_pop
    reg = gui.dlg.cb_region.currentText()
    dsct = gui.dlg.cb_district.currentText()

    pop_df = pd.read_csv(root_dir)

    pop_project = ""

    if yr in config.yr_selectable:
        if dsct == "Entire Region":
            try:
                pop_project = pop_df.loc[pop_df["reg_prop"] == reg, f"{yr}_dif"].sum()
            except IndexError:
                pass
            proj_pop_field.setText(str(pop_project))
        else:  # district specific
            try:
                pop_project = pop_df.loc[pop_df["dist"] == dsct, f"{yr}_dif"].values[0]
            except IndexError:
                pass
            proj_pop_field.setText(str(pop_project))
    else:  # custom year
        pass


def disable_gpkg(gui):
    shp = gui.dlg.check_shp
    gpk = gui.dlg.check_gpk

    shp_state = shp.checkState()

    if shp_state == 0:
        gpk.setCheckable(True)
    elif shp_state == 2:
        gpk.setCheckable(False)
    else:
        pass


def disable_shp(gui):
    shp = gui.dlg.check_shp
    gpk = gui.dlg.check_gpk

    gpk_state = gpk.checkState()

    if gpk_state == 0:
        shp.setCheckable(True)
    elif gpk_state == 2:
        shp.setCheckable(False)
    else:
        pass


def disable_study_area(gui):
    dsct_txt = gui.dlg.cb_district.currentText()

    if dsct_txt == "N/A":
        gui.dlg.file_area.setEnabled(True)
    else:
        gui.dlg.file_area.setEnabled(False)
