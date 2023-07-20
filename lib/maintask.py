import os
from pathlib import Path

from qgis.utils import iface
from qgis.core import *
from qgis.core import QgsTask, QgsProcessingParameterVectorDestination, QgsVectorLayer

from PyQt5.QtGui import *
from PyQt5.QtCore import *

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon

import rasterio
from rasterio import mask
from rasterstats import zonal_stats

import pylusat
from pylusat import zonal, geotools

import logging

from . import config
from . import lucis
from . import ga

from ..gui.future_land_use_simulator_dialog import FutureLandUseSimulatorDialog

from ..gui.textlogger import QPlainTextEditLogger


class FutureLandUseSimulatorTask(QgsTask):
    """Future Land Use Simulation"""

    FUTURE_SCENARIO = "FUTURE_SCENARIO"

    def __init__(self, description, gui):
        super().__init__(description, QgsTask.CanCancel)
        self.gui = gui

        self.logTextBox = QPlainTextEditLogger(self.gui.qpt_edit)
        self.logTextBox.setFormatter(
            logging.Formatter("%(asctime)s - %(message)s", "%Y-%m-%d %H:%M:%S")
        )
        logging.getLogger().addHandler(self.logTextBox)
        logging.getLogger().setLevel(logging.INFO)

    def run(self):
        try:
            prog = self.gui.progressBar

            prog.setMinimum(0)
            prog.setMaximum(100)

            temp_csv = self.gui.file_slopes.filePath()
            temp_df = pd.read_csv(temp_csv)
            temp_border_str = temp_df["TEMP_PER_BORDER"].values[0]
            temp_ag_str = temp_df["TEMP_PER_AG"].values[0]
            temp_con_str = temp_df["TEMP_PER_CON"].values[0]
            temp_urb_str = temp_df["TEMP_PER_URB"].values[0]
            temp_water_str = temp_df["TEMP_PER_WATER"].values[0]
            temp_herb_str = temp_df["TEMP_PER_HERB"].values[0]
            temp_constraint = self.gui.le_temp.text()
            prog.setValue(3)
            logging.info("Reading temperature slopes")

            config.TEMP_PER_BORDER = float(temp_border_str)
            config.TEMP_PER_AG = float(temp_ag_str)
            config.TEMP_PER_CON = float(temp_con_str)
            config.TEMP_PER_URB = float(temp_urb_str)
            config.TEMP_PER_WATER = float(temp_water_str)
            config.TEMP_PER_HERB = float(temp_herb_str)
            config.TEMP_TARGET = float(temp_constraint)

            current_lc_path = self.gui.file_lc.filePath()
            ag_suit_path = self.gui.file_ag_suit.filePath()
            con_suit_path = self.gui.file_con_suit.filePath()
            urb_suit_path = self.gui.file_urb_suit.filePath()
            prog.setValue(15)

            dist_name = self.gui.cb_district.currentText()
            proj_pop_txt = self.gui.le_proj_pop.text()
            pop_rast_path = self.gui.file_pop.filePath()
            region_name = self.gui.cb_region.currentText()
            logging.info("Grabbing files")

            urb_ex = self.gui.check_urb_ex.checkState()
            ag_ex = self.gui.check_ag_ex.checkState()
            con_ex = self.gui.check_con_ex.checkState()
            herb_ex = self.gui.check_herb_ex.checkState()
            excess_lu_list = []
            if urb_ex == 2:
                excess_lu_list.append(3)
            if ag_ex == 2:
                excess_lu_list.append(1)
            if con_ex == 2:
                excess_lu_list.append(2)
            if herb_ex == 2:
                excess_lu_list.append(5)
            excess_lu = tuple(excess_lu_list)
            config.LU_EXCESS = excess_lu

            urb_insuf = self.gui.check_urb_insuf.checkState()
            ag_insuf = self.gui.check_ag_insuf.checkState()
            con_insuf = self.gui.check_con_insuf.checkState()
            herb_insuf = self.gui.check_herb_insuf.checkState()
            insuf_lu_list = []
            if urb_insuf == 2:
                insuf_lu_list.append(3)
            if ag_insuf == 2:
                insuf_lu_list.append(1)
            if con_insuf == 2:
                insuf_lu_list.append(2)
            if herb_insuf == 2:
                insuf_lu_list.append(5)
            insuf_lu = tuple(insuf_lu_list)
            config.LU_INSUFFICIENT = insuf_lu
            logging.info("Evaluating excess and insufficient land covers")

            study_area_path = ""
            logging.info("Entering study area analysis")

            if region_name == "N/A":
                study_area_path = self.gui.file_area.filePath()
                filename, file_extension = os.path.splitext(study_area_path)
                if file_extension == ".shp":
                    gdf = gpd.read_file(study_area_path)
                    gdf.drop(
                        columns=[
                            "right",
                            "left",
                            "bottom",
                            "top",
                            "LC_majorit",
                            "ag",
                            "con",
                            "urb",
                        ],
                        axis=1,
                        inplace=True,
                    )
                else:
                    gdf = gpd.read_file(study_area_path, layer="kumasi_area")
                    gdf.drop(
                        columns=[
                            "right_",
                            "left_",
                            "bottom",
                            "top",
                            "LC_majorit",
                            "ag",
                            "con",
                            "urb",
                        ],
                        axis=1,
                        inplace=True,
                    )

                mgrs_gdf = gpd.read_file(study_area_path)
                config.base_df = mgrs_gdf

                row_dif_val = (
                    config.base_df["right"].max() - config.base_df["right"].min()
                )
                row_div_val = row_dif_val / 250 + 1
                n_row = int(row_div_val)

                col_dif_val = (
                    config.base_df["bottom"].max() - config.base_df["bottom"].min()
                )
                col_div_val = col_dif_val / 250 + 1
                n_col = int(col_div_val)

                if self.isCanceled():
                    self.cancel()
                    return False
                else:
                    pass

                logging.info("Creating bounding box around study area")

                mgrs_gdf = pylusat.zonal.zonal_stats_raster(
                    mgrs_gdf,
                    current_lc_path,
                    stats="majority",
                    stats_prefix="lc",
                    nodata=0,
                )
                logging.info("Calculating majority land cover per cell")
                mgrs_gdf = pylusat.zonal.zonal_stats_raster(
                    mgrs_gdf, pop_rast_path, stats="sum", stats_prefix="pop", nodata=0
                )
                logging.info("Calculating sum of population per cell")
                pop_current_val = mgrs_gdf["pop_sum"].sum()
                mgrs_gdf["pop_sum"].fillna(0, inplace=True)
                mgrs_gdf["lc_majority"].fillna(0, inplace=True)

                config.temp_slope = pd.Series(
                    (
                        config.TEMP_PER_AG,
                        config.TEMP_PER_CON,
                        config.TEMP_PER_URB,
                        config.TEMP_PER_WATER,
                        config.TEMP_PER_HERB,
                    ),
                    index=np.arange(1, 6),
                )

                self.gui.file_area.setEnabled(True)
            elif dist_name == "Entire Region":
                current_lc_path = self.gui.file_lc.filePath()
                lc_raster_data = rasterio.open(current_lc_path)

                ghana_shp = (Path(__file__).parents[1]).joinpath('data/Ghana.shp')
                region_gdf = gpd.read_file(ghana_shp)
                reg_name_str = str(config.reg_dict[region_name])
                reg_sel_gdf = region_gdf.loc[(region_gdf["REGION"] == reg_name_str)]
                # create tuple to grab bounds of selected district
                boundary_tup = reg_sel_gdf.total_bounds
                minx = boundary_tup[0]
                miny = boundary_tup[1]
                maxx = boundary_tup[2]
                maxy = boundary_tup[3]
                # create bounding box with district bondaries
                poly = Polygon([(minx, miny), (minx, maxy), (maxx, maxy), (maxx, miny)])
                out_img_arr, out_transform_aff = rasterio.mask.mask(
                    dataset=lc_raster_data, shapes=[poly], crop=True
                )
                out_img_arr = out_img_arr[0]
                if self.isCanceled():
                    self.cancel()
                    return False
                else:
                    pass

                gdf = pylusat.geotools.gridify(
                    reg_sel_gdf, cell_x=250, cell_y=250, n_cols=None, n_rows=None
                )
                logging.info("Creating bounding box around study area")
                zonal_output = zonal_stats(
                    vectors=gdf,
                    raster=out_img_arr,
                    nodata=0,
                    affine=out_transform_aff,
                    stats="majority",
                    all_touched=True,
                )
                zone = pd.DataFrame(zonal_output)

                gdf = gdf.join(zone, lsuffix="_l", rsuffix="_r")
                gdf.rename(
                    columns={"majority": "lc_majority"},
                    inplace=True,
                )

                gdf = pylusat.zonal.zonal_stats_raster(
                    gdf, pop_rast_path, stats="sum", stats_prefix="pop", nodata=0
                )
                logging.info("Calculating sum of population per cell")

                pop_current_val = gdf["pop_sum"].sum()
                gdf["pop_sum"].fillna(0, inplace=True)
                gdf["lc_majority"].fillna(0, inplace=True)
                mgrs_gdf = gdf

                boundary_tup = mgrs_gdf.total_bounds
                row_dif_val = boundary_tup[3] - boundary_tup[1]
                row_div_val = row_dif_val / 250
                n_row = int(row_div_val)

                col_dif_val = boundary_tup[2] - boundary_tup[0]
                col_div_val = col_dif_val / 250
                n_col = int(col_div_val)
                if self.isCanceled():
                    self.cancel()
                    return False
                else:
                    pass

                config.temp_slope = pd.Series(
                    (
                        config.TEMP_PER_BORDER,
                        config.TEMP_PER_AG,
                        config.TEMP_PER_CON,
                        config.TEMP_PER_URB,
                        config.TEMP_PER_WATER,
                        config.TEMP_PER_HERB,
                    ),
                    index=np.arange(0, 6),
                )

                self.gui.file_area.setEnabled(False)
            else:
                current_lc_path = self.gui.file_lc.filePath()
                lc_raster_data = rasterio.open(current_lc_path)

                ghana_shp = (Path(__file__).parents[1]).joinpath('data/Ghana.shp')
                dist_gdf = gpd.read_file(ghana_shp)
                dist_name_str = str(config.dist_dict[dist_name])
                dist_sel_gdf = dist_gdf.loc[(dist_gdf["DISTRICT"] == dist_name_str)]
                prog.setValue(26)

                # create tuple to grab bounds of selected district
                boundary_tup = dist_sel_gdf.total_bounds
                minx = boundary_tup[0]
                miny = boundary_tup[1]
                maxx = boundary_tup[2]
                maxy = boundary_tup[3]

                # create bounding box with district bondaries
                poly = Polygon([(minx, miny), (minx, maxy), (maxx, maxy), (maxx, miny)])
                prog.setValue(28)
                out_img_arr, out_transform_aff = rasterio.mask.mask(
                    dataset=lc_raster_data, shapes=[poly], crop=True
                )
                out_img_arr = out_img_arr[0]

                if self.isCanceled():
                    self.cancel()
                    return False
                else:
                    pass
                gdf = pylusat.geotools.gridify(
                    dist_sel_gdf, cell_x=250, cell_y=250, n_cols=None, n_rows=None
                )
                zonal_output = zonal_stats(
                    vectors=gdf,
                    raster=out_img_arr,
                    nodata=0,
                    affine=out_transform_aff,
                    stats="majority",
                    all_touched=True,
                )
                zone = pd.DataFrame(zonal_output)
                gdf = gdf.join(zone, lsuffix="_l", rsuffix="_r")
                gdf.rename(
                    columns={"majority": "lc_majority"},
                    inplace=True,
                )
                logging.info("Calculating majority land cover per cell")

                gdf = pylusat.zonal.zonal_stats_raster(
                    gdf, pop_rast_path, stats="sum", stats_prefix="pop", nodata=0
                )
                logging.info("Calculating sum of population per cell")

                pop_current_val = gdf["pop_sum"].sum()
                gdf["pop_sum"].fillna(0, inplace=True)
                gdf["lc_majority"].fillna(0, inplace=True)

                mgrs_gdf = gdf

                boundary_tup = mgrs_gdf.total_bounds
                row_dif_val = boundary_tup[3] - boundary_tup[1]
                row_div_val = row_dif_val / 250
                n_row = int(row_div_val)

                col_dif_val = boundary_tup[2] - boundary_tup[0]
                col_div_val = col_dif_val / 250
                n_col = int(col_div_val)

                if self.isCanceled():
                    self.cancel()
                    return False
                else:
                    pass

                config.temp_slope = pd.Series(
                    (
                        config.TEMP_PER_BORDER,
                        config.TEMP_PER_AG,
                        config.TEMP_PER_CON,
                        config.TEMP_PER_URB,
                        config.TEMP_PER_WATER,
                        config.TEMP_PER_HERB,
                    ),
                    index=np.arange(0, 6),
                )

                self.gui.file_area.setEnabled(False)
                prog.setValue(37)

            mgrsbb = study_area_path
            prog.setValue(40)

            mgrs_gdf = pylusat.zonal.zonal_stats_raster(
                mgrs_gdf, ag_suit_path, stats="mean", stats_prefix="ag", nodata=0
            )
            mgrs_gdf = pylusat.zonal.zonal_stats_raster(
                mgrs_gdf, con_suit_path, stats="mean", stats_prefix="con", nodata=0
            )
            mgrs_gdf = pylusat.zonal.zonal_stats_raster(
                mgrs_gdf, urb_suit_path, stats="mean", stats_prefix="urb", nodata=0
            )
            logging.info("Calculating average suitability per cell")

            prog.setValue(45)

            mgrs_gdf.rename(
                columns={
                    "lc_majority": "current_lu",
                    "con_mean": "con_suit",
                    "urb_mean": "urb_suit",
                    "ag_mean": "ag_suit",
                },
                inplace=True,
            )

            config.base_df = mgrs_gdf

            mask = mgrs_gdf["current_lu"].isin([1, 2, 3, 5])
            prog.setValue(55)

            mgrs_gdf["ag_pref"] = pd.cut(
                mgrs_gdf["ag_suit"],
                bins=[1, 1 + 8 / 3, 1 + 16 / 3, 9],
                labels=[1, 2, 3],
                include_lowest=True,
            ).astype(int)
            mgrs_gdf["con_pref"] = pd.cut(
                mgrs_gdf["con_suit"],
                bins=[1, 1 + 8 / 3, 1 + 16 / 3, 9],
                labels=[1, 2, 3],
                include_lowest=True,
            ).astype(int)
            mgrs_gdf["urb_pref"] = pd.cut(
                mgrs_gdf["urb_suit"],
                bins=[1, 1 + 8 / 3, 1 + 16 / 3, 9],
                labels=[1, 2, 3],
                include_lowest=True,
            ).astype(int)
            mgrs_gdf.loc[mask, "urban_conflict"] = mgrs_gdf[mask].apply(
                lambda x: lucis.LandUseConflict(
                    x["ag_pref"], x["con_pref"], x["urb_pref"]
                ).urban_conflict(),
                axis=1,
            )
            logging.info("Masking rasters for specific land cover classes")

            prog.setValue(60)

            mgrs_gdf.loc[(mgrs_gdf["current_lu"] == 3)]

            if self.isCanceled():
                self.cancel()
                return False
            else:
                pass

            config.GRID_R = n_row
            config.GRID_C = n_col
            proj_pop_txt = self.gui.le_proj_pop.text()
            config.PPL_PER_URB = 514
            config.PPL_CURRENT = pop_current_val
            config.PPL_GROWTH = int(float(proj_pop_txt)) - int(float(pop_current_val))
            config.min_conf = 0
            config.max_conf = config.GRID_C * config.GRID_R

            # --- Column names for the input DataFrame ---
            config.AG_SUIT = "ag_suit"  # agriculture suitability
            config.CON_SUIT = "con_suit"  # conservation suitability
            config.URB_SUIT = "urb_suit"  # urban suitability
            config.CURRENT_LU = "current_lu"  # current land use
            prog.setValue(75)

            # --- GA parameters ---
            cxpb_str = self.gui.le_cxpb.text()
            config.CXPB = float(cxpb_str)
            mut_str = self.gui.le_mutpb.text()
            config.MUTPB = float(mut_str)
            gen_str = self.gui.le_gen.text()
            config.GEN_NUM = float(gen_str)
            pop_str = self.gui.le_init_pop.text()
            config.POP_SIZE = int(pop_str)
            logging.info("Reading GA parameters")

            # --- Initial sampling probability for each land use ---

            prog.setValue(90)
            ag_low = float(self.gui.le_ag_low.text())
            ag_med = float(self.gui.le_ag_med.text())
            ag_high = float(self.gui.le_ag_high.text())
            con_low = float(self.gui.le_con_low.text())
            con_med = float(self.gui.le_con_med.text())
            con_high = float(self.gui.le_con_high.text())
            urb_low = float(self.gui.le_urb_low.text())
            urb_med = float(self.gui.le_urb_med.text())
            urb_high = float(self.gui.le_urb_high.text())
            config.INIT_PROB_AG = {1: ag_low, 2: ag_med, 3: ag_high}
            config.INIT_PROB_CON = {1: con_low, 2: con_med, 3: con_high}
            config.INIT_PROB_URB = {1: urb_low, 2: urb_med, 3: urb_high}
            sample_str = self.gui.le_init_sample.text()
            config.INIT_SAMPLE_SIZE = int(
                sample_str
            )  # number of cells changing its land use

            if self.isCanceled():
                self.cancel()
                return False
            else:
                pass

            self.gui.scrollArea.setWidgetResizable(True)

            dist_text = self.gui.cb_district.currentText()
            dist_text_str = str(dist_text)

            ga_result = ga.main()

            new_lu = self.gui.le_out_field.text()

            if self.isCanceled():
                self.cancel()
                return False
            else:
                pass

            config.base_df[new_lu] = ""
            config.base_df[new_lu] = ga_result[0]

            prog.setValue(100)

            return True

        except Exception as e:
            self.exception = e
            return False

    def finished(self, result):
        if result:
            project = QgsProject.instance()
            shp = self.gui.check_shp
            gpk = self.gui.check_gpk

            shp_state = shp.checkState()
            gpk_state = gpk.checkState()

            if self.isCanceled():
                self.cancel()
                return False
            else:
                pass

            if shp_state == 2:
                output_file = self.gui.output.filePath()
                config.base_df.to_file(output_file, driver="ESRI Shapefile")
                lyr = QgsVectorLayer(output_file, self.FUTURE_SCENARIO, "ogr")
                project.addMapLayer(lyr)
            elif gpk_state == 2:
                output_file = self.gui.output.filePath()
                config.base_df.to_file(output_file, driver="GPKG")
                iface.addVectorLayer(output_file, self.FUTURE_SCENARIO, "ogr")
            else:
                pass
        else:
            pass

    def cancel(self):
        QgsMessageLog.logMessage(
            "{name} was canceled".format(name=self.description()), level=Qgis.Info
        )
        super().cancel()
