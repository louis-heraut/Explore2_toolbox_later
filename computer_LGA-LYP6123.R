# Copyright 2021-2023 Louis Héraut (louis.heraut@inrae.fr)*1,
#                     Éric Sauquet (eric.sauquet@inrae.fr)*1
#
# *1   INRAE, France
#
# This file is part of Explore2 R toolbox.
#
# Explore2 R toolbox is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# Explore2 R toolbox is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Explore2 R toolbox.
# If not, see <https://www.gnu.org/licenses/>.


# Main script that regroups all command lines needed to interact with
# this toolbox. Choose your parameters before executing all the script
# (RStudio : Ctrl+Alt+R) or line by line.


#  ___         __                         _    _                
# |_ _| _ _   / _| ___  _ _  _ __   __ _ | |_ (_) ___  _ _   ___
#  | | | ' \ |  _|/ _ \| '_|| '  \ / _` ||  _|| |/ _ \| ' \ (_-<
# |___||_||_||_|  \___/|_|  |_|_|_|\__,_| \__||_|\___/|_||_|/__/ _____
# If you want to contact the author of the code you need to contact
# first Louis Héraut who is the main developer. If it is not possible,
# Éric Sauquet is the main referent at INRAE to contact.
#
# Louis Héraut : <louis.heraut@inrae.fr>
# Éric Sauquet : <eric.sauquet@inrae.fr>
#
# See the 'README.md' file for more information about the utilisation
# of this toolbox.


#   ___                          _             
#  / __| ___  _ __   _ __  _  _ | |_  ___  _ _ 
# | (__ / _ \| '  \ | '_ \| || ||  _|/ -_)| '_|
#  \___|\___/|_|_|_|| .__/ \_,_| \__|\___||_|  
## 1. INFO ________ |_| ______________________________________________
# Work path
computer_work_path = '/home/lheraut/Documents/INRAE/projects/Explore2_project/Explore2_toolbox_later'

# Library path for package dev
dev_lib_path = '/home/lheraut/Documents/INRAE/project/'

## 2. INPUT DIRECTORIES ______________________________________________
# archive_data_path = "/media/lheraut/Explore2"
archive_data_path = "/home/lheraut/Bureau"

### 2.1 Hydro _____________________________________________________
hydro_data_dirpath = file.path(archive_data_path, "")
projections_selection_file = "chaines_simulations_Explore2.csv"
stations_selection_file = "stations_selection.csv"
variables_hydro_selection_file = "variables_selection.csv"

### 2.2 Climate ___________________________________________________
climate_data_dirpath = file.path(computer_work_path, "climate_data")
climate_data_file = "climatZH.csv"

### 2.3 Recharge ___________________________________________________
recharge_data_dirpath = file.path(archive_data_path, "")


### 2.3. Resources ___________________________________________________
resources_path = file.path(computer_work_path, 'resources')
logo_dir = 'logo'
icon_dir = 'icon/SVG'
font_dir = 'fonts'

### 2.4. Shapefile ________________________________________________
computer_shp_path =
    '/home/lheraut/Documents/INRAE/data/map'
# Path to the shapefile for france contour from 'computer_data_path' 
france_shp_path = 'france/gadm36_FRA_0.shp'
# Path to the shapefile for basin shape from 'computer_data_path' 
bassinHydro_shp_path = 'bassinHydro/bassinHydro.shp'
# Path to the shapefile for sub-basin shape from 'computer_data_path' 
regionHydro_shp_path = 'regionHydro/regionHydro.shp'
# Path to the shapefile for sub-basin shape from 'computer_data_path' 
secteurHydro_shp_path = 'secteurHydro/secteurHydro.shp'
# Path to the shapefile for station basins shape from 'computer_data_path' 
entiteHydro_shp_path = c('entiteHydro/BV_4207_stations.shp',
                         'entiteHydro/3BVs_FRANCE_L2E_2018.shp')
river_shp_path = 'coursEau/CoursEau_FXX.shp'
# piezo
entitePiezo_shp_path = "entitePiezo_niveau1_extension/entitePiezo_niveau1_extension.shp"

## 3. OUTPUT DIRECTORIES _____________________________________________
### 3.0. Info ________________________________________________________
today = format(Sys.Date(), "%Y_%m_%d")
now = format(Sys.time(), "%H_%M_%S")

### 3.1. Results _____________________________________________________
resdir = file.path(computer_work_path, 'results')
today_resdir = file.path(resdir, today)
now_resdir = file.path(today_resdir, now)

### 3.2. Figures  ____________________________________________________
figdir = file.path(computer_work_path, 'figures')
today_figdir = file.path(figdir, today)
now_figdir = file.path(today_figdir, now)

### 3.3. Tmp  ________________________________________________________
tmpdir = "tmp"
