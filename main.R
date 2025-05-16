lib_path = "./"

verbose =
    TRUE
    # FALSE

MPI = NULL

Colors_of_HM = c(
    "CTRIP"="#A88D72", #marron
    "EROS"="#CECD8D", #vert clair
    "GRSD"="#619C6C", #vert foncé
    "J2000"="#74AEB9", #bleu clair
    "MORDOR-SD"="#D8714E", #orange
    "MORDOR-TS"="#AE473E", #rouge
    "ORCHIDEE"="#EFA59D", #rose
    "SIM2"="#475E6A", #bleu foncé
    "SMASH"="#F6BA62", #mimosa

    "AquiFR"="#AF3FA5", #violet
    "EROS Bretagne"="#CECD8D", #vert clair
    "MONA"="#F5D80E" #jaune
)

Colors_of_storylines =
    c("HadGEM2-ES|historical-rcp85|ALADIN63|ADAMONT"="#569A71", #vert
      "CNRM-CM5|historical-rcp85|ALADIN63|ADAMONT"="#EECC66", #jaune
      "EC-EARTH|historical-rcp85|HadREM3-GA7|ADAMONT"="#E09B2F", #orange
      "HadGEM2-ES|historical-rcp85|CCLM4-8-17|ADAMONT"="#791F5D" #violet
      )
Colors_light_of_storylines = # 60% lighter
    c("HadGEM2-ES|historical-rcp85|ALADIN63|ADAMONT"="#BAD8C6", #vert
      "CNRM-CM5|historical-rcp85|ALADIN63|ADAMONT"="#F8EBC2", #jaune
      "EC-EARTH|historical-rcp85|HadREM3-GA7|ADAMONT"="#F3D7AC", #orange
      "HadGEM2-ES|historical-rcp85|CCLM4-8-17|ADAMONT"="#E9A9D5" #violet
      )

storylines =
    c("HadGEM2-ES|historical-rcp85|ALADIN63|ADAMONT"="Réchauffement marqué et augmentation des précipitations", #vert
      "CNRM-CM5|historical-rcp85|ALADIN63|ADAMONT"="Changements futurs relativement peu marqués", #jaune
      "EC-EARTH|historical-rcp85|HadREM3-GA7|ADAMONT"="Fort réchauffement et fort assèchement en été (et en annuel)", #orange
      "HadGEM2-ES|historical-rcp85|CCLM4-8-17|ADAMONT"="Fort réchauffement et forts contrastes saisonniers en précipitations" #violet
      )

logo_info = list(
    "EX2"=c(file='LogoExplore2.png', y=0.4, height=0.7, width=0.5)
)

n_projections_by_code = 4

# If the hydrological network needs to be plot
river_selection =
    # NULL
    c("la Durance", "la Marne", "la Vienne", "le Loir", "la Loire",
      "l'Oise", "la Seine", "le Lot", "l'Adour", "le Rhône",
      "la Moselle", "l'Aisne", "la Garonne", "le Tarn", "le Doubs",
      "la Dordogne", "la Charente", "le Cher", "la Saône", "l'Allier",
      "Fleuve la Loire", "la Meuse", "la Sarthe", "la Somme",
      "l'Isère", "la Vilaine", "l'Aude", "l'Yonne")
river_selection = paste0("^", river_selection, "$")
river_length =
    # NULL
    30000
# 300000

# Tolerance of the simplification algorithm for shapefile in sf
toleranceRel =
    1000 # normal map
    # 9000 # mini map





#  ___        _  _    _        _  _            _    _            
# |_ _| _ _  (_)| |_ (_) __ _ | |(_) ___ __ _ | |_ (_) ___  _ _  
#  | | | ' \ | ||  _|| |/ _` || || |(_-</ _` ||  _|| |/ _ \| ' \ 
# |___||_||_||_| \__||_|\__,_||_||_|/__/\__,_| \__||_|\___/|_||_| ____
##### /!\ Do not touch if you are not aware #####
## 0. LIBRARIES ______________________________________________________
# Computer
computer = Sys.info()["nodename"]
print(paste0("Computer ", computer))
computer_file_list = list.files(path=lib_path,
                                pattern="computer[_].*[.]R")
computer_list = gsub("(computer[_])|([.]R)", "", computer_file_list)
computer_file = computer_file_list[sapply(computer_list,
                                          grepl, computer)]
computer_path = file.path(lib_path, computer_file)
print(paste0("So reading file ", computer_path))
source(computer_path, encoding='UTF-8')

# Sets working directory
setwd(computer_work_path)

library(dplyr)

devtools::load_all("../../dataSHEEP_project/dataSHEEP/")
dev_path = "../../SHEEPfold_project/SHEEPfold/__SHEEP__"
list_path = list.files(dev_path, pattern='*.R$', full.names=TRUE, recursive=TRUE)
for (path in list_path) {
    source(path, encoding='UTF-8')
}

add_path = function (x) {
    x = c(x, file.path(resources_path, logo_dir, x["file"]))
    names(x)[length(x)] = "path"
    return (x)
}
logo_info = lapply(logo_info, add_path)
icon_path = file.path(resources_path, icon_dir)


Stations_path = file.path(archive_data_path,
                          stations_selection_file)
Stations = ASHE::read_tibble(Stations_path)
Stations = filter(Stations, n_rcp85 >=4)

Projections_path = file.path(archive_data_path,
                             projections_selection_file)
Projections = ASHE::read_tibble(Projections_path)
Projections = filter(Projections,
                     EXP == "historical-rcp85" &
                     grepl("ADAMONT", BC))

Variables_hydro_path = file.path(archive_data_path,
                                 variables_hydro_selection_file)
Variables_hydro = ASHE::read_tibble(Variables_hydro_path)


dataEX_serie = NULL
metaEX_serie = NULL

dataEX_criteria_climate =
    ASHE::read_tibble(file.path(climate_data_dirpath, climate_data_file))
dataEX_criteria_climate$EXP = "historical-rcp85" 
dataEX_criteria_climate$BC = "ADAMONT"
dataEX_criteria_climate$HM = NA

dataEX_criteria_climate$variable =
    paste0(dataEX_criteria_climate$Variable,
           "_", gsub("seas[-]", "",
                     dataEX_criteria_climate$Saison))
dataEX_criteria_climate =
    dataEX_criteria_climate %>%
    select(-Variable, -Saison, -surface) %>%
    rename(SH=ZH) %>%
    rename(GWL=NivRechauf) %>%
    relocate(variable, .after=RCM) %>%
    relocate(EXP, .before=GCM) %>%
    relocate(BC, .before=GCM)

Ok = dataEX_criteria_climate$RCM == "SMHI-RCA4"
dataEX_criteria_climate$RCM[Ok] = "RCA4"
Ok = grepl("REMO", dataEX_criteria_climate$RCM)
dataEX_criteria_climate$RCM[Ok] = "REMO"

dataEX_criteria_climate = tidyr::pivot_wider(dataEX_criteria_climate,
                                             values_from=delta,
                                             names_from=variable,
                                             names_prefix="delta_")

metaEX_criteria_climate =
    tibble(variable=c("RR_DJF", "RR_JJA", "TMm_DJF", "TMm_JJA"),
           name=c("Précipitations hivernales",
                  "Précipitations estivales",
                  "Température moyenne hivernale",
                  "Température moyenne estivale"))

dataEX_criteria = dataEX_criteria_climate #full_join

dataEX_criteria$climateChain = paste(dataEX_criteria$GCM,
                                     dataEX_criteria$EXP,
                                     dataEX_criteria$RCM,
                                     dataEX_criteria$BC, sep="|")
dataEX_criteria$Chain = paste(dataEX_criteria$climateChain,
                              dataEX_criteria$HM, sep="|")


metaEX_criteria = bind_rows(metaEX_criteria_climate)

# stop()
if (!exists("Shapefiles")) {
    ASHE::post("### Loading shapefiles")

    Shapefiles = load_shapefile(
        computer_shp_path, Code=NULL,
        france_shp_path=france_shp_path,
        bassinHydro_shp_path=bassinHydro_shp_path,
        regionHydro_shp_path=regionHydro_shp_path,
        secteurHydro_shp_path=secteurHydro_shp_path,
        river_shp_path=river_shp_path,
        river_selection=river_selection,
        river_length=river_length,
        toleranceRel=toleranceRel)
}


sheet_projection_secteur(
    Stations,
    dataEX_serie,
    metaEX_serie,
    dataEX_criteria,
    metaEX_criteria,
    Colors=Colors_of_storylines,
    Colors_light=Colors_light_of_storylines,
    Names=storylines,
    historical=historical,
    icon_path=icon_path,
    logo_info=logo_info,
    Shapefiles=Shapefiles,
    figdir=figdir,
    Pages=NULL,
    verbose=subverbose)


# confiance au dessus de 80%


### /!\ ?
# PAS de Z
# RCM : REMO2015 et REMO2009 -> REMO2009
# RCM : SMHI-RCA4 -> RCA4
