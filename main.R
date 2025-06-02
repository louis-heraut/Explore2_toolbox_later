lib_path = "./"

verbose =
    TRUE
# FALSE

subverbose =
    TRUE
# FALSE

MPI = NULL

logo_info = list(
    "Explore2"=c(file='LogoExplore2.png'),
    "TRACC"=c(file='la_france_s_adapte.png')
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
toleranceRel_normal = 1000
toleranceRel_mini = 9000



# grep("aterial", systemfonts::system_fonts()$family, value=TRUE)


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
library(ggplot2)
library(latex2exp)


# if (!exists("font")) {
#     library(extrafont)
#     font_import(paths=file.path(resources_path, "fonts"),
#                 prompt=FALSE)
#     # loadfonts()
#     loadfonts(device="pdf")
#     font = TRUE
# }



devtools::load_all("../../dataSHEEP_project/dataSHEEP/")
dev_path = "../../SHEEPfold_project/SHEEPfold/__SHEEP__"
list_path = list.files(dev_path, pattern='*.R$', full.names=TRUE, recursive=TRUE)
for (path in list_path) {
    source(path, encoding='UTF-8')
}
assign_colors(refCOL="TRACC")
# dataSHEEP::load_fonts()

# showtext::showtext_auto()
# showtext::showtext_opts(dpi = 300)
# sysfonts::font_add_google("Lato", "Lato")

# library(showtext)
# font_add("Lato",
#          regular="./resources/fonts/Lato/Lato-Regular.ttf",
#          bold="./resources/fonts/Lato/Lato-Bold.ttf")
# font_add("Raleway",
#          regular="./resources/fonts/Raleway/Raleway-Regular.ttf",
#          bold="./resources/fonts/Raleway/Raleway-Bold.ttf")
# showtext::showtext_auto()




add_path = function (x) {
    x = c(x, file.path(resources_path, logo_dir, x["file"]))
    names(x)[length(x)] = "path"
    return (x)
}
logo_info = lapply(logo_info, add_path)

icons_dirpath = file.path(resources_path, icons_dir)
icons_paths = list.files(icons_dirpath,
                         pattern="[.]svg",
                         recursive=TRUE, full.names=TRUE)
icons = lapply(icons_paths, svgparser::read_svg)
names(icons) = gsub(".svg", "", basename(icons_paths))

Stations_path = file.path(archive_data_path,
                          stations_selection_file)
Stations = ASHE::read_tibble(Stations_path)
Stations = filter(Stations, n_rcp85 >=4)


Secteurs_path = file.path(archive_data_path,
                          secteurs_selection_file)
Secteurs = ASHE::read_tibble(Secteurs_path)


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
        toleranceRel=toleranceRel_mini)
}


WL = list(
    "GWL-30"=c(GWL=3,
               RWL=4,
               GWLfull="GWL-3.0",
               RWLfull="RWL-4.0",
               GWLclean="GWL-30",
               RWLclean="RWL-40",
               color="#AE1C27")
    # "GWL-20"=c(GWL=2,
               # RWL=2.7,
               # GWLfull="GWL-2.0",
               # RWLfull="RWL-2.7",
               # GWLclean="GWL-20",
               # RWLclean="RWL-27",
               # color="#F47216")
)

NarraTRACC = list(
    "A"=c(name="Argousier",
          name_short="A",
          description="Débits réduits et étiages sévères",
          climateChain="HadGEM2-ES|historical-rcp85|ALADIN63|ADAMONT",
          Chain="XXX",
          color="#E66912",
          color_light="#f7c39e"),
    
    "G"=c(name="Genévrier",
          name_short="G",
          description="Débits en légère hausse et crues plus intenses",
          climateChain="IPSL-CM5A-MR|historical-rcp85|HIRHAM5|ADAMONT",
          Chain="YYY",
          color="#0f063b",
          color_light="#765def"),
    
    "E"=c(name="Érable",
          name_short="E",
          description="Intensification des extrêmes",
          climateChain="MPI-ESM-LR|historical-rcp85|CCLM4-8-17|ADAMONT",
          Chain="ZZZ",
          color="#870000",
          color_light="#ff6969"),
    
    "C"=c(name="Cèdre",
          name_short="C",
          description="évolutions modérées",
          climateChain="NorESM1-M|historical-rcp85|REMO|ADAMONT",
          Chain="AAA",
          color="#016367",
          color_light="#5ef7fd")
)

sheet_projection_secteur(
    Stations,
    Secteurs,
    dataEX_serie,
    metaEX_serie,
    dataEX_criteria,
    metaEX_criteria,
    WL=WL,
    NarraTRACC=NarraTRACC,
    icons=icons,
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

# capitalize_first <- function(s) {
  # paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
# }

