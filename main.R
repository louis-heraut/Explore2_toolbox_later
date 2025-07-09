lib_path = "./"

verbose =
    TRUE
# FALSE

subverbose =
    # TRUE
    FALSE

to_do = c(
    # "compute_delta"
    # "concatenate_delta"
    # "filter_concatenate_delta"
    # "reshape_filter_concatenate_delta"
    "plot"
)


MPI =
    ""
    # "file"
    # "secteur"

path_to_load =
    "/home/lheraut/Documents/INRAE/projects/Explore2_project/Explore2_toolbox_later/results/2025_06_19"


GWL = c("GWL-15", "GWL-20", "GWL-30")
period_reference_TRACC = c("1991-01-01", "2020-12-31")

NarraTRACC_order = c("X1", "X2", "X3",
                     "E1", "E2", "E3",
                     "C1", "C2",
                     "M1", "M2")

logo_info = list(
    "Explore2"=c(file='LogoExplore2.png'),
    "TRACC"=c(file='la_france_s_adapte.png')
)

n_projections_by_code = 4

# If the hydrological network needs to be plot
river_length = 30000

river_selection_mini =
    # NULL
    c("la Durance", "la Marne", "la Vienne", "le Loir", "la Loire",
      "l'Oise", "la Seine", "le Lot", "l'Adour", "le Rhône",
      "la Moselle", "l'Aisne", "la Garonne", "le Tarn", "le Doubs",
      "la Dordogne", "la Charente", "le Cher", "la Saône", "l'Allier",
      "Fleuve la Loire", "la Meuse", "la Sarthe", "la Somme",
      "l'Isère", "la Vilaine", "l'Aude", "l'Yonne")
river_selection_mini = paste0("^", river_selection_mini, "$")
river_length_mini = 30000


# Tolerance of the simplification algorithm for shapefile in sf
toleranceRel = 1000
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

post = function (x, ...) {
    if (verbose) {
        if (MPI != "") {
            print(paste0(formatC(as.character(rank), width = 3, 
                                 flag = " "), "/", size - 1, " > ", x), ...)
        }
        else {
            print(x, ...)
        }
    }
}

if (MPI != "") {
    library(Rmpi)
    rank = mpi.comm.rank(comm=0)
    size = mpi.comm.size(comm=0)
    if (size > 1) {
        if (rank == 0) {
            Rrank_sample = sample(0:(size-1))
            for (root in 1:(size-1)) {
                Rmpi::mpi.send(as.integer(Rrank_sample[root+1]),
                               type=1, dest=root,
                               tag=1, comm=0)
            }
            Rrank = Rrank_sample[1]
        } else {
            Rrank = Rmpi::mpi.recv(as.integer(0),
                                   type=1,
                                   source=0,
                                   tag=1, comm=0)
        }
    } else {
        Rrank = 0
    }
    post(paste0("Random rank attributed : ", Rrank))
    
} else {
    rank = 0
    size = 1
    Rrank = 0
}


# if (!exists("font")) {
#     library(extrafont)
#     font_import(paths=file.path(resources_path, "fonts"),
#                 prompt=FALSE)
#     # loadfonts()
#     loadfonts(device="pdf")
#     font = TRUE
# }


if ("compute_delta" %in% to_do) {

    variable_to_compute = c(
        "^QA[_]",
        "^QMA[_]",
        "^QSA[_]",
        "^VCN10[_]",
        "^QJXA[_]"
    )
    
    pivot_year_TRACC_path = file.path(archive_data_path,
                                      archive_metadata_dir, 
                                      pivot_year_TRACC_file)
    pivot_year_TRACC = ASHE::read_tibble(pivot_year_TRACC_path)

    data_path = file.path(
        archive_data_path,
        "hydrological-projections_indicateurs",
        "hydrological-projections_yearly-variables_by-chain_fst")

    outdir = "delta"
    
    Paths = list.files(data_path, pattern=".fst",
                       full.names=TRUE, recursive=TRUE)

    Files = basename(Paths)
    is_SAFRAN = grepl("SAFRAN", Files)
    is_RCP85 = grepl("rcp85", Files)
    is_ADAMONT = grepl("ADAMONT", Files)

    variable_to_compute_pattern =
        paste0("(", paste0(variable_to_compute, collapse=")|("), ")")
    is_variable = grepl(variable_to_compute_pattern, Files)
    
    Paths = Paths[!is_SAFRAN &
                  is_RCP85 &
                  is_ADAMONT &
                  is_variable]
    
    nPaths = length(Paths)
    nGWL = length(GWL)

    # stop()

    modify_filename = function (x, variable, gwl) {
        info = unlist(strsplit(x, "_"))
        if (info[2] == "yr") {
            info[1] = variable
        } else {
            info[1:2] = unlist(strsplit(variable,
                                        "_"))
        }
        filename = c(info[1:2], gwl, info[3:length(info)])
        paste0(filename, collapse="_")
    }

    format_dataEX_gwl = function (dataEX_gwl, variable, gwl) {
        dataEX_gwl$GWL = gwl
        dataEX_gwl = relocate(dataEX_gwl,
                              code, .after=HM)
        dataEX_gwl = relocate(dataEX_gwl,
                              GCM, .after=EXP)
        dataEX_gwl = relocate(dataEX_gwl,
                              GWL, .before=EXP)
        return (dataEX_gwl)
    }

    get_deltaEX_gwl = function (deltaEX_gwl, variable, gwl) {
        deltaEX_gwl$delta =
            (deltaEX_gwl$futur - deltaEX_gwl$historical) /
            deltaEX_gwl$historical * 100            
        deltaEX_gwl = rename(deltaEX_gwl,
                             !!variable:=delta)
        deltaEX_gwl = select(deltaEX_gwl, -historical, -futur)
        deltaEX_gwl = format_dataEX_gwl(deltaEX_gwl, variable, gwl)
        return (deltaEX_gwl)
    }



    if (MPI == "file") {
        start = ceiling(seq(1,  nPaths,
                            by=(nPaths/size)))
        if (any(diff(start) == 0)) {
            start = 1:nPaths
            end = start
        } else {
            end = c(start[-1]-1, nPaths)
        }
        if (rank == 0) {
            post(paste0(paste0("rank ", 0:(size-1), " get ",
                                     end-start+1, " files"),
                              collapse="    "))
        }
        if (Rrank+1 > nPaths) {
            Paths = NULL
            Rmpi::mpi.send(as.integer(1), type=1,
                           dest=0, tag=1, comm=0)
            post(paste0("End signal from rank ", rank))
        } else {
            Paths = Paths[start[Rrank+1]:end[Rrank+1]]
        }
    } else {
        Paths = Paths
    } 
    nPaths = length(Paths)
    
    post(paste0("All ", nPaths, " paths: ",
                paste0(basename(Paths), collapse=" | ")))
    
    for (i in 1:nPaths) {
        path = Paths[i]
        post(paste0(i, "/", nPaths,
                          " paths -> ",
                          round(i/nPaths*100, 1), "%"))
        
        dataEX = ASHE::read_tibble(path)
        
        for (j in 1:nGWL) {
            gwl = GWL[j]
            
            climateChain_regexp =
                paste0(".*",
                       paste0(unlist(strsplit(basename(path),
                                              "_"))[3:6],
                              collapse=".*"), ".*")
            climateChain_regexp = gsub("historical-", "",
                                       climateChain_regexp)
            climateChain_regexp = gsub("[-]", "[-]",
                                       climateChain_regexp)

            pivot_year_TRACC_climateChain =
                dplyr::filter(pivot_year_TRACC,
                          grepl(climateChain_regexp, climateChain))

            if (nrow(pivot_year_TRACC_climateChain) == 0) {
                stop(paste0("no pivot year for ", path))
            }
            if (nrow(pivot_year_TRACC_climateChain) > 1) {
                stop(paste0("more than one pivot year for ", path))
            }

            start_gwl =
                pivot_year_TRACC_climateChain[[paste0("start_", gwl)]]
            end_gwl =
                pivot_year_TRACC_climateChain[[paste0("end_", gwl)]]
            period_futur_TRACC = c(start_gwl, end_gwl)
            
            variable = paste0(unlist(strsplit(basename(path),
                                              "_"))[1:2],
                              collapse="_")
            variable = gsub("_yr", "", variable)

            # stop()
            
            if (grepl("QMA[_]", variable)) {
                mean_variable = paste0("mean", variable)
                
                dataEX_gwl = summarise(group_by(
                    filter(dataEX,
                           period_futur_TRACC[1] <= date &
                           date <= period_futur_TRACC[2]),
                    code, GCM, EXP, RCM, BC, HM),
                    !!mean_variable:=mean(get(variable),
                                     na.rm=TRUE),
                    .groups="drop")

                dataEX_gwl = format_dataEX_gwl(dataEX_gwl, mean_variable, gwl)
                
                outfile =
                    modify_filename(basename(path),
                                    variable=mean_variable,
                                    gwl=gwl)
                outdirpath = gsub(".*[/]", "", dirname(path))
                outpath = file.path(today_resdir, outdir,
                                    outdirpath, outfile)
                ASHE::write_tibble(dataEX_gwl, outpath)
            }
            
            if (variable == "VCN10") {
                returnPeriod = 5
                waterType = "low"
            }
            if (variable == "QJXA") {
                returnPeriod = 10
                waterType = "high"
            }

            if (variable %in% c("VCN10", "QJXA")) {
                delta_rp_variable = paste0("delta", variable,
                                           "-", returnPeriod)
                deltaEX_gwl =
                    full_join(
                        summarise(group_by(
                            filter(dataEX,
                                   period_reference_TRACC[1] <= date &
                                   date <= period_reference_TRACC[2]),
                            code, GCM, EXP, RCM, BC, HM),
                            historical=CARD::get_Xn(get(variable),
                                                    returnPeriod=
                                                        returnPeriod,
                                                    waterType=
                                                        waterType),
                            .groups="drop"),
                        
                        summarise(group_by(
                            filter(dataEX,
                                   period_futur_TRACC[1] <= date &
                                   date <= period_futur_TRACC[2]),
                            code, GCM, EXP, RCM, BC, HM),
                            futur=CARD::get_Xn(get(variable),
                                               returnPeriod=
                                                   returnPeriod,
                                               waterType=
                                                   waterType),
                            .groups="drop"),
                        by=c("code", "GCM", "EXP", "RCM", "BC", "HM"))

                deltaEX_gwl = get_deltaEX_gwl(deltaEX_gwl,
                                              variable=delta_rp_variable,
                                              gwl=gwl)

                outfile =
                    modify_filename(basename(path),
                                    variable=
                                        delta_rp_variable,
                                    gwl=gwl)
                outdirpath = gsub(".*[/]", "", dirname(path))
                outpath = file.path(today_resdir, outdir,
                                    outdirpath, outfile)
                ASHE::write_tibble(deltaEX_gwl, outpath)
                
            }

            delta_variable = paste0("delta", variable)
            deltaEX_gwl =
                full_join(
                    summarise(group_by(
                        filter(dataEX,
                               period_reference_TRACC[1] <= date &
                               date <= period_reference_TRACC[2]),
                        code, GCM, EXP, RCM, BC, HM),
                        historical=mean(get(variable),
                                        na.rm=TRUE),
                        .groups="drop"),
                    
                    summarise(group_by(
                        filter(dataEX,
                               period_futur_TRACC[1] <= date &
                               date <= period_futur_TRACC[2]),
                        code, GCM, EXP, RCM, BC, HM),
                        futur=mean(get(variable),
                                   na.rm=TRUE),
                        .groups="drop"),
                    by=c("code", "GCM", "EXP", "RCM", "BC", "HM"))

            deltaEX_gwl = get_deltaEX_gwl(deltaEX_gwl,
                                          variable=delta_variable,
                                          gwl=gwl)

            outfile =
                modify_filename(basename(path),
                                variable=delta_variable,
                                gwl=gwl)
            outdirpath = gsub(".*[/]", "", dirname(path))
            outpath = file.path(today_resdir, outdir,
                                outdirpath, outfile)
            ASHE::write_tibble(deltaEX_gwl, outpath)

            
            # stop()
        }
    }
}


get_var = function (x) {
    paste0(unlist(strsplit(x, "_"))[1:2], collapse="_")
}


if ("concatenate_delta" %in% to_do) {
    outdir = "concatenated_delta"

    Paths = list.files(file.path(path_to_load, "delta"),
                       pattern=".fst",
                       recursive=TRUE, full.names=TRUE)

    Variables_ALL = sapply(basename(Paths), get_var, USE.NAMES=FALSE)
    Variables = unique(Variables_ALL)
    nVariables = length(Variables)

    if (MPI == "file") {
        start = ceiling(seq(1,  nVariables,
                            by=(nVariables/size)))
        if (any(diff(start) == 0)) {
            start = 1:nVariables
            end = start
        } else {
            end = c(start[-1]-1, nVariables)
        }
        if (rank == 0) {
            post(paste0(paste0("rank ", 0:(size-1), " get ",
                               end-start+1, " files"),
                        collapse="    "))
        }
        if (Rrank+1 > nVariables) {
            Variables = NULL
            Rmpi::mpi.send(as.integer(1), type=1,
                           dest=0, tag=1, comm=0)
            post(paste0("End signal from rank ", rank))
        } else {
            Variables = Variables[start[Rrank+1]:end[Rrank+1]]
        }
    }
    nVariables = length(Variables)

    post(paste0("All ", nVariables, " variables: ",
                paste0(Variables, collapse=" | ")))
    
    for (i in 1:nVariables) {
        variable = Variables[i]
        post(paste0("* ", i, "/", nVariables, " -> ",
                    round(i/nVariables, 1)*100, "%"))
        Paths_var = Paths[grepl(variable, basename(Paths))]
        nPaths_var = length(Paths_var)
        deltaEX = dplyr::tibble()

        for (j in 1:nPaths_var) {
            if (j %% 10 == 0) {
                post(paste0("** ", j, "/", nPaths_var, " -> ",
                            round(j/nPaths_var, 1)*100, "%"))
            }
            path_var = Paths_var[j]
            deltaEX_tmp = ASHE::read_tibble(path_var)
            deltaEX = dplyr::bind_rows(deltaEX, deltaEX_tmp)
        }

        outfile =
            paste0(variable,
                   "_GWL-all_historical-rcp85_all_all_ADAMONT_all.fst")
        outpath = file.path(today_resdir,
                            outdir, outfile)
        ASHE::write_tibble(deltaEX, outpath)
    }
}


if ("filter_concatenate_delta" %in% to_do) {

    outdir = "filtered_concatenated_delta"

    Paths = list.files(file.path(path_to_load, "concatenated_delta"),
                       pattern=".fst",
                       recursive=TRUE, full.names=TRUE)
    
    Stations_path = file.path(archive_data_path,
                              archive_metadata_dir, 
                              stations_selection_file)
    Stations = ASHE::read_tibble(Stations_path)
    Stations = filter(Stations, n_rcp85 >=4)
    Code_selection = dplyr::filter(Stations, n_rcp85 >=4)$code
    
    code_Chain_outliers_path = file.path(archive_data_path,
                                         archive_metadata_dir, 
                                         code_Chain_outliers_file)
    code_Chain_outliers = ASHE::read_tibble(code_Chain_outliers_path)

    
    for (path in Paths) {
        deltaEX = ASHE::read_tibble(path)
        
        deltaEX = dplyr::filter(deltaEX, code %in% Code_selection)

        deltaEX = tidyr::unite(deltaEX,
                               "code_Chain",
                               "code",
                               "EXP", "GCM", ,
                               "RCM", "BC",
                               "HM", sep="_",
                               remove=FALSE)
        deltaEX = dplyr::filter(deltaEX,
                                !(code_Chain %in%
                                  code_Chain_outliers$code_Chain))
        deltaEX = dplyr::select(deltaEX, -code_Chain)
        
        outpath = file.path(today_resdir,
                            outdir,
                            basename(path))
        outpath = gsub("_all.fst", "_all_filtered.fst", outpath)
        ASHE::write_tibble(deltaEX, outpath)
    }
}


Months = c("jan", "feb", "mar", "apr",
           "may", "jun", "jul", "aug",
           "sep", "oct", "nov", "dec")
GWL_year = c(2030, 2050, 2100)
# names(GWL_year) = GWL

if ("reshape_filter_concatenate_delta" %in% to_do) {

    outdir_criteria = "reshaped_filterd_concatenated_delta"
    outdir_serie = "reshaped-serie_filterd_concatenated_delta"

    Paths = list.files(file.path(path_to_load,
                                 "filtered_concatenated_delta"),
                       pattern=".fst",
                       recursive=TRUE, full.names=TRUE)


    Variables_ALL = sapply(basename(Paths), get_var, USE.NAMES=FALSE)
    Variables = unique(Variables_ALL)
    Variables = Variables[!grepl("QMA", Variables)]
    Variables = c("meanQMA_month", Variables, "deltaQMA_month")
    
    nVariables = length(Variables)

    if (MPI == "file") {
        start = ceiling(seq(1, nVariables,
                            by=(nVariables/size)))
        if (any(diff(start) == 0)) {
            start = 1:nVariables
            end = start
        } else {
            end = c(start[-1]-1, nVariables)
        }
        if (rank == 0) {
            post(paste0(paste0("rank ", 0:(size-1), " get ",
                               end-start+1, " files"),
                        collapse="    "))
        }
        if (Rrank+1 > nVariables) {
            Variables = NULL
            Rmpi::mpi.send(as.integer(1), type=1,
                           dest=0, tag=1, comm=0)
            post(paste0("End signal from rank ", rank))
        } else {
            Variables = Variables[start[Rrank+1]:end[Rrank+1]]
        }
    }
    nVariables = length(Variables)

    post(paste0("All ", nVariables, " variables: ",
                paste0(Variables, collapse=" | ")))
    
    for (i in 1:nVariables) {
        variable = Variables[i]
        post(paste0("* ", i, "/", nVariables, " -> ",
                    round(i/nVariables, 1)*100, "%"))

        if (variable == "deltaQMA_month") {
            Paths_var = Paths[grepl("deltaQMA", Paths)]
        }
        if (variable == "meanQMA_month") {
            Paths_var = Paths[grepl("meanQMA", Paths)]
        }
                
        if (variable == "deltaQMA_month" |


variable == "meanQMA_month") {
            deltaEX = dplyr::tibble()
            for (path_var in Paths_var) {                
                deltaEX_tmp = ASHE::read_tibble(path_var)
                
                var = get_var(basename(path_var))
                month = gsub(".*[_]", "", var)
                month_id = formatC(which(Months == month), width=2,
                                   flag="0")
                deltaEX_tmp$date = GWL_year[sapply(deltaEX_tmp$GWL,
                                                   match, table=GWL)]
                deltaEX_tmp$date = as.Date(paste0(deltaEX_tmp$date,
                                                  "-", month_id, "-01"))
                deltaEX_tmp = dplyr::relocate(deltaEX_tmp, date,
                                              .after=code)
                names(deltaEX_tmp)[grepl("QMA",
                                         names(deltaEX_tmp))] =
                    gsub("[_]month", "", variable)
                deltaEX = dplyr::bind_rows(deltaEX, deltaEX_tmp)
            }
            outdir = outdir_serie
        } else {
            path_var = Paths[grepl(variable, Paths)]
            deltaEX = ASHE::read_tibble(path_var)
            outdir = outdir_criteria
        }
        
        SH = unique(substr(deltaEX$code, 1, 2))
        nSH = length(SH)

        for (j in 1:nSH) {
            if (j %% 10 == 0) {
                post(paste0("** ", j, "/", nSH, " -> ",
                            round(j/nSH, 1)*100, "%"))
            }
            sh = SH[j]
            deltaEX_sh = dplyr::filter(deltaEX, substr(code, 1, 2) == sh)
            outfile = paste0(
                variable,
                "_GWL-all_historical-rcp85_all_all_ADAMONT_all_filtered_",
                sh, ".fst")
            outpath = file.path(today_resdir,
                                outdir, sh,
                                outfile)
            ASHE::write_tibble(deltaEX_sh, outpath)
        }
    }
}






if ("plot" %in% to_do) {
    library(ggplot2)
    library(latex2exp)
    
    library(ggtext)
    library(gridtext)
    library(Cairo)
    options(bitmapType="cairo")
    # Cairo::Cairo.options(antialias="subpixel")
    
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
                              archive_metadata_dir, 
                              stations_selection_file)
    Stations = ASHE::read_tibble(Stations_path)
    Stations = filter(Stations, n_rcp85 >=4)


    Secteurs_path = file.path(archive_data_path,
                              archive_metadata_dir,
                              secteurs_selection_file)
    Secteurs = ASHE::read_tibble(Secteurs_path)
    Secteurs$id = 1:nrow(Secteurs)

    Projections_path = file.path(archive_data_path,
                                 archive_metadata_dir,
                                 projections_selection_file)
    Projections = ASHE::read_tibble(Projections_path)
    Projections = filter(Projections,
                         EXP == "historical-rcp85" &
                         grepl("ADAMONT", BC))

    Variables_hydro_path = file.path(archive_data_path,
                                     archive_metadata_dir,
                                     variables_hydro_selection_file)
    Variables_hydro = ASHE::read_tibble(Variables_hydro_path)


    dataEX_criteria_climate_secteur_path =
        file.path(climate_data_dirpath,
                  climate_data_file)
    dataEX_criteria_climate_secteur =
        ASHE::read_tibble(dataEX_criteria_climate_secteur_path)

    dataEX_criteria_climate_secteur =
        tidyr::unite(dataEX_criteria_climate_secteur,
                     climateChain,
                     EXP, GCM,
                     RCM, BC,
                     sep="_",
                     remove=FALSE)

    dataEX_criteria_climate_secteur =
        dplyr::relocate(dataEX_criteria_climate_secteur,
                        climateChain, .after=BC)

    dataEX_criteria_recharge_MESO_path =
        list.files(file.path(recharge_data_dirpath,
                             recharge_criteria_MESO_dir),
                   pattern=".csv", full.names=TRUE)
    dataEX_criteria_recharge_MESO =
        ASHE::read_tibble(dataEX_criteria_recharge_MESO_path)
    
    dataEX_criteria_recharge_MESO =
        tidyr::unite(dataEX_criteria_recharge_MESO,
                     climateChain,
                     EXP, GCM,
                     RCM, BC,
                     sep="_",
                     remove=FALSE)

    
    dataEX_criteria_recharge_secteur_path =
        list.files(file.path(recharge_data_dirpath,
                             recharge_criteria_secteur_dir),
                   pattern=".csv", full.names=TRUE)
    dataEX_criteria_recharge_secteur =
        ASHE::read_tibble(dataEX_criteria_recharge_secteur_path)
    
    dataEX_criteria_recharge_secteur =
        tidyr::unite(dataEX_criteria_recharge_secteur,
                     climateChain,
                     EXP, GCM,
                     RCM, BC,
                     sep="_",
                     remove=FALSE)

    if (!exists("Shapefiles") | !exists("Shapefiles_mini")) {
        post("### Loading shapefiles")

        Shapefiles = load_shapefile(
            computer_shp_path, Code=NULL,
            europe_shp_path=europe_shp_path,
            france_shp_path=france_shp_path,
            bassinHydro_shp_path=bassinHydro_shp_path,
            regionHydro_shp_path=regionHydro_shp_path,
            secteurHydro_shp_path=secteurHydro_shp_path,
            MESO_shp_path=MESO_shp_path,
            river_shp_path=river_shp_path,
            river_selection=NULL,
            river_length=river_length,
            toleranceRel=toleranceRel)
        
        Shapefiles_mini = load_shapefile(
            computer_shp_path, Code=NULL,
            france_shp_path=france_shp_path,
            bassinHydro_shp_path=bassinHydro_shp_path,
            regionHydro_shp_path=regionHydro_shp_path,
            secteurHydro_shp_path=secteurHydro_shp_path,
            river_shp_path=river_shp_path,
            river_selection=river_selection_mini,
            river_length=river_length_mini,
            toleranceRel=toleranceRel_mini)
    }


    WL = list(
        # "GWL-30"=c(GWL=3,
                   # RWL=4,
                   # GWLfull="GWL-3.0",
                   # RWLfull="RWL-4.0",
                   # GWLclean="GWL-30",
                   # RWLclean="RWL-40",
                   # color="#AE1C27"),
        "GWL-20"=c(GWL=2,
                   RWL=2.7,
                   GWLfull="GWL-2.0",
                   RWLfull="RWL-2.7",
                   GWLclean="GWL-20",
                   RWLclean="RWL-27",
                   color="#F47216")
    )


    ####### /!\ climateChain pas le meme ordre que Chain
    # NarraTRACC = list(
    #     "A"=c(name="Argousier",
    #           name_short="A",
    #           description="Débits réduits et étiages sévères",
    #           climateChain="HadGEM2-ES_historical-rcp85_ALADIN63_ADAMONT",
    #           Chain="historical-rcp85_HadGEM2-ES_ALADIN63_ADAMONT_SMASH",
    #           color="#E66912",
    #           color_light="#f7c39e"),
        
    #     "G"=c(name="Genévrier",
    #           name_short="G",
    #           description="Débits en légère hausse et crues plus intenses",
    #           climateChain="IPSL-CM5A-MR_historical-rcp85_HIRHAM5_ADAMONT",
    #           Chain="historical-rcp85_IPSL-CM5A-MR_HIRHAM5_ADAMONT_SMASH",
    #           color="#0f063b",
    #           color_light="#765def"),
        
    #     "E"=c(name="Érable",
    #           name_short="E",
    #           description="Intensification des extrêmes",
    #           climateChain="MPI-ESM-LR_historical-rcp85_CCLM4-8-17_ADAMONT",
    #           Chain="historical-rcp85_MPI-ESM-LR_CCLM4-8-17_ADAMONT_SMASH",
    #           color="#870000",
    #           color_light="#ff6969"),
        
    #     "C"=c(name="Cèdre",
    #           name_short="C",
    #           description="Évolutions modérées",
    #           climateChain="NorESM1-M_historical-rcp85_REMO_ADAMONT",
    #           Chain="historical-rcp85_NorESM1-M_REMO_ADAMONT_SMASH",
    #           color="#016367",
    #           color_light="#5ef7fd")
    # )

    NarraTRACC_selection_Paths =
        list.files(file.path(archive_data_path,
                             archive_metadata_dir),
                   pattern="narraTRACC",
                   full.names=TRUE)
    NarraTRACC_selection = lapply(NarraTRACC_selection_Paths,
                                  ASHE::read_tibble)
    names(NarraTRACC_selection) =
        gsub("[_].*", "", basename(NarraTRACC_selection_Paths))
    

    SH = unique(substr(Stations$code, 1, 2))

    ###
    # SH = SH[grepl("(O)|(M)", SH)]
    # SH = c("K2", "M0", "Q0")
    SH = "W0"
    ###
    nSH = length(SH) 

    if (MPI == "secteur") {
        start = ceiling(seq(1, nSH,
                            by=(nSH/size)))
        if (any(diff(start) == 0)) {
            start = 1:nSH
            end = start
        } else {
            end = c(start[-1]-1, nSH)
        }
        if (rank == 0) {
            post(paste0(paste0("rank ", 0:(size-1), " get ",
                               end-start+1, " files"),
                        collapse="    "))
        }
        if (Rrank+1 > nSH) {
            SH = NULL
            Rmpi::mpi.send(as.integer(1), type=1,
                           dest=0, tag=1, comm=0)
            post(paste0("End signal from rank ", rank))
        } else {
            SH = SH[start[Rrank+1]:end[Rrank+1]]
        }
    }
    nSH = length(SH)

    post(paste0("All ", nSH, " secteurs: ",
                paste0(SH, collapse=" | ")))

    
    for (i in 1:nSH) {
        sh = SH[i]
        post(paste0(i, "/", nSH, " so ", round(i/nSH*100, 1),
                    "% done -> ", sh))

        secteur = Secteurs[Secteurs$id_secteur == sh,]
        Stations_sh = filter(Stations, substr(code_hydro2, 1, 2) == sh)

        dataEX_criteria_hydro_dirpath =
            file.path(hydro_data_dirpath, hydro_criteria_dir,sh)
        dataEX_criteria_hydro_paths =
            list.files(dataEX_criteria_hydro_dirpath,
                       full.names=TRUE)
        dataEX_criteria_hydro_list = lapply(dataEX_criteria_hydro_paths,
                                            ASHE::read_tibble)

        dataEX_criteria_hydro =
            purrr::reduce(dataEX_criteria_hydro_list, full_join,
                          by=c("code", "GWL", "EXP",
                               "GCM", "RCM", "BC", "HM"))
        dataEX_criteria_hydro$SH = substr(dataEX_criteria_hydro$code,
                                          1, 2)
        dataEX_criteria_hydro = dplyr::relocate(dataEX_criteria_hydro,
                                                SH, .before=code)

        dataEX_serie_hydro_dirpath = file.path(hydro_data_dirpath,
                                               hydro_serie_dir, sh)
        dataEX_serie_hydro_paths = list.files(dataEX_serie_hydro_dirpath,
                                              full.names=TRUE)
        dataEX_serie_hydro = lapply(dataEX_serie_hydro_paths,
                                    ASHE::read_tibble)
        names(dataEX_serie_hydro) =
            gsub("[_].*", "", basename(dataEX_serie_hydro_paths))

        for (j in 1:length(dataEX_serie_hydro)) {
            dataEX_serie_hydro[[j]]$SH =
                substr(dataEX_serie_hydro[[j]]$code, 1, 2)
            dataEX_serie_hydro[[j]] =
                dplyr::relocate(dataEX_serie_hydro[[j]],
                                SH, .before=code)
        }
        
        sheet_projection_secteur(
            Stations_sh,
            secteur,
            dataEX_criteria_climate_secteur,
            dataEX_criteria_hydro,
            dataEX_serie_hydro,
            dataEX_criteria_recharge_MESO,
            dataEX_criteria_recharge_secteur,
            # metaEX_criteria_climate,
            # metaEX_criteria_hydro,
            # metaEX_serie_hydro,
            WL=WL,
            NarraTRACC_selection=NarraTRACC_selection,
            NarraTRACC_order=NarraTRACC_order,
            icons=icons,
            logo_info=logo_info,
            Shapefiles=Shapefiles,
            Shapefiles_mini=Shapefiles_mini,
            figdir=figdir,
            Pages=NULL,
            is_MPI=MPI!="",
            verbose=subverbose)
    }

    # confiance au dessus de 80%


    ### /!\ ?
    # PAS de Z
    # RCM : REMO2015 et REMO2009 -> REMO2009
    # RCM : SMHI-RCA4 -> RCA4

    # capitalize_first <- function(s) {
    # paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
    # }
}



if (MPI != "") {
    Sys.sleep(5)
    mpi.finalize()
}







Stations = ASHE::read_tibble("Selection_points_simulation_20240219.csv")
Stations$XL93 = as.numeric(gsub(",", ".", Stations$XL93))
Stations$YL93 = as.numeric(gsub(",", ".", Stations$YL93))

SH = unique(substr(Stations$Code, 1, 2))
# SH = SH[grepl("H", SH)]
nSH = length(SH) 

for (sh in SH) {
    Stations_sh = filter(Stations, substr(Code, 1, 2) == sh)
    secteurHydro_shp =
        Shapefiles$secteurHydro[Shapefiles$secteurHydro$CdSecteurH ==
                                sh,]
    
    plot =
        ggplot() + theme_void() +
        ggtitle(sh) +
        geom_sf(data=Shapefiles$secteurHydro) +
        geom_sf(data=secteurHydro_shp, fill=IPCCgold) +
        geom_point(data=Stations_sh, aes(x=XL93, y=YL93), size=0.5)

    ggsave(plot=plot,
           filename=paste0(sh, ".pdf"),
           path=file.path(figdir, "SH"),
           height=10, width=10, units="cm",
           device=cairo_pdf)
    
}

secteurHydro_centroid=sf::st_centroid(Shapefiles$secteurHydro)

plot =
    ggplot() + theme_void() +
    ggtitle("localisation secteur") +
    geom_sf(data=Shapefiles$secteurHydro) +
    geom_sf_text(data=secteurHydro_centroid,
                 aes(label=CdSecteurH), size=3)

ggsave(plot=plot,
       filename="_SH_location_.pdf",
       path=file.path(figdir, "SH"),
       height=10, width=10, units="cm",
       device=cairo_pdf)


Paths = list.files(file.path(figdir, "SH"),
                   pattern=".pdf",
                   full.names=TRUE)
output_file = file.path(file.path(figdir, "SH"), "_SH_.pdf")
qpdf::pdf_combine(input=Paths, output=output_file)
