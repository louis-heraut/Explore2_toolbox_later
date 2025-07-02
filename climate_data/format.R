
library(dplyr)

Files = list.files(pattern="NivRechauff")


Variables = c("RR", "TMm")
Variables_nc = c("prtotAdjust", "tasAdjust")

RWL_eric = c("2.7", "4")
GWL = c("GWL-20", "GWL-30")


data_list = list()

for (file in Files) {
    data_tmp = ASHE::read_tibble(file, sep=";")
    names(data_tmp)[1:3] = c("ncfile", "SH", "delta")
    
    data_tmp = select(data_tmp, ncfile, SH, delta)
    
    data_tmp = tidyr::separate(data_tmp,
                           ncfile,
                           c("variable", "xx", "EXP",
                             "GCM", "RCM", "BC", "yy", "zz", "season"),
                           sep="_")

    data_tmp = select(data_tmp, -xx, -yy, -zz)

    data_tmp$season = gsub(".nc", "", data_tmp$season)
    data_tmp$EXP = paste0("historical-", data_tmp$EXP)

    data_tmp$variable = Variables[match(data_tmp$variable, Variables_nc)]
    data_tmp$variable = paste0(data_tmp$variable, "_", data_tmp$season)
    data_tmp = select(data_tmp, -season)
    
    rwl = gsub("[_].*", "" ,gsub(".*NivRechauff", "", file))
    gwl = GWL[RWL_eric == rwl]
    
    data_tmp$GWL = gwl
    data_tmp = tidyr::pivot_wider(data_tmp,
                                  values_from=delta,
                                  names_from=variable,
                                  names_prefix="delta")

    data_tmp = relocate(data_tmp, GWL, .before=EXP)
    data_tmp = relocate(data_tmp, SH, .after=BC)

    data_tmp$GCM = gsub("CNRM-CM5-LR", "CNRM-CM5", data_tmp$GCM)
    data_tmp$RCM = gsub("HadREM3-GA7-05", "HadREM3-GA7", data_tmp$RCM)

    data_list = append(data_list, list(data_tmp))
}

data_RR = bind_rows(data_list[[1]], data_list[[2]])
data_TMm = bind_rows(data_list[[3]], data_list[[4]])

data = full_join(data_RR, data_TMm, by=c("GWL", "EXP", "GCM",
                                         "RCM", "BC", "SH"))
ASHE::write_tibble(data, "deltaRR_deltaTMm_season_GWL.csv")

