
library(dplyr)

Files = list.files(pattern="^GWL.*[.]csv", full.names=TRUE)


GCM_antoine = c("CNRM",
                "EARTH",
                "Had2",
                "IPSL",
                "MPI",
                "Nor")
GCM = c(
    "CNRM-CM5",
    "EC-EARTH",
    "HadGEM2-ES",
    "IPSL-CM5A-MR",
    "MPI-ESM-LR",  
    "NorESM1-M"  
)

RCM_antoine = c(
    "ALAD",
    "Had3",
    "RAC",
    "RCA4",
    "CCLM4",
    "Reg6",
    "HIR5",
    "REMO",
    "WRF"
)
RCM = c(
    "ALADIN63",
    "HadREM3-GA7",
    "RACMO22E",
    "RCA4",
    "CCLM4-8-17",
    "RegCM4-6",
    "HIRHAM5",
    "REMO",
    "WRF381P"
)    


data = tibble()

for (file in Files) {
    data_tmp = ASHE::read_tibble(file, sep=";")
    data_tmp = select(data_tmp, -Med_models)

    data_tmp = rename(data_tmp, SH=CdSecteurH)
    data_tmp = rename(data_tmp, code_MESO=CdEuMasseD)

    cols = names(data_tmp)[3:ncol(data_tmp)]
    data_tmp = filter(data_tmp, !if_all(all_of(cols), is.na))
    data_tmp = tidyr::pivot_longer(data_tmp, cols=cols,
                               names_to="Chain",
                               values_to="deltaRecharge")

    data_tmp = tidyr::separate(data_tmp, col=Chain,
                           into=c("GCM", "RCM"))

    gwl = gsub("[_].*", "", basename(file))
    data_tmp$GWL = gwl
    data_tmp$EXP = "historical-rcp85"
    data_tmp$BC = "ADAMONT"
    data_tmp$HM = "RECHARGE"
    
    data_tmp$GCM = GCM[match(data_tmp$GCM, GCM_antoine)]
    data_tmp$RCM = RCM[match(data_tmp$RCM, RCM_antoine)]

    data_tmp = relocate(data_tmp, code_MESO, .before=deltaRecharge)
    data_tmp = relocate(data_tmp, GWL, .before=GCM)
    data_tmp = relocate(data_tmp, EXP, .before=GCM)
    data_tmp = relocate(data_tmp, BC, .before=code_MESO)
    data_tmp = relocate(data_tmp, HM, .before=code_MESO)
    data_tmp = relocate(data_tmp, SH, .before=code_MESO)
    
    data = bind_rows(data, data_tmp)
}

outfile =
    "deltaRecharge_yr_GWL-all_historical-rcp85_all_all_ADAMONT_all_on-MESO.csv"

ASHE::write_tibble(data,
                   outfile)

