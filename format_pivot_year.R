
# pivot_year_TRACC[["start_GWL-15"]] = pivot_year_TRACC[["GWL-15"]] - lubridate::years(10)
# pivot_year_TRACC[["end_GWL-15"]] = pivot_year_TRACC[["GWL-15"]] + lubridate::years(9)
# pivot_year_TRACC[["start_GWL-20"]] = pivot_year_TRACC[["GWL-20"]] - lubridate::years(10)
# pivot_year_TRACC[["end_GWL-20"]] = pivot_year_TRACC[["GWL-20"]] + lubridate::years(9)
# pivot_year_TRACC[["start_GWL-30"]] = pivot_year_TRACC[["GWL-30"]] - lubridate::years(10)
# pivot_year_TRACC[["end_GWL-30"]] = pivot_year_TRACC[["GWL-30"]] + lubridate::years(9)
# ASHE::write_tibble(pivot_year_TRACC, pivot_year_TRACC_path)


# pivot_year_TRACC = tidyr::unite(pivot_year_TRACC, "climateChain", "EXP", "GCM", "RCM", "BC", sep="_", remove=FALSE)    
# pivot_year_TRACC = tidyr::unite(pivot_year_TRACC, "regexp", "EXP", "GCM", "RCM", "BC", sep=".*", remove=FALSE)
# pivot_year_TRACC = tidyr::unite(pivot_year_TRACC, "regexp_simulation_DRIAS", "GCM", "EXP", "RCM", "BC", sep=".*", remove=FALSE)
# pivot_year_TRACC = tidyr::unite(pivot_year_TRACC, "regexp_indicateur_DRIAS", "BC", "EXP", "GCM", "RCM", sep=".*", remove=FALSE)

# pivot_year_TRACC$regexp = paste0(".*", pivot_year_TRACC$regexp, ".*")
# pivot_year_TRACC$regexp_simulation_DRIAS = paste0(".*", pivot_year_TRACC$regexp_simulation_DRIAS, ".*")
# pivot_year_TRACC$regexp_indicateur_DRIAS = paste0(".*", pivot_year_TRACC$regexp_indicateur_DRIAS, ".*")

# pivot_year_TRACC$regexp = gsub("[-]", "[-]", pivot_year_TRACC$regexp)
# pivot_year_TRACC$regexp_simulation_DRIAS = gsub("[-]", "[-]", pivot_year_TRACC$regexp_simulation_DRIAS)
# pivot_year_TRACC$regexp_indicateur_DRIAS = gsub("[-]", "[-]", pivot_year_TRACC$regexp_indicateur_DRIAS)

