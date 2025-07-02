
library(dplyr)

Files = list.files(pattern=".csv")


blaise = c(
    "S1"="#E66912",
    "H1"="#016367",
    "C1"="#16085C",
    "E2"="#D15252",
    "C2"="#3114C4",
    "S2"="#D9905D",
    "A1"="#3D038F",
    "S3"="#E0AD89",
    "E1"="#870000"
)
lou = c(
    "S1"="#D04435",
    "S2"="#D88A1C",
    "S3"="#FABF41",
    "C1"="#072327",
    "H1"="#1C3C2F",
    "C2"="#2977A8",
    "E1"="#653b2a",
    "E2"="#b76f52",
    "A1"="#5F4B66"
)

blaise = toupper(blaise)
lou = toupper(lou)
blaise = blaise[sort(names(blaise))]
lou = lou[sort(names(lou))]




for (file in Files) {
    meta = ASHE::read_tibble(file)
    meta_name = names(meta)
    meta_name = gsub("1", "_1", meta_name)
    meta_name = gsub("2", "_2", meta_name)
    meta_name = gsub("3", "_3", meta_name)
    meta_name = gsub("4", "_4", meta_name)
    names(meta) = meta_name

    meta = mutate(meta,
                  across(.cols=where(is.character),
                         .fns=~gsub("[|]", "_", .x)))
    meta = mutate(meta,
                  across(.cols=contains("Chain"),
                         .fns=~paste0("historical-rcp85_", .x)))
    meta = mutate(meta,
                  across(.cols=contains("color"),
                         .fns=~lou[match(toupper(.x), blaise)]))
    
    ASHE::write_tibble(meta, file)
}


