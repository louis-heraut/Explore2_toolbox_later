
library(dplyr)

Files = list.files(pattern=".csv")


blaise = c(
    C1="#603f8b",
    C2="#a16ae8",
    E1="#993404",
    E2="#d95f0e",
    E3="#fe9929",
    M1="#1f992d",
    M2="#94c973",
    X1="#252525",
    X2="#575656",
    X3="#7a7979" 
)
lou = c(
    C1="#0e464e", # bleu foncé
    C2="#2977A8", # bleu clair
    E1="#D04435", # rouge
    E2="#D88A1C", # orange
    E3="#FABF41", # jaune
    M1="#1C3C2F", # vert foncé
    M2="#52AE89", # vert clair
    X1="#57445d", # violet foncé
    X2="#715979", # violet moyen
    X3="#C1B2C7"  # violet clair
)

blaise = toupper(blaise)
lou = toupper(lou)
blaise = blaise[order(names(blaise))]
lou = lou[order(names(lou))]


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


