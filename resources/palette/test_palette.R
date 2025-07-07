



meta = ASHE::read_tibble("families.csv", sep=";")

sous_famille = unique(meta[["sous-famille_description"]])
names(sous_famille) = unique(meta[["sous-famille"]])


Narratif_blaise = unique(meta$narratif_couleur)
names(Narratif_blaise) = unique(meta$narratif_id)
Narratif_blaise = Narratif_blaise[order(names(Narratif_blaise))]


# > sous_famille
# E : "Intensification des étiages" 
# X : "Intensification des évènements extrêmes" 
# M : "Changements modérés" 
# C : "Intensification des crues" 

Narratif_blaise = c(
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


Palette = c(
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

dataSHEEP::test_palette(Palette, colorStep=length(Palette))

Narratif = c(
    C1="#194765", # "#0e464e", # bleu foncé
    C2="#2977A8", # bleu clair
    E1="#D04435", # rouge
    E2="#D88A1C", # orange
    E3="#FABF41", # jaune
    M1="#37765c", #"#1C3C2F", # vert foncé
    M2="#52AE89", # vert clair
    X1="#57445d", # violet foncé
    X2="#977c9f", # violet moyen
    X3="#C1B2C7" # violet clair
)

dataSHEEP::test_palette(Narratif, colorStep=length(Narratif),
                        outname="Narratif")
