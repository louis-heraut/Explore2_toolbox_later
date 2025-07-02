



# 1: "Forte diminution généralisée des débits" nommé "S1" ROUGE
# 6: "Diminution généralisée des débits" nommé "S2" ORANGE
# 8: "Changements futurs relativement peu marqués" nommé "S3" JAUNE

# 3: "Très forte intensification des crues et diminution des débits d'étiage" nommé "C1" BLEU
# 5: "Intensification des crues et légère diminution des débits d'étiage" nommé "C2" BLEU CLAIR

# 9: "Très forte diminution des débits d'étiage et intensification des crues" nommé "E1" MARRON
# 4: "Diminution des débits d'étiage et légère intensification des crues" nommé "E2" MARRON CLAIR

# 2: "Augmentation généralisée des débits" nommé "H1" VERT

# 7: "Baisse drastique des débits d'étiage" nommé "A1" VIOLET




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
