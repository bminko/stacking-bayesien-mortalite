# Données HMD requises

Ce dépôt ne redistribue pas les données de la Human Mortality Database.
Chaque utilisateur doit les télécharger sous son propre compte et respecter
les conditions d'utilisation de la HMD.

Pour la Belgique, enregistrer à la racine du dépôt :

| Fichier local | Table HMD |
|---|---|
| `death.txt` | `Deaths_1x1` |
| `exposure.txt` | `Exposures_1x1` |

Les fichiers doivent rester dans le format texte HMD original : deux lignes
d'en-tête descriptives, puis les colonnes `Year`, `Age`, `Female`, `Male` et
`Total`. Le pipeline utilise la population `Total`, les âges 50–90 et les
années 1970–2024. Il arrête le traitement en présence de cellules manquantes
ou dupliquées, d'expositions non positives ou de décès fractionnaires ; il
n'effectue ni imputation ni arrondi silencieux.

Les observations 1970–2015 servent à l'ajustement final et au calibrage de
la simulation. Les années 2016–2024 restent hors apprentissage et servent
uniquement au test empirique depuis l'origine fixe 2015.
