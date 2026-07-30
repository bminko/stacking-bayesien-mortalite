# Données de mortalité

L’étude de simulation est autonome. Les fichiers décrits ici sont requis
uniquement pour reproduire l’étude empirique belge.

## Source

Human Mortality Database. Max Planck Institute for Demographic Research
(Allemagne), University of California, Berkeley (États-Unis) et Institut
national d’études démographiques (France), données belges consultées le
13 juillet 2026 : <https://www.mortality.org/>.

Le projet utilise les fichiers annuels par âge simple :

- `death.txt` : décès ;
- `exposure.txt` : expositions au risque.

Les en-têtes des fichiers archivés indiquaient « Belgium, Deaths (period 1x1) »
et « Belgium, Exposure to risk (period 1x1) », dernière modification le
21 octobre 2025.

## Installation locale

Après authentification et téléchargement auprès de la HMD, créer le dossier
`data/raw` puis y placer :

```text
data/raw/death.txt
data/raw/exposure.txt
```

Le pipeline sélectionne la catégorie `Total`, femmes et hommes réunis, les âges
50 à 90 ans et les années 1970 à 2024. Il vérifie les années, les âges, les
valeurs manquantes, les expositions et les décès fractionnaires avant tout
ajustement.

Pour utiliser un autre emplacement, copier `config/paths_example.R` vers
`config/paths_local.R`, puis modifier uniquement cette copie locale.

## Redistribution

Les données brutes ne sont pas incluses dans ce dépôt. Chaque utilisateur doit
les obtenir directement auprès de la Human Mortality Database et respecter ses
conditions d’utilisation.
