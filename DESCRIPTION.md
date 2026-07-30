# Description technique

## Objet

Le projet compare cinq experts de mortalité et trois formes de stacking :
global, contextuel non régularisé et contextuel hiérarchique. L’étude empirique
utilise les données belges de la Human Mortality Database. L’étude de simulation
évalue la récupération de surfaces de poids connues et la performance
prédictive dans quatre scénarios.

## Composants

- `config/` centralise les profils, les chemins et les graines.
- `R/` contient les fonctions de données, d’estimation, de prévision, de
  scoring, de stacking, de simulation, de diagnostic et de rapport.
- `stan/` contient les cinq modèles de mortalité et le méta-modèle
  hiérarchique.
- `scripts/` contient les points d’entrée exécutables.
- `tests/` contient les contrôles déterministes du pipeline.
- `results/` contient uniquement les sorties synthétiques partageables.
- `figures/` contient les figures PDF vectorielles.
- `docs/` contient les rapports techniques.

## Décisions de publication

Les données HMD brutes, les chaînes Stan, les caches, la bibliothèque R locale,
les journaux temporaires et les grands objets intermédiaires sont exclus. Le
fichier de poids cellule par cellule des 120 répétitions est également exclu en
raison de sa taille ; les résumés nécessaires à la vérification sont conservés.

Les scripts utilisent uniquement des chemins relatifs. Les chemins particuliers
à une machine peuvent être définis dans `config/paths_local.R`, fichier ignoré
par Git.
