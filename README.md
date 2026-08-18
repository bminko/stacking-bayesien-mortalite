# Contextual Bayesian mortality stacking

Code minimal associé au mémoire *Agrégation bayésienne contextuelle de
modèles de mortalité : une approche par stacking hiérarchique selon l'âge et
l'horizon* (UCLouvain, 2026).

Le dépôt contient uniquement le code nécessaire pour reproduire l'analyse :
cinq modèles de mortalité bayésiens, la validation Leave-Future-Out (LFO),
les méthodes d'agrégation, l'évaluation hors échantillon, les quantités
actuarielles et l'étude de simulation calibrée sur des données réelles. Les
figures du mémoire, rapports intermédiaires, tirages MCMC et bibliothèques R
locales ne sont pas versionnés.

## Méthode reproduite

- Belgique, population `Total`, âges 50 à 90 ans ;
- modèles LC, RH, APC, CBD et M6 sous vraisemblance binomiale négative NB2 ;
- ajustement final sur 1970–2015 et test hors échantillon sur 2016–2024 ;
- validation LFO aux origines 2000–2014 et horizons 1 à 10 ;
- stacking global, stacking contextuel non régularisé et stacking contextuel
  hiérarchique entièrement estimé sous Stan ;
- calculs actuariels à partir de la moyenne pondérée des forces de mortalité,
  avec le même indice de tirage pour toute la trajectoire ;
- simulation HMD-calibrée : quatre scénarios, 30 répétitions par scénario,
  sans approximation de Laplace.

Le critère de validité Stan de l'étude de simulation reste fixé à
`R-hat <= 1.05`, ESS bulk et tail `>= 400`, aucune divergence, aucune
atteinte de la profondeur maximale, aucune chaîne bloquée et E-BFMI
`>= 0.30`. Dans les résultats définitifs du mémoire, 113 répétitions sur 120
satisfont ces critères ; les minima parmi les répétitions valides sont
409,29 pour l'ESS bulk et 401,86 pour l'ESS tail.

## Données

Les données HMD ne sont pas redistribuées. Après avoir obtenu l'autorisation
d'accès à la [Human Mortality Database](https://www.mortality.org/), placer à
la racine du dépôt :

```text
death.txt
exposure.txt
```

Il s'agit des tables belges `Deaths_1x1` et `Exposures_1x1`, conservées dans
leur format HMD original. Les colonnes attendues sont `Year`, `Age`, `Female`,
`Male` et `Total`. Voir [data/README.md](data/README.md).

## Installation

Prérequis : R >= 4.3, une chaîne de compilation C++ et soit RStan, soit
CmdStanR. Depuis la racine du dépôt :

```sh
Rscript scripts/00_setup.R
```

Le script installe dans `R/library/<version>` les dépendances manquantes et,
si nécessaire, CmdStan.

## Reproduction

Vérification courte de l'installation :

```sh
Rscript scripts/00_run_pipeline.R --profile=smoke
Rscript scripts/07_simulation_study.R --stage=pilot --profile=smoke
```

Analyse empirique complète :

```sh
Rscript scripts/00_run_pipeline.R --profile=full
```

Étude de simulation définitive :

```sh
Rscript scripts/07_simulation_study.R --stage=reference --profile=full
Rscript scripts/07_simulation_study.R --stage=main --profile=full
Rscript scripts/07_simulation_study.R --stage=finalize --profile=full
```

Les ajustements terminés sont mis en cache. Une commande relancée reprend
donc les calculs sans réestimer les répétitions déjà sauvegardées. Définir
`MEMOIRE_FORCE=1` uniquement pour invalider volontairement les caches.

## Arborescence

```text
config/   Paramètres, périodes, graines et seuils MCMC
R/        Fonctions de données, modèles, prévisions et agrégation
scripts/  Six étapes empiriques et étude de simulation
stan/     Cinq modèles de mortalité et méta-modèle hiérarchique
tests/    Tests numériques sans données HMD
data/     Instructions d'accès aux données
```

Les sorties sont créées sous `data/processed/<profil>/` et
`results/<profil>/`. Elles sont volontairement exclues de Git.

## Tests

Après avoir placé les deux fichiers HMD à la racine, exécuter :

```sh
Rscript tests/run_tests.R
Rscript tests/test_simulation_protocol.R
```

Ces tests contrôlent notamment le format des données, les transformations de
poids, l'agrégation actuarielle et le DGP HMD-calibré. Le profil `smoke`
vérifie ensuite l'orchestration complète avec un coût réduit ; il ne remplace
pas l'analyse définitive.
