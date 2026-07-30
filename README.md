# Stacking bayésien contextuel de modèles de mortalité

Ce dépôt accompagne le mémoire consacré à l’agrégation bayésienne contextuelle
de modèles de mortalité. Il contient le code de l’étude empirique et de l’étude
de simulation, les programmes Stan, les diagnostics synthétiques, les figures
finales et les rapports techniques.

La contribution principale est un stacking contextuel hiérarchique dont les
poids dépendent de l’âge et de l’horizon. Tous les ajustements hiérarchiques
rapportés ici sont réalisés sous Stan. Aucune approximation de Laplace n’est
utilisée dans les résultats finaux.

## Étude de simulation

La grille comprend les âges 50 à 90 et les horizons 1 à 10, soit 410 cellules
par répétition. Les noms LC, RH, APC, CBD et M6 désignent dans cette étude cinq
experts synthétiques, et non cinq modèles complets de mortalité réestimés à
chaque répétition.

Les quatre scénarios sont :

- `constant_weights` : poids vrais constants et nettement différenciés ;
- `horizon_only` : poids vrais variant uniquement avec l’horizon ;
- `age_horizon` : surface de poids variant avec l’âge et l’horizon ;
- `low_information` : poids constants et diffus, expositions plus faibles,
  surdispersion plus forte et experts plus proches.

L’étude principale comporte 30 répétitions par scénario, soit 120 répétitions.
Les 120 ajustements finaux sont valides. Dix ont nécessité une relance ciblée ;
aucun échec définitif, aucune divergence et aucune atteinte finale de la
profondeur maximale n’ont été observés.

## Environnement vérifié

Les simulations finales ont été exécutées sous Windows 11 64 bits avec :

- R 4.5.3 ;
- `rstan` 2.32.7 ;
- Stan 2.32.2 ;
- `posterior` 1.7.0 ;
- `bridgesampling` 1.2-1.

L’environnement détaillé est consigné dans
[`sessionInfo.txt`](sessionInfo.txt) et
[`DEPENDENCIES.md`](DEPENDENCIES.md). Ces deux fichiers constituent
l’équivalent reproductible disponible du verrou de dépendances : aucun
`renv.lock` original n’a été archivé au moment des calculs, donc un verrou
artificiel n’a pas été reconstruit après coup.

Un compilateur C++ compatible avec R et Stan est nécessaire. Le projet accepte
RStan lorsqu’il est installé ; sinon, `scripts/00_setup.R` installe CmdStanR et
CmdStan.

## Installation

Depuis la racine du dépôt :

```powershell
Rscript scripts/00_setup.R
```

Le script crée une bibliothèque locale dans `R/library`, vérifie les
dépendances et prépare un moteur Stan. Ce dossier local est ignoré par Git.

Un contrôle rapide du code ne nécessitant pas les données HMD s’exécute avec :

```powershell
Rscript tests/test_basic_pipeline.R
```

## Données HMD

L’étude de simulation est autonome et ne requiert pas les données HMD.
L’étude empirique utilise les fichiers annuels belges de décès et d’exposition
par âge simple, catégorie `Total`, femmes et hommes réunis.

Les données brutes ne sont pas redistribuées. Après téléchargement auprès de
la Human Mortality Database, placer les fichiers ainsi :

```text
data/raw/death.txt
data/raw/exposure.txt
```

Les instructions détaillées et la référence officielle se trouvent dans
[`data/README.md`](data/README.md). Un emplacement différent peut être indiqué
avec une copie locale non suivie de
[`config/paths_example.R`](config/paths_example.R).

## Exécuter la simulation

Les commandes suivantes doivent être lancées depuis la racine du dépôt.

### Une répétition de test

Ce test vérifie l’orchestration avec une configuration volontairement courte.
Il ne sert pas à interpréter la convergence finale.

```powershell
Rscript scripts/07_simulation_study.R --stage=pilot --profile=smoke --scenario=constant_weights --repetition=1
```

### Une répétition avec la configuration principale

```powershell
Rscript scripts/07_simulation_study.R --stage=main --profile=full --scenario=age_horizon --repetition=1
```

### Un scénario complet

```powershell
Rscript scripts/07_simulation_study.R --stage=main --profile=full --scenario=age_horizon
```

### Étude complète

```powershell
Rscript scripts/07_simulation_study.R --stage=pilot --profile=full
Rscript scripts/07_simulation_study.R --stage=main --profile=full
Rscript scripts/07_simulation_study.R --stage=validation --profile=full
Rscript scripts/07_simulation_study.R --stage=finalize --profile=full
Rscript scripts/13_generate_simulation_technical_report.R
```

Le pilote et l’étude principale utilisent des cohortes de graines distinctes.
Les résultats du pilote ne sont jamais agrégés aux 120 répétitions
principales.

## Configuration MCMC et relances

La configuration principale de la simulation utilise quatre chaînes, 500
itérations de chauffe et 500 itérations conservées par chaîne,
`adapt_delta = 0.98`, `max_treedepth = 12` et une métrique dense.

Un ajustement qui ne satisfait pas les diagnostics est réestimé avec quatre
chaînes, 1 000 itérations de chauffe et 1 000 itérations conservées,
`adapt_delta = 0.99`. En cas d’atteinte de la profondeur maximale, une dernière
configuration fixe `max_treedepth = 15`. Seule la première configuration
valide est conservée.

Les dix contrôles de robustesse réestiment exactement les mêmes données
simulées avec la configuration complète. Ils comparent les surfaces de poids,
les coefficients, les paramètres de régularisation, la RMSE, le LogS, le CRPS,
les classements et les diagnostics MCMC.

## Temps de calcul

La somme des temps Stan de l’étude principale est de 11,28 heures. Cette mesure
n’inclut pas la préparation des données, les dix contrôles complets ni la
production des rapports. La machine utilisée est un Lenovo 82BH équipé d’un
Intel Core i5-1135G7 à 2,40 GHz, de quatre cœurs physiques, de huit processeurs
logiques et de 16 Go de mémoire vive. Les répétitions sont séquentielles ;
les quatre chaînes d’un ajustement utilisent les quatre cœurs physiques.

## Étude empirique

Le pipeline empirique suit l’ordre :

```text
01_data_preprocessing.R
02_fit_individual_models.R
03_lfo_validation_loop.R
04_fit_aggregation_methods.R
05_evaluate_test_metrics.R
06_actuarial_quantities.R
08_prepare_chapter3_results.R
```

Le script `scripts/00_run_pipeline.R` orchestre ces étapes. Le profil `smoke`
sert uniquement à vérifier l’installation et les chemins ; les résultats du
mémoire proviennent du profil `full`.

## Organisation des sorties

```text
results/summaries/     résultats agrégés et performances
results/diagnostics/   diagnostics des 120 répétitions et des relances
results/validation/    sélection et comparaison des dix contrôles
figures/main/          figures vectorielles principales
figures/validation/    figures des contrôles numériques
docs/                  rapports PDF
```

Les grands objets RDS, les chaînes Stan, les caches de compilation et le fichier
`weights_by_repetition.csv` d’environ 80 Mo ne sont pas distribués. Ils sont
recréés par le pipeline.

## Graines

La graine maîtresse est `26052026`. Les graines dérivées dépendent de la
version du protocole, de la cohorte, du scénario, de la répétition et de
l’étape. La règle exacte est implémentée dans `stable_seed()` et
`simulation_repetition_seed()`. Les éléments de construction sont documentés
dans [`config/seeds.csv`](config/seeds.csv).

## Licence et citation

Le code est distribué sous licence MIT. Les données HMD restent soumises aux
conditions de la Human Mortality Database et ne sont pas couvertes par cette
licence.

Citation recommandée :

> Minko, B. A. J. (2026). *Agrégation bayésienne contextuelle de modèles de
> mortalité : code et résultats reproductibles*. Dépôt GitHub, version 1.0.0.

Les mêmes métadonnées sont fournies dans [`CITATION.cff`](CITATION.cff), afin
que GitHub et les gestionnaires bibliographiques puissent proposer directement
la référence du logiciel. L’URL définitive du dépôt doit être ajoutée ici et
dans le mémoire après sa création ; aucune URL n’a été inventée.
