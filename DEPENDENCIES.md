# Dépendances reproductibles

Le calcul final de simulation a été réalisé sous R 4.5.3, Windows 11 64 bits,
avec RStan 2.32.7 et Stan 2.32.2.

## Dépendances directes

| Paquet | Version vérifiée | Rôle |
|---|---:|---|
| rstan | 2.32.7 | compilation et échantillonnage Stan |
| posterior | 1.7.0 | extraction et diagnostics des tirages |
| bridgesampling | 1.2-1 | calculs BMA de l’étude empirique |
| ggplot2 | 4.0.3 | certaines sorties graphiques |

Le code utilise aussi les packages de base `stats`, `utils`, `graphics`,
`grDevices`, `parallel` et `tools`.

## Installation

Exécuter :

```powershell
Rscript scripts/00_setup.R
```

Le script privilégie RStan lorsqu’il est déjà disponible. Sinon, il installe
CmdStanR et CmdStan. Les résultats archivés ont été produits avec RStan.

## Pourquoi il n’y a pas de renv.lock

Aucun fichier `renv.lock` n’a été archivé au moment du calcul final. Le
reconstituer après coup donnerait une fausse impression d’exactitude pour les
dépendances transitives. Ce fichier, `sessionInfo.txt` et les journaux
d’installation constituent donc la solution équivalente vérifiable. Une
nouvelle archive pourra créer un verrou `renv` avant tout recalcul futur.
