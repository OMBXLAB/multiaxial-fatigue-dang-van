# Évaluation de la fatigue multiaxiale — critère de Dang Van

Ce projet est un démonstrateur Python consacré à l’évaluation d’un historique
de contraintes multiaxiales périodiques avec le critère de fatigue à grand
nombre de cycles de Dang Van. Il utilise exclusivement des données synthétiques
et ne contient aucune information industrielle ou confidentielle.

## Objectif

Pour un historique discrétisé du tenseur des contraintes, le programme :

1. reconstruit le tenseur symétrique à partir des six composantes de Voigt ;
2. calcule la contrainte hydrostatique $p(t)=\mathrm{tr}(\boldsymbol\sigma)/3$ ;
3. discrétise les normales aux plans et les directions de cisaillement ;
4. détermine le cisaillement alterné résolu pour chaque orientation ;
5. recherche l’indicateur de Dang Van maximal et l’orientation critique ;
6. calcule un coefficient de sécurité en fatigue.

Le critère mis en œuvre s’écrit :

$$\max_t\left[\tau_a(t)+a\,p(t)\right]\leq b$$

Les paramètres sont calibrés à partir des limites d’endurance en flexion
alternée et en torsion alternée :

$$a=3\frac{\tau_{-1}}{\sigma_{-1}}-\frac{3}{2},\qquad b=\tau_{-1}$$

## Organisation du projet


```text
.
├── src/dang_van.py       # critère, orientations et génération des données
├── tests/test_dang_van.py
├── exemple.py             # démonstration et création du graphique
├── requirements.txt
└── LICENSE
```

## Installation et exécution

```bash
python -m venv .venv
pip install -r requirements.txt
python exemple.py
python -m unittest discover -s tests -v
```

La démonstration applique un chargement traction–torsion déphasé, représentant
un trajet de chargement non proportionnel simple. Elle affiche les paramètres
du critère, l’indicateur maximal et le coefficient de sécurité, puis génère le
fichier `resultats/chargement_synthetique.png`.

## Convention d’entrée

L’historique des contraintes est un tableau NumPy de dimensions `N × 6`,
exprimé en MPa et organisé de la manière suivante :

```text
[sigma_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_xz]
```

Exemple minimal :

```python
from dang_van import MateriauDangVan, evaluer_dang_van

materiau = MateriauDangVan(
    limite_flexion_alternee=220.0,
    limite_torsion_alternee=130.0,
)

resultat = evaluer_dang_van(historique_contraintes, materiau)
print(resultat.coefficient_securite)
```

## Hypothèses et limites

- les contraintes sont exprimées en MPa et décrivent un cycle périodique stabilisé ;
- les limites d’endurance doivent correspondre au matériau et au domaine de durée de vie étudiés ;
- la recherche des orientations est discrétisée et nécessite une étude de convergence ;
- le cisaillement alterné résolu est évalué par retrait de la valeur médiane ;
- l’état de surface, l’effet d’échelle, la température, les défauts et la dispersion statistique nécessitent un traitement spécifique ;
- ce démonstrateur a une vocation pédagogique et ne permet pas de certifier un composant critique.

## Auteur

**Oumar Mbengue** — Génie mécanique, calcul par éléments finis et programmation scientifique.

## English summary

This educational Python project evaluates a periodic multiaxial stress history
with the Dang Van high-cycle fatigue criterion. It reconstructs stress tensors,
calculates hydrostatic stress, searches discretised material-plane orientations
and returns the maximum criterion value and a fatigue safety factor.

The repository uses synthetic data only. It is intended for learning and
portfolio review, not for certifying a safety-critical component.

Released under the MIT License.
