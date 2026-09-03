"""Implémentation pédagogique d'une vérification macroscopique de Dang Van.

Les composantes du tenseur sont ordonnées comme suit :
[sigma_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_xz], en MPa.
"""

from dataclasses import dataclass
from itertools import product

import numpy as np


@dataclass(frozen=True)
class MateriauDangVan:
    """Limites d'endurance utilisées dans la droite tau + a*p <= b."""

    limite_flexion_alternee: float
    limite_torsion_alternee: float

    def __post_init__(self) -> None:
        if self.limite_flexion_alternee <= 0 or self.limite_torsion_alternee <= 0:
            raise ValueError("Les limites d'endurance doivent être strictement positives.")

    @property
    def a(self) -> float:
        """Pente de la droite de Dang Van."""
        return 3.0 * self.limite_torsion_alternee / self.limite_flexion_alternee - 1.5

    @property
    def b(self) -> float:
        """Ordonnée à l'origine de la droite de Dang Van, en MPa."""
        return self.limite_torsion_alternee


@dataclass(frozen=True)
class ResultatDangVan:
    """Résultat de la recherche du plan critique."""

    contrainte_equivalente: float
    coefficient_securite: float
    normale_critique: np.ndarray
    direction_critique: np.ndarray
    indice_temps_critique: int


def voigt_vers_tenseur(historique_contraintes: np.ndarray) -> np.ndarray:
    """Convertit un historique (N, 6) en N tenseurs symétriques 3 x 3."""
    historique_contraintes = np.asarray(historique_contraintes, dtype=float)
    if historique_contraintes.ndim != 2 or historique_contraintes.shape[1] != 6:
        raise ValueError("L'historique des contraintes doit avoir la forme (N, 6).")

    tenseurs = np.zeros((len(historique_contraintes), 3, 3), dtype=float)
    tenseurs[:, 0, 0], tenseurs[:, 1, 1], tenseurs[:, 2, 2] = historique_contraintes[:, :3].T
    tenseurs[:, 0, 1] = tenseurs[:, 1, 0] = historique_contraintes[:, 3]
    tenseurs[:, 1, 2] = tenseurs[:, 2, 1] = historique_contraintes[:, 4]
    tenseurs[:, 0, 2] = tenseurs[:, 2, 0] = historique_contraintes[:, 5]
    return tenseurs


def contrainte_hydrostatique(tenseurs: np.ndarray) -> np.ndarray:
    """Calcule p(t) = trace(sigma(t))/3, en MPa."""
    return np.trace(tenseurs, axis1=1, axis2=2) / 3.0


def _orientations(nb_theta: int, nb_phi: int, nb_psi: int):
    """Génère les normales de plans et les directions tangentielles."""
    for theta, phi in product(
        np.linspace(0.0, np.pi / 2.0, nb_theta, endpoint=True),
        np.linspace(0.0, 2.0 * np.pi, nb_phi, endpoint=False),
    ):
        normale = np.array([
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ])
        base_1 = np.array([
            np.cos(theta) * np.cos(phi),
            np.cos(theta) * np.sin(phi),
            -np.sin(theta),
        ])
        base_2 = np.array([-np.sin(phi), np.cos(phi), 0.0])
        for psi in np.linspace(0.0, np.pi, nb_psi, endpoint=False):
            yield normale, np.cos(psi) * base_1 + np.sin(psi) * base_2


def evaluer_dang_van(
    historique_contraintes: np.ndarray,
    materiau: MateriauDangVan,
    nb_theta: int = 10,
    nb_phi: int = 20,
    nb_psi: int = 18,
) -> ResultatDangVan:
    """Recherche le plan et la direction les plus sollicités sur un cycle.

    La contrainte de cisaillement alternée est obtenue en retirant la valeur
    médiane de l'historique de cisaillement résolu. Le maximum temporel de
    |tau_a(t)| + a*p(t) est ensuite recherché sur les orientations discrétisées.
    """
    if min(nb_theta, nb_phi, nb_psi) < 2:
        raise ValueError("Les discrétisations d'orientation doivent être au moins égales à 2.")

    tenseurs = voigt_vers_tenseur(historique_contraintes)
    pression_hydrostatique = contrainte_hydrostatique(tenseurs)
    valeur_maximale = -np.inf
    meilleure_normale = meilleure_direction = None
    meilleur_indice_temps = 0

    for normale, direction in _orientations(nb_theta, nb_phi, nb_psi):
        cisaillement_resolu = np.einsum("i,tij,j->t", direction, tenseurs, normale)
        valeur_mediane = 0.5 * (cisaillement_resolu.max() + cisaillement_resolu.min())
        cisaillement_alterne = np.abs(cisaillement_resolu - valeur_mediane)
        indicateur = cisaillement_alterne + materiau.a * pression_hydrostatique
        indice_temps = int(np.argmax(indicateur))
        if indicateur[indice_temps] > valeur_maximale:
            valeur_maximale = float(indicateur[indice_temps])
            meilleure_normale = normale.copy()
            meilleure_direction = direction.copy()
            meilleur_indice_temps = indice_temps

    coefficient_securite = np.inf if valeur_maximale <= 0 else materiau.b / valeur_maximale
    return ResultatDangVan(
        valeur_maximale, float(coefficient_securite), meilleure_normale,
        meilleure_direction, meilleur_indice_temps
    )


def creer_historique_non_proportionnel(nb_points: int = 361) -> tuple[np.ndarray, np.ndarray]:
    """Crée un cycle générique traction-torsion déphasé, exprimé en MPa."""
    phase = np.linspace(0.0, 2.0 * np.pi, nb_points)
    contraintes = np.zeros((nb_points, 6))
    contraintes[:, 0] = 80.0 + 140.0 * np.sin(phase)
    contraintes[:, 3] = 75.0 * np.sin(phase + np.pi / 2.0)
    return phase, contraintes
