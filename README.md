# Multiaxial Fatigue Assessment — Dang Van Criterion

Educational Python demonstrator for evaluating a periodic multiaxial stress
history with the Dang Van high-cycle fatigue criterion. The repository uses
only synthetic data and contains no industrial or confidential information.

## Engineering purpose

For a sampled stress tensor history, the program:

1. reconstructs the symmetric stress tensor from six Voigt components;
2. calculates the hydrostatic stress \(p(t)=\mathrm{tr}(\boldsymbol\sigma)/3\);
3. samples material-plane normals and in-plane shear directions;
4. determines the alternating resolved shear stress for every orientation;
5. finds the maximum Dang Van indicator and its critical orientation;
6. returns a fatigue safety factor.

The implemented criterion is

$$\max_t\left[\tau_a(t)+a\,p(t)\right]\leq b$$

with the parameters calibrated from fully reversed bending and torsion limits:

$$a=3\frac{\tau_{-1}}{\sigma_{-1}}-\frac{3}{2},\qquad b=\tau_{-1}$$

## Repository structure

```text
.
├── src/dang_van.py       # criterion, orientation search and data generator
├── tests/test_dang_van.py
├── example.py            # executable demonstration and chart generation
├── requirements.txt
└── LICENSE
```

## Quick start

```bash
python -m venv .venv
```

Activate the environment, then install and run:

```bash
pip install -r requirements.txt
python example.py
python -m unittest discover -s tests -v
```

The demonstration applies an out-of-phase tension–torsion history, which gives
a simple non-proportional loading path. It prints the criterion parameters,
maximum indicator and safety factor, then generates
`results/synthetic_loading.png`.

## Input convention

The stress history is an `N × 6` NumPy array expressed in MPa:

```text
[sigma_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_xz]
```

Minimal use:

```python
from dang_van import DangVanMaterial, evaluate_dang_van

material = DangVanMaterial(
    reversed_bending_limit=220.0,
    reversed_torsion_limit=130.0,
)
result = evaluate_dang_van(stress_history, material)
print(result.safety_factor)
```

## Assumptions and limitations

- stresses are supplied in MPa and represent one stabilized periodic cycle;
- fatigue limits must be consistent with the chosen material and life domain;
- plane and direction searches are discretized, so convergence must be checked;
- the alternating resolved shear is evaluated with a midrange-removal approach;
- mean stress, surface condition, size, temperature, defects and statistical
  scatter require application-specific treatment;
- this demonstrator is intended for learning and portfolio review, not for
  certifying a safety-critical component.

## Résumé en français

Ce projet illustre le post-traitement d'un historique de tenseurs de contraintes
par le critère de fatigue multiaxiale de Dang Van. Il recherche une orientation
critique, calcule l'indicateur maximal et fournit un coefficient de sécurité.
Les données sont entièrement fictives. Le code explicite les hypothèses et les
limites afin de rester vérifiable et adaptable.

## Author

**Oumar Mbengue** — Mechanical engineering, finite element analysis and
scientific computing.

Released under the MIT License.
