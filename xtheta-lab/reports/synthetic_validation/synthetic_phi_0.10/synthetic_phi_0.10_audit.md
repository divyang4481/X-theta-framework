# X-Theta Audit Report: synthetic_phi_0.10

**Scientific Warning:** Phi_eff is an effective phenomenological parameter only. Without gravitational path, altitude, curvature, or spacetime-baseline metadata, this is not evidence of spacetime-induced X-Theta holonomy.

## Claim Level

- **Level:** simulation

## Summary Results

- **Valid Bell Trial Count:** 10000
- **CHSH S-statistic:** 2.807060 ± 0.028489 (Standard Error)
- **95% Bootstrap Confidence Interval:** [2.753761, 2.858892]
- **Bootstrap Samples:** 100

## Effective X-Theta Fit

- **Effective Phase ($\Phi_{eff}$):** 0.087193 rad
- **Effective Anisotropy ($R_{\Theta, eff}$):** 0.060206
- **Fit Status:** Success
- **Fit Warning:** Phi_eff is an effective phenomenological parameter only. Without gravitational path, altitude, curvature, or spacetime-baseline metadata, this is not evidence of spacetime-induced X-Theta holonomy.

## Setting Expectations and Counts

|   setting_pair |   count |   expectation |
|---------------:|--------:|--------------:|
|             00 |    2521 |      0.699326 |
|             01 |    2480 |      0.715323 |
|             10 |    2448 |      0.70915  |
|             11 |    2551 |     -0.683261 |

## Schema Validation (First Chunk)

- Missing columns: ['trial_id', 'timestamp', 'source_file']
- Alice unique settings: [0, 1]
- Bob unique settings: [1, 0]
