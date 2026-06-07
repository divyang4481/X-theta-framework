# X-Theta Audit Report: synthetic_phi_0.30

**Scientific Warning:** Phi_eff is an effective phenomenological parameter only. Without gravitational path, altitude, curvature, or spacetime-baseline metadata, this is not evidence of spacetime-induced X-Theta holonomy.

## Claim Level

- **Level:** simulation

## Summary Results

- **Valid Bell Trial Count:** 10000
- **CHSH S-statistic:** 2.587723 ± 0.030404 (Standard Error)
- **95% Bootstrap Confidence Interval:** [2.535843, 2.642170]
- **Bootstrap Samples:** 100

## Effective X-Theta Fit

- **Effective Phase ($\Phi_{eff}$):** 0.303799 rad
- **Effective Anisotropy ($R_{\Theta, eff}$):** 0.651845
- **Fit Status:** Success
- **Fit Warning:** Phi_eff is an effective phenomenological parameter only. Without gravitational path, altitude, curvature, or spacetime-baseline metadata, this is not evidence of spacetime-induced X-Theta holonomy.

## Setting Expectations and Counts

|   setting_pair |   count |   expectation |
|---------------:|--------:|--------------:|
|             00 |    2521 |      0.699326 |
|             01 |    2480 |      0.715323 |
|             10 |    2448 |      0.597222 |
|             11 |    2551 |     -0.575853 |

## Schema Validation (First Chunk)

- Missing columns: ['trial_id', 'timestamp', 'source_file']
- Alice unique settings: [0, 1]
- Bob unique settings: [1, 0]
