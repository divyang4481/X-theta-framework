# Open Bell/CHSH Dataset Survey

This report documents candidate datasets for validating the Bell/CHSH processing pipeline and fitting effective X-Theta parameters.

**Mandatory Scientific Warning:**
The fitted Phi_eff is an effective phenomenological parameter only. Without gravitational path, altitude, curvature, or spacetime-baseline metadata, this is not evidence of spacetime-induced X-Theta holonomy.

---

## 1. Weihs / Zeilinger 1998 photon Bell-test data
- **dataset_name:** Weihs 1998
- **experiment_type:** Photon polarization (PDC source)
- **data_url:** Locally available in repo (`Bell_test_with_weihs_data/`)
- **paper_url:** https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.81.5039
- **license_or_access_note:** Open research data preserved in repository artifacts.
- **available_format:** `.npz` (Numpy compressed archive)
- **event_level_available:** yes
- **columns_or_schema_if_known:** `settings` (0-3 representing Alice/Bob setting pairs), `products` (±1 outcome product), `phases` (likely analyzer settings).
- **expected_CHSH_or_published_result:** S ≈ 2.73 ± 0.02
- **recommended_use:** Primary validation for large streaming data processing.
- **limitations:** Data is already pre-processed into 'mined' setting pairs and products in some files.
- **classification:** confirmed_downloadable

---

## 2. Hensen et al. 2015 Delft loophole-free Bell-test data
- **dataset_name:** Hensen 2015 Delft
- **experiment_type:** Entangled electron spins (NV centers in diamond)
- **data_url:** Locally available in repo (`bell_open_data.txt`); also via NIST: https://www.nist.gov/pml/applied-physics-division/bell-test-research-software-and-data/repository-bell-test-research-3
- **paper_url:** https://www.nature.com/articles/nature15759
- **license_or_access_note:** Publicly released for the NIST Bell Test project.
- **available_format:** Plain text (space/comma separated)
- **event_level_available:** yes
- **columns_or_schema_if_known:** `timestamp`, `alice_setting`, `bob_setting`, `alice_outcome`, `bob_outcome` (mapped from raw indexes).
- **expected_CHSH_or_published_result:** S ≈ 2.42 ± 0.20 (based on 245 trials)
- **recommended_use:** Validation for 'loophole-free' statistical workflows and raw event mapping.
- **limitations:** Small number of events compared to photon tests.
- **classification:** confirmed_downloadable

---

## 3. BIG Bell Test 2018 public datasets
- **dataset_name:** BIG Bell Test 2018
- **experiment_type:** Multiple worldwide labs (photons, ions, etc.) with human-generated settings.
- **data_url:** NIST Repository: https://www.nist.gov/pml/applied-physics-division/bell-test-research-software-and-data/repository-bell-test-research-3
- **paper_url:** https://www.nature.com/articles/s41586-018-0085-3
- **license_or_access_note:** Publicly available research data.
- **available_format:** `.dat.zip`, `.csv` via NIST servers.
- **event_level_available:** yes
- **columns_or_schema_if_known:** Standard Bell schema with human choices as settings.
- **expected_CHSH_or_published_result:** Varies by lab; generally strong violation of S > 2.
- **recommended_use:** Testing multi-lab adapter support and high-volume data validation.
- **limitations:** Requires separate download from NIST servers.
- **classification:** confirmed_downloadable
