# Data Sources and Provenance

This file tracks where each piece of data or method comes from, what license applies,
and what paper or resource backs it. It is meant to stay current as the project grows.

## External APIs

- **PubChem**
  - Use: molecule search and compound metadata.
  - License: consult PubChem's current terms before commercial or large-scale use.
  - Notes: results are cached with a finite TTL. Input is normalized before lookup.

- **RCSB PDB**
  - Use: protein structure fetch and metadata.
  - License: consult RCSB's current usage policy.
  - Notes: PDB ID format is validated before download. Structures are cached with a finite TTL.

## Quantum

- **Qiskit and Qiskit Aer**
  - Use: Hamiltonian simulation, statevector estimation, noisy estimation, and VQE loop.
  - Notes: exact estimator import paths are confirmed against the pinned version before use.
  - Hamiltonians for the demo molecules are precomputed offline into JSON and loaded by the API.

- **Scientific reference point**
  - H2 in STO-3G near 0.74 angstrom is about -1.137 Ha. This is used as a sanity check for the exact solver and VQE loop.

## Frontend vendor libraries

- **3Dmol.js 2.5.5**
  - Use: the 3D structure viewers (Drug-Target Studio and protein panel).
  - License: BSD-3-Clause. The license file ships next to the code as
    `frontend/static/vendor/3Dmol-LICENSE-2.5.5.txt`.
  - Version: 2.5.5, the newest upstream release at the time of vendoring
    (published 2026-05-22; releases page: https://github.com/3dmol/3Dmol.js/releases/tag/2.5.5).
  - Source: the pinned npm tarball `https://registry.npmjs.org/3dmol/-/3dmol-2.5.5.tgz`,
    file `package/build/3Dmol-min.js`, copied unmodified to
    `frontend/static/vendor/3Dmol-min-2.5.5.js` (537,792 bytes).
  - Why vendored: the viewer previously loaded the library from the third-party CDN at
    3dmol.csb.pitt.edu, so a CDN hiccup showed "3Dmol.js did not load" until a hard
    refresh. The app now serves the pinned file itself.

## Molecules and properties

- **RDKit**
  - Use: SMILES normalization and validation, molecular descriptors, Lipinski rule-of-five computation.
  - Notes: descriptors are computed by us from the molecule, not fetched from an external "truth" service.

## Drug-target complexes

- Docking validation and interaction tables are served from curated artifact files.
- Every curated entry has an explicit status: `computed` or `placeholder`.
- Placeholder entries use `null` values and must not present invented RMSD or interaction data.
- Docking-validation demos should use non-covalent complexes where possible.

### The 12 curated entries

Every title, ligand code, resolution and source organism below was read from the RCSB entry
metadata (`https://data.rcsb.org/rest/v1/core/entry/{id}`,
`.../nonpolymer_entity/{id}/{entity}`, `.../polymer_entity/{id}/{entity}`) in this change, except
where a row says the original entry was left as shipped. The ligand is the drug-like non-polymer
group: water, ions, sugars, cofactors, buffers/additives and standard amino acids are excluded
and at least 10 heavy atoms are required (same threshold as the redocking pipeline). Where an
entry deposits several copies, the first deposited copy is the one analysed and the payload
lists the others.

| PDB | Target and drug | Ligand | Resolution | Tier | Why this entry | Source |
| --- | --- | --- | --- | --- | --- | --- |
| 1HSG | HIV-1 protease + indinavir | MK1 | 2.0 A (RCSB; the deposited title says 1.9 A) | curated | original viral example | As shipped in the manifest; not re-verified here |
| 1VRT | HIV-1 reverse transcriptase + nevirapine | NVP | 2.2 A | curated | original allosteric example | As shipped; not re-verified here |
| 2HU4 | Influenza A N1 neuraminidase + oseltamivir carboxylate | G39 | 2.5 A | curated | original viral example | As shipped; not re-verified here |
| 1O86 | Human ACE + lisinopril | LPR | 2.0 A | curated | original human enzyme example | As shipped; not re-verified here |
| 1M17 | Human EGFR kinase domain + erlotinib | AQ4 | 2.6 A | curated | human receptor tyrosine kinase; erlotinib is the 4-anilinoquinazoline in the deposit | Stamos, Sliwkowski, Eigenbrot, *J Biol Chem* 2002, DOI 10.1074/jbc.M207135200, PMID 12196540 |
| 2HYY | Human Abl kinase domain + imatinib | STI | 2.4 A | curated | human kinase; imatinib (STI571) is the drug in the deposit | Cowan-Jacob et al., *Acta Crystallogr D* 2007, DOI 10.1107/S0907444906047287, PMID 17164530 |
| 3ERT | Human estrogen receptor alpha + 4-hydroxytamoxifen | OHT | 1.9 A | curated | human nuclear receptor; antagonist mechanism | Shiau et al., *Cell* 1998, DOI 10.1016/S0092-8674(00)81717-1, PMID 9875847 |
| 4EY7 | Human acetylcholinesterase + donepezil | E20 | 2.35 A | curated | human enzyme, and a structure of the human rather than Torpedo enzyme | Cheung et al., *J Med Chem* 2012, DOI 10.1021/jm300871x, PMID 23035744 |
| 2H42 | Human phosphodiesterase-5 + sildenafil | VIA | 2.3 A | curated | human enzyme; the only new entry with a literature-named residue (Gly659) | Wang et al., *J Biol Chem* 2006, DOI 10.1074/jbc.M512527200, PMID 16735511 |
| 3HS4 | Human carbonic anhydrase II + acetazolamide | AZM | 1.1 A | curated | human zinc enzyme; highest-resolution entry in the set | Sippel et al., *Acta Crystallogr F* 2009, DOI 10.1107/S1744309109036665, PMID 19851004 |
| 1HWK | Human HMG-CoA reductase + atorvastatin | 117 | 2.22 A | curated | human enzyme; statin mechanism | Istvan & Deisenhofer, *Science* 2001, DOI 10.1126/science.1059344, PMID 11349148 |
| 6W63 | SARS-CoV-2 main protease + X77 | X77 | 2.1 A | contacts only | viral; the inhibitor is experimental, **not an approved drug**, and the entry has no published citation | RCSB citation is "To be published" (Mesecar), no DOI |

Tier meaning: `curated` entries carry a mechanism narrative whose statements are labelled
LITERATURE and backed by the citation in the table; `contacts_only` entries carry measured
geometry and no narrative. **No catalytic-residue list was written for any new entry from
memory.** Only 2H42 has a residue named by its cited abstract (Gly659); the others carry
`catalytic_residues_note` saying that no literature-verified residue list is curated.

### Candidates that were rejected

- **1IEP** (Abl kinase + imatinib): the deposited protein is *Mus musculus* Abl, not human, so
  it does not belong in a human-target set. 2HYY (human ABL1 + imatinib) replaces it.
- **3LN1** (COX-2 + celecoxib): the deposited protein is *Mus musculus* COX-2. No human COX-2
  celecoxib entry was added, because the available candidates checked (6COX) deposit a different
  inhibitor (SC-558), also from mouse.
- **6LU7** (SARS-CoV-2 Mpro + N3): the N3 inhibitor is deposited as a **polymer entity** (a
  peptide chain), and the entry has no non-polymer ligand at all. The project's ligand rule is
  metadata/HETATM based and finds zero candidates, so 6LU7 cannot be rendered by this pipeline.
  6W63 (same protein, a non-polymer inhibitor) is the viral replacement.

## Offline redocking validation (SPEC.md P3-3, Part A)

Part A ships the pipeline that would redock three known complexes; **it runs no docking and
produced no result file**. Part B is the Colab run that writes
`app/data/docking/redock_validation.json`. Everything below is what the pipeline reads, and
why; the sources are checked against the live services, so the field names and the ligand
codes in this section are the ones the entries actually return.

- **RCSB PDB entries (structures and metadata).** 1VRT, 2HU4 and 1HSG, chosen from SPEC.md
  section 7's rule that redocking uses **non-covalent** complexes. Which residue is the ligand
  is read from the metadata rather than assumed:
  - `https://files.rcsb.org/download/{id}.pdb` — the deposited entry (coordinates for the
    ligand, the crystal poses the RMSD is measured against, and the receptor chains).
  - `https://data.rcsb.org/rest/v1/core/entry/{id}` — `non_polymer_entity_ids`, the list the
    ligand rule starts from.
  - `https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{id}/{entity_id}` — each entity's
    `pdbx_entity_nonpoly.comp_id` and name.
  - `https://data.rcsb.org/rest/v1/core/chemcomp/{code}` — the wwPDB chemical component
    definition: `rcsb_chem_comp_info.atom_count_heavy` (the number the ligand rule filters on)
    and `rcsb_chem_comp_descriptor.SMILES_stereo` (the bond orders, formal charges and
    stereochemistry the deposited ligand is rebuilt with).
  - Licence and terms: consult RCSB's current usage policy and the wwPDB CCD terms before
    commercial or large-scale use. Structures are not cached in this repository.
  - Ligands as returned by that metadata (sanity check for a human reader, never hardcoded in
    the script): 1VRT nevirapine (`NVP`, 20 heavy atoms) in the RT p66 chain A, with a
    magnesium ion (`MG`, 1 heavy atom) present and excluded by the rule; 1HSG indinavir
    (`MK1`, 45 heavy atoms) in chain B, with both chains A and B kept because the site sits at
    the dimer interface; 2HU4 an oseltamivir-type neuraminidase inhibitor (`G39`, 20 heavy
    atoms) deposited in eight copies (A-H), of which the first deposited instance is used.
- **AutoDock Vina** (`vina==1.2.7`, Python bindings). Use: the docking search and its scores.
  Cite: Trott & Olson, *AutoDock Vina: improving the speed and accuracy of docking with a new
  scoring function, efficient optimization and multithreading*, J. Comput. Chem. 31 (2010)
  455-461, and the AutoDock Vina repository (ccsb-scripps/AutoDock-Vina) for the version's
  terms. The box is centred on the crystal ligand and one fixed setting is used for every
  target.
- **meeko** (`meeko==0.8.0`, forlilab). Use: receptor and ligand PDBQT preparation, including
  the missing receptor hydrogens. Terms: the meeko repository's licence.
- **RDKit**. Use: bond-order assignment from the RCSB ligand definition
  (`AllChem.AssignBondOrdersFromTemplate`), hydrogen addition, stereochemistry from the
  deposited coordinates, and the heavy-atom symmetry-aware RMSD without re-alignment
  (`rdMolAlign.CalcRMS`). Terms: the RDKit licence (BSD-3-Clause).
- **What is ours, not fetched.** The ligand rule, the chain-selection cutoff (5 A), the box
  (22 A cube on the crystal ligand centroid), and the RMSD method are project choices and are
  recorded in every result file. The docked poses are the one thing in this pipeline that is a
  **prediction rather than a measurement**, and the file labels them that way.

## AI and property prediction

- Optional module. Built last.
- If present, trained models are produced by us and stored with their preprocessing pipeline and manifest.
- The API loads only model files it produced itself.
- If no model is loaded, the endpoint returns an honest "not trained" or "heuristic" response.

## Licensing and attribution

- Keep this file updated whenever a new dataset, API, or method is added.
- Record the license, source URL, and the paper/resource behind each method.
- Do not paste external terms here from memory if they may have changed; verify current terms before relying on them.
