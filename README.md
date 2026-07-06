# Developability-Aware Optimization of an AI-Designed PD-L1 Antibody

## Project Overview

This project builds a developability-aware optimization pipeline for an AI-designed PD-L1 antibody/binder candidate.

The workflow starts from a final PD-L1-binding antibody candidate generated in a prior computational design project and asks the next biologics-discovery question:

> Can an AI-designed antibody candidate be triaged and improved for developability while preserving the predicted PD-L1 binding mode?

This repository demonstrates a practical antibody-engineering workflow that connects computational antibody design, sequence liability analysis, rational mutation proposal, structural re-evaluation, and final candidate ranking.

---

## Motivation

AI-based antibody and protein design can generate candidate binders, but therapeutic antibody discovery requires more than predicted binding. Candidates must also be evaluated for developability risks such as sequence liabilities, aggregation-prone motifs, high hydrophobicity, unfavorable charge properties, chemical instability, and loss of binding mode after optimization.

This project was designed to show the next step after computational antibody generation:

1. Start with an AI-designed PD-L1 antibody candidate.
2. Annotate CDR and framework regions.
3. Identify sequence liabilities.
4. Propose rational developability-improving mutations.
5. Re-evaluate whether variants preserve the predicted PD-L1 binding mode.
6. Rank optimized variants using both developability and structural-retention criteria.

---

## Workflow

```text
AI-designed PD-L1 antibody candidate
        ↓
Antibody numbering and CDR annotation
        ↓
Sequence liability analysis
        ↓
Developability feature calculation
        ↓
Rational mutation proposal
        ↓
Variant developability re-scoring
        ↓
Boltz structure re-evaluation
        ↓
Epitope/contact retention scoring
        ↓
Final ranked variants
```

---

## Key Features

This pipeline performs:

* IMGT, Kabat, and Chothia antibody numbering
* CDR and framework region annotation
* Sequence liability screening
* Developability feature calculation
* Rational mutation proposal
* Variant developability re-scoring
* Boltz-based structure re-evaluation
* PD-L1 epitope/contact retention scoring
* Final integrated variant ranking

---

## Repository Structure

```text
.
├── scripts/
│   ├── 01_number_antibody_regions.py
│   ├── 02_scan_sequence_liabilities.py
│   ├── 03_calculate_developability_features.py
│   ├── 04_propose_rational_variants.py
│   ├── 05_calculate_variant_developability_features.py
│   ├── 06_make_summary_figures.py
│   ├── 07_prepare_structure_evaluation_inputs.py
│   ├── 08_parse_boltz_outputs.py
│   ├── 09_score_epitope_retention.py
│   ├── 10_final_variant_ranking.py
│   └── 11_make_final_figures.py
│
├── data/
│   ├── input/
│   ├── intermediate/
│   ├── structure_eval/
│   └── output/
│
├── figures/
├── docs/
└── README.md
```

---

## Input Data

The main input files are:

```text
data/input/final_pdl1_candidate.fasta
data/input/pdl1_target.fasta
data/input/reference_wt_complex.pdb
```

The reference complex represents the parent AI-designed PD-L1 binder before developability optimization. In this workflow, “WT” refers to the original parent candidate sequence, not a naturally occurring wild-type antibody.

---

## 1. Antibody Numbering and CDR Annotation

The parent antibody candidate was numbered using multiple antibody numbering schemes:

* IMGT
* Kabat
* Chothia

CDR and framework regions were extracted to support downstream liability mapping and mutation design.

Output files:

```text
data/intermediate/numbered_sequences.csv
data/intermediate/cdr_annotation.csv
```

---

## 2. Sequence Liability Analysis

The candidate was screened for common early developability risk motifs, including:

* N-linked glycosylation motifs
* Asn deamidation-prone motifs
* Asp isomerization motifs
* oxidation-prone Met and Trp residues
* cysteines requiring disulfide/free-Cys review
* polybasic patches
* hydrophobic patches

Each liability was mapped to the antibody region where it occurred, allowing CDR-localized liabilities to be prioritized for structural review.

Output files:

```text
data/intermediate/liability_report.csv
data/intermediate/liability_summary.csv
```

---

## 3. Developability Feature Calculation

The pipeline calculated simple sequence-level and region-level developability descriptors, including:

* sequence length
* molecular weight
* pI
* approximate net charge at pH 7.4
* whole-sequence hydrophobic fraction
* CDR hydrophobic fraction
* CDR aromatic fraction
* acidic/basic residue fractions
* liability counts
* CDR-localized liability counts

Output files:

```text
data/intermediate/developability_features.csv
data/intermediate/region_developability_features.csv
```

---

## 4. Rational Variant Proposal

Rational mutations were proposed using conservative antibody-engineering rules.

Examples of mutation logic:

```text
N-linked glycosylation motif  → N→Q
Asn deamidation motif         → N→Q
Asp isomerization motif       → D→E
Met oxidation outside CDR     → M→L
Polybasic patch outside CDR   → K/R→Q
```

Hydrophobic patches and cysteines were flagged but not automatically mutated from sequence alone, because they require structural context to distinguish true developability liabilities from buried stabilizing residues or antigen-contacting paratope residues.

Output files:

```text
data/intermediate/proposed_mutations.csv
data/intermediate/proposed_variants.csv
data/intermediate/proposed_variants.fasta
```

---

## 5. Pre-Structure Variant Ranking

Proposed variants were re-scored using the same developability pipeline applied to the parent candidate. Each variant was compared against the parent for:

* total sequence-liability burden
* CDR-localized liabilities
* pI
* net charge
* CDR hydrophobicity
* CDR aromaticity
* simple developability risk score

Output files:

```text
data/intermediate/variant_liability_report.csv
data/intermediate/variant_developability_features.csv
data/intermediate/variant_ranking_pre_structure.csv
data/output/top_pre_structure_variants.csv
```

Key figures:

![Pre-structure risk score](figures/pre_structure_risk_score_by_variant.png)

![Pre-structure optimization score](figures/pre_structure_optimization_score_by_variant.png)

![Liability count by variant](figures/liability_count_by_variant.png)

![CDR hydrophobicity by variant](figures/cdr_hydrophobicity_by_variant.png)

![Liability type summary](figures/liability_type_summary.png)

---

## 6. Boltz Structure Re-Evaluation

Top-ranked variants from the pre-structure screen were prepared for Boltz-based structure prediction.

The structure-evaluation inputs were generated as:

```text
data/structure_eval/boltz_inputs/
data/structure_eval/rf2_inputs/
data/structure_eval/proteinmpnn_inputs/
```

Boltz-predicted structures were parsed and linked back to each variant.

Output files:

```text
data/structure_eval/top_variants.fasta
data/structure_eval/top_variants_metadata.csv
data/structure_eval/boltz_structure_summary.csv
data/structure_eval/predicted_structures/
```

Raw Boltz output directories are excluded from the GitHub repository because they can be large. Representative structures and parsed summary tables are included instead.

---

## 7. Epitope and Binding-Mode Retention Scoring

Each optimized variant was compared against the parent PD-L1/binder reference complex.

The structure-based scoring included:

* PD-L1 epitope contact residues
* binder contact residues
* interface contact-pair overlap
* epitope contact retention
* binder contact retention
* binder RMSD after PD-L1 alignment
* integrated epitope retention score

Output files:

```text
data/structure_eval/interface_contact_summary.csv
data/structure_eval/epitope_retention_scores.csv
```

This step evaluates whether developability-improving mutations preserve the predicted PD-L1 binding mode.

---

## 8. Final Integrated Ranking

The final ranking combines:

```text
developability improvement
+ structure retention
+ interface preservation
+ mutation burden
```

Final output files:

```text
data/output/final_ranked_variants.csv
data/output/final_ranked_variants_readme_summary.csv
data/output/final_recommendation.md
```

The final ranking is intended as an in silico triage framework. It is not experimental evidence of binding, stability, or therapeutic developability.

Key final figures:

![Final integrated score](figures/final_integrated_score_by_variant.png)

![Final score components](figures/final_score_components_top_variants.png)

![Developability vs structure retention](figures/developability_vs_structure_retention.png)

![Final decision counts](figures/final_decision_counts.png)

---

## Main Result

The final output identifies variants that improve sequence-level developability metrics while preserving the predicted PD-L1 binding mode.

See:

```text
data/output/final_ranked_variants.csv
data/output/final_recommendation.md
```

For a compact summary table:

```text
data/output/final_ranked_variants_readme_summary.csv
```

---

## How to Run the Pipeline

The pipeline was developed in a Linux/RunPod environment. Local execution is recommended through Linux or WSL Ubuntu rather than native Windows because AbNumber/ANARCI depends on HMMER.

Activate the sequence-analysis environment:

```bash
conda activate abdev
```

Run scripts in order:

```bash
python scripts/01_number_antibody_regions.py
python scripts/02_scan_sequence_liabilities.py
python scripts/03_calculate_developability_features.py
python scripts/04_propose_rational_variants.py
python scripts/05_calculate_variant_developability_features.py
python scripts/06_make_summary_figures.py
python scripts/07_prepare_structure_evaluation_inputs.py
python scripts/08_parse_boltz_outputs.py
python scripts/09_score_epitope_retention.py
python scripts/10_final_variant_ranking.py
python scripts/11_make_final_figures.py
```

Boltz predictions were run separately using the Boltz environment and the YAML files in:

```text
data/structure_eval/boltz_inputs/
```

---

## Environment

Main sequence-analysis environment:

```text
environment_abdev.yml
```

Boltz structure-prediction environment:

```text
environment_boltz.yml
```

Core Python packages used include:

* pandas
* matplotlib
* Biopython
* AbNumber
* Boltz

---

## Important Limitations

This project is an in silico portfolio workflow and should not be interpreted as experimental validation.

Key limitations:

1. Sequence-liability rules are heuristic.
2. pI and net charge calculations are approximate.
3. Sequence-only hydrophobic patch analysis lacks solvent-exposure information.
4. Predicted structures may not reproduce the true binding pose.
5. Epitope/contact retention depends on the quality of the reference complex.
6. Boltz-based structural results are used for triage, not definitive validation.
7. No experimental expression, purification, stability, aggregation, binding, or functional data are included.

The purpose of this workflow is candidate triage and portfolio demonstration, not therapeutic nomination.

---

## Recommended Experimental Validation

Top-ranked variants would require experimental validation by:

* transient mammalian expression
* affinity purification
* analytical SEC for aggregation/purity
* DSF or nanoDSF thermal stability
* SPR/BLI binding to recombinant PD-L1
* flow cytometry binding to PD-L1-expressing cells
* specificity and polyspecificity assessment

A realistic experimental follow-up would compare the parent candidate and top-ranked variants side by side for expression, aggregation, thermal stability, binding kinetics, and cell-surface PD-L1 binding.

---

## Portfolio Message

This project demonstrates the next practical step after AI antibody design:

> I can not only generate an AI-designed antibody candidate, but also evaluate whether it has a path toward developable biologic-like properties.

The workflow integrates computational design, antibody-specific annotation, developability triage, rational mutation design, structural re-evaluation, and final candidate ranking.

