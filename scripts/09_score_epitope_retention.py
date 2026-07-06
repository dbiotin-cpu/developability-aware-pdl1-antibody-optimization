from pathlib import Path
import csv
import math
import re
import numpy as np

from Bio.PDB import PDBParser, MMCIFParser, Superimposer


BOLTZ_SUMMARY = Path("data/structure_eval/boltz_structure_summary.csv")

REFERENCE_PDB = Path("data/input/reference_wt_complex.pdb")
REFERENCE_CIF = Path("data/input/reference_wt_complex.cif")

OUT_INTERFACE_SUMMARY = Path("data/structure_eval/interface_contact_summary.csv")
OUT_RETENTION_SCORES = Path("data/structure_eval/epitope_retention_scores.csv")

TARGET_CHAIN_ID = "A"   # PD-L1 chain from Boltz YAML
BINDER_CHAIN_ID = "H"   # Antibody/binder chain from Boltz YAML

CONTACT_CUTOFF_A = 5.0


def read_csv(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")

    with open(path, "r") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_structure(path):
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Structure file not found: {path}")

    if path.suffix.lower() == ".pdb":
        parser = PDBParser(QUIET=True)
    elif path.suffix.lower() == ".cif":
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported structure format: {path}")

    return parser.get_structure(path.stem, str(path))


def get_first_model(structure):
    return next(structure.get_models())


def get_chain(model, preferred_id, fallback_index=None):
    chains = list(model.get_chains())

    for chain in chains:
        if chain.id == preferred_id:
            return chain

    if fallback_index is not None and len(chains) > fallback_index:
        print(
            f"WARNING: Preferred chain {preferred_id} not found. "
            f"Using fallback chain {chains[fallback_index].id}."
        )
        return chains[fallback_index]

    available = [chain.id for chain in chains]
    raise ValueError(f"Chain {preferred_id} not found. Available chains: {available}")


def is_standard_residue(residue):
    hetflag, resseq, icode = residue.id
    return hetflag == " "


def residue_key(residue):
    hetflag, resseq, icode = residue.id
    return f"{resseq}{icode.strip()}"


def residue_index(residue):
    hetflag, resseq, icode = residue.id
    return int(resseq)


def get_heavy_atom_coords(residue):
    coords = []

    for atom in residue.get_atoms():
        if atom.element == "H":
            continue
        coords.append(atom.coord)

    if not coords:
        return np.empty((0, 3))

    return np.array(coords)


def get_chain_residues(chain):
    return [res for res in chain.get_residues() if is_standard_residue(res)]


def calculate_interface_contacts(target_chain, binder_chain, cutoff=CONTACT_CUTOFF_A):
    target_residues = get_chain_residues(target_chain)
    binder_residues = get_chain_residues(binder_chain)

    contact_pairs = set()
    target_contact_residues = set()
    binder_contact_residues = set()

    for t_res in target_residues:
        t_coords = get_heavy_atom_coords(t_res)
        if t_coords.size == 0:
            continue

        for b_res in binder_residues:
            b_coords = get_heavy_atom_coords(b_res)
            if b_coords.size == 0:
                continue

            diff = t_coords[:, None, :] - b_coords[None, :, :]
            distances = np.sqrt(np.sum(diff * diff, axis=2))

            if np.any(distances <= cutoff):
                t_key = residue_key(t_res)
                b_key = residue_key(b_res)

                contact_pairs.add((t_key, b_key))
                target_contact_residues.add(t_key)
                binder_contact_residues.add(b_key)

    return {
        "contact_pairs": contact_pairs,
        "target_contact_residues": target_contact_residues,
        "binder_contact_residues": binder_contact_residues,
    }


def get_ca_atoms(chain):
    atoms = []

    for residue in get_chain_residues(chain):
        if "CA" in residue:
            atoms.append(residue["CA"])

    return atoms


def calculate_binder_rmsd_after_target_alignment(ref_structure, mob_structure):
    """
    Align mobile structure onto reference using target chain CA atoms.
    Then calculate binder CA RMSD.

    This approximates binding-mode preservation.
    """
    ref_model = get_first_model(ref_structure)
    mob_model = get_first_model(mob_structure)

    ref_target = get_chain(ref_model, TARGET_CHAIN_ID, fallback_index=0)
    mob_target = get_chain(mob_model, TARGET_CHAIN_ID, fallback_index=0)

    ref_binder = get_chain(ref_model, BINDER_CHAIN_ID, fallback_index=1)
    mob_binder = get_chain(mob_model, BINDER_CHAIN_ID, fallback_index=1)

    ref_target_ca = get_ca_atoms(ref_target)
    mob_target_ca = get_ca_atoms(mob_target)

    n_target = min(len(ref_target_ca), len(mob_target_ca))

    if n_target < 3:
        return None

    sup = Superimposer()
    sup.set_atoms(ref_target_ca[:n_target], mob_target_ca[:n_target])

    # Apply transformation to the whole mobile structure.
    sup.apply(list(mob_structure.get_atoms()))

    ref_binder_ca = get_ca_atoms(ref_binder)
    mob_binder_ca = get_ca_atoms(mob_binder)

    n_binder = min(len(ref_binder_ca), len(mob_binder_ca))

    if n_binder < 3:
        return None

    squared = []

    for ref_atom, mob_atom in zip(ref_binder_ca[:n_binder], mob_binder_ca[:n_binder]):
        diff = ref_atom.coord - mob_atom.coord
        squared.append(np.sum(diff * diff))

    rmsd = math.sqrt(sum(squared) / len(squared))

    return round(rmsd, 3)


def retention_fraction(ref_set, var_set):
    if not ref_set:
        return ""

    return round(len(ref_set & var_set) / len(ref_set), 3)


def jaccard_index(set_a, set_b):
    if not set_a and not set_b:
        return ""

    union = set_a | set_b

    if not union:
        return ""

    return round(len(set_a & set_b) / len(union), 3)


def rmsd_component(rmsd):
    if rmsd == "" or rmsd is None:
        return ""

    try:
        rmsd = float(rmsd)
    except ValueError:
        return ""

    # 0 Å = 1.0, 10 Å or worse = 0.0
    return round(max(0.0, 1.0 - rmsd / 10.0), 3)


def calculate_retention_score(row):
    """
    Higher is better.

    Main idea:
    - epitope retention is most important
    - contact-pair preservation matters
    - binder contact retention matters
    - binder RMSD after PD-L1 alignment captures gross binding-mode preservation
    """
    required = [
        "epitope_contact_retention",
        "contact_pair_jaccard_vs_reference",
        "binder_contact_retention",
    ]

    for key in required:
        if row.get(key, "") == "":
            return "", "Reference unavailable"

    epitope = float(row["epitope_contact_retention"])
    pair_jaccard = float(row["contact_pair_jaccard_vs_reference"])
    binder = float(row["binder_contact_retention"])

    rmsd_comp = rmsd_component(row.get("binder_ca_rmsd_after_target_alignment_A", ""))

    if rmsd_comp == "":
        rmsd_comp = 0.5

    score = (
        0.45 * epitope
        + 0.25 * pair_jaccard
        + 0.20 * binder
        + 0.10 * float(rmsd_comp)
    )

    score = round(score, 3)

    if score >= 0.75:
        decision = "Binding mode likely preserved"
    elif score >= 0.50:
        decision = "Partially preserved; review structure"
    else:
        decision = "Low retention; deprioritize"

    return score, decision


def find_reference_structure():
    if REFERENCE_PDB.exists():
        return REFERENCE_PDB

    if REFERENCE_CIF.exists():
        return REFERENCE_CIF

    return None


def analyze_structure(path):
    structure = load_structure(path)
    model = get_first_model(structure)

    target_chain = get_chain(model, TARGET_CHAIN_ID, fallback_index=0)
    binder_chain = get_chain(model, BINDER_CHAIN_ID, fallback_index=1)

    contacts = calculate_interface_contacts(
        target_chain=target_chain,
        binder_chain=binder_chain,
        cutoff=CONTACT_CUTOFF_A,
    )

    return structure, contacts


def main():
    summary_rows = read_csv(BOLTZ_SUMMARY)

    reference_path = find_reference_structure()

    reference_structure = None
    reference_contacts = None

    if reference_path:
        print(f"Using reference structure: {reference_path}")
        reference_structure, reference_contacts = analyze_structure(reference_path)
        print(f"Reference epitope residues: {len(reference_contacts['target_contact_residues'])}")
        print(f"Reference binder contact residues: {len(reference_contacts['binder_contact_residues'])}")
        print(f"Reference contact pairs: {len(reference_contacts['contact_pairs'])}")
    else:
        print("WARNING: No WT/reference complex found.")
        print("Expected one of:")
        print(f"  {REFERENCE_PDB}")
        print(f"  {REFERENCE_CIF}")
        print("The script will calculate interface metrics, but not WT-relative retention.")

    interface_rows = []
    retention_rows = []

    for meta in summary_rows:
        variant_id = meta["variant_id"]
        structure_path = meta.get("representative_structure_copied", "")

        if not structure_path:
            print(f"Skipping {variant_id}: no representative structure.")
            continue

        print(f"\nAnalyzing {variant_id}")
        print(f"  Structure: {structure_path}")

        try:
            variant_structure, variant_contacts = analyze_structure(structure_path)
        except Exception as e:
            print(f"  ERROR: {e}")
            continue

        target_contacts = variant_contacts["target_contact_residues"]
        binder_contacts = variant_contacts["binder_contact_residues"]
        contact_pairs = variant_contacts["contact_pairs"]

        interface_rows.append(
            {
                "variant_id": variant_id,
                "structure_file": structure_path,
                "target_chain_id": TARGET_CHAIN_ID,
                "binder_chain_id": BINDER_CHAIN_ID,
                "contact_cutoff_A": CONTACT_CUTOFF_A,
                "target_contact_residue_count": len(target_contacts),
                "binder_contact_residue_count": len(binder_contacts),
                "contact_pair_count": len(contact_pairs),
                "target_contact_residues": ";".join(sorted(target_contacts, key=lambda x: int(re.sub(r'[^0-9]', '', x) or 0))),
                "binder_contact_residues": ";".join(sorted(binder_contacts, key=lambda x: int(re.sub(r'[^0-9]', '', x) or 0))),
            }
        )

        if reference_contacts:
            epitope_retention = retention_fraction(
                reference_contacts["target_contact_residues"],
                target_contacts,
            )

            binder_retention = retention_fraction(
                reference_contacts["binder_contact_residues"],
                binder_contacts,
            )

            pair_jaccard = jaccard_index(
                reference_contacts["contact_pairs"],
                contact_pairs,
            )

            interface_ratio = ""
            if reference_contacts["contact_pairs"]:
                interface_ratio = round(
                    len(contact_pairs) / len(reference_contacts["contact_pairs"]),
                    3,
                )

            try:
                # Reload variant structure because RMSD alignment modifies coordinates.
                variant_structure_for_rmsd = load_structure(structure_path)
                binder_rmsd = calculate_binder_rmsd_after_target_alignment(
                    ref_structure=reference_structure,
                    mob_structure=variant_structure_for_rmsd,
                )
            except Exception as e:
                print(f"  RMSD calculation failed: {e}")
                binder_rmsd = ""

            retention_row = {
                "variant_id": variant_id,
                "mutations": meta.get("mutations", ""),
                "num_mutations": meta.get("num_mutations", ""),
                "pre_structure_risk_score": meta.get("pre_structure_risk_score", ""),
                "pre_structure_optimization_score": meta.get("pre_structure_optimization_score", ""),
                "structure_file": structure_path,
                "epitope_contact_retention": epitope_retention,
                "binder_contact_retention": binder_retention,
                "contact_pair_jaccard_vs_reference": pair_jaccard,
                "interface_contact_count_ratio_vs_reference": interface_ratio,
                "binder_ca_rmsd_after_target_alignment_A": binder_rmsd if binder_rmsd is not None else "",
                "target_contact_residue_count": len(target_contacts),
                "binder_contact_residue_count": len(binder_contacts),
                "contact_pair_count": len(contact_pairs),
                "reference_target_contact_residue_count": len(reference_contacts["target_contact_residues"]),
                "reference_binder_contact_residue_count": len(reference_contacts["binder_contact_residues"]),
                "reference_contact_pair_count": len(reference_contacts["contact_pairs"]),
            }

            score, decision = calculate_retention_score(retention_row)

            retention_row["epitope_retention_score"] = score
            retention_row["structure_retention_decision"] = decision

            retention_rows.append(retention_row)

            print(f"  Epitope retention: {epitope_retention}")
            print(f"  Contact-pair Jaccard: {pair_jaccard}")
            print(f"  Binder RMSD after PD-L1 alignment: {binder_rmsd}")
            print(f"  Retention score: {score}")
            print(f"  Decision: {decision}")

        else:
            retention_rows.append(
                {
                    "variant_id": variant_id,
                    "mutations": meta.get("mutations", ""),
                    "num_mutations": meta.get("num_mutations", ""),
                    "pre_structure_risk_score": meta.get("pre_structure_risk_score", ""),
                    "pre_structure_optimization_score": meta.get("pre_structure_optimization_score", ""),
                    "structure_file": structure_path,
                    "epitope_contact_retention": "",
                    "binder_contact_retention": "",
                    "contact_pair_jaccard_vs_reference": "",
                    "interface_contact_count_ratio_vs_reference": "",
                    "binder_ca_rmsd_after_target_alignment_A": "",
                    "target_contact_residue_count": len(target_contacts),
                    "binder_contact_residue_count": len(binder_contacts),
                    "contact_pair_count": len(contact_pairs),
                    "reference_target_contact_residue_count": "",
                    "reference_binder_contact_residue_count": "",
                    "reference_contact_pair_count": "",
                    "epitope_retention_score": "",
                    "structure_retention_decision": "Reference unavailable",
                }
            )

            print(f"  Target contact residues: {len(target_contacts)}")
            print(f"  Binder contact residues: {len(binder_contacts)}")
            print(f"  Contact pairs: {len(contact_pairs)}")

    interface_fields = [
        "variant_id",
        "structure_file",
        "target_chain_id",
        "binder_chain_id",
        "contact_cutoff_A",
        "target_contact_residue_count",
        "binder_contact_residue_count",
        "contact_pair_count",
        "target_contact_residues",
        "binder_contact_residues",
    ]

    retention_fields = [
        "variant_id",
        "mutations",
        "num_mutations",
        "pre_structure_risk_score",
        "pre_structure_optimization_score",
        "structure_file",
        "epitope_contact_retention",
        "binder_contact_retention",
        "contact_pair_jaccard_vs_reference",
        "interface_contact_count_ratio_vs_reference",
        "binder_ca_rmsd_after_target_alignment_A",
        "target_contact_residue_count",
        "binder_contact_residue_count",
        "contact_pair_count",
        "reference_target_contact_residue_count",
        "reference_binder_contact_residue_count",
        "reference_contact_pair_count",
        "epitope_retention_score",
        "structure_retention_decision",
    ]

    write_csv(OUT_INTERFACE_SUMMARY, interface_rows, interface_fields)
    write_csv(OUT_RETENTION_SCORES, retention_rows, retention_fields)

    print("\nSaved:")
    print(f"  {OUT_INTERFACE_SUMMARY}")
    print(f"  {OUT_RETENTION_SCORES}")


if __name__ == "__main__":
    main()
