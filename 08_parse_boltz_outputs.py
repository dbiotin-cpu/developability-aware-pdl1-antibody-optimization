from pathlib import Path
import csv
import json
import shutil
import re


BOLTZ_OUT_DIR = Path("data/structure_eval/boltz_outputs/top_variants")
TOP_METADATA = Path("data/structure_eval/top_variants_metadata.csv")

OUT_PREDICTED_STRUCTURES = Path("data/structure_eval/predicted_structures")
OUT_SUMMARY = Path("data/structure_eval/boltz_structure_summary.csv")
OUT_UNMATCHED = Path("data/structure_eval/boltz_unmatched_files.txt")


STRUCTURE_SUFFIXES = {".cif", ".pdb"}
CONFIDENCE_SUFFIXES = {".json", ".npz"}


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


def sanitize_filename(name):
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
    return name.strip("_")


def normalize_for_matching(text):
    text = str(text).lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    text = re.sub(r"_+", "_", text)
    return text.strip("_")


def collect_files(root):
    if not root.exists():
        raise FileNotFoundError(f"Boltz output directory not found: {root}")

    files = [p for p in root.rglob("*") if p.is_file()]
    return files


def associate_files_to_variant(files, variant_id):
    """
    Match files by variant_id or sanitized variant_id in path.
    This is robust to Boltz creating nested output directories.
    """
    variant_norm = normalize_for_matching(variant_id)
    variant_safe_norm = normalize_for_matching(sanitize_filename(variant_id))

    matched = []

    for path in files:
        path_norm = normalize_for_matching(str(path))

        if variant_norm in path_norm or variant_safe_norm in path_norm:
            matched.append(path)

    return matched


def structure_priority(path):
    """
    Prefer CIF/PDB files that look like model outputs.
    Lower score is better.
    """
    name = path.name.lower()
    score = 100

    if path.suffix.lower() == ".cif":
        score -= 10
    if "model" in name:
        score -= 10
    if "rank" in name:
        score -= 5
    if "prediction" in str(path).lower():
        score -= 3
    if "aligned" in name:
        score += 5

    return score


def choose_representative_structure(structure_files):
    if not structure_files:
        return None

    structure_files = sorted(structure_files, key=lambda p: (structure_priority(p), len(str(p))))
    return structure_files[0]


def flatten_json(obj, prefix=""):
    """
    Flatten JSON recursively and keep simple scalar values.
    """
    items = {}

    if isinstance(obj, dict):
        for k, v in obj.items():
            new_prefix = f"{prefix}.{k}" if prefix else str(k)
            items.update(flatten_json(v, new_prefix))

    elif isinstance(obj, list):
        # Keep short numeric lists summarized.
        numeric_values = [x for x in obj if isinstance(x, (int, float))]
        if numeric_values and len(numeric_values) == len(obj):
            items[f"{prefix}.mean"] = sum(numeric_values) / len(numeric_values)
            items[f"{prefix}.min"] = min(numeric_values)
            items[f"{prefix}.max"] = max(numeric_values)
        else:
            for i, v in enumerate(obj[:10]):
                new_prefix = f"{prefix}[{i}]"
                items.update(flatten_json(v, new_prefix))

    else:
        if isinstance(obj, (int, float, str, bool)) or obj is None:
            items[prefix] = obj

    return items


def extract_json_metrics(json_files):
    """
    Boltz output schema can vary by version.
    This extracts useful-looking scalar metrics from JSON files.
    """
    extracted = {}

    useful_patterns = [
        "confidence",
        "plddt",
        "ptm",
        "iptm",
        "pae",
        "score",
        "prob",
        "clash",
        "ranking",
    ]

    for path in json_files:
        try:
            with open(path, "r") as f:
                data = json.load(f)
        except Exception:
            continue

        flat = flatten_json(data)

        for key, value in flat.items():
            key_lower = key.lower()

            if not any(pattern in key_lower for pattern in useful_patterns):
                continue

            if isinstance(value, (int, float)):
                metric_name = f"{path.stem}:{key}"
                extracted[metric_name] = round(float(value), 6)
            elif isinstance(value, str):
                # Keep short string labels if relevant.
                if len(value) <= 80:
                    metric_name = f"{path.stem}:{key}"
                    extracted[metric_name] = value

    return extracted


def compact_metrics(metrics, max_items=12):
    if not metrics:
        return ""

    preferred_order = [
        "confidence",
        "iptm",
        "ptm",
        "plddt",
        "score",
        "pae",
        "prob",
        "clash",
    ]

    def metric_sort_key(item):
        key, _value = item
        key_lower = key.lower()

        for i, pattern in enumerate(preferred_order):
            if pattern in key_lower:
                return (i, key_lower)

        return (999, key_lower)

    items = sorted(metrics.items(), key=metric_sort_key)

    chunks = []
    for key, value in items[:max_items]:
        chunks.append(f"{key}={value}")

    return "; ".join(chunks)


def main():
    metadata_rows = read_csv(TOP_METADATA)
    all_files = collect_files(BOLTZ_OUT_DIR)

    OUT_PREDICTED_STRUCTURES.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    all_matched_files = set()

    for meta in metadata_rows:
        variant_id = meta["variant_id"]

        matched = associate_files_to_variant(all_files, variant_id)
        all_matched_files.update(matched)

        structure_files = [
            p for p in matched
            if p.suffix.lower() in STRUCTURE_SUFFIXES
        ]

        json_files = [
            p for p in matched
            if p.suffix.lower() == ".json"
        ]

        npz_files = [
            p for p in matched
            if p.suffix.lower() == ".npz"
        ]

        representative_structure = choose_representative_structure(structure_files)

        copied_structure_path = ""

        if representative_structure:
            safe_id = sanitize_filename(variant_id)
            out_name = f"{safe_id}_boltz{representative_structure.suffix.lower()}"
            copied_structure = OUT_PREDICTED_STRUCTURES / out_name

            shutil.copy2(representative_structure, copied_structure)
            copied_structure_path = str(copied_structure)

        metrics = extract_json_metrics(json_files)
        compact_metric_text = compact_metrics(metrics)

        summary_rows.append(
            {
                "variant_id": variant_id,
                "variant_type": meta.get("variant_type", ""),
                "mutations": meta.get("mutations", ""),
                "num_mutations": meta.get("num_mutations", ""),
                "pre_structure_risk_score": meta.get("risk_score", ""),
                "pre_structure_optimization_score": meta.get("optimization_score_pre_structure", ""),
                "boltz_matched_file_count": len(matched),
                "boltz_structure_file_count": len(structure_files),
                "boltz_json_file_count": len(json_files),
                "boltz_npz_file_count": len(npz_files),
                "representative_structure_source": str(representative_structure) if representative_structure else "",
                "representative_structure_copied": copied_structure_path,
                "available_boltz_metrics": compact_metric_text,
            }
        )

    unmatched_files = [p for p in all_files if p not in all_matched_files]

    with open(OUT_UNMATCHED, "w") as f:
        for path in unmatched_files:
            f.write(str(path) + "\n")

    fieldnames = [
        "variant_id",
        "variant_type",
        "mutations",
        "num_mutations",
        "pre_structure_risk_score",
        "pre_structure_optimization_score",
        "boltz_matched_file_count",
        "boltz_structure_file_count",
        "boltz_json_file_count",
        "boltz_npz_file_count",
        "representative_structure_source",
        "representative_structure_copied",
        "available_boltz_metrics",
    ]

    write_csv(OUT_SUMMARY, summary_rows, fieldnames)

    print("Saved:")
    print(f"  {OUT_SUMMARY}")
    print(f"  {OUT_PREDICTED_STRUCTURES}")
    print(f"  {OUT_UNMATCHED}")

    print("\nParsed variants:")
    for row in summary_rows:
        print(
            f"  {row['variant_id']} | "
            f"structures={row['boltz_structure_file_count']} | "
            f"json={row['boltz_json_file_count']}"
        )


if __name__ == "__main__":
    main()