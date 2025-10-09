import argparse
from pathlib import Path
from typing import Dict, List

import pandas as pd

PLANT_CONFIG: Dict[str, Dict[str, List[str]]] = {
    "indigo": {
        "IndigoRoots": [
            "Indigo_RootsS4",
            "Indigo_RootsS5",
            "Indigo_RootsS6",
        ],
        "IndigoNods": [
            "Indigo_nodS4",
            "Indigo_nodS5",
            "Indigo_nodS6",
        ],
    }
}


def summarise_group(group_name: str, samples: List[str], salmon_dir: Path, output_dir: Path) -> Dict[str, float]:
    summary_out = output_dir / f"{group_name}_salmon_triplicate_tpm.tsv"
    merged_df = None
    length_lookup: Dict[str, float] = {}

    for sample in samples:
        quant_path = salmon_dir / sample / "quant.sf"
        if not quant_path.exists():
            raise FileNotFoundError(f"Cannot find Salmon quant file: {quant_path}")

        quant_df = pd.read_csv(
            quant_path,
            sep="\t",
            usecols=["Name", "Length", "NumReads", "TPM"],
        )

        length_lookup.update(dict(zip(quant_df["Name"], quant_df["Length"])))

        sample_df = quant_df[["Name", "NumReads", "TPM"]].rename(
            columns={
                "NumReads": f"{sample}_count",
                "TPM": f"{sample}_tpm",
            }
        )

        if merged_df is None:
            merged_df = sample_df
        else:
            merged_df = merged_df.merge(sample_df, on="Name", how="outer")


    merged_df = merged_df.fillna(0.0)
    count_cols = [col for col in merged_df.columns if col.endswith("_count")]
    tpm_cols = [col for col in merged_df.columns if col.endswith("_tpm")]

    merged_df[count_cols] = merged_df[count_cols].astype(float)
    merged_df[tpm_cols] = merged_df[tpm_cols].astype(float)

    merged_df["Average_count"] = merged_df[count_cols].mean(axis=1)
    merged_df["Average_tpm"] = merged_df[tpm_cols].mean(axis=1)

    ordered_cols = ["Name"] + count_cols + tpm_cols + ["Average_count", "Average_tpm"]
    merged_df[ordered_cols].rename(columns={"Name": "Geneid"}).to_csv(
        summary_out, sep="\t", index=False
    )

    return length_lookup


def write_quants_csv(config: Dict[str, List[str]], salmon_dir: Path, output_dir: Path) -> Path:
    quant_rows = []

    for group_name, samples in config.items():
        for sample in samples:
            quant_path = (salmon_dir / sample / "quant.sf").resolve()
            quant_rows.append(
                {
                    "Sample": sample,
                    "quant_file": str(quant_path),
                }
            )

    quant_rows.sort(key=lambda row: (0 if "Roots" in row["Sample"] else 1, row["Sample"]))

    quants_path = output_dir / "Quants.csv"
    pd.DataFrame(quant_rows).to_csv(quants_path, index=False)
    return quants_path


def write_tx2gene(length_lookup: Dict[str, float], output_dir: Path) -> Path:
    tx2gene_records = []
    for transcript_id in sorted(length_lookup.keys()):
        gene_id = transcript_id
        tx2gene_records.append((transcript_id, gene_id))

    tx2gene_path = output_dir / "geneID.csv"
    pd.DataFrame(tx2gene_records).to_csv(tx2gene_path, index=False, header=False)
    return tx2gene_path


def write_length_table(length_lookup: Dict[str, float], output_dir: Path, plant: str) -> Path:
    length_records = sorted(length_lookup.items())
    length_path = output_dir / f"{plant}_length.txt"
    pd.DataFrame(length_records, columns=["Geneid", "Length"]).to_csv(
        length_path, sep="\t", index=False
    )
    return length_path


def merge_salmon(salmon_dir: Path, output_dir: Path, plant: str) -> None:
    
    output_dir.mkdir(parents=True, exist_ok=True)

    plant_config = PLANT_CONFIG[plant]
    combined_lengths: Dict[str, float] = {}

    for group_name, samples in plant_config.items():
        group_lengths = summarise_group(group_name, samples, salmon_dir, output_dir)
        combined_lengths.update(group_lengths)

    write_quants_csv(plant_config, salmon_dir, output_dir)
    write_tx2gene(combined_lengths, output_dir)
    write_length_table(combined_lengths, output_dir, plant)


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge Salmon results for replicates")
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Folder containing per-sample Salmon output directories",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Destination folder for merged tables and metadata files",
    )
    parser.add_argument(
        "-p",
        "--plant",
        required=True,
        help="Plant identifier (choices: {})".format(", ".join(PLANT_CONFIG.keys())),
    )
    args = parser.parse_args()

    salmon_dir = Path(args.input).expanduser().resolve()
    output_dir = Path(args.output).expanduser().resolve()

    merge_salmon(salmon_dir, output_dir, args.plant.lower())


if __name__ == "__main__":
    main()
