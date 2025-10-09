#!/usr/bin/env python3
"""
Update the fold_tree Snakemake configuration file
Usage:
    python update_foldtree_config.py path/to/config_vars.yaml
"""

import sys
from pathlib import Path


def update_lines(text: str, replacements: dict[str, str]) -> str:
    lines = text.splitlines()
    updated = []
    for line in lines:
        stripped = line.strip()
        replaced = False
        for key, value in replacements.items():
            prefix = f"{key}:"
            if stripped.startswith(prefix):
                indent = line[: line.index(prefix)]
                updated.append(f"{indent}{prefix} {value}")
                replaced = True
                break
        if not replaced:
            updated.append(line)
    return "\n".join(updated) + "\n"


def main() -> None:
    
    config_path = Path(sys.argv[1])

    content = config_path.read_text()

    #Add and update project-specific defaults (could be also other parameters)
    overrides = {
        "filter": "False",
        "foldseek_path": "foldseek",
        "foldseek_cores": "4",
        "custom_structs": "True",
        "clean_folder": "False",
    }
    updated = update_lines(content, overrides)
    config_path.write_text(updated)


if __name__ == "__main__":
    main()
