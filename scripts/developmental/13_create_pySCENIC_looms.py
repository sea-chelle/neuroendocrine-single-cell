from pathlib import Path

import loompy
import numpy as np
from scipy.io import mmread


# project directory: 05_analysis
PROJECT_DIR = Path(__file__).resolve().parents[2]

# directory containing the exported R inputs and loom outputs
BASE_DIR = PROJECT_DIR / "results" / "pySCENIC"

GENE_FILE = BASE_DIR / "pyscenic_genes.txt"

MATRICES = [
    {
        "name": "inference",
        "matrix": BASE_DIR / "pyscenic_inference_counts.mtx",
        "cells": BASE_DIR / "pyscenic_inference_cells.txt",
        "output": BASE_DIR / "pyscenic_inference.loom",
    },
    {
        "name": "scoring",
        "matrix": BASE_DIR / "pyscenic_scoring_counts.mtx",
        "cells": BASE_DIR / "pyscenic_scoring_cells.txt",
        "output": BASE_DIR / "pyscenic_scoring.loom",
    },
]

CHUNK_SIZE = 1000


def read_lines(path):
    """Read a text file while preserving line order."""
    with path.open("r", encoding="utf-8") as handle:
        return np.asarray(
            [line.rstrip("\n\r") for line in handle],
            dtype=str,
        )


def create_loom(matrix_path, cell_path, gene_path, output_path, name):
    """Convert a genes-by-cells Matrix Market file into a loom file."""

    print(f"\nLoading {name} matrix:")
    print(f"  {matrix_path}")

    genes = read_lines(gene_path)
    cells = read_lines(cell_path)

    # Matrix Market files were exported as genes × cells
    matrix = mmread(matrix_path).tocsc().astype(np.float32)

    expected_shape = (len(genes), len(cells))

    if matrix.shape != expected_shape:
        raise ValueError(
            f"{name} matrix shape is {matrix.shape}, "
            f"but gene/cell files imply {expected_shape}"
        )

    if len(np.unique(genes)) != len(genes):
        raise ValueError(f"Duplicated genes detected in {name} input")

    if len(np.unique(cells)) != len(cells):
        raise ValueError(f"Duplicated cells detected in {name} input")

    print(f"  Genes: {len(genes):,}")
    print(f"  Cells: {len(cells):,}")
    print(f"  Non-zero entries: {matrix.nnz:,}")

    # Remove an incomplete output from a previous attempt
    if output_path.exists():
        output_path.unlink()

    row_attributes = {
        "Gene": genes,
    }

    # Write manageable blocks of cells
    with loompy.new(str(output_path)) as loom:
        for start in range(0, matrix.shape[1], CHUNK_SIZE):
            end = min(
                start + CHUNK_SIZE,
                matrix.shape[1],
            )

            sparse_block = matrix[:, start:end]

            dense_block = sparse_block.toarray().astype(
                np.float32,
                copy=False,
            )

            column_attributes = {
                "CellID": cells[start:end],
                "nGene": np.asarray(
                    (sparse_block > 0).sum(axis=0)
                ).ravel(),
                "nUMI": np.asarray(
                    sparse_block.sum(axis=0)
                ).ravel(),
            }

            if start == 0:
                loom.add_columns(
                    dense_block,
                    column_attributes,
                    row_attrs=row_attributes,
                )
            else:
                loom.add_columns(
                    dense_block,
                    column_attributes,
                )

            print(
                f"  Written cells {start + 1:,}–{end:,} "
                f"of {matrix.shape[1]:,}"
            )

    # Verify the completed loom
    with loompy.connect(
        str(output_path),
        mode="r",
        validate=True,
    ) as loom:
        if loom.shape != expected_shape:
            raise ValueError(
                f"Completed {name} loom has shape {loom.shape}; "
                f"expected {expected_shape}"
            )

        if not np.array_equal(
            np.asarray(loom.ra["Gene"]).astype(str),
            genes,
        ):
            raise ValueError(
                f"Gene order changed in {name} loom"
            )

        if not np.array_equal(
            np.asarray(loom.ca["CellID"]).astype(str),
            cells,
        ):
            raise ValueError(
                f"Cell order changed in {name} loom"
            )

    print(f"Completed: {output_path}")
    print(
        f"File size: "
        f"{output_path.stat().st_size / (1024 ** 3):.2f} GB"
    )


def main():
    required_files = [GENE_FILE]

    for specification in MATRICES:
        required_files.extend(
            [
                specification["matrix"],
                specification["cells"],
            ]
        )

    missing_files = [
        str(path)
        for path in required_files
        if not path.exists()
    ]

    if missing_files:
        raise FileNotFoundError(
            "Missing required files:\n" +
            "\n".join(missing_files)
        )

    for specification in MATRICES:
        create_loom(
            matrix_path=specification["matrix"],
            cell_path=specification["cells"],
            gene_path=GENE_FILE,
            output_path=specification["output"],
            name=specification["name"],
        )

    print("\nBoth pySCENIC loom files were created successfully.")


if __name__ == "__main__":
    main()