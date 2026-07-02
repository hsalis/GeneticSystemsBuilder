# Genetic Systems Builder AI Usage Guide

This file documents how an AI coding agent should use `GeneticSystemsCalculator_Build.py` and the packaged `gsb` command to design Golden Gate oligo pools from sequence specification files.

## Purpose

`GeneticSystemsCalculator_Build.py` builds Golden Gate assembly oligo pools from input DNA sequences. It:

- Loads sequence specifications from CSV, Excel, FASTA, GenBank, or JSON.
- Splits target sequences into Golden Gate fragments.
- Chooses assembly overhangs using bundled ligation-frequency data.
- Designs PCR primers for each assembly well.
- Exports master build spreadsheets, oligo order sheets, unique oligo sheets, and primer plate maps.

The package exposes a console command:

```bash
gsb
```

`gsb ...` is the installable-package equivalent of:

```bash
python3 GeneticSystemsCalculator_Build.py ...
```

When working from a source checkout before installation, either command form is acceptable.

## Important Runtime Requirements

The code is pure Python, but it depends on scientific and bioinformatics packages listed in `pyproject.toml`, including `numpy`, `pandas`, `biopython`, `numba`, `scipy`, `networkx`, `rocksdict`, `openpyxl`, and `ViennaRNA`.

The most common runtime failure is:

```text
ModuleNotFoundError: No module named 'RNA'
```

That means the ViennaRNA Python binding is not installed or not visible in the active Python environment.

## Basic Command Shape

```bash
gsb INPUT_FILE -o OUTPUT_PREFIX [options]
```

Required:

- `INPUT_FILE`: CSV, Excel, FASTA, GenBank, or JSON sequence specification.
- `-o, --output-prefix`: prefix used for all generated output files.

Useful options:

- `--file-format auto|csv|excel|fasta|genbank|json`
- `--sheet-name SHEET_NAME_OR_INDEX`
- `--default-is-circular 0|1`
- `--assembly-enzyme BsaI_HFv2|BbsI_HF|BsmBI_v2|Esp3I|SapI`
- `--level-two-assembly-enzyme BbsI_HF`
- `--maximum-oligo-length 300`
- `--primer-sequence-length 18`
- `--padding-end 5p|3p`
- `--parallel --workers N`
- `--verbose`

Supported assembly enzymes are:

- `BsaI_HFv2`
- `BbsI_HF`
- `BsmBI_v2`
- `Esp3I`
- `SapI`

Common aliases such as `BsaI`, `BbsI`, `BsmBI`, and lowercase variants are normalized by the code.

## Input File Schema

For CSV and Excel input, column names are normalized to lowercase. Required columns:

- `name`
- `sequence` or `nucleotide_sequence`

Optional columns:

- `is_circular`: accepted values include `1`, `0`, `true`, `false`, `yes`, `no`.
- `combinatorial_set`: integer grouping value that determines which DNA assemblies take place in the same well. See the detailed section below.
- `assembly_enzyme` or `enzyme`: row-specific assembly enzyme.
- `no_overhangs_position_range_list`: Python/JSON-like list of `[begin, end]` ranges where overhangs should not be placed, for example `"[[100, 180], [420, 480]]"`.

Sequences must contain only `A`, `C`, `G`, and `T`. Whitespace is stripped and sequences are uppercased.

For FASTA and GenBank, record IDs become names and record sequences become DNA sequences. Use `--default-is-circular` and enzyme options to supply missing metadata.

For JSON input, the file must contain a list of objects with at least:

```json
[
  {
    "name": "example construct",
    "sequence": "ACGTACGTACGT",
    "is_circular": 0
  }
]
```

`nucleotide_sequence` can be used instead of `sequence`.

## Bundled Example Files

The repo contains two example input files:

- `examples/structural_proteins.csv`
- `examples/vector_sequences.xlsx`

`examples/structural_proteins.csv` has columns:

```text
name,is_circular,sequence
```

It contains linear protein-coding DNA constructs such as silk, resilin, keratin, aspein, and shematrin examples.

`examples/vector_sequences.xlsx` has one sheet, `Sheet1`, with columns:

```text
name,is_circular,sequence
```

It contains long circular vector/plasmid-like sequence designs.

## The `combinatorial_set` Column

`combinatorial_set` controls physical well assignment. It tells Genetic Systems Builder which DNA assemblies should be performed together in the same assembly well.

Use this column when multiple input sequences are alternative variants that belong to the same combinatorial library. For example, if three plasmid designs differ only by promoter or coding sequence choice and you want their assembly reactions to occur in the same physical well, give all three rows the same `combinatorial_set` value.

Example CSV fragment:

```csv
name,is_circular,combinatorial_set,sequence
PromoterA-GeneX-Terminator1,0,0,ACGT...
PromoterB-GeneX-Terminator1,0,0,ACGT...
PromoterC-GeneX-Terminator1,0,0,ACGT...
PromoterA-GeneY-Terminator2,0,1,ACGT...
PromoterB-GeneY-Terminator2,0,1,ACGT...
```

In this example:

- Rows with `combinatorial_set = 0` are assigned to the same well.
- Rows with `combinatorial_set = 1` are assigned to another well.
- The actual numeric values are labels; they do not need to be consecutive, but using small integers makes files easier to audit.

If `combinatorial_set` is omitted, the loader assigns each row its own default value based on row order. That means each sequence is treated as its own independent assembly well unless you explicitly provide shared values.

### What Happens Internally

During `GoldenGatePool.addMultipleAssemblies`, each sequence spec is normalized and its `combinatorial_set` is mapped to an internal compact integer. The original value is used as a group label, not as a literal plate position.

For regular one-level assemblies:

- The first successful assembly in a combinatorial set is assigned the next available well number.
- Later successful assemblies with the same combinatorial set reuse that same `well_number`.
- Exported overview/input sheets report that shared `well_number`.
- The internal `well_to_assembly_names` mapping records all assembly names assigned to that well.

For long assemblies:

- Long inputs are split into level-two fragments.
- Assemblies with the same combinatorial set share the same series of wells by long-fragment number.
- Fragment 0 for each variant maps to the first well for that combinatorial set, fragment 1 maps to the next corresponding well, and so on.

The builder also uses `combinatorial_set` to improve overhang placement for one-pot combinatorial assemblies. For sequences in the same set, it compares the grouped sequences and identifies variable or nonmatching regions. Those regions are added to `no_overhangs_position_range_list`, which prevents selecting Golden Gate overhangs where the sequence differs among variants. This helps keep shared assembly logic compatible across the variants in the same well.

### Practical Rules

Use the same `combinatorial_set` only when the rows are intended to be assembled together in one reaction/well.

Grouped sequences should generally represent variants of the same overall construct architecture. They should have compatible assembly boundaries, compatible primer design assumptions, and the same intended assembly enzyme unless the workflow has been deliberately designed otherwise.

Be cautious when grouping unrelated constructs. If unrelated sequences share a combinatorial set, the builder will place them in the same well, compare their sequence differences, and may produce confusing or infeasible overhang constraints.

For plate-map exports, primers are recorded per well. When multiple assemblies share a well, the primer plate map uses the first assembly assigned to that well as the representative primer pair. Therefore, combinatorial variants in the same well should be designed so a shared primer strategy is biologically and experimentally valid.

### When to Use Separate Values

Use different `combinatorial_set` values when:

- Assemblies should be carried out in separate wells.
- Constructs are unrelated.
- Variants require different outer primers.
- Variants require different reaction conditions.
- You want a one-to-one mapping of input sequence to assembly well.

### When to Reuse the Same Value

Reuse the same `combinatorial_set` value when:

- Rows are combinatorial variants meant to be pooled or assembled together.
- The variants share the same construct architecture.
- You want the outputs to show those assemblies as belonging to the same physical well.
- Variable regions should be considered when avoiding overhang placement.

## Example Uses

### 1. Show Help

From a source checkout:

```bash
python3 GeneticSystemsCalculator_Build.py --help
```

After package installation:

```bash
gsb --help
```

### 2. Build Oligos From the Structural Proteins CSV

```bash
gsb examples/structural_proteins.csv \
  --file-format csv \
  --assembly-enzyme BsaI_HFv2 \
  --level-two-assembly-enzyme BbsI_HF \
  -o outputs/structural_proteins_bsaI
```

This reads the CSV examples, treats each row's `is_circular` value as authoritative, uses BsaI-HFv2 for normal assembly, and uses BbsI-HF as the level-two enzyme for long assemblies.

Expected output files include:

- `outputs/structural_proteins_bsaI.xlsx`
- `outputs/structural_proteins_bsaI_oligos.xlsx`
- `outputs/structural_proteins_bsaI_oligos_unique.xlsx`
- `outputs/structural_proteins_bsaI_primer_plates_5p_1.xls`
- `outputs/structural_proteins_bsaI_primer_plates_3p_1.xls`

Additional primer plate map files may be produced if the build spans multiple plates.

### 3. Build Circular Vectors From the XLSX Example

```bash
gsb examples/vector_sequences.xlsx \
  --file-format excel \
  --sheet-name Sheet1 \
  --assembly-enzyme BsaI_HFv2 \
  --level-two-assembly-enzyme BbsI_HF \
  -o outputs/vector_sequences_bsaI
```

The workbook already includes `is_circular` values. If a workbook lacks that column, add:

```bash
--default-is-circular 1
```

### 4. Use an Alternate Assembly Enzyme

```bash
gsb examples/structural_proteins.csv \
  --file-format csv \
  --assembly-enzyme Esp3I \
  --level-two-assembly-enzyme BbsI_HF \
  -o outputs/structural_proteins_esp3i
```

Use this when BsaI sites or overhang constraints make BsaI-HFv2 undesirable.

### 5. Increase Oligo Length and Run in Parallel

```bash
gsb examples/vector_sequences.xlsx \
  --file-format excel \
  --sheet-name 0 \
  --maximum-oligo-length 350 \
  --parallel \
  --workers 8 \
  -o outputs/vector_sequences_parallel_350nt
```

`--parallel` uses Python multiprocessing. `--workers` must be a positive integer. Parallel mode can be faster for large inputs but may make debugging harder.

### 6. Tighten Primer Design Constraints

```bash
gsb examples/structural_proteins.csv \
  --file-format csv \
  --pcr-tm-target 62 \
  --primer-sequence-length 20 \
  --tolerance-delta-tm 1.0 \
  --tolerance-dg-folding -2.0 \
  --tolerance-dg-homodimer -8.0 \
  --tolerance-primer-dimer-folding-energy -5.0 \
  -o outputs/structural_proteins_primer_tuned
```

Use this if the default primer design settings fail for difficult sequences or if a user needs a different PCR temperature target.

### 7. Export Excel Files From a Previously Generated Build JSON

If another workflow has already produced a GSC build JSON containing `oligo_pool_specification`, regenerate only the Excel outputs:

```bash
gsb --export-build-json prior_build.json \
  -o outputs/reexported_build
```

This does not rerun assembly design. It only writes Excel artifacts.

## Output Workbook Contents

The master workbook `${OUTPUT_PREFIX}.xlsx` contains sheets including:

- `inputs`: normalized input records.
- `overview`: one row per successful assembly with enzyme, overhangs, oligo count, and predicted ligation fidelity.
- `oligos`: all designed oligos.
- `unique_oligos`: deduplicated oligo order list.
- `primers`: primer sequences and thermodynamic annotations.
- `assemblies`: fragment positions, overhangs, cut sites, padding, and oligo sequences.
- `material costs`: estimated cost summary.
- `screening results`: present only when screening is enabled and results exist.

The separate `${OUTPUT_PREFIX}_oligos.xlsx` and `${OUTPUT_PREFIX}_oligos_unique.xlsx` files are simplified oligo order sheets.

Primer plate maps are written as `${OUTPUT_PREFIX}_primer_plates_5p_N.xls` and `${OUTPUT_PREFIX}_primer_plates_3p_N.xls`, where `N` is the plate number.

## AI Agent Workflow Recommendations

1. Confirm dependencies before running a full build:

```bash
python3 GeneticSystemsCalculator_Build.py --help
```

2. Inspect input headers before invoking the build:

```bash
python3 - <<'PY'
import pandas as pd
print(pd.read_csv("examples/structural_proteins.csv").columns.tolist())
print(pd.read_excel("examples/vector_sequences.xlsx", sheet_name=0).columns.tolist())
PY
```

3. Use an output directory prefix, not just a bare filename, to keep generated spreadsheets organized:

```bash
mkdir -p outputs
gsb examples/structural_proteins.csv -o outputs/structural_proteins_demo
```

4. Start with serial mode unless the input is large or the user explicitly requests speed:

```bash
gsb examples/structural_proteins.csv -o outputs/debug_run --verbose
```

5. Use `--parallel --workers N` only after a serial run succeeds.

6. If a build reports no successful assemblies, retry with:

- A different `--assembly-enzyme`.
- A larger `--maximum-oligo-length`.
- Less strict primer tolerances.
- Adjusted long-assembly fragment thresholds for very long constructs.

## Common Failure Modes

### Missing Input Columns

Error:

```text
Missing required column: 'name'
Missing required column: 'sequence' or 'nucleotide_sequence'
```

Fix: rename spreadsheet columns to `name` and `sequence`, or provide JSON/FASTA/GenBank input.

### Invalid DNA Characters

Error:

```text
Sequence '...' contains non-ACGT characters
```

Fix: remove ambiguous bases such as `N`, `R`, `Y`, gaps, spaces, or annotations. This builder requires concrete DNA sequences.

### Unknown Enzyme

Error:

```text
Unsupported assembly enzyme
```

Fix: use one of `BsaI_HFv2`, `BbsI_HF`, `BsmBI_v2`, `Esp3I`, or `SapI`.

### ViennaRNA Missing

Error:

```text
ModuleNotFoundError: No module named 'RNA'
```

Fix: install ViennaRNA into the same Python environment used to run `gsb`.

### Output Directory Missing

If `-o outputs/run_name` is used and `outputs/` does not exist, spreadsheet writing may fail.

Fix:

```bash
mkdir -p outputs
```

## Programmatic Use

The main callable API is available from Python:

```python
from GeneticSystemsCalculator_Build import (
    GoldenGatePool,
    loadSequencesFromFile,
    build_arg_parser,
)
```

For most AI-agent workflows, prefer invoking the CLI through `gsb` or `python3 GeneticSystemsCalculator_Build.py`. The CLI performs input normalization, validation, build execution, finalization, and export in the expected order.

## Quick Smoke-Test Commands

Use these when validating a fresh checkout or package install:

```bash
python3 -m py_compile GeneticSystemsCalculator_Build.py
python3 GeneticSystemsCalculator_Build.py --help
gsb --help
```

Full example build:

```bash
mkdir -p outputs
gsb examples/structural_proteins.csv \
  --file-format csv \
  --assembly-enzyme BsaI_HFv2 \
  --level-two-assembly-enzyme BbsI_HF \
  -o outputs/structural_proteins_smoke
```
