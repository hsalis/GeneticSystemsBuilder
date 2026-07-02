# Genetic Systems Builder

Genetic Systems Builder (GSB) designs oligopools and PCR primers for massively parallel DNA assembly, enabling the construction of hundreds to thousands of complex plasmids with up 24-fold lower material costs. The algorithm facilitates the "build" step of an integrated, high-throughput design-build-test pipeline for engineering genetic systems. The pipeline was stress-tested on complex genetic systems (structural protein expression systems), including repetitive spider silks, biocements, reflectins, and talins.

GSB takes desired DNA assemblies as input, selects high-fidelity Type IIS overhangs already present in the target sequences, designs oligopool oligonucleotides with cut sites and primer-binding regions, designs well-specific PCR primers, and exports plate-ready spreadsheets for synthesis and assembly.

## Citation and Authorship

This software repository is accompanied by a publication:

**A Low-Cost, High-Throughput Design-Build-Test Pipeline for Engineering Genetic Systems: Stress Testing with Complex Structural Proteins**

Authors:

- Harry E. Adamson
- James R. McLellan
- Khushank Singhal
- Melik C. Demirel
- Howard M. Salis

Affiliations listed in the manuscript:

- Department of Chemical Engineering, The Pennsylvania State University, University Park, PA 16802
- Department of Engineering Science and Mechanics, The Pennsylvania State University, University Park, PA 16802
- Huck Institutes of Life Sciences and Materials Research Institute, The Pennsylvania State University, University Park, PA 16802
- Department of Chemical & Environmental Engineering, University of California Riverside, Riverside, CA 92521

Correspondence:

- Howard M. Salis, howard.salis@ucr.edu

Preprint:

- bioRxiv version: v1
- Posted: June 8, 2026
- DOI: [10.64898/2026.06.08.729977](https://doi.org/10.64898/2026.06.08.729977)
- Preprint URL: <https://www.biorxiv.org/content/10.64898/2026.06.08.729977v1>

Recommended citation:

> Adamson, H. E.; McLellan, J. R.; Singhal, K.; Demirel, M. C.; Salis, H. M. *A Low-Cost, High-Throughput Design-Build-Test Pipeline for Engineering Genetic Systems: Stress Testing with Complex Structural Proteins.* bioRxiv 2026.06.08.729977, version 1. doi: [10.64898/2026.06.08.729977](https://doi.org/10.64898/2026.06.08.729977).

## What the Genetic Systems Builder (GSB) Does for You

The Genetic Systems Builder is part of an integrated pipeline that combines computational genetic systems design, low-cost many-plasmid DNA assembly from oligopools, nanopore-sequencing-based validation, and label-free biosensor measurements. The pipeline was demonstrated by building 240 complex plasmids and lowering material costs by up to 24-fold compared with conventional commercial synthesis workflows.

The GSB will automatically:
- Identify optimal Golden Gate overhangs within inputted DNA sequences.
- Take into account low probability events (e.g. shifted cuts by Type IIS enzymes).
- Design oligopools and PCR primers for high-throughput, low-cost DNA assembly
- Export oligo lists and primer plate maps, ready for ordering

The GSB algorithm supports:

- Linear assemblies into an added vector fragment.
- Circular assemblies without an added vector fragment.
- One-pot combinatorial assemblies, where multiple DNA assemblies occur in the same well.
- Automated splitting into multi-level assemblies.
- Deterministic and multi-core parallel operation

### A step-by-step protocol is provided to carry out the PCR and Golden Gate assembly using a high-throughput plate-based format (suitable for either robots or multi-channel pipettors). 

## Installation

From this source checkout:

```bash
python3 -m pip install .
```

After the package is published on PyPI:

```bash
python3 -m pip install gsb-dna
```

After installation, the package provides a command-line executable:

```bash
gsb --help
```

`gsb` is an alias for running the main build script:

```bash
python3 GeneticSystemsCalculator_Build.py
```

Important dependency note: GSB uses the ViennaRNA Python binding, imported as `RNA`. If `gsb --help` fails with `ModuleNotFoundError: No module named 'RNA'`, install ViennaRNA into the same Python environment.

### Installing ViennaRNA

GSB imports ViennaRNA through its Python module name:

```python
import RNA
```

Install ViennaRNA in the same Python environment where you install and run GSB.

Using `pip`:

```bash
python3 -m pip install ViennaRNA
```

Using `conda`:

```bash
conda install -c bioconda viennarna
```

If you are using a virtual environment, activate it first:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install ViennaRNA
python3 -m pip install .
```

Verify that the binding is available:

```bash
python3 -c "import RNA; print(RNA.__version__ if hasattr(RNA, '__version__') else 'ViennaRNA import OK')"
```

Then verify GSB:

```bash
gsb --help
```

## Input Files

GSB accepts CSV, Excel, FASTA, GenBank, and JSON inputs. For CSV and Excel files, column names are normalized to lowercase.

Required columns:

- `name`: human-readable assembly name.
- `sequence` or `nucleotide_sequence`: DNA sequence containing only `A`, `C`, `G`, and `T`.

Optional columns:

- `is_circular`: `1` for circular assembly, `0` for linear assembly.
- `combinatorial_set`: grouping value that determines which DNA assemblies occur in the same well.
- `assembly_enzyme` or `enzyme`: row-specific Type IIS enzyme.
- `no_overhangs_position_range_list`: list of sequence coordinate ranges where overhangs should not be placed, for example `[[100, 180], [420, 480]]`.

The repository includes example files:

- `examples/structural_proteins.csv`
- `examples/vector_sequences.xlsx`

Both examples use the simple column layout:

```text
name,is_circular,sequence
```

## The `combinatorial_set` Column

`combinatorial_set` controls physical well assignment. Rows with the same `combinatorial_set` value are treated as belonging to the same combinatorial assembly group and are assigned to the same assembly well when possible.

Use the same value when variants should be assembled together in one reaction. Use different values when assemblies should be physically separated.

Example:

```csv
name,is_circular,combinatorial_set,sequence
PromoterA-GeneX-Terminator1,0,0,ACGT...
PromoterB-GeneX-Terminator1,0,0,ACGT...
PromoterC-GeneX-Terminator1,0,0,ACGT...
PromoterA-GeneY-Terminator2,0,1,ACGT...
PromoterB-GeneY-Terminator2,0,1,ACGT...
```

Rows in set `0` are assigned to one well. Rows in set `1` are assigned to another well. If `combinatorial_set` is omitted, each row is assigned its own default set based on row order.

For grouped sequences, GSB also compares sequences in the same set to identify variable regions. Those variable regions are avoided when selecting Golden Gate overhang positions, helping preserve compatibility across one-pot combinatorial variants.

## Running GSB

Basic command:

```bash
gsb INPUT_FILE -o OUTPUT_PREFIX
```

The output prefix determines the generated spreadsheet names.

### Example: Structural Protein Assemblies From CSV

```bash
mkdir -p outputs
gsb examples/structural_proteins.csv \
  --file-format csv \
  --assembly-enzyme BsaI_HFv2 \
  --level-two-assembly-enzyme BbsI_HF \
  -o outputs/structural_proteins_bsaI
```

### Example: Circular Vector Assemblies From Excel

```bash
mkdir -p outputs
gsb examples/vector_sequences.xlsx \
  --file-format excel \
  --sheet-name Sheet1 \
  --assembly-enzyme BsaI_HFv2 \
  --level-two-assembly-enzyme BbsI_HF \
  -o outputs/vector_sequences_bsaI
```

### Example: Parallel Build

```bash
mkdir -p outputs
gsb examples/vector_sequences.xlsx \
  --file-format excel \
  --sheet-name 0 \
  --maximum-oligo-length 350 \
  --parallel \
  --workers 8 \
  -o outputs/vector_sequences_parallel_350nt
```

Start with serial mode when debugging. Use `--parallel --workers N` after a serial build succeeds.

## Output Files

For `-o outputs/example_run`, GSB writes:

- `outputs/example_run.xlsx`: master workbook.
- `outputs/example_run_oligos.xlsx`: oligo order workbook.
- `outputs/example_run_oligos_unique.xlsx`: deduplicated oligo order workbook.
- `outputs/example_run_primer_plates_5p_1.xls`: 5' primer plate map.
- `outputs/example_run_primer_plates_3p_1.xls`: 3' primer plate map.

Additional primer plate-map files may be created for larger builds.

The master workbook includes sheets such as:

- `inputs`: normalized input records.
- `overview`: assembly names, enzymes, overhangs, well numbers, oligo counts, and predicted ligation fidelity.
- `oligos`: all designed oligos.
- `unique_oligos`: deduplicated oligo order list.
- `primers`: primer sequences and thermodynamic annotations.
- `assemblies`: fragment positions, overhangs, cut sites, padding, and final oligo sequences.
- `material costs`: estimated material costs.
- `screening results`: present when screening is enabled and results exist.

## Overview of the Assembly Protocol

The wet-lab workflow described in `GeneticSystemsBuilder_Protocol_v1_1.xlsx` has eight major steps:

1. Design oligopool and primer sequences using GSB.
2. Purchase the oligopool, primers, enzymes, master mix, ligase, buffers, and purification kits.
3. Prepare an oligopool template stock solution.
4. Amplify each well-specific DNA fragment library by PCR.
5. Purify PCR product in each well.
6. Carry out Golden Gate assembly.
7. Transform assembly product into chemically competent cells.
8. Purify assembled plasmids.

The manuscript emphasizes two practical points: keep PCR cycle number low enough to preserve uniform fragment amplification, and heat-soak PCR product before Golden Gate assembly to unfold intramolecular DNA structures.

## Step-by-Step Protocol Using the Excel File

Use the supplementary Excel protocol as the authoritative wet-lab worksheet. It contains editable input cells for concentrations, fragment counts, and volumes, along with per-well and per-plate calculations.

### Step 1. Design the Oligopool and Primers

Open `GeneticSystemsBuilder_Protocol_v1_1.xlsx` and begin with Step 1.

For each desired assembly, define:

1. DNA nucleotide sequence.
2. Type IIS assembly enzyme: BsaI, BbsI, BsmBI, Esp3I, or SapI.
3. Assembly mode:
   - Mode `0`: linear assembly into an added vector fragment.
   - Mode `1`: circular assembly without an added vector fragment.

Set global design specifications:

- Maximum oligo length: default `300 nt`.
- Primer target melting temperature: default `60 degC`.
- Primer length: default `18 nt`.
- Padding placement: default `5'` end.

Run GSB with matching values. For example:

```bash
gsb input_sequences.xlsx \
  --maximum-oligo-length 300 \
  --pcr-tm-target 60 \
  --primer-sequence-length 18 \
  --padding-end 5p \
  -o outputs/my_build
```

Use the exported oligo workbook to order the oligopool and the primer plate-map files to order or arrange primers.

### Step 2. Purchase Oligopool, Primers, and Reagents

The protocol notes example suppliers:

- Oligopools: Twist Bioscience, IDT, or GenScript.
- Primers: IDT or equivalent provider.
- Reagents: NEB or equivalent provider.
- DNA purification kits: Zymo or equivalent provider.

Example reagent list from the protocol:

- NEBNext Ultra II Q5 Master Mix (2X), NEB M0544X.
- Nuclease-free water, NEB B1500.
- T4 DNA Ligase Buffer, 10X, NEB B0202S.
- T4 DNA Ligase, NEB M0202M.
- BsaI-HFv2, NEB R3733.
- BbsI-HF, NEB R3539L.
- BsmBI-v2, NEB R0739L.
- Esp3I, NEB R0734L.
- SapI, NEB R0569L.
- ZR-96 DNA Clean and Concentrator-5, Zymo D4024.
- ZymoPURE 96 Plasmid Miniprep Kit, Zymo D4215.

The protocol assumes primer wet formulation in TE buffer. The example primer concentration is `10 uM`; the worksheet notes a usable range of `2.5 to 10 uM`.

### Step 3. Prepare the Oligopool Template Stock

Prepare a `0.5 ng/uL` oligopool template solution.

Example from the workbook:

- Dry oligopool template: `250 ng`.
- Nuclease-free water: `500 uL`.

Mix thoroughly and use this stock as the template in the PCR master mix.

### Step 4. Amplify the Fragment Library by PCR

The protocol uses a `25 uL` PCR reaction per well.

Prepare PCR master mix per well:

- NEBNext Ultra II Q5 Master Mix (2X): `12.5 uL`.
- Oligopool template solution: `1.0 uL`.
- Nuclease-free water: `9.0 uL`.

Aliquot `22.5 uL` PCR master mix into each PCR plate well.

Add well-specific primers:

- Forward primer: `1.25 uL`.
- Reverse primer: `1.25 uL`.

Thermocycler program:

1. Initial denaturation: `98 degC`, `120 s`.
2. Denaturation: `98 degC`, `20 s`.
3. Annealing: `61 degC`, `30 s`.
4. Elongation: `72 degC`, `10 s`.
5. Repeat denaturation, annealing, and elongation for `26 cycles`.
6. Final elongation: `72 degC`, `60 s`.
7. Hold: `4 degC`.

The worksheet notes that PCR cycles can be reduced when more oligopool template is added, for example `25 cycles` with `1 ng`, `24 cycles` with `2 ng`, or `23 cycles` with `4 ng`.

### Step 5. Purify DNA Fragment Libraries

Purify PCR products in each well.

Protocol steps:

1. Add `100 uL` DNA Binding Buffer to each well.
2. Transfer each well to a Zymo-Spin I-96 Plate.
3. Centrifuge for `10 min`; discard flow-through.
4. Add `300 uL` DNA Wash Buffer.
5. Centrifuge for `10 min`; discard flow-through.
6. Add another `300 uL` DNA Wash Buffer.
7. Centrifuge for `15 min`; discard flow-through.
8. Add `16 uL` DNA Elution Buffer.
9. Centrifuge for `5 min` into a PCR Product Elution Plate.
10. Measure DNA concentration for selected wells, for example by Nanodrop.

The workbook notes that DNA concentration is typically proportional to the number of fragments in the fragment library.

### Step 6. Carry Out Golden Gate Assembly

The protocol uses a `25 uL` assembly reaction.

First prepare the vector fragment or plasmid stock if using linear assembly into a vector. Example workbook inputs:

- Vector fragment/plasmid concentration: `30 ng/uL`.
- Vector fragment/plasmid length: `1500 bp`.

Before assembly, heat-soak PCR product:

1. Incubate PCR Product Elution Plate at `60 degC` for `30 min`.
2. Slowly cool to `16 degC` using a ramp slower than `1 degC/s`.

This heat-soak and slow-cool step is intended to unfold DNA structures and reduce undesired intramolecular structure formation.

Prepare assembly master mix per well:

- T4 DNA Ligase Buffer: `2.5 uL`.
- T4 DNA Ligase: `0.5 uL`.
- Assembly enzyme: `1.5 uL` for BsaI or BbsI to add 30 units.
- Vector fragment/plasmid: calculated by the worksheet when using a common vector.
- Nuclease-free water: bring master mix volume to `15 uL`.

For BsmBI, Esp3I, or SapI, the worksheet notes that `3 uL` enzyme is typically needed to add 30 units.

For circular assemblies, do not add vector. Add nuclease-free water to a final master mix volume of `15 uL`.

Aliquot `15 uL` assembly master mix to each well, then transfer PCR product:

- PCR product elution: default `10 uL` per well.
- The worksheet notes this can be adjusted from `8 to 16 uL` depending on fragment count.

Golden Gate thermocycler program:

1. Digestion: `37 degC`, `5 min`.
2. Ligation: `16 degC`, `5 min`.
3. Repeat digestion and ligation for `60 cycles`.
4. Hold: `4 degC`.

The worksheet notes that the number of cycles can be changed from `16` to `60`, with more cycles improving assembly efficiency with diminishing returns.

Optional: purify the assembly product by repeating the PCR-product purification workflow if increased transformation efficiency is needed.

### Step 7. Transform Assembly Product

Use chemically competent cells in plate format.

Protocol steps:

1. Heat-soak assembly product or purified assembly product at `60 degC` for `5 min`.
2. Add `10 uL` assembly product to a fresh plate containing chemically competent cells.
3. Mix gently.
4. Incubate on ice for `30 min`.
5. Heat shock at `42 degC` for `30 s`.
6. Add recovery medium.
7. Recover for `90 min` with shaking at `300+ RPM`, using `28 to 37 degC` as appropriate.
8. Transfer cultures to a 96-well deep-well plate containing `1 mL` LB medium with selective antibiotic per well.
9. Seal with an air-permeable cover.
10. Incubate shaking at `300+ RPM` for `18 to 24 h` at the desired growth temperature.

### Step 8. Purify Assembled Plasmids

The protocol points users to the manufacturer instructions for the ZymoPURE 96 Plasmid Miniprep Kit:

```text
https://www.zymoresearch.com/products/zymopure-96-plasmid-miniprep-kit
```

Use either centrifuge or vacuum manifold workflows according to the kit instructions and lab equipment.

## Practical Notes From the Manuscript

- The workflow is designed for 96-well plate operation and can be paired with multichannel pipetting or robotic liquid handlers.
- PCR product concentration tended to increase with the number of amplified fragments, though yields plateaued after roughly 20 fragments.
- More fragments generally reduced assembly efficiency; the assembly efficiency decreases by about 1.85% per DNA fragment.
- The material per-plasmid costs are lowest when many assemblies are carried out at the same using a single oligopool, called molecular economy of scale.

## Troubleshooting

No successful assemblies:

- Try a different `--assembly-enzyme`.
- Increase `--maximum-oligo-length`.
- Relax primer design tolerances.
- Check for forbidden or repeated sequence regions.
- Confirm sequences contain only `A`, `C`, `G`, and `T`.

Missing output directory:

```bash
mkdir -p outputs
```

Missing ViennaRNA:

```text
ModuleNotFoundError: No module named 'RNA'
```

Install ViennaRNA in the same Python environment used to run `gsb`; see [Installing ViennaRNA](#installing-viennarna).

Unexpected well grouping:

- Inspect the `combinatorial_set` column.
- Rows with the same value are intentionally assigned to the same assembly well.
- Remove the column or give each row a unique value for one assembly per well.

## Quick Test

```bash
python3 -m py_compile GeneticSystemsCalculator_Build.py
python3 GeneticSystemsCalculator_Build.py --help
gsb --help
```

Example design run:

```bash
mkdir -p outputs
gsb examples/structural_proteins.csv \
  --file-format csv \
  --assembly-enzyme BsaI_HFv2 \
  --level-two-assembly-enzyme BbsI_HF \
  -o outputs/structural_proteins_smoke
```
