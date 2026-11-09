# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

FEC (Finding ECCDNA) detects extrachromosomal circular DNA from long-read (PacBio HiFi / ONT) alignments. Input is a sorted, indexed BAM plus the reference FASTA (with `.fai`); output is a set of eccDNA cycle paths and their consensus sequences. It is a flat collection of Python scripts with no package, build system, test suite, or linter config.

## Running

```sh
# Full pipeline
./fclr_main.py -b <sorted_indexed.bam> -r <ref.fa> [-o fec_output] [-t 4] [-ms INT] [-d FLOAT]
./fclr_main.py -h          # all tuning flags (junction support, depth filtering, consensus options)

# Consensus step only, re-run on an existing summary_bigcycle.txt
./fclr_consensus.py -s fec_output/summary_bigcycle.txt -b <bam> -r <ref.fa> -o consensus_out [--ids 0,3]
```

There are no tests; verify changes by running the pipeline on a BAM. Run from the repo directory (or with it on `PATH`) because modules import each other by bare name. Python deps: pysam, pandas, numpy, scipy (`scipy.optimize.milp`, imported lazily in `fclr_cycle.py`), matplotlib (Agg backend). External binary: `bsalign` must be on `PATH` unless `--no_consensus` is passed — `fclr_main.py` checks this up front.

`generate_read/` holds a standalone eccDNA read simulator (`sim_ecread.py`) and a SLURM wrapper; it is not imported by the pipeline.

## Pipeline architecture

`fclr_main.py` is the orchestrator: it parses args and calls each stage in order, passing in-memory results forward. Modules are pulled in with `from fclr_x import *`, so function names share one namespace. Stages and their modules:

1. `fclr_minisup.py` — estimate genome-wide depth and min-support by sampling windows (skips centromere/telomere-adjacent windows).
2. `fclr_amplified.py` — scan for high-depth (amplified) regions with a multiprocessing pool; `merge_regions` merges overlaps. Also returns depth `gaps` used later.
3. `fclr_bpmq.py` — extract split-read breakpoints from SA tags.
   `fclr_refine.py` — 3a split regions at bridged depth dips that carry a junction; 3b extend region edges to nearby anchored breakpoints; 4 refine boundaries by clustering high-MQ breakpoints.
5. `fclr_junctiongraph.py` — build a graph whose nodes are stranded region fragments and whose edges are junctions supported by split reads.
6. `fclr_cycled.py` — DFS cycle enumeration (bounded by `--max_cycle_len`), cycle normalization, merging overlapping cycles into events, writing coordinate files.
7. `fclr_filterc.py` — drop cycles whose regions fail read-quality checks.
8. `fclr_cycle.py` (largest module) — `merge_small_cycles` combines small cycles into "big cycles": groups cycles by shared regions, splits groups by a depth cutoff (fold-gap), estimates per-region copy number from depth, solves appearance counts with an integer program (`milp`), reconstructs an Eulerian-style walk over junction uses, checks junction support, and emits low-depth variants and weak 2-fragment pairs to separate files. `fclr_plot.py` then draws each big cycle.
9. `fclr_consensus.py` — for each big cycle, cut the circular path into windows, gather read pieces from the BAM per window, run `bsalign poa` per window (threaded), and stitch into a consensus FASTA. Also runnable standalone (above); `add_consensus_arguments` shares its flags with `fclr_main.py`.

`fclr_login.py` is the shared utility module despite its name: logger setup (writes `<output_dir>/login_note.txt` and stderr), `tic`/`toc` timing, and the junction MQ rules `junction_mq_passes` / `partner_side_allowed`. Those two rules implement the asymmetric mapping-quality policy (one side MQ≥30, partner side only `--junction_partner_mq`, gated by `--partner_mq_min_region_len`) and are used by refine, graph, and big-cycle stages — keep them consistent when changing junction filtering.

## Data conventions

- Region nodes are strings like `chr14:40103961-40230096` with a strand suffix in cycle paths; `parse_node`/`parse_region`/`strandless_node` helpers are duplicated across `fclr_cycle.py`, `fclr_plot.py`, and `fclr_consensus.py` — a format change must be updated in all of them.
- Intermediate files in the output dir are prefixed by stage number and sample name (BAM basename), e.g. `03_<sample>_breakpoints.txt`, `04_<sample>_refined_regions.txt`, `06_<sample>_raw_cycles.txt`. Later stages sometimes re-read these (e.g. `fclr_cycle.py` locates the breakpoint table via glob), so renaming them breaks downstream lookup.
- Final outputs: `summary.txt`, `summary_bigcycle.txt` (+ `_variants`, `_weak_pairs`), `bigcycle_plots/`, `consensus/`, and `summary_bigcycle.fa`.

The README describes the output as VCF and the log as `login_not.txt`; both are out of date relative to the code.
