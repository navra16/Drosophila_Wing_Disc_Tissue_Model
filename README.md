# Drosophila Wing Disc Eversion Tissue Model

A 3D, multi-layer, GPU-accelerated computational model of Drosophila wing disc eversion. Starting from a triangulated dome representing the wing pouch at wandering third instar (wL3), the simulation applies a region-specific active strain field to drive the morphological changes that occur over the wL3 ? 4 hours after puparium formation (4hAPF) window.

The mesh is composed of stacked layers (apical, body, basal) connected by linear, area, and bending springs, plus a per-prism volume-conservation term. The strain field formulation extends the active shape programming framework of Fuhrmann et al. (Science Advances, 2024) — which applied strain only to the apical surface of a thin shell — to a fully three-dimensional multi-layer mesh, with a separate dorsoventral (DV) boundary stripe receiving its own lambda values.

## What this code does

Given an input mesh and a schedule of lambda parameters defining the strain field, the simulation:

1. Initializes the mesh and computes a per-vertex orthonormal basis (e_R, e_phi, e_h) using the geometry described in Fuhrmann et al. 2024.
2. Constructs a spontaneous-strain tensor ? at each vertex using lambda values that differ between the inDV (DV stripe) and outDV regions.
3. Projects ? onto each edge to compute new rest lengths.
4. Relaxes the mesh to mechanical equilibrium under linear-spring, volume-conservation, and (optionally) bending and area forces.
5. Repeats over a 3-stage schedule (wL3 ? 0hAPF ? 2hAPF ? 4hAPF) with intermediate stage targets pulled from Fuhrmann's `input_lambda_df.csv`.
6. Pulls volume targets at each stage from biological measurements in Liu, O'Connell, Wall & Carthew (2024, *eLife* 12:RP91572).

VTK frames are written every substep so the trajectory can be visualized and analyzed.

## Code layout

| Area | Files |
|---|---|
| Simulation driver | `System.cu`, `System.h` |
| Strain field & basis vectors | `StrainTensor.cu`, `StrainTensor.h` |
| Springs & forces | `LinearSprings.cu/.h`, `VolumeComp.cu/.h`, `VolumeSprings.cu/.h` |
| Time stepping & relaxation | `NodeAdvance.cu/.h`, `gradientRelax.cu/.h` |
| I/O | `Storage.cpp/.h` (VTK output naming), `XMLParser.cu` |
| Build | `Makefile`, `SBATCH_try_this_one_if_the_other_does_not_work.sh`, `run.sh` |
| Initial mesh | `Data_Structures/Data_Structure.xml` (built via MATLAB scripts) |
| Post-processing | `postprocess_eversion.py` |

The general flow of simulation steps is implemented in `void System::solveSystem()` in `System.cu`. For specific force terms (linear spring, bending, area), see the corresponding `.cu` and `.h` files. Animation/data output names are configured in `Storage.cpp` and `Storage.h`. Job titles and time-step settings are in the SBATCH script. The initial data structure (mesh) lives in `Data_Structures/Data_Structure.xml`.

## Overall simulation flow

1. **Initialization** of global parameters and data structures (XML mesh load, layer flags, prism construction, edge classification).
2. **Initial relaxation**: a predetermined number of relaxation steps to bring the mesh to quasi-steady state with no applied strain.
3. **Main loop** — for each stage (wL3 ? 0hAPF, ? 2hAPF, ? 4hAPF):
   - Update lambda parameters for the current stage from the pre-computed `lam_per_stage` table.
   - Pre-compute the edge target rest lengths by applying that stage's lambdas to the wL3 reference geometry.
   - Within the stage, ramp `edge_rest_length` linearly from the previous stage's target to the current stage's target across `Nsteps` substeps.
   - Linearly ramp the volume-spring target from the start-of-stage to end-of-stage Liu/Carthew value within the same substeps.
   - At each substep, run an energy-minimization relaxation until forces converge.
   - Write a VTK frame.

## Running the simulation on UCR HPCC

1. Make sure you have every file and folder in this repository in your HPCC account (clone the repo: `git clone https://github.com/navra16/Drosophila_Wing_Disc_Eversion_Tissue_Model`).
2. `cd` into the folder containing `System.cu`, `Makefile`, and the SBATCH script.
3. `module load cuda; module load cmake`
4. `make` (or `make -j N` where N is 2-12, only if you are in an interactive GPU session — see UCR HPCC docs).
5. `sbatch -p gpu --gres=gpu:1 --time=01:00:00 SBATCH_try_this_one_if_the_original_does_not_work.sh`

## Running the simulation on Notre Dame CRC

1. Make sure you have every file and folder in this repository on your CRC account.
2. `cd` into the folder containing `System.cu`, `Makefile`, and `run.sh`.
3. `module load cuda cmake`
4. `make`
5. `qsub run.sh`

## Post-processing

Once a simulation has finished (or while it is still running), use the included Python script `postprocess_eversion.py` to generate three diagnostic figures:

1. **Cross-sections** of the apical and basal surfaces along and across the DV midline, directly comparable to Ecad:GFP imaging cross-sections.
2. **Curvature profiles** along each cross-section, computed via local quadratic fits, separately for apical and basal surfaces.
3. **Volume trajectory** showing simulated volume vs. the Liu/Carthew 2024 biological target schedule, with the wildtype pouch measurements overlaid as reference points.

### Setting up Python on the HPCC (one-time)

The system Python on UCR HPCC is in a read-only shared miniconda installation, so you need to install the three Python dependencies into your user space:

```bash
python3 -m pip install --user numpy matplotlib scipy
```

Verify with:

```bash
python3 -c "import numpy, matplotlib, scipy; print('all three imported OK')"
```

If `pip install` fails with permission errors, try `PYTHONUSERBASE=$HOME/.local python3 -m pip install --user numpy matplotlib scipy` instead.

### Running the post-processing script

The script auto-discovers VTK frames in a directory and parses the simulation log for volume data. From the directory containing the script:

```bash
python3 postprocess_eversion.py \
    --vtk-dir /path/to/your/vtk/output/directory/ \
    --log     /path/to/your/run.log \
    --num-layers 5 \
    --out-dir ./figures
```

Flags:
- `--vtk-dir` — directory containing the `_NNNNN.vtk` frame files (script picks up all matching files automatically and uses the latest by default)
- `--log` — captured stdout from the simulation run, containing the `V_target = ... | V_curr = ...` lines
- `--num-layers` — match what your mesh was built with (5 for the test mesh, 12 for the production res4 mesh)
- `--out-dir` — output directory for the three PNGs (created if it doesn't exist)
- `--frame-index N` — optional, plot a specific frame index instead of the last one
- `--vtk-file path.vtk` — optional, use a single named VTK file instead of `--vtk-dir`
- `--band-along F`, `--band-across F` — optional, in micrometers, controls the slice thickness for cross-sections (default 3.0)

If you only have the log (no VTK files yet, simulation still running), pass just `--log` and the script will produce the volume plot only. If you only have a VTK file (no log), pass just `--vtk-file` or `--vtk-dir` and you'll get cross-sections and curvature only.

To copy the figures back to your local machine for viewing:

```bash
# from your local terminal (not the cluster)
scp 'username@cluster.hpcc.ucr.edu:/path/to/figures/*.png' .
```

## References

Fuhrmann, J. F., Krishna, A., Paijmans, J., Duclut, C., Cwikla, G., Eaton, S., Popovic, M., Jülicher, F., Modes, C. D., & Dye, N. A. (2024). *Active shape programming drives Drosophila wing disc eversion.* Science Advances. <https://www.science.org/doi/10.1126/sciadv.adp0860>

Liu, Y., O'Connell, J. M., Wall, M. E., & Carthew, R. W. (2024). *Robust and rapid larva-to-pupa body shape transformation in Drosophila.* eLife 12:RP91572. <https://doi.org/10.7554/eLife.91572.3>

Sui, L., Alt, S., Weigert, M., Dye, N., Eaton, S., Jug, F., Myers, E. W., Jülicher, F., Salbreux, G., & Dahmann, C. (2018). *Differential lateral and basal tension drive folding of Drosophila wing discs through two distinct mechanisms.* Nature Communications 9:4620. <https://doi.org/10.1038/s41467-018-06497-3>

Tozluoglu, M., Maraspini, R., Sharma, M., Honigmann, A., & Mao, Y. (2020). *Epithelial organ shape is generated by patterned actomyosin contractility and maintained by the extracellular matrix.* bioRxiv. <https://doi.org/10.1101/2020.01.22.915272>

Nematbakhsh, A., Sun, W., Brodskiy, P. A., Amiri, A., Narciso, C., Xu, Z., Zartman, J. J., & Alber, M. (2024). *Balancing competing effects of tissue growth and cytoskeletal regulation during Drosophila wing disc development.* Nature Communications 15:2477. <https://doi.org/10.1038/s41467-024-46698-7>

## Acknowledgments

This implementation builds on previous work in the lab on tissue mechanics simulation. The strain field formulation is adapted from Fuhrmann et al. 2024, with extensions for multi-layer 3D meshes.