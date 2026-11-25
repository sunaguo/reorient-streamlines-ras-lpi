<h1 align="center">fix RAS/LPI labeling</h1>

<p align="center">
  Fixes inconsistent RAS/LPI endpoint labels of left and right MDLFang, MDLFspl, Uncinate and Aslant ROIs by comparing ROI centroids along each tract's main axis and swapping file names when needed.
</p>

<p align="center">
  <a href="https://github.com/brainlife/abcd-spec"><img src="https://img.shields.io/badge/ABCD_Spec-v1.1-green.svg" alt="ABCD Spec v1.1"></a>
  <a href="https://doi.org/10.25663/brainlife.app.795"><img src="https://img.shields.io/badge/Brainlife-brainlife.app.795-blue.svg" alt="Run on brainlife.io"></a>
</p>

## At a glance

<!-- brainlife-app:begin (generated from the App manifest; edit the manifest, not this block) -->

> **Fixes inconsistent RAS/LPI endpoint labels of left and right MDLFang, MDLFspl, Uncinate and Aslant ROIs by comparing ROI centroids along each tract's main axis and swapping file names when needed.**

For each listed tract, loads the left/right RAS and LPI endpoint ROIs, checks that left and right are in consistent coordinates, and swaps the RAS and LPI file names whenever the RAS ROI is not on the expected side along the tract's main orientation (z for MDLF and Aslant, y for Uncinate). Use it to make endpoint labels consistent across subjects before tract profile or endpoint analyses. It needs an ROIs dataset containing left/right RAS and LPI ROIs for all four tracts.

| | |
|---|---|
| **Inputs** | `neuro/rois` |
| **Outputs** | `neuro/rois` |
| **Container** | `brainlife/pyafq:dev` |
| **Resources** | 4 CPUs · 28 GB memory · 10 h walltime |
| **Run it** | [doi.org/10.25663/brainlife.app.795](https://doi.org/10.25663/brainlife.app.795) |

#### Inputs

| Input | Datatype | Required | Description |
|---|---|---|---|
| `rois` | `neuro/rois` | yes | ROIs directory with left and right RAS and LPI endpoint ROIs per tract. Must contain left/right *RAS* and *LPI* NIfTI ROIs for MDLFang, MDLFspl, Uncinate and Aslant. |

#### Outputs

| Output | Datatype | Description |
|---|---|---|
| `rois_relabeled` | `neuro/rois` | Same ROIs with the RAS/LPI labels corrected where they were swapped. |

<!-- brainlife-app:end -->

## Running the App

### On brainlife.io

Open the App at [doi.org/10.25663/brainlife.app.795](https://doi.org/10.25663/brainlife.app.795), choose `neuro/rois` from your project and press **Execute**. The results are saved to your project as `neuro/rois`.

### Locally

You need [Singularity/Apptainer](https://apptainer.org/) (the App runs in `brainlife/pyafq:dev`), [`jq`](https://jqlang.github.io/jq/).

1. Clone this repository.
2. Create `config.json` with the paths to your inputs and the parameters:

   ```json
   {
       "rois": "/path/to/rois/rois"
   }
   ```

3. Run `./main`.

## Troubleshooting

- **something failed. checked logs**: Check that all eight ROIs per tract exist in the input.

## Citing

If you use this App, please cite it, the software it builds on and brainlife.io:

- fix RAS/LPI labeling. brainlife.io App. [doi:10.25663/brainlife.app.795](https://doi.org/10.25663/brainlife.app.795)
- Kruper et al. Evaluating the Reliability of Human Brain White Matter Tractometry. *Aperture Neuro* 2021 (2021). [doi:10.52294/e6198273-b8e3-4b63-babb-6e6b0da10669](https://doi.org/10.52294/e6198273-b8e3-4b63-babb-6e6b0da10669)
- Hayashi S, Caron BA, Heinsfeld AS, et al. brainlife.io: a decentralized and open-source cloud platform to support neuroscience research. *Nature Methods* 21:809–813 (2024). [doi:10.1038/s41592-024-02237-2](https://doi.org/10.1038/s41592-024-02237-2)
- Avesani P, McPherson B, Hayashi S, et al. The open diffusion data derivatives, brain data upcycling via integrated publishing of derivatives and reproducible open cloud services. *Scientific Data* 6:69 (2019). [doi:10.1038/s41597-019-0073-y](https://doi.org/10.1038/s41597-019-0073-y)

## Authors

Maintained by [@sunaguo](https://github.com/sunaguo).
