# Thousands of taxa, a handful of states: Arctic protist communities follow a few environmental drivers.

**Draft submission package** with source code and data for this publication.

**Authors:** Alexey Ryabov<sup>* </sup>, Rahel Vortmeyer-Kley, Bernd Blasius, Thilo Gross, Rebecca McPherson, Sinhué Torres-Valdés, Katja Metfies<sup>*</sup>.

Please contact Alexey Ryabov (ICBM, Uni Oldenburg) or Katja Metfies (AWI) if you want to use the data.

The MATLAB source covers Fig. 2–3 and 5–7 and Fig. S1–S11. Fig. 1 and 4 are outside this package.

## Folder layout

- `source/`: figure scripts, project configuration, and required local helpers.
- `data/`: the input files read by the original `dt_loadData`, plus retained
  station coordinates. `source/dt_loadData.m` retains its filtering and
  normalization. The selected settings are `DataID = 'HG_20092023'` and
  `JobID = 'ASVNoMetazoaJ'`.
- `results/`: generated output location. It is excluded by `.gitignore`; a fresh clone starts with an empty `results/`.

The four retained data files are `HGSamples.mat`, `Samples_imp.xlsx`,
`Phyla_roles_and_sizes.xlsx`, and `Stations_coorsd.xlsx`. The first three
are read by `dt_loadData`; station coordinates are retained but are not
used by the current scripts.
`HGSamples.mat` and `Samples_imp.xlsx` are prepared upstream data, not raw
sequencing reads or unimputed environmental observations. Their upstream
preparation is outside the scope of this package.

## Requirements

- MATLAB R2025b was used to check this package. A static dependency scan also
  identified Statistics and Machine Learning Toolbox, Bioinformatics Toolbox,
  Financial Toolbox, and Parallel Computing Toolbox. The exact product set
  needed for each figure has not yet been verified by a clean run.
- [Diffusion map](https://www.mathworks.com/matlabcentral/fileexchange/180223-diffusion-map)
  is an external dependency. Install it separately; its functions are not
  included in this repository.
- [Violin Plots for Matlab](https://github.com/bastibe/Violinplot-Matlab)
  must be installed separately for `Fram_0273_Violins_CommComposition.m`.
  The installation must provide both `violinplot.m` and `Violin.m`.
- `source/lib/` contains the other functions found by the dependency scan.
  `bcdistcorr.m` is distributed with its original license in
  `source/lib/LICENSE-bcdistcorr.txt`.

Run scripts with `source/` as the MATLAB current folder. Initialization adds
`source/lib/` to the MATLAB path. `DT_Params.m` resolves the input files
relative to this package's `data/` folder. Running `proj_Initialize.m`
performs the original filtering and writes `tSamplesNorm.xlsx` and
`tStPredictors.xlsx` under `results/HG_20092023/ASVNoMetazoaJ/`.

## Figure scripts

| Figure | Script(s) |
| --- | --- |
| 2 | `Fram_0260_DMClusterPlot.m` |
| 3a–c | `Fram_0310_SVM_Type_Init_LonMonth.m`, `Fram_0322_SVM_Type_Train_LonMonth.m`, `Fram_0330_SVM_Type_PD_2D_LonMonth.m` |
| 3d–f | `Fram_0350_ML_Init_TempDepDoy.m`, `Fram_0352_ML_Train_TempDepDoy.m`, `Fram_0354_ML_Plot_TempDepDoy.m` |
| 3g, 5 | `Fram_0273_Violins_CommComposition.m` |
| 6 | `Fram_0296_CoreContribFinalClean.m` |
| 7 | `Fram_0347_Chl_LongTempIceSeasonUpper30m.m` |
| S1, S2 | `Fram_0100_DM.m` |
| S3 | `Fram_0105_DM_Cluster_Test.m` |
| S4a,b | `Fram_0310_SVM_Type_Init_LonMonth.m`, `Fram_0322_SVM_Type_Train_LonMonth.m` (`PredictorSetID = 'LonMonthDepth'`) |
| S4c,d | `Fram_0350_ML_Init_TempDepDoy.m`, `Fram_0352_ML_Train_TempDepDoy.m` |
| S5 | `Fram_0310_SVM_Type_Init_5drivers.m`, `Fram_0322_SVM_Type_Train_5drivers.m`, `Fram_0332_SVM_Type_PD_1D_5drivers.m` (`PredictorSetID = 'TempPARDepthNitrateIce'`) |
| S6 | `Fram_006_Env.m` |
| S7 | `Fram_0273_Violins_CommComposition.m` |
| S8 | `Fram_0250_Cluster.m` |
| S9 | `Fram_0275_CommCoreSpec.m` |
| S10, S11 | `Fram_0296_CoreContribFinalClean.m` |

The diffusion map from `Fram_0100_DM.m` is an input to clustering; the cluster output is an
input to several later scripts. `Fram_0275_CommCoreSpec.m` creates the core
species data used by `Fram_0296_CoreContribFinalClean.m`.

For a fresh clone, run `Fram_0100_DM.m` before `Fram_0250_Cluster.m`.
The latter creates the cluster labels needed by the later figure and model
scripts. Run `Fram_0275_CommCoreSpec.m` before
`Fram_0296_CoreContribFinalClean.m`.

Train each SVM before running its partial dependence plot script. Partial
dependence values are not included in the repository: they must be calculated
by running the plot scripts. All three partial dependence scripts set
`calcPlot = true`, so a fresh clone calculates the values rather than
loading absent cached files. The calculated values are saved under
`results/` for later reuse.

The LonMonth and 5drivers chains use separate initialization and training
scripts; each plot script loads the model trained with its own predictor set.
