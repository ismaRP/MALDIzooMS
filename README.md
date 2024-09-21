
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MALDIzooMS

<!-- badges: start -->

<!-- badges: end -->

A package to manage large datasets of ZooMS MALDI data, change format,
collect from platemap coordinates. Preprocess the spectra and align
markers.

## Installation

Installing MALDIzooMS automatically installs all the packages it depends
on. However, MALDIZooMS depends on some Bioconductor packages. So before
installing it, we first need to tell R to also look on Bioconductor
(besided CRAN) repositories for the dependencies:

``` r
# Install Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.19")
# Tell R to also look for packages on Bioconductor repos
options(repos=BiocManager::repositories())
```

Install the development version from [GitHub](https://github.com/):

``` r
install.packages('devtools')
devtools::install_github("ismaRP/MALDIzooMS")
```

## Example

First we might want to collect all sample triplicates from different
plates. This will also change the names of the samples from map
coordinates to a sample name.

``` r
library(Spectra)
library(MsCoreUtils)
library(MALDIzooMS)


platemap = read_csv('path/to/platemap.csv')

platemap_1 = collect_triplicates(
  platemap,
  basepath = 'path/to/data', ext = "txt",
  outfolder = 'path/to/triplicates/folder',
  keep_dupl = T, keep_incomplete_tripl = T, dry_run=F)

platemap_1 = platemap_1 %>% filter(exists)

change_format_chunks(
  spectra_names = paste0(platemap_1$sample_name, "_", platemap_1$replicate),
  indir = 'path/to/triplicates/folder', in_fmt = 'txt', readf = "fread",
  outpath = 'path/to/formated/data', writef = "mzML",
  mc.cores = 8, nchunks = 20)
```

MALDIzooMS uses a lot of Spectra and MsCoreUtils packages functionality
for a lot of the preprocessing. They provide the basic infrastructure to
handle large datasets and also functions for preprocessing.

``` r

mzml_files = file.path(mzml_folder, denisova_metadata$file)
ncores = 6
sps_mzr = suppressMessages(
  Spectra(mzml_files, source = MsBackendMzR(), centroided = FALSE,
          BPPARAM = MulticoreParam(workers = ncores)))

transform_sqrt = function(x, ...) return(sqrt(x[,2]))
sps_mzr = addProcessing(sps_mzr, transform_sqrt)

smooth_sg_hws = 8
iterations = 50
snr = 2
k = 10
pickpeaks_hws = 20L
threshold = 0
# For deisotoping
tolerance = 0.4
ppm = 50

coefs = coefSG(smooth_sg_hws, k)
sps_mzr = addProcessing(sps_mzr, MsCoreUtils::smooth, coefs=coefs)
sps_mzr = addProcessing(
      sps_mzr, MALDIzooMS::baseline_correction, int_index = 'intensity',
      keep_bl = FALSE, substract_index = 'intensity', in_place = TRUE,
      method = 'SNIP', iterations = iterations, decreasing = TRUE)

sps_mzr = pickPeaks(sps_mzr, halfWindowSize=pickpeaks_hws, snr=snr, k=k,
                    descending=TRUE)
centroided(sps_mzr) = TRUE
```

Up until now we’ve just added preprocessing step to a processing queue
that needs to be applied. The next function will load data in chunks and
apply the preprocessing in parallel using `ncores` cores. We can also
tell if and where to write the resulting preprocessed and centroided
spectra.

After this point, we can have all the peaks data loaded together and
apply the rest of functions sequentially, instead of adding them to a
queue like before. Internally, we are using some useful MALDIquant
functions for aligning or filtering.

``` r
# We can write all spectra in separate files
processed_files = file.path('path/to/processed', platemap_1$sample_name)
sps_mzr = apply_preprocess(
  sps_mzr, write_data = TRUE,
  file=processed_files)
# Or in one big file
sps_mzr = apply_preprocess(
  sps_mzr, write_data = TRUE,
  file='path/to/preprocessed.mzML')

sps_df_deiso = deisotopeSpectra(sps_mzr, tolerance=tolerance, ppm=ppm)

sps_df_aligned = align_peaks(
  sps_df_deiso, tolerance = 0.005, minFreq = 0.9,
  labels=NULL)

sps_df_binned = bin_peaks(
  sps_df_aligned, method='strict', tolerance=0.02)
sps_df_binned = bin_peaks(
  sps_df_binned, method='relaxed', tolerance=0.02)

# First peaks that only appear in 1 of the replicates
sps_df_filtered = filter_peaks(
  sps_df_binned, minFreq=0.35, labels=platemap_1$sample_name)

sps_df_merged = combineSpectra(
  sps_df_filtered, f=platemap_1$sample_name)

int_mat = intensity_matrix(sps_df_merged, 'sample_name')

int_mat[is.na(int_mat)] = 0
```
