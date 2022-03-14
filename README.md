# mixlasso

TBA...


# Troubleshooting

## Using all cores on HPC clusters

When running this package without parallelization, it might occupy all cores on HPC clusters. This is because an R dependency `mlegp` with setting BLAS/LAPACK uses all available cores by default. [One solution](https://stackoverflow.com/questions/30791550/limit-number-of-threads-in-numpy) is to set BLAS/OpenMP environment variables to control allocation:

```bash
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1
```
