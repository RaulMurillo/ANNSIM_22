# ANNSIM_22

Code for reproducing ANNSIM 2022 results.  
For more information, refer to the original paper: https://scs.org/wp-content/uploads/2022/07/39_Paper_THE-EFFECTS-OF-NUMERICAL-PRECISION-IN-SCIENTIFIC-APPLICATIONS.pdf

## Quick start
Sources in `physic_probs/` make use of [Universal Numbers Library](https://github.com/stillwater-sc/universal), while the models in `DNN/` are implemented on [Qtorch+](https://github.com/minhhn2910/QPyTorch).
You would need to install those libraries to use the sources in this repo.

Universal library is necessary to compile the projects.  
Use the script `update_cmake.sh` with option `-p <universal directory>` for quick setup.

## Citation 
If you find this repo useful, please cite our [paper](https://scs.org/wp-content/uploads/2022/07/39_Paper_THE-EFFECTS-OF-NUMERICAL-PRECISION-IN-SCIENTIFIC-APPLICATIONS.pdf) listed below.

```bib
@inproceedings{Murillo2022Effects,
    title   = {The Effects of Numerical Precision in Scientific Applications},
    author  = {Murillo, Raul and {Del Barrio}, Alberto A. and Botella, Guillermo},
    booktitle={2022 Annual Modeling and Simulation Conference (ANNSIM)},
    pages   = {},
    year    = {2022},
}
```

## Credits
Some of the original sources are taken from 
https://people.sc.fsu.edu/~jburkardt/
and 
https://people.math.sc.edu/Burkardt/