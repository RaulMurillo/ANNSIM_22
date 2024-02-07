# ANNSIM_22

Code for reproducing ANNSIM 2022 results.  
For more information, refer to the original paper: [The Effects of Numerical Precision In Scientific Applications](https://doi.org/10.23919/ANNSIM55834.2022.9859379)

## Quick start
Sources in `physic_probs/` make use of [Universal Numbers Library](https://github.com/stillwater-sc/universal), while the models in `DNN/` are implemented on [Qtorch+](https://github.com/minhhn2910/QPyTorch).
You would need to install those libraries to use the sources in this repo.

Universal library is necessary to compile the projects.  
Use the script `update_cmake.sh` with option `-p <universal directory>` for quick setup.

## Citation 
If you find this repo useful, please cite our [paper](https://doi.org/10.23919/ANNSIM55834.2022.9859379) listed below.

```bib
@inproceedings{Murillo2022Effects,
  author={Murillo, Raul and Barrio, Alberto A. Del and Botella, Guillermo},
  booktitle={2022 Annual Modeling and Simulation Conference (ANNSIM)}, 
  title={The Effects of Numerical Precision In Scientific Applications}, 
  year={2022},
  volume={},
  number={},
  pages={152--163},
  doi={10.23919/ANNSIM55834.2022.9859379}
}
```

## Credits
Some of the original sources are taken from 
https://people.sc.fsu.edu/~jburkardt/
and 
https://people.math.sc.edu/Burkardt/
