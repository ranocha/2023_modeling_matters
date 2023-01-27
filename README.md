# Modeling still matters: a surprising instance of catastrophic floating point errors in mathematical biology and numerical methods for ODEs

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/TODO.svg)](https://doi.org/TODO)

This repository contains information and code to reproduce the results presented in the
article
```bibtex
@online{reisch2023modeling,
  title={Modeling still matters: a surprising instance of catastrophic floating 
         point errors in mathematical biology and numerical methods for {ODEs}},
  author={Reisch, Cordula and Ranocha, Hendrik},
  year={2023},
  month={TODO},
  eprint={TODO},
  eprinttype={arxiv},
  eprintclass={TODO}
}
```

If you find these results useful, please cite the article mentioned above. If you
use the implementations provided here, please **also** cite this repository as
```bibtex
@misc{reisch2023modelingRepro,
  title={Reproducibility repository for
         {M}odeling still matters: a surprising instance of catastrophic floating 
         point errors in mathematical biology and numerical methods for {ODEs}},
  author={Reisch, Cordula and Ranocha, Hendrik},
  year={2023},
  month={TODO},
  howpublished={\url{https://github.com/ranocha/2022\_modeling\_matters}},
  doi={TODO}
}
```

## Abstract

TODO


## Numerical experiments

To reproduce the numerical experiments presented in this article, you need 
to install [Julia](https://julialang.org/). The numerical experiments presented 
in this article were performed using Julia v1.8.3.

Then, you need to start Julia in this directory and execute the following commands
in the REPL

```julia
import Pkg; Pkg.activate("."); Pkg.instantiate(); import Pluto; Pluto.run()
```

You can combine this with starting Julia from the command line as follows:

```bash
julia -e 'import Pkg; Pkg.activate("."); Pkg.instantiate(); import Pluto; Pluto.run()'
```

Then, the web server of [Pluto.jl](https://github.com/fonsp/Pluto.jl) should start
and open a browser window for you. There, you need to select the file `notebook.jl`.
Then, Pluto.jl should load the file and you should start to see some text. The
setup will take up to several minutes since Julia needs to install all dependencies
and execute the code. When everything is finished, the notebook is ready for interactive
use and exploration. In particular, you will be able to change several parameters
interactively and Julia will make sure that everything is sufficiently fast to be
updated on the fly.

A static preview is also available [online](https://ranocha.de/2022_modeling_matters/).


## Authors

- [Cordula Reisch](https://www.tu-braunschweig.de/ipde/personal/creisch) (TU Braunschweig, Germany)
- [Hendrik Ranocha](https://ranocha.de) (University of Hamburg, Germany)


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
