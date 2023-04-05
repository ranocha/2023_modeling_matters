# Modeling still matters: a surprising instance of catastrophic floating point errors in mathematical biology and numerical methods for ODEs

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7801250.svg)](https://doi.org/10.5281/zenodo.7801250)

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
  month={04},
  howpublished={\url{https://github.com/ranocha/2023\_modeling\_matters}},
  doi={10.5281/zenodo.7801250}
}
```

## Abstract

We guide the reader on a journey through mathematical modeling and 
numerical analysis, emphasizing the crucial interplay of both disciplines.
Targeting undergraduate students with basic knowledge in dynamical
systems and numerical methods for ordinary differential equations,
we explore a model from mathematical biology where numerical methods
fail badly due to catastrophic floating point errors. We analyze the
reasons for this behavior by studying the steady states of the model
and use the theory of invariants to develop an alternative model that
is suited for numerical simulations. Our story intends 
to motivate combining analytical and numerical knowledge, even in 
cases where the world looks fine at first sight. We have set up an
online repository containing an interactive notebook with all
numerical experiments to make this study fully reproducible and
useful for classroom teaching.


## Numerical experiments

To reproduce the numerical experiments presented in this article, you need 
to install [Julia](https://julialang.org/). The numerical experiments presented 
in this article were performed using Julia v1.8.3.

First, you need to download this repository, e.g., by cloning it with `git`
or by downloading an archive via the GitHub interface.

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
