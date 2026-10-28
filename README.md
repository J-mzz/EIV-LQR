# Errors-in-Variables Data-Driven Receding-Horizon LQR

MATLAB code for the paper

> J. Zheng and M. Sznaier, "Errors-in-variables data-driven receding-horizon LQR," *IEEE Control Systems Letters*, vol. 10, pp. 1483–1488, 2026.
> [doi:10.1109/LCSYS.2026.3709150](https://doi.org/10.1109/LCSYS.2026.3709150)

The method designs a state-feedback controller for the input-constrained LQR problem from noisy data. Both the state and the input measurements contain errors (errors in variables), and the system has process noise. The controller is robust over every linear model consistent with the data and the noise bounds. The theorem of alternatives eliminates the measurement-error variables, and sum-of-squares relaxations based on Scherer's or Putinar's Positivstellensatz give tractable semidefinite programs.

## Requirements

- MATLAB (tested with R2024a)
- [YALMIP](https://yalmip.github.io/)
- [MOSEK](https://www.mosek.com/) (free academic licenses are available); tested with MOSEK 9.1

## Usage

To draw the figures of the paper from the saved results in `data/`, run:

```matlab
plot_results('example1')   % Fig. 1
plot_results('example2')   % Fig. 2
```

To recompute the results, set `method` at the top of a script to `'proposed'` or `'unconstrained'` and run it:

```matlab
example1_linf          % Example 1 (Fig. 1)
example2_l2            % Example 2 (Fig. 2)
compare_relaxations    % comparison of the Scherer and Putinar relaxations
```

The scripts save their results in `results/`. After running both methods of an example, `plot_results('example1', 'results')` draws the figure from your runs.

| Script | Paper | Setting |
|---|---|---|
| `example1_linf.m` | Example 1, Fig. 1 | l_inf noise: eps_x = eps_u = 0.02, eps_w = 0.1; input constraint U = 0.5 I |
| `example2_l2.m` | Example 2, Fig. 2 | pointwise l_2 noise: eps_x = eps_u = 0.03, eps_w = 0.1; U = 0.3 I |
| `compare_relaxations.m` | Example 1, comparison of the relaxations | Scherer with r = 1 and Putinar with r = 2 on the same data (N = 19 transitions, only the state measurements are noisy), five receding-horizon steps |

Both examples use the open-loop unstable system of the paper and one noisy trajectory with N = 19 transitions. They compare the proposed input-constrained EIV controller (`'proposed'`) with the EIV controller without the input constraint (`'unconstrained'`).

## Repository layout

```
example1_linf.m, example2_l2.m   Examples 1 and 2
compare_relaxations.m            Scherer vs Putinar relaxation
plot_results.m                   Figs. 1 and 2 from saved results
data/                            results reported in the paper
figures/                         the figures of the paper
utils/
  EIVLQR_linf.m, psatz_linf.m    SDP and certificate for l_inf noise, (10)
  EIVLQR_l2.m, psatz_l2.m        SDP and certificate for pointwise l_2 noise, (11)
  EIVLQR_xnoise_scherer.m, EIVLQR_xnoise_putinar.m, psatz_xnoise_*.m
                                 Scherer and Putinar versions when only the
                                 state measurements are noisy (comparison)
  generate_sys.m, getl2noise.m   performance matrices and l_2 noise
  third_party/                   external code, see utils/third_party/README.md
```

## Citation

If you use this code, please cite:

```bibtex
@article{zheng2026errors,
  title={Errors-in-variables data-driven receding-horizon {LQR}},
  author={Zheng, Jian and Sznaier, Mario},
  journal={IEEE Control Systems Letters},
  volume={10},
  pages={1483--1488},
  year={2026},
  publisher={IEEE},
  doi={10.1109/LCSYS.2026.3709150}
}
```

## License

The code in this repository is released under the [MIT License](LICENSE), except the files in `utils/third_party/`, which keep their original licenses.
