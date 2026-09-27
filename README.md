# When does sequence order matter?

**An information-theoretic taxonomy of biological representation tasks and the limits of local encoding**

Zhengyi Chen<sup>1</sup> and Ruqing Chen<sup>2</sup>

<sup>1</sup> Guangxi Key Laboratory of Drug Discovery and Optimization, Guangxi Engineering Research Center for Pharmaceutical Molecular Screening and Druggability Evaluation, School of Pharmacy, Guilin Medical University, Guilin 541199, P. R. China (chenzhengyi@glmc.edu.cn)
<sup>2</sup> GUT Geoservice Inc., Montreal, Québec, Canada (ruqing@hotmail.com)

## About

An encoder is *local of order k* if it depends on a sequence only through its multiset of k-mers. Such an encoder is constant on the *spectral fibre* F<sub>k</sub>(S) of every sequence S, and the fibre entropy B<sub>k</sub>(S) = log<sub>2</sub> |F<sub>k</sub>(S)| is the number of bits of sequence identity it cannot see. This closing paper of an eight-part series turns that observation into an exact, proteome-wide accounting of sequence order information.

- **A three-term decomposition.** For every sequence of length L over an alphabet of size q and every order k, the specification cost L log<sub>2</sub> q splits exactly into the composition constraint L log<sub>2</sub> q − B<sub>1</sub>(S), the resolved local-syntax information B<sub>1</sub>(S) − B<sub>k</sub>(S), and the unresolved global-order entropy B<sub>k</sub>(S). Evaluated on the 20,431 reviewed human proteins for 15 orders from k = 1 to 50: composition fixes a median of 0.3947 bits per residue; tripeptides resolve 83.34% of the permutation entropy B<sub>1</sub>(S) on average and pentapeptides 99.75%, yet 99.64% of proteins still have B<sub>3</sub>(S) > 0 and 4.18% (854 proteins) have B<sub>10</sub>(S) > 0.
- **Six classes, three regimes.** Catalytic enzymes (n = 3,449), receptors, channels and transporters (n = 1,816) and other globular proteins (n = 14,153) have median spectral resolution κ = 5 and 95th percentile κ ≤ 8. The unresolved tail at k > 10 is concentrated in C2H2 zinc-finger regulators (n = 747, 49.13% with κ > 10), structural tandem-repeat proteins (n = 199, 39.70%) and hominoid segmental-duplication arrays (n = 67, 74.63%, 95th percentile κ = 257).
- **A decision rule for encoder selection.** The paper synthesises the seven companion studies (convolutional invariance theorems, ClinVar copy ambiguity, C2H2 zinc-finger reorderings, FibreBench, symmetry quotients) into four verifiable checks for choosing a sequence representation.

The fibre sizes come from the companion paper Z. Chen and R. Chen, *Separation without resolution: the spectral fibre geometry of the human proteome and the blind spots of local sequence encoders*, Zenodo (2026), [doi:10.5281/zenodo.22959211](https://doi.org/10.5281/zenodo.22959211) ([repository](https://github.com/Ruqing1963/separation-resolution)). The mathematical background is the monograph R. Chen, *Algebraic Genomics: The Combinatorics of Sequence Reconstruction and Representation*, Zenodo (2026), [doi:10.5281/zenodo.22943487](https://doi.org/10.5281/zenodo.22943487).

## Repository layout

```
order-taxonomy/
├── paper/      order_taxonomy.tex, order_taxonomy.pdf   the paper
├── code/       order_taxonomy.py                         Tables 1 and 2 and Figure 1
├── data/       fibre_results.json                        N_k for every k and κ, per protein
│                                                         (copied from separation-resolution)
├── results/    order_taxonomy_output.txt                 transcript of order_taxonomy.py
│               order_taxonomy_summary.json               the same numbers, machine-readable
├── figures/    fig_taxonomy.pdf                          Figure 1
├── README.md
└── LICENSE (code, MIT), LICENSE-CC-BY-4.0 (text, figure and derived data)
```

## Reproducing the results

Run the script from the `code` folder:

```
cd code
python order_taxonomy.py          # a few seconds; writes ../results/*.txt, ../results/*.json and ../figures/fig_taxonomy.pdf
```

and compile the paper from the `paper` folder:

```
cd paper
pdflatex order_taxonomy
pdflatex order_taxonomy
```

`order_taxonomy.py` needs Python 3.10 or later with `numpy` and `matplotlib`. It reads `data/fibre_results.json` (falling back to `../separation-resolution/results/fibre_results.json` if the companion repository is checked out alongside), assigns each protein to one of six classes from its gene symbol and UniProt family annotation, and evaluates the decomposition exactly. Every number in Tables 1 and 2 of the paper is taken from `results/order_taxonomy_output.txt`.

## How to cite

```bibtex
@misc{ChenChen2026OrderTaxonomy,
  author    = {Chen, Zhengyi and Chen, Ruqing},
  title     = {When does sequence order matter? An information-theoretic taxonomy of
               biological representation tasks and the limits of local encoding},
  year      = {2026},
  publisher = {Zenodo}
}
```

## License

Text, figure and derived data: [CC BY 4.0](LICENSE-CC-BY-4.0). Code: [MIT](LICENSE). The underlying UniProt data are distributed by UniProt under CC BY 4.0.
