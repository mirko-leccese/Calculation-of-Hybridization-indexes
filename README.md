# Calculation-of-Hybridization-indexes

## Description
This is a simple Python code that allows computing the hybridization indexes of C atoms of a $\pi$ conjugated molecule. Namely, for any C atom with hybrid orbitals:

$$| \psi_i \rangle = |s\rangle + \sqrt{\tau_i} |p_\rangle$$

it computes the hybridization indexes $\tau_i$ (and corresponding $s$ weights, $s_i=1/(1+\tau_i)$ ), which are entirely determined by the geometrical structure in the hypothesis of non-bent bonds.  Details of the calculation and of the code implementation are given in the Jupyter Notebook ```hybrids_notebook.ipynb```. The latter also contains a brief overview on the fractional hybridization theory. 



