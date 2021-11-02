# Understanding X-ray Absorption Spectra by Means of Descriptors and Machine Learning Algorithms

### A. A. Guda, S. A. Guda, A. Martini, A.N. Kravtsova, A. Algasov, A. Bugaev, S. P. Kubrin, L. V. Guda, P. Šot, J. A. van Bokhoven, C. Copéret, A. V. Soldatov

The notebook contains the code for using machine learning algorithms and establishing the relationship between intuitive descriptors of spectra, such as edge position, intensities, positions and curvatures of minima and maxima on the one side, and those of the local atomic and electronic structure which are the coordination numbers, bond distances and angles, and oxidation state on the other. This approach allows overcoming the problem of the systematic difference between theoretical and experimental spectra. Furthermore, the numerical relations can be expressed in analytical formulas providing a simple and fast tool to extract structural parameters based on spectral shape. The methodology was successfully applied to experimental data of the multicomponent Fe:SiO2 system and reference iron compounds, demonstrating the high prediction quality for both theoretical validation sets and experimental data.

The code is based on [PyFitIt project](https://github.com/gudasergey/pyFitIt)

### Repository contains:
- exp subfolder with experimental spectra
- generated subfolder with data generated during calculations and used afterwards
- results subfolder
- samples subfolder with theoretical spectra database
- xyz subfolder with .xyz files used in molecule constructor and samples generation
- exp_true_values.csv file with known parameters for experimental spectra database
- instructions.docx file with install instructions
- paper_calculations.ipynb - this notebook
- Fe_project.py - project file with different settings

### Install instructions

Anaconda with Python 3.7.6 is needed for running the code.

Install additional packages with the command:

```pip install jupytext ipykernel ipywidgets tqdm scipy numba cycler statsmodels lmfit matplotlib numpy pandas parsy notebook nbformat pulp scikit_learn seaborn```

Open file paper_calculations.ipynb in Jupyter notebook and run it
