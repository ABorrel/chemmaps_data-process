# update for production 12-15-2022

## update python -> 3.9
conda update python
conda -V

## update compdesc
pip uninstall  CompDesc
pip install -i https://test.pypi.org/simple/ CompDesc

## check rdkit version
python
import rdkit
rdkit.__version__

