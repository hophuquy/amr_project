# Set environment variables before running any other code.
import os
os.environ['TF_FORCE_UNIFIED_MEMORY'] = '1'
os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '4.0'

#@title 1. Install third-party software

#@markdown Please execute this cell by pressing the _Play_ button
#@markdown on the left to download and import third-party software
#@markdown in this Colab notebook. (See the [acknowledgements](https://docs.anaconda.com/miniconda/) in our readme.)

#@markdown **Note**: This installs the software on the Colab
#@markdown notebook in the cloud and not on your computer.

from IPython.utils import io
import os
import subprocess
import tqdm.notebook

TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'

try:
  with tqdm.notebook.tqdm(total=100, bar_format=TQDM_BAR_FORMAT) as pbar:
    with io.capture_output() as captured:

      %shell pip install condacolab
      import condacolab
      condacolab.install()
      pbar.update(6)

      %shell apt-get update -qq && apt-get upgrade -qq -y
      %shell apt-get install -qq build-essential gawk
      %shell apt autoremove -qq && apt autoclean -qq
      %shell conda update --all --quiet
      pbar.update(6)

      %shell mamba install bioconda::kraken2
      pbar.update(6)
      %shell mamba install bioconda::bracken
      pbar.update(6)
      # mamba install conda-forge::aria2
      # mamba install bioconda::entrez-direct
      %shell mamba install bioconda::bioconductor-ancombc
      %shell mamba install bioconda::bioconductor-enhancedvolcano
      %shell mamba install bioconda::samtools
      pbar.update(6)
      %shell mamba install bioconda::bowtie2
      pbar.update(6)
      %shell mamba install bioconda::blast
      pbar.update(6)
      # mamba install bioconda::sra-tools
      %shell mamba install bioconda::prodigal
      pbar.update(6)
      %shell mamba install bioconda::ncbi-amrfinderplus
      pbar.update(6)
      %shell mamba install bioconda::fastp
      pbar.update(6)
      %shell mamba install bioconda::fusioncatcher-seqtk
      pbar.update(6)
      %shell mamba install bioconda::bbmap
      pbar.update(6)
      %shell mamba install bioconda::subread
      pbar.update(6)
      %shell pip install bioinfokit
      pbar.update(6)
      %shell pip install scipy
      pbar.update(6)
      %shell pip install ete3 PyQt5
      %shell pip install scikit-posthocs
      pbar.update(6)
      %shell pip install pandas numpy scipy scikit-bio scikit-learn
      pbar.update(2)
      %shell pip install biopython
      %shell pip install kraken-biom
      pbar.update(2)


except subprocess.CalledProcessError:
  print(captured)
  raise
