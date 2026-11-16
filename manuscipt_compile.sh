cd ~/Desktop/DPS_PRL_PUBLICATION/two-loop-dps-prl

source .venv/bin/activate

make -C paper2/manuscript_prl pdf PY=python

make -C paper2/manuscript_prl words PY=python

# To generate arXiv tar 
# python paper2/manuscript_prl/make_arxiv.py --check
