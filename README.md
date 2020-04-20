# pipeline
Data reduction pipeline for OARPAF data and much more 

curl -L https://raw.githubusercontent.com/pyenv/pyenv-installer/master/bin/pyenv-installer | bash
pyenv install 3.7

Add to ~/.basrc: 

export PATH="/home/indy/.pyenv/bin:$PATH"                                                                                                                          
eval "$(pyenv init -)"                                                                                                                                             
eval "$(pyenv virtualenv-init -)"                                                                                                                                  

source ~/.basrc

pyenv install 3.7.7
pyenv global 3.7.7

python --version #check if 3.7.7

python -m pip install --upgrade pip
python -m pip install --upgrade ipython
pip install --upgrade  astropy numpy scipy matplotlib pyds9 scikit-image astroplan
