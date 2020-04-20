# pipeline
Data reduction pipeline for OARPAF data and much more 

    curl -L https://raw.githubusercontent.com/pyenv/pyenv-installer/master/bin/pyenv-installer | bash

Add to `~/.basrc`: 

    export PATH="/home/indy/.pyenv/bin:$PATH"                                                                                                                          
    eval "$(pyenv init -)"                                                                                                                                             
    eval "$(pyenv virtualenv-init -)"
   
Reload `~/.basrc`:
   
    source ~/.basrc

Install specific Python version

    pyenv install 3.6.8
    pyenv global 3.6.8

Check if alias is ok:

    python --version 

Upgrade pip and ipyhton

    pip install --upgrade pip
    pip install --upgrade ipython
    
Install packages
    
    pip install --upgrade  astropy numpy scipy matplotlib pyds9 scikit-image astroplan
    

