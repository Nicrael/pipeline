# pipeline
Data reduction pipeline for OARPAF data and much more 

    curl -L https://raw.githubusercontent.com/pyenv/pyenv-installer/master/bin/pyenv-installer | bash

Add to `~/.bashrc`: 

    export PATH="/home/indy/.pyenv/bin:$PATH"                                                                                                                          
    eval "$(pyenv init -)"                                                                                                                                             
    eval "$(pyenv virtualenv-init -)"
   
Reload `~/.bashrc`:
   
    source ~/.bashrc

Install prerequisites (as suggested here: https://github.com/pyenv/pyenv/wiki/common-build-problems)

    sudo apt-get install -y make build-essential libssl-dev zlib1g-dev libbz2-dev \
    libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev libncursesw5-dev \
    xz-utils tk-dev libffi-dev liblzma-dev python-openssl git

Install specific Python version

    pyenv install 3.6.8
    pyenv global 3.6.8

Check if alias is ok (if not, reload `~/.basrc` or open new terminal):

    python --version 

Upgrade pip and ipyhton

    pip install --upgrade pip
    pip install --upgrade ipython
    
Install packages
    
    pip install --upgrade  astropy numpy scipy matplotlib pyds9 scikit-image astroplan
    

