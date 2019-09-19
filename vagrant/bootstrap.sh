#!/usr/bin/env bash
export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get install -y python-pip python3-pip
pip install numpy biopython
pip3 install numpy biopython

# create Trinity output, DATA, scripts folder
sudo -S -u vagrant mkdir -p /vagrant/output
sudo -S -u vagrant ln -s /vagrant/output output

sudo -S -u vagrant mkdir -p /vagrant/DATA
sudo -S -u vagrant ln -s /vagrant/DATA DATA

sudo -S -u vagrant mkdir -p /vagrant/scripts
sudo -S -u vagrant ln -s /vagrant/scripts scripts

# Download Trinity singularity file
# echo "Download Trinity-v2.8.4.simg from github"
# sudo -S -u vagrant wget -q https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.8.4/Trinity-v2.8.4.simg

# Install beakerX
apt-get install -y openjdk-8-jdk maven
apt-get install -y nodejs npm
pip3 install beakerx ipywidgets pandas py4j jupyterlab
#pip3 install jupyterlab
beakerx install
#beakerx install --lab

# Install javascript kernel for jupyter
npm install -g ijavascript
/usr/local/bin/ijsinstall --install=global

# Move node dependencies in /usr/local/lib/node_modules/
mv /root/.npm/* /usr/local/lib/node_modules/
rm -f /root/.npm

#jupyter labextension install bqplot
#jupyter labextension install beakerx-jupyterlab
#jupyter labextension install @jupyter-widgets/jupyterlab-manager

# Install igv Jupyter Extension
pip3 install igv-jupyter

jupyter serverextension enable --py igv --system
jupyter nbextension install --py igv --system
jupyter nbextension enable --py igv --system

# plotly
pip3 install plotly==4.1.0

# tidyverse dependencies
apt-get install -y libssl-dev libcurl4-openssl-dev libxml2-dev pandoc

# Install R
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
apt-get update
apt-get install -y r-base

echo '
install.packages("tidyverse", repos="http://cran.us.r-project.org")
install.packages("IRkernel", repos="http://cran.us.r-project.org")
install.packages("plotly", repos="http://cran.us.r-project.org")
install.packages("devtools", repos="http://cran.us.r-project.org")
IRkernel::installspec(user = FALSE)
' | R -q --vanilla

# Change the ownership to vagrant for the convenience to install user packages
chown -R vagrant:vagrant /usr/local/lib

# Install beakerx examples
mkdir examples
cd examples/
git init
git remote add -f origin https://github.com/twosigma/beakerx.git
git config core.sparseCheckout true
echo "doc/" >> .git/info/sparse-checkout
git pull origin master
mv doc/* .
rmdir doc
git remote rm origin
cd ..
chown -R vagrant:vagrant examples
chown vagrant:vagrant .wget-hsts

echo "Provision task finished."
