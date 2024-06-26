# install astral-pro
git clone https://github.com/chaoszhang/ASTER $CONDA_PREFIX/lib/ASTER
cd $CONDA_PREFIX/lib/ASTER
make astral-pro
cp bin/astral-pro $CONDA_PREFIX/bin/
cd -

# install Ranger-DTL
# wget https://compbio.mit.edu/ranger-dtl/ranger-dtl-linux.tar.gz -O $CONDA_PREFIX/lib/RangerDTL.tar.gz
# tar xvzf $CONDA_PREFIX/lib/RangerDTL.tar.gz -C $CONDA_PREFIX/bin/ #--strip-components=2 ./bin/foldmason

# 404 error. let's wait and see if they fix it, for now use one already installed.
