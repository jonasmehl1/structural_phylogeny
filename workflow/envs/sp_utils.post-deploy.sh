# install foldmason
# find other binaries here: https://mmseqs.com/foldmason/
mkdir -p $CONDA_PREFIX/lib/
wget https://mmseqs.com/foldmason/foldmason-linux-avx2.tar.gz -O $CONDA_PREFIX/lib/foldmason.tar.gz

tar xvzf $CONDA_PREFIX/lib/foldmason.tar.gz -C $CONDA_PREFIX/bin/ --strip-components=2 ./bin/foldmason