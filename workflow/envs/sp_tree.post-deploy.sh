# install astral-pro
git clone https://github.com/chaoszhang/ASTER $CONDA_PREFIX/lib/ASTER
cd $CONDA_PREFIX/lib/ASTER
make astral-pro
cp bin/astral-pro $CONDA_PREFIX/bin/
cd -

# install Ranger-DTL
wget https://compbio.engr.uconn.edu/wp-content/uploads/sites/2447/2019/08/Linux.zip -O $CONDA_PREFIX/lib/RangerDTL.zip
unzip $CONDA_PREFIX/lib/RangerDTL.zip -d $CONDA_PREFIX/lib/RangerDTL
chmod +x $CONDA_PREFIX/lib/RangerDTL/Linux/CorePrograms/*
mv $CONDA_PREFIX/lib/RangerDTL/Linux/CorePrograms/* $CONDA_PREFIX/bin/

