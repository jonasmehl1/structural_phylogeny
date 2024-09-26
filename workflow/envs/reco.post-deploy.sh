# install Notung
wget https://amberjack.compbio.cs.cmu.edu/Notung/Notung-2.9.1.5.zip -O $CONDA_PREFIX/lib/Notung-2.9.1.5.zip
unzip $CONDA_PREFIX/lib/Notung-2.9.1.5.zip  -d $CONDA_PREFIX/lib/

echo -e '#!/usr/bin/env bash\njava -jar $CONDA_PREFIX/lib/Notung-2.9.1.5.jar $@' > $CONDA_PREFIX/bin/notung
chmod +x $CONDA_PREFIX/bin/notung
