echo "List modules:"
module list

autoreconf --version

echo *************************
echo *************************

autoreconf -i
echo *************************
echo *************************
./configure --prefix=/usr/workspace/pham20/codes/libxc/libxc-7.0.0/INSTALL_DIR
make
make install
