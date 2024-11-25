#!/bin/bash
set -e -u -x

function repair_wheel {
    wheel="$1"
    if ! auditwheel show "$wheel"; then
        echo "Skipping non-platform wheel $wheel"
    else
        auditwheel repair "$wheel" --plat "$PLAT" -w /io/wheelhouse/
    fi
}

# Install a system package required by our library
# yum install -y wget
# wget https://archives.boost.io/release/1.75.0/source/boost_1_75_0.tar.gz
cp io/boost_1_75_0.tar.gz .
tar -xzf boost_1_*
cd boost_1_75_0
./bootstrap.sh
./b2 install --with-system --with-iostreams
cd ../..

# tbb
# https://github.com/oneapi-src/oneTBB/archive/refs/tags/v2021.13.0.zip
cp io/v2021.13.0.zip .
unzip v2021.13.0.zip
cd oneTBB-2021.13.0
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DTBB_TEST=OFF ..
cmake --build .
cmake --install .
cd ../..

# blosc
# https://github.com/Blosc/c-blosc/archive/refs/tags/v1.21.5.tar.gz
cp io/v1.21.5.tar.gz .
tar -xzf v1.21.5.tar.gz
cd c-blosc-1.21.5
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
cmake --install .
cd ../..

export MAKEFLAGS="-j$(nproc)"
# Compile wheels

for pypy_ver in $PYPY_VERS; do
  py_prefix=("/opt/python/pp${pypy_ver/./}-"*)
  # "$py_prefix/bin/pip" install -r /io/dev-requirements.txt
  "$py_prefix/bin/pip" wheel /io/ --no-deps -w wheelhouse/
done

for py_ver in $PY_VERS; do
  py_prefix=("/opt/python/cp${py_ver/./}-"*)
  # "$py_prefix/bin/pip" install -r /io/dev-requirements.txt
  "$py_prefix/bin/pip" wheel /io/ --no-deps -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    repair_wheel "$whl"
done

# # Install packages and test
# for PYBIN in /opt/python/*/bin/; do
#     "${PYBIN}/pip" install python-manylinux-demo --no-index -f /io/wheelhouse
#     (cd "$HOME"; "${PYBIN}/nosetests" pymanylinuxdemo)
# done
#
