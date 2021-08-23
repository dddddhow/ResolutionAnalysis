set -x

#g++ main_laplace.cpp -lmkl_rt
icpc main_generateseisprofile.cpp -DARMA_ALLOW_FAKE_GCC -lmkl_rt

./a.out

rm a.out

#./show.sh
