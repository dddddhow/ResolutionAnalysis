set -x
icpc main_generateseisprofile.cpp  -DARMA_ALLOW_FAKE_GCC -lmkl_rt
./a.out
rm a.out
