set -x
icpc main.cpp  -std=c++11  -DARMA_ALLOW_FAKE_GCC -lpowershen -lopencv_core -lopencv_imgproc -lopencv_highgui -lopencv_imgcodecs -lmkl_rt

./a.out

rm a.out
