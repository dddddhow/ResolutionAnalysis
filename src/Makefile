divider="---------------------------------------------------------------------"
InstallDir=/
BinDir =./bin
LibDir = ./lib
IncDirList =./inc
ObjDir =./obj
BackDir=/home/ss/WINDOWS/ResolutionAnalysis/bck/
Date=$(shell date +%F_time%H)
#--------------------------------------#
Projn = a_out
#--------------------------------------#
#CXX=g++
CXX=icpc
ifeq ($(CXX),g++)
	flag_warn = -Wno-unused-result
endif
CPP_CMPL=mpicxx -cxx=${CXX}
C_CMPL  =mpicc
F_CMPL  =gfortran
LINKER  =$(CPP_CMPL)

flag_profiler=
flag_debug=
flag_armadillo=-DARMA_DONT_USE_WRAPPER -DARMA_ALLOW_FAKE_GCC

#flag_profiler= -pg     # use this if gprof
#flag_debug= -DDEBUG

C_CMPL_FLAG   = -fopenmp -O3
CPP_CMPL_FLAG = $(flag_debug) -std=c++11 -fopenmp -fmax-errors=2 -O3 \
				$(flag_profiler) $(flag_warn) \
				$(flag_armadillo)


F_CMPL_FLAG   = -cpp -O3
# list of libs in use (with no -l)
LIB_LIST	  = powershen mkl_rt opencv_core opencv_imgproc opencv_highgui opencv_imgcodecs
# list of libs in use (with -l)
LIB_FLAG      = $(addprefix -L,${LibDir}) $(addprefix -l,${LIB_LIST})
# list of inc
INC_FLAG      = $(addprefix -I,${IncDirList})
# finally, the LINKER_FLAG!
LINKER_FLAG   = -std=c++11 -fopenmp -fmax-errors=2 -O3 \
				$(flag_profiler) ${LIB_FLAG}


CppSrc =$(wildcard *.cpp)
F90Src =$(wildcard *.f90)
CppObj0=$(patsubst %.cpp,%.o,${CppSrc})
F90Obj0=$(patsubst %.f90,%.o,${F90Src})
CppObj =$(addprefix $(ObjDir)/,$(CppObj0))
F90Obj =$(addprefix $(ObjDir)/,$(F90Obj0))
#ProjLib=./lib_local/libpowerdeconv.a
OBJ    =$(CppObj) $(F90Obj)
EXE    =$(BinDir)/$(Projn).e
BACK   =$(BackDir)/$(Projn)_back_$(Date).tgz
Back_Content=$(shell basename $(shell pwd))
ex_list=$(Back_Content)/result *.[oae] *.sw? *log .git
ex_stuff=$(addprefix --exclude=,$(ex_list))

ifeq ($(JOBS),)
  JOBS := $(shell grep -c ^processor /proc/cpuinfo 2>/dev/null)
  ifeq ($(JOBS),)
    JOBS := 1
  endif
endif



#JOBS = 1

all:
	mkdir -p ${ObjDir} ${BinDir};
	make $(EXE) -j ${JOBS};

$(EXE):$(ProjLib) $(OBJ)
	@echo $(divider)
	@echo '			Linking ---> ' $(EXE)
	@echo $(divider)
	$(LINKER) $^ -o $(EXE)  $(LINKER_FLAG) -fopenmp
$(CppObj):$(ObjDir)/%.o:%.cpp
	@echo $(divider)
	@echo '			Compiling '$<'  --> *.o'
	@echo $(divider)
	$(CPP_CMPL) -c  $(INC_FLAG) $< -o $@ $(CPP_CMPL_FLAG) -fopenmp
$(F90Obj):$(ObjDir)/%.o:%.f90
	@echo $(divider)
	@echo '			Compiling '$<'  --> *.o'
	@echo $(divider)
	$(F_CMPL)  -c  $< -o $@ $(F_CMPL_FLAG)
#$(ProjLib):
	#cd lib_local/src && make

backup:
	mkdir -p $(BackDir)
	cd .. && tar $(ex_stuff) -czvf $(BACK) $(Back_Content)
	ls -hl $(BACK)

clean:
	@#CREAT EMPTY files
	@touch $(EXE)
	@touch $(OBJ)
	rm  $(EXE) $(OBJ)
	#cd lib_local/src && make clean
	@pwd; echo $(Date)

relink:
	@touch $(EXE) && rm ${EXE}
	make

b:backup
c:clean

install: ${EXE}
	cp ${EXE} ${InstallDir}
