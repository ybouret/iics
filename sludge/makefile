YROOT  := $(HOME)/program/yocto4/sdk
CREATE := $(YROOT)/share/yocto/create.sh

all:
	@echo $(YROOT)

clean:
	@echo "-- removing .dat" && rm -f *.dat
	@echo "-- removing .curve" && rm -f *.curve
	@echo "-- removing .vtk" && rm -f *.vtk
	@echo "-- removing .ppm" && rm -f *.ppm


veryclean: clean
	@echo "-- removing out of sources builds" && cd forge && touch targets && ( cat targets | xargs rm -rf ) && rm -f targets

gnu:
	@bash $(CREATE) src gnu ${BUILD_TYPE}

intel:
	@bash $(CREATE) src intel ${BUILD_TYPE}

clang:
	@bash $(CREATE) src clang ${BUILD_TYPE}

xcode:
	@bash $(CREATE) src xcode ${BUILD_TYPE}
	
vs9:
	@bash $(CREATE) src vs9 ${BUILD_TYPE}

vs10:
	@bash $(CREATE) src vs10 ${BUILD_TYPE}

