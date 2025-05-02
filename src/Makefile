CXXFLAGS  = -Wall -std=c++11


all: Non_ecological_clean

Non_ecological_clean: multimodeFSD_multispecies_nonecological_v2.cpp
	$(CXX) $(CXXFLAGS) $? -o $@

.PHONY : clean

clean :
	-rm Non_ecological_clean
