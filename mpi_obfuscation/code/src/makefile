# To adapt to your settings
INCLUDE_DIRS 		= -I/opt/local/include 

LIB_DIRS 			= -L/usr/bin/ -L/usr/local/lib
CXX 				= g++
PARALLEL 			= yes

################################################################################################
# Targets
CLT_DRIVER_TARGET 	= clt_driver


# flags
CXXFLAGS			= -Wall -Ofast -g -std=c++11 -funroll-loops $(INCLUDE_DIRS)
LDLIBS 				= $(LIB_DIRS) -lgmp -lgmpxx -ltcmalloc_and_profiler -dbg


ifeq ($(PARALLEL),yes)
  CXXFLAGS+= -fopenmp
endif

# Files 
CLT_DRIVER_SRCS 		= clt_driver.cpp clt_mlm.cpp utils.cpp
CLT_DRIVER_OBJS 		= $(CLT_DRIVER_SRCS:%.cpp=%.o)
HDRS 				= clt_mlm.h utils.h gperftools-2.4/src/gperftools/profiler.h

# All & targets
all: $(CLT_DRIVER_TARGET)

$(CLT_DRIVER_TARGET): $(CLT_DRIVER_OBJS) 
		$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) -o $@ $(CLT_DRIVER_OBJS) $(LDLIBS) 


# Dependencies of objs
$(CLT_DRIVER_OBJS): $(CLT_CDRIVER_SRCS) $(HDRS) makefile


# cleaning
clean:
	$(RM) -r $(CLT_DRIVER_TARGET) $(CLT_DRIVER_TARGET).dSYM $(CLT_DRIVER_OBJS) 

