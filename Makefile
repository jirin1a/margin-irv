PROGRAM = marginirv

RM = rm -rf
OBJDIR = obj

CPLEX=/Applications/CPLEX_Studio_Community221/
CPLEXLIB=$(CPLEX)/cplex/lib/x86-64_osx/static_pic
CONCERTLIB=$(CPLEX)/concert/lib/x86-64_osx/static_pic/


BASEDIRS = \
	-I. \
	-I$(CPLEX)/cplex/include \
	-I$(CPLEX)/concert/include 


INCLUDEDIRS = $(BASEDIRS)

CXX = g++
LD =
SUFFIX = o

CXXFLAGS = -O0 -Wall -pedantic -g $(INCLUDEDIRS) -m64 -fPIC \
	-fexceptions -DNEBUG -DIL_STD -Wno-long-long \
	-Wno-attributes -Wno-ignored-attributes -fpermissive -Wno-sign-compare


#LDFLAGS =  -lboost_system  -lboost_filesystem \
#	-L$(CPLEXLIB) -lilocplex -lcplex \
#	-L$(CONCERTLIB) -lconcert -m64 -lm -pthread -lrt -ldl
LDFLAGS =  -lboost_system  -lboost_filesystem \
	-L$(CPLEXLIB) -lilocplex -lcplex \
	-L$(CONCERTLIB) -lconcert -m64 -lm -pthread -lSystem.B -ldl


RENAME = -o

CXXSOURCES = \
	marginirv.cpp \
	sim_irv.cpp \
	model.cpp \
	tree_irv.cpp \
	irv_distance.cpp \
	nonmono_irv_distance.cpp
	
CXXOBJECTS = $(patsubst %.cpp, $(OBJDIR)/%.$(SUFFIX), $(CXXSOURCES))

all : clean $(PROGRAM)

$(PROGRAM) : $(CXXOBJECTS)
	$(CXX) -o ${@} $(CXXOBJECTS) $(LD) $(LDFLAGS) 

$(OBJDIR)/%.$(SUFFIX) : %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(RENAME) $(@D)/$(@F) -c $(<)

clean:
	$(RM) $(CXXOBJECTS) $(PROGRAM) $(OBJDIR)


