cmake_minimum_required(VERSION 3.0)

project(quantiNemo)

SET(GCC_COVERAGE_COMPILE_FLAGS "-Iinclude")
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

# files to be compiled
set(SRCS
src/bandmat.cpp
src/mtrand.cpp
src/statservices.cpp
src/basicsimulation.cpp
src/myexcept.cpp
src/submat.cpp
src/cholesky.cpp
src/newmat1.cpp
src/svd.cpp
src/coal_deme.cpp
src/newmat2.cpp
src/tequation.cpp
src/filehandler.cpp
src/newmat3.cpp
src/tgenome.cpp
src/fileservices.cpp
src/newmat4.cpp
src/tlocus.cpp
src/functions.cpp
src/newmat5.cpp
src/tmatrix.cpp
src/genealogy.cpp
src/newmat6.cpp
src/tpatchfitness.cpp
src/hholder.cpp
src/newmat7.cpp
src/tree.cpp
src/indfactory.cpp
src/newmat8.cpp
src/treplicate.cpp
src/tindividual.cpp
src/newmat9.cpp
src/tselection.cpp
src/lce_breed.cpp
src/newmatex.cpp
src/tselectiontrait.cpp
src/lce_coalescence_base.cpp 
src/newmatnl.cpp
src/tselectiontype.cpp
src/lce_coalescence.cpp
src/newmatrm.cpp
src/tsim_manager.cpp
src/lce_disperse.cpp
src/node.cpp
src/tsimulation.cpp
src/lce_extinction.cpp
src/param.cpp
src/tstat_db.cpp
src/lce_misc.cpp
src/patch.cpp
src/tstring.cpp
src/lce_regulation.cpp
src/sort.cpp
src/ttneutral.cpp
src/main.cpp
src/stathandlerbase.cpp 
src/ttquanti.cpp
src/tmetapop.cpp
src/stathandler.cpp 
src/ttrait.cpp
src/metapop_sh.cpp
src/stat_rec_base.cpp
)

set(HEADERS
include/basicsimulation.h
include/lce_misc.h
include/precisio.h
include/treplicate.h
include/coal_deme.h
include/lce_regulation.h
include/randomC11.h
include/tselection.h
include/controlw.h
include/lifecycleevent.h
include/random.h
include/tselectiontrait.h
include/filehandler.h
include/tmetapop.h
include/service.h
include/tselectiontype.h
include/fileservices.h
include/metapop_sh.h
include/simcomponent.h
include/tsim_manager.h
include/functions.h
include/mtrand.h
include/stathandlerbase.h 
include/tsimulation.h
include/genealogy.h
include/myexcept.h
include/stathandler.h
include/tstat_db.h
include/handler.h
include/newmatap.h
include/stat_rec_base.h
include/tstring.h
include/include.h
include/newmat.h
include/statservices.h
include/ttneutral.h
include/indfactory.h
include/newmatio.h
include/tarray.h
include/ttquanti.h
include/tindividual.h
include/newmatnl.h
include/tequation.h
include/ttrait.h
include/lce_breed.h
include/newmatrc.h
include/tgenome.h
include/ttree.h
include/lce_coalescence_base.h 
include/newmatrm.h
include/tlocus.h
include/types.h
include/lce_coalescence.h
include/node.h
include/tmatrix.h
include/version.h
include/lce_disperse.h
include/param.h
include/tpatchfitness.h
include/lce_extinction.h
include/patch.h
include/tree.h
)

#What to compile
add_executable(quantinemo ${SRCS} ${HEADERS})
