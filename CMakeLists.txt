cmake_minimum_required(VERSION 3.3)
project(MSFragmentDeisotoper)

include(FindXercesC FindZLIB)
find_package(XercesC REQUIRED)
find_package(ZLIB REQUIRED)


include_directories( ${XercesC_INCLUDE_DIR} ${ZLIB_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES AveragineModel.cpp FragmentIsotopeCalculator.cpp FragmentIsotopeCalculator.h libmercury++.cpp libmercury++.h Residue.cpp Residue.h Peptide.cpp Peptide.h Molecule.cpp Molecule.h MolecularFormula.cpp MolecularFormula.h SubatomicParticle.cpp SubatomicParticle.h Element.cpp Element.h)
set(MSFRAGMENTDEISOTOPER_SOURCE_FILES ${SOURCE_FILES} main.cpp SpectralMultiplexerMain.cpp)
set(SPECTRALMULTIPLEXER_SOURCE_FILES ${SOURCE_FILES} SpectralMultiplexerMain.cpp AveragineModel.cpp AveragineModel.h)
set(GLOBALISOTOPICDISTRIBUTIONSTATS_SOURCE_FILES ${SOURCE_FILES} GlobalIsotopicDistributionStats.cpp b_ion.cpp b_ion.h y_ion.cpp y_ion.h)
set(MODELEVALUATOR_SOURCE_FILES ${SOURCE_FILES} ModelEvaluator.cpp Base64.h Base64.cpp TensorSplineModel.cpp TensorSplineModel.h FragmentIsotopeApproximator.cpp FragmentIsotopeApproximator.h kseq.h b_ion.h b_ion.cpp y_ion.h y_ion.cpp Peptide.h Peptide.cpp Histogram.cpp Histogram.h)

#add_executable(MSFragmentDeisotoper ${MSFRAGMENTDEISOTOPER_SOURCE_FILES})
#add_executable(SpectralMultiplexer ${SPECTRALMULTIPLEXER_SOURCE_FILES})
add_executable(GlobalIsotopicDistributionStats ${GLOBALISOTOPICDISTRIBUTIONSTATS_SOURCE_FILES})
add_executable(ModelEvaluator ${MODELEVALUATOR_SOURCE_FILES})

target_link_libraries (ModelEvaluator ${XercesC_LIBRARY} ${ZLIB_LIBRARIES})
