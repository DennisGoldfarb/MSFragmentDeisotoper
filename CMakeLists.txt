#cmake_minimum_required(VERSION 3.3)
project(MSFragmentDeisotoper)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES FragmentIsotopeCalculator.cpp FragmentIsotopeCalculator.h libmercury++.cpp libmercury++.h Residue.cpp Residue.h Peptide.cpp Peptide.h Molecule.cpp Molecule.h MolecularFormula.cpp MolecularFormula.h SubatomicParticle.cpp SubatomicParticle.h Element.cpp Element.h)
set(MSFRAGMENTDEISOTOPER_SOURCE_FILES ${SOURCE_FILES} main.cpp SpectralMultiplexerMain.cpp)
set(SPECTRALMULTIPLEXER_SOURCE_FILES ${SOURCE_FILES} SpectralMultiplexerMain.cpp AveragineModel.cpp AveragineModel.h)
set(GLOBALISOTOPICDISTRIBUTIONSTATS_SOURCE_FILES ${SOURCE_FILES} GlobalIsotopicDistributionStats.cpp b_ion.cpp b_ion.h y_ion.cpp y_ion.h)

add_executable(MSFragmentDeisotoper ${MSFRAGMENTDEISOTOPER_SOURCE_FILES})
add_executable(SpectralMultiplexer ${SPECTRALMULTIPLEXER_SOURCE_FILES})
add_executable(GlobalIsotopicDistributionStats ${GLOBALISOTOPICDISTRIBUTIONSTATS_SOURCE_FILES})