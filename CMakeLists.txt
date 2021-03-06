cmake_minimum_required(VERSION 3.3)

# Use the default gcc/g++ version
find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

# require at least gcc 5.0
if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
    message(WARNING "Require at least gcc-5.0")
endif()

project(InSeqt)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -w")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/bin)

# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")  


# include directory in find path where all dependency modules exist
include_directories(./)

# InSeqt binaries
set(SOURCE_FILES_BASICSTATS ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc ryggrad/src/visual/Color.cc ryggrad/src/visual/Eps.cc ryggrad/src/visual/Histogram.cc src/BasicStats.cc) 
set(SOURCE_FILES_BASICSTATS2HTML ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc  ryggrad/src/base/StringUtil.cc ryggrad/src/util/mutil.cc src/BasicStats2HTML.cc)
set(SOURCE_FILES_ESTIMATEERRORS  ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc  ryggrad/src/base/StringUtil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc ryggrad/src/visual/Color.cc ryggrad/src/visual/Eps.cc ryggrad/src/visual/Histogram.cc src/EstimateErrors.cc) 

set(SOURCE_FILES_FINDOVERLAPCANDS ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/util/mutil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/aligns/KmerAlignCore.cc src/FindOverlapCands.cc)

set(SOURCE_FILES_PULLREADS ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/util/mutil.cc ryggrad/src/general/DNAVector.cc src/PullReads.cc)

set(SOURCE_FILES_INSEQT ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc ryggrad/src/visual/Color.cc ryggrad/src/visual/Eps.cc ryggrad/src/visual/Histogram.cc src/InSeqt.cc) 
set(SOURCE_FILES_LAPSTATS  ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc  ryggrad/src/base/StringUtil.cc ryggrad/src/util/mutil.cc ryggrad/src/visual/Color.cc ryggrad/src/visual/Eps.cc ryggrad/src/visual/Histogram.cc src/LapStats.cc)  
set(SOURCE_FILES_MAKEBASICREPORT  ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc  ryggrad/src/base/StringUtil.cc ryggrad/src/util/mutil.cc ryggrad/src/visual/Color.cc ryggrad/src/visual/Eps.cc ryggrad/src/visual/Histogram.cc src/Database.cc src/HTMLTemp.cc src/MakeBasicReport.cc)
set(SOURCE_FILES_MAKEOVERLAPREPORT  ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc  ryggrad/src/base/StringUtil.cc ryggrad/src/util/mutil.cc ryggrad/src/visual/Color.cc ryggrad/src/visual/Eps.cc ryggrad/src/visual/Histogram.cc src/Database.cc src/HTMLTemp.cc src/MakeOverlapReport.cc)
set(SOURCE_FILES_PARSEALIGNS ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc  ryggrad/src/base/StringUtil.cc ryggrad/src/util/mutil.cc src/ParseAligns.cc)
set(SOURCE_FILES_PLOTGCCONTENT ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc  ryggrad/src/base/StringUtil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc ryggrad/src/visual/Color.cc ryggrad/src/visual/Eps.cc ryggrad/src/visual/Histogram.cc src/PlotGCContent.cc)  
set(SOURCE_FILES_RUNFINDOVERLAPCANDS ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc  ryggrad/src/base/StringUtil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc src/RunFindOverlapCands.cc) 
set(SOURCE_FILES_SCATTERPLOT ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc  ryggrad/src/base/StringUtil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc ryggrad/src/visual/Color.cc ryggrad/src/visual/Compounds.h ryggrad/src/visual/Eps.cc ryggrad/src/visual/ScatterPlot.cc)  
set(SOURCE_FILES_MAKESITEMAP ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc  ryggrad/src/base/StringUtil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc src/MakeSiteMap.cc)  
set(SOURCE_FILES_FINDSITELAPS ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc  ryggrad/src/base/StringUtil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc src/FindSiteLaps.cc)  
set(SOURCE_FILES_POOLBYSITES ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc  ryggrad/src/base/StringUtil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc src/PoolBySites.cc)  
set(SOURCE_FILES_KMERPROFILE ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/util/mutil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/aligns/KmerAlignCore.cc src/KmerProfile.cc)
set(SOURCE_FILES_SCOREKMERS ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/util/mutil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/aligns/KmerAlignCore.cc src/ScoreKmers.cc)

set(SOURCE_FILES_SIMULATEREADDATA ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/util/mutil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/aligns/KmerAlignCore.cc src/SimulateReadData.cc)

set(SOURCE_FILES_BUILDKMERPROFILE ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/util/mutil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/aligns/KmerAlignCore.cc src/BuildKmerProfile.cc)
set(SOURCE_FILES_TESTKMERPROFILE ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/util/mutil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/aligns/KmerAlignCore.cc src/TestKmerProfile.cc)
set(SOURCE_FILES_ASSEMBLESITELAPS ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc  ryggrad/src/base/StringUtil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc src/AssembleSiteLaps.cc) 
set(SOURCE_FILES_REPLACEFASTAHEAD ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc  ryggrad/src/base/StringUtil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc src/ReplaceFastaHead.cc) 


add_executable(BasicStats           ${SOURCE_FILES_BASICSTATS})
add_executable(BasicStats2HTML      ${SOURCE_FILES_BASICSTATS2HTML})
add_executable(EstimateErrors       ${SOURCE_FILES_ESTIMATEERRORS})
add_executable(FindOverlapCands     ${SOURCE_FILES_FINDOVERLAPCANDS})
add_executable(InSeqt               ${SOURCE_FILES_INSEQT})
add_executable(LapStats             ${SOURCE_FILES_LAPSTATS})
add_executable(MakeBasicReport      ${SOURCE_FILES_MAKEBASICREPORT})
add_executable(MakeOverlapReport    ${SOURCE_FILES_MAKEOVERLAPREPORT})
add_executable(ParseAligns          ${SOURCE_FILES_PARSEALIGNS})
add_executable(PlotGCContent        ${SOURCE_FILES_PLOTGCCONTENT})
add_executable(RunFindOverlapCands  ${SOURCE_FILES_RUNFINDOVERLAPCANDS})
add_executable(ScatterPlot          ${SOURCE_FILES_SCATTERPLOT}) 
add_executable(PullReads            ${SOURCE_FILES_PULLREADS}) 
add_executable(MakeSiteMap          ${SOURCE_FILES_MAKESITEMAP}) 
add_executable(FindSiteLaps         ${SOURCE_FILES_FINDSITELAPS}) 
add_executable(KmerProfile          ${SOURCE_FILES_KMERPROFILE}) 
add_executable(ScoreKmers           ${SOURCE_FILES_SCOREKMERS}) 
add_executable(PoolBySites          ${SOURCE_FILES_POOLBYSITES}) 
add_executable(SimulateReadData     ${SOURCE_FILES_SIMULATEREADDATA})
add_executable(BuildKmerProfile     ${SOURCE_FILES_BUILDKMERPROFILE})
add_executable(TestKmerProfile      ${SOURCE_FILES_TESTKMERPROFILE})
add_executable(AssembleSiteLaps     ${SOURCE_FILES_ASSEMBLESITELAPS})
add_executable(ReplaceFastaHead     ${SOURCE_FILES_REPLACEFASTAHEAD})

