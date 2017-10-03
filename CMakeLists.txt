cmake_minimum_required(VERSION 3.7)
project(signalAlign_c)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -g -Wall --pedantic -fopenmp -funroll-loops -lm")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -DNDEBUG")

set(LASTZ_SOURCE
        externalTools/lastz-distrib-1.03.54/src/align_diffs.c
        externalTools/lastz-distrib-1.03.54/src/align_diffs.h
        externalTools/lastz-distrib-1.03.54/src/axt.c
        externalTools/lastz-distrib-1.03.54/src/axt.h
        externalTools/lastz-distrib-1.03.54/src/build_options.h
        externalTools/lastz-distrib-1.03.54/src/capsule.c
        externalTools/lastz-distrib-1.03.54/src/capsule.h
        externalTools/lastz-distrib-1.03.54/src/chain.c
        externalTools/lastz-distrib-1.03.54/src/chain.h
        externalTools/lastz-distrib-1.03.54/src/cigar.c
        externalTools/lastz-distrib-1.03.54/src/cigar.h
        externalTools/lastz-distrib-1.03.54/src/continuity_dist.c
        externalTools/lastz-distrib-1.03.54/src/continuity_dist.h
        externalTools/lastz-distrib-1.03.54/src/coverage_dist.c
        externalTools/lastz-distrib-1.03.54/src/coverage_dist.h
        externalTools/lastz-distrib-1.03.54/src/diag_hash.c
        externalTools/lastz-distrib-1.03.54/src/diag_hash.h
        externalTools/lastz-distrib-1.03.54/src/dna_utilities.c
        externalTools/lastz-distrib-1.03.54/src/dna_utilities.h
        externalTools/lastz-distrib-1.03.54/src/edit_script.c
        externalTools/lastz-distrib-1.03.54/src/edit_script.h
        externalTools/lastz-distrib-1.03.54/src/gapped_extend.c
        externalTools/lastz-distrib-1.03.54/src/gapped_extend.h
        externalTools/lastz-distrib-1.03.54/src/genpaf.c
        externalTools/lastz-distrib-1.03.54/src/genpaf.h
        externalTools/lastz-distrib-1.03.54/src/gfa.c
        externalTools/lastz-distrib-1.03.54/src/gfa.h
        externalTools/lastz-distrib-1.03.54/src/identity_dist.c
        externalTools/lastz-distrib-1.03.54/src/identity_dist.h
        externalTools/lastz-distrib-1.03.54/src/infer_scores.c
        externalTools/lastz-distrib-1.03.54/src/infer_scores.h
        externalTools/lastz-distrib-1.03.54/src/lastz.c
        externalTools/lastz-distrib-1.03.54/src/lastz.h
        externalTools/lastz-distrib-1.03.54/src/lav.c
        externalTools/lastz-distrib-1.03.54/src/lav.h
        externalTools/lastz-distrib-1.03.54/src/maf.c
        externalTools/lastz-distrib-1.03.54/src/maf.h
        externalTools/lastz-distrib-1.03.54/src/masking.c
        externalTools/lastz-distrib-1.03.54/src/masking.h
        externalTools/lastz-distrib-1.03.54/src/output.c
        externalTools/lastz-distrib-1.03.54/src/output.h
        externalTools/lastz-distrib-1.03.54/src/pos_table.c
        externalTools/lastz-distrib-1.03.54/src/pos_table.h
        externalTools/lastz-distrib-1.03.54/src/quantum.c
        externalTools/lastz-distrib-1.03.54/src/quantum.h
        externalTools/lastz-distrib-1.03.54/src/sam.c
        externalTools/lastz-distrib-1.03.54/src/sam.h
        externalTools/lastz-distrib-1.03.54/src/seed_search.c
        externalTools/lastz-distrib-1.03.54/src/seed_search.h
        externalTools/lastz-distrib-1.03.54/src/seeds.c
        externalTools/lastz-distrib-1.03.54/src/seeds.h
        externalTools/lastz-distrib-1.03.54/src/segment.c
        externalTools/lastz-distrib-1.03.54/src/segment.h
        externalTools/lastz-distrib-1.03.54/src/sequences.c
        externalTools/lastz-distrib-1.03.54/src/sequences.h
        externalTools/lastz-distrib-1.03.54/src/text_align.c
        externalTools/lastz-distrib-1.03.54/src/text_align.h
        externalTools/lastz-distrib-1.03.54/src/tweener.c
        externalTools/lastz-distrib-1.03.54/src/tweener.h
        externalTools/lastz-distrib-1.03.54/src/utilities.c
        externalTools/lastz-distrib-1.03.54/src/utilities.h
        externalTools/lastz-distrib-1.03.54/src/version.h
        )


set(SONLIB_HEADERS
        sonLib/C/inc
        sonLib/C/impl
        sonLib/C/tests
        sonLib/externalTools/cutest
        sonLib/externalTools/quicktree_1.1/include
        )

set(SIGNALALIGN_HEADERS
        inc
#        inc/continuousHmm.h
#        inc/discreteHmm.h
#        inc/hdp.h
#        inc/hdp_math_utils.h
#        inc/multipleAligner.h
#        inc/nanopore.h
#        inc/nanopore_hdp.h
#        inc/pairwiseAligner.h
#        inc/ranlib.h
#        inc/rnglib.h
#        inc/signalMachineUtils.h
#        inc/stateMachine.h
        )


#include_directories(${LASTZ_SOURCE})
include_directories(externalTools/lastz-distrib-1.03.54/src)
include_directories(${SONLIB_HEADERS})
include_directories(${SIGNALALIGN_HEADERS})

set(SONLIB_SOURCE_FILES
#        sonLib/C/impl/lz4.h
#        sonLib/C/impl/lz4hc.h
#        sonLib/C/impl/sonLibGlobalsInternal.h
#        sonLib/C/impl/sonLibKVDatabasePrivate.h
#        sonLib/C/impl/sonLibListPrivate.h
        sonLib/C/impl/avl.c
        sonLib/C/impl/bioioC.c
        sonLib/C/impl/commonC.c
        sonLib/C/impl/fastCMaths.c
        sonLib/C/impl/hashTableC.c
        sonLib/C/impl/hashTableC_itr.c
        sonLib/C/impl/lz4.c
        sonLib/C/impl/lz4hc.c
        sonLib/C/impl/pairwiseAlignment.c
        sonLib/C/impl/sonLibAlign.c
        sonLib/C/impl/sonLibCache.c
        sonLib/C/impl/sonLibCommon.c
        sonLib/C/impl/sonLibCompression.c
        sonLib/C/impl/sonLibConnectivity.c
        sonLib/C/impl/sonLibEulerTour.c
        sonLib/C/impl/sonLibExcept.c
        sonLib/C/impl/sonLibFile.c
        sonLib/C/impl/sonLibHash.c
        sonLib/C/impl/sonLibKVDatabase.c
        sonLib/C/impl/sonLibKVDatabase_BigRecordFile.c
#        sonLib/C/impl/sonLibKVDatabase_KyotoTycoon.cpp
        sonLib/C/impl/sonLibKVDatabase_MySql.c
        sonLib/C/impl/sonLibKVDatabase_TokyoCabinet.c
        sonLib/C/impl/sonLibKVDatabaseConf.c
        sonLib/C/impl/sonLibList.c
        sonLib/C/impl/sonLibNaiveConnectivity.c
        sonLib/C/impl/sonLibRandom.c
        sonLib/C/impl/sonLibSet.c
        sonLib/C/impl/sonLibSortedSet.c
        sonLib/C/impl/sonLibString.c
        sonLib/C/impl/sonLibTreap.c
        sonLib/C/impl/sonLibTree.c
        sonLib/C/impl/sonLibTuples.c
        sonLib/C/impl/stGraph.c
        sonLib/C/impl/stMatrix.c
        sonLib/C/impl/stPhylogeny.c
        sonLib/C/impl/stPosetAlignment.c
        sonLib/C/impl/stSafeC.c
        sonLib/C/impl/stThreadPool.c

### Headers
#        sonLib/C/inc/avl.h
#        sonLib/C/inc/bioioC.h
#        sonLib/C/inc/commonC.h
#        sonLib/C/inc/fastCMaths.h
#        sonLib/C/inc/hashTableC.h
#        sonLib/C/inc/hashTableC_itr.h
#        sonLib/C/inc/hashTablePrivateC.h
#        sonLib/C/inc/pairwiseAlignment.h
#        sonLib/C/inc/sonLib.h
#        sonLib/C/inc/sonLibAlign.h
#        sonLib/C/inc/sonLibCache.h
#        sonLib/C/inc/sonLibCommon.h
#        sonLib/C/inc/sonLibCompression.h
#        sonLib/C/inc/sonLibConnectivity.h
#        sonLib/C/inc/sonLibEulerTour.h
#        sonLib/C/inc/sonLibExcept.h
#        sonLib/C/inc/sonLibFile.h
#        sonLib/C/inc/sonLibHash.h
#        sonLib/C/inc/sonLibKVDatabase.h
#        sonLib/C/inc/sonLibKVDatabaseConf.h
#        sonLib/C/inc/sonLibList.h
#        sonLib/C/inc/sonLibNaiveConnectivity.h
#        sonLib/C/inc/sonLibRandom.h
#        sonLib/C/inc/sonLibSet.h
#        sonLib/C/inc/sonLibSortedSet.h
#        sonLib/C/inc/sonLibString.h
#        sonLib/C/inc/sonLibTreap.h
#        sonLib/C/inc/sonLibTree.h
#        sonLib/C/inc/sonLibTuples.h
#        sonLib/C/inc/sonLibTypes.h
#        sonLib/C/inc/stGraph.h
#        sonLib/C/inc/stMatrix.h
#        sonLib/C/inc/stPhylogeny.h
#        sonLib/C/inc/stPosetAlignment.h
#        sonLib/C/inc/stSafeC.h
#        sonLib/C/inc/stSpimapLayer.h
#        sonLib/C/inc/stThreadPool.h

### Tests
#        sonLib/C/tests/allTests.c
#        sonLib/C/tests/cigarsTest.c
#        sonLib/C/tests/dynamicConnectivity_performance.c
#        sonLib/C/tests/fastaCTest.c
#        sonLib/C/tests/kt_connect_test.cpp
#        sonLib/C/tests/kvDatabaseTest.c
#        sonLib/C/tests/kvDatabaseTestCommon.c
#        sonLib/C/tests/kvDatabaseTestCommon.h
#        sonLib/C/tests/sonLibAlignTest.c
#        sonLib/C/tests/sonLibCacheTest.c
#        sonLib/C/tests/sonLibCommonTest.c
#        sonLib/C/tests/sonLibCompressionTest.c
#        sonLib/C/tests/sonLibConnectivityTests.c
#        sonLib/C/tests/sonLibDoubleTuplesTest.c
#        sonLib/C/tests/sonLibEdgeContainerTests.c
#        sonLib/C/tests/sonLibEulerTest.c
#        sonLib/C/tests/sonLibExceptTest.c
#        sonLib/C/tests/sonLibFileTest.c
#        sonLib/C/tests/sonLibGlobalsTest.h
#        sonLib/C/tests/sonLibGraphTest.c
#        sonLib/C/tests/sonLibHashTest.c
#        sonLib/C/tests/sonLibIntTuplesTest.c
#        sonLib/C/tests/sonLibListTest.c
#        sonLib/C/tests/sonLibPosetAlignmentTest.c
#        sonLib/C/tests/sonLibRandomTest.c
#        sonLib/C/tests/sonLibSetTest.c
#        sonLib/C/tests/sonLibSortedSetTest.c
#        sonLib/C/tests/sonLibStringTest.c
#        sonLib/C/tests/sonLibTreapTest.c
#        sonLib/C/tests/sonLibTreeTest.c
#        sonLib/C/tests/sonLibTuplesTest.c
#        sonLib/C/tests/stMatrixTest.c
#        sonLib/C/tests/stPhylogenyTest.c
#        sonLib/C/tests/stThreadPoolTest.c
        )


set(CUTEST_SOURCE_FILES
### CuTest
#        sonLib/externalTools/cutest/CuTest.h
#        sonLib/externalTools/cutest/AllTests.c
        sonLib/externalTools/cutest/CuTest.c
        sonLib/externalTools/cutest/CuTestTest.c
        )

set(QUICKTREE_SOURCE_FILES
### Quicktree headers
#        sonLib/externalTools/quicktree_1.1/include/align.h
#        sonLib/externalTools/quicktree_1.1/include/aln2dist.h
#        sonLib/externalTools/quicktree_1.1/include/aln2tree.h
#        sonLib/externalTools/quicktree_1.1/include/buildtree.h
#        sonLib/externalTools/quicktree_1.1/include/cluster.h
#        sonLib/externalTools/quicktree_1.1/include/dist2tree.h
#        sonLib/externalTools/quicktree_1.1/include/distancemat.h
#        sonLib/externalTools/quicktree_1.1/include/options.h
#        sonLib/externalTools/quicktree_1.1/include/sequence.h
#        sonLib/externalTools/quicktree_1.1/include/tree.h
#        sonLib/externalTools/quicktree_1.1/include/util.h
### these have a 'main' function
#        sonLib/externalTools/quicktree_1.1/src/aln2dist.c
#        sonLib/externalTools/quicktree_1.1/src/aln2tree.c
#        sonLib/externalTools/quicktree_1.1/src/dist2tree.c
#        sonLib/externalTools/quicktree_1.1/src/distancemat_merops.c
#        sonLib/externalTools/quicktree_1.1/src/quicktree.c
### actual source files
        sonLib/externalTools/quicktree_1.1/src/align.c
        sonLib/externalTools/quicktree_1.1/src/buildtree.c
        sonLib/externalTools/quicktree_1.1/src/cluster.c
        sonLib/externalTools/quicktree_1.1/src/distancemat.c
        sonLib/externalTools/quicktree_1.1/src/options.c
        sonLib/externalTools/quicktree_1.1/src/sequence.c
        sonLib/externalTools/quicktree_1.1/src/tree.c
        sonLib/externalTools/quicktree_1.1/src/util.c
        )

set(SOURCE_FILES
        impl/continuousHmm.c
        impl/discreteHmm.c
        impl/hdp.c
        impl/hdp_math_utils.c
        impl/multipleAligner.c
        impl/nanopore.c
        impl/nanopore_hdp.c
        impl/pairwiseAligner.c
        impl/ranlib.c
        impl/rnglib.c
        impl/signalMachineUtils.c
        impl/stateMachine.c
### tests
#        tests/allTests.c
#        tests/hdpTests.c
#        tests/nanoporeHdpTests.c
#        tests/randomSequences.c
#        tests/randomSequences.h
#        tests/signalPairwiseAlignerTest.c
#        tests/stateMachineTests.c
#        tests/variableOrderPairwiseAlignerTests.c
## our executables
#        buildHdpUtil.c
#        compareDistributions.c
#        estimateNanoporeParams.c
#        signalMachine.c
        )


find_library(ZLIB NAMES z)
if(${ZLIB} STREQUAL "ZLIB-NOTFOUND")
    message(FATAL_ERROR "Couldn't find the 'z' library")
endif()



set(ALL_SRC
        ${QUICKTREE_SOURCE_FILES}
        ${CUTEST_SOURCE_FILES}
        ${SONLIB_SOURCE_FILES}
        ${SOURCE_FILES})

add_executable(signalAlign ${ALL_SRC} signalMachine.c)
target_link_libraries(signalAlign m z)

add_executable(buildHdpUtil ${ALL_SRC} buildHdpUtil.c)
target_link_libraries(buildHdpUtil m z)

add_executable(compareDistributions ${ALL_SRC} compareDistributions.c)
target_link_libraries(compareDistributions m z)

add_executable(estimateNanoporeParams ${ALL_SRC} estimateNanoporeParams.c)
target_link_libraries(estimateNanoporeParams m z)
