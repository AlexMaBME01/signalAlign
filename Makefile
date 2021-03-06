rootPath = ./
include ./include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

signalAlignDependencies =  ${basicLibsDependencies}
signalAlignLib = ${basicLibs}

#coverage = -lgcov -coverage

test_directory = ${rootPath}/src/tests/

htsLib = -L././htslib -lhts

all : sL bD hs python-utils ${libPath}/signalAlignLib.a ${signalAlignBin}/signalAlignLibTests index_fasta \
	  ${signalAlignBin}/compareDistributions \
	  ${signalAlignBin}/signalMachine ${signalAlignBin}/runSignalAlign \
	  ${signalAlignBin}/signalAlignLib.py ${signalAlignBin}/variantCallingLib.py ${signalAlignBin}/alignmentAnalysisLib.py \
	  ${signalAlignBin}/buildHdpUtil ${signalAlignBin}/trainModels ${signalAlignBin}/hdp_pipeline ${signalAlignBin}/test_SignalAlign.py \
	  externals nanoporeParams python_setup  \
	  #${signalAlignBin}/zayante ${signalAlignBin}/bonnyDoon \
	  #${signalAlignBin}/empire ${signalAlignBin}/jamison \

python-utils :
#	echo "NOT PYPORE MAN"
	cd python_utils && python setup.py install


index_fasta : hs ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -I${htsLibRootPath} -o ${signalAlignBin}/index_fasta index_fasta.c ${libPath}/signalAlignLib.a ${signalAlignLib} ${htsLib}

# -I${htsLibPath}  -I${htsLibRootPath}
#_curl_easy_init
core : sL bD ${libPath}/signalAlignLib.a ${signalAlignBin}/signalAlignLibTests ${signalAlignBin}/signalMachine

install: all pip_install

clean_light:
	if [ -d ${signalAlignBin} ]; then rm -r ${signalAlignBin}; fi
	rm -f ${libPath}/signalAlignLib.a

clean :
	if [ -d ${signalAlignBin} ]; then rm -r ${signalAlignBin}; fi
	#rm -r ${signalAlignBin}
	rm -f ${libPath}/signalAlignLib.a
	cd externalTools && make clean

python_setup :
	python setup.py install

pip_install :
	pip install -e .

signalAlignLib : ${libPath}/signalAlignLib.a

sL :
	cd sonLib && make

bD :
	mkdir -v -p ${rootPath}bin

externals :
	cd externalTools && make all

test :
	cd ${signalAlignBin} && ./test_SignalAlign.py
	cd ${binPath} && ./sonLibTests
	cd python_utils && pytest

${signalAlignBin}/compareDistributions : compareDistributions.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags}  -I inc -I${libPath} -o ${signalAlignBin}/compareDistributions compareDistributions.c ${libPath}/signalAlignLib.a ${signalAlignLib}

${signalAlignBin}/signalAlignLibTests : ${libTests} tests/*.h ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags}  -I inc -I${libPath} -I${htsLibRootPath} -I${htsLibPath} -Wno-error -o ${signalAlignBin}/signalAlignLibTests ${libTests} ${libPath}/signalAlignLib.a ${signalAlignLib}  ${htsLib}

${signalAlignBin}/signalMachine : signalMachine.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags}  -I inc -I${libPath} -I${htsLibRootPath} -I${htsLibPath} -o ${signalAlignBin}/signalMachine signalMachine.c ${libPath}/signalAlignLib.a ${signalAlignLib}  ${htsLib}

nanoporeParams : estimateNanoporeParams.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags}  -I inc -I${libPath} -o ${signalAlignBin}/estimateNanoporeParams estimateNanoporeParams.c ${libPath}/signalAlignLib.a ${signalAlignLib}
	cp ${rootPath}src/signalalign/scripts/nanoporeParamRunner.py ${signalAlignBin}/nanoporeParamRunner
	chmod +x ${signalAlignBin}/nanoporeParamRunner

${signalAlignBin}/buildHdpUtil : buildHdpUtil.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags}   -I inc -I${libPath} -o ${signalAlignBin}/buildHdpUtil buildHdpUtil.c ${libPath}/signalAlignLib.a ${signalAlignLib}

${signalAlignBin}/runSignalAlign : ${rootPath}src/signalalign/scripts/runSignalAlign.py
	cp ${rootPath}src/signalalign/scripts/runSignalAlign.py ${signalAlignBin}/runSignalAlign
	chmod +x ${signalAlignBin}/runSignalAlign

${signalAlignBin}/trainModels : ${rootPath}src/signalalign/scripts/trainModels.py
	cp ${rootPath}src/signalalign/scripts/trainModels.py ${signalAlignBin}/trainModels
	chmod +x ${signalAlignBin}/trainModels

${signalAlignBin}/hdp_pipeline : ${rootPath}src/signalalign/scripts/hdp_pipeline.py
	cp ${rootPath}src/signalalign/scripts/hdp_pipeline.py ${signalAlignBin}/hdp_pipeline
	chmod +x ${signalAlignBin}/hdp_pipeline

${signalAlignBin}/test_SignalAlign.py : ${rootPath}src/signalalign/tests/test_SignalAlign.py
	cp ${rootPath}src/signalalign/tests/test_SignalAlign.py ${signalAlignBin}/test_SignalAlign.py
	cp ${rootPath}src/signalalign/tests/event_detection_test.py ${signalAlignBin}/event_detection_test.py
	cp ${rootPath}src/signalalign/tests/mea_algorithm_test.py ${signalAlignBin}/mea_algorithm_test.py
	cp ${rootPath}src/signalalign/tests/fast5_test.py ${signalAlignBin}/fast5_test.py

	chmod +x ${signalAlignBin}/test_SignalAlign.py
	chmod +x ${signalAlignBin}/fast5_test.py
	chmod +x ${signalAlignBin}/mea_algorithm_test.py
	chmod +x ${signalAlignBin}/event_detection_test.py

${signalAlignBin}/zayante : ${rootPath}src/signalalign/scripts/zayante.py
	cp ${rootPath}src/signalalign/scripts/zayante.py ${signalAlignBin}/zayante
	chmod +x ${signalAlignBin}/zayante

${signalAlignBin}/bonnyDoon : ${rootPath}src/signalalign/scripts/bonnyDoon.py
	cp ${rootPath}src/signalalign/scripts/bonnyDoon.py ${signalAlignBin}/bonnyDoon
	chmod +x ${signalAlignBin}/bonnyDoon

${signalAlignBin}/empire : ${rootPath}src/signalalign/scripts/empire.py
	cp ${rootPath}src/signalalign/scripts/empire.py ${signalAlignBin}/empire
	chmod +x ${signalAlignBin}/empire

${signalAlignBin}/jamison : ${rootPath}src/signalalign/scripts/jamison.py
	cp ${rootPath}src/signalalign/scripts/jamison.py ${signalAlignBin}/jamison
	chmod +x ${signalAlignBin}/jamison

${signalAlignBin}/signalAlignLib.py : ${rootPath}src/signalalign/scripts/signalAlignLib.py
	cp ${rootPath}src/signalalign/scripts/signalAlignLib.py ${signalAlignBin}/signalAlignLib.py

${signalAlignBin}/variantCallingLib.py : ${rootPath}src/signalalign/scripts/variantCallingLib.py
	cp ${rootPath}src/signalalign/scripts/variantCallingLib.py ${signalAlignBin}/variantCallingLib.py

${signalAlignBin}/alignmentAnalysisLib.py : ${rootPath}src/signalalign/scripts/alignmentAnalysisLib.py
	cp ${rootPath}src/signalalign/scripts/alignmentAnalysisLib.py ${signalAlignBin}/alignmentAnalysisLib.py

${libPath}/signalAlignLib.a : ${libSources} ${libHeaders} ${stBarDependencies}
	${cxx} ${cflags} -I inc -I ${libPath}/ -I ${htsLibRootPath} -I ${htsLibPath}  ${htsLib} -c ${libSources}
	ar rc signalAlignLib.a *.o
	ranlib signalAlignLib.a
	rm *.o
	mv signalAlignLib.a ${libPath}/
	cp ${libHeaders} ${libPath}/


hs :
	cd htslib && make
