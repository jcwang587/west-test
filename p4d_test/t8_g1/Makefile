# Makefile

include ../maketest.inc

# Hard-coded GPU count for this test
export NP=1
export NT=8
export PARA_PREFIX=mpirun -np ${NP}
export PARA_PREFIX_QE=mpirun -np ${NP}

default: title description all

all: prepare_inputs run

prepare_inputs:
	bash prepare_inputs.sh

run: pw wstat

pw:
	${PARA_PREFIX_QE} ${BINDIR}/pw.x -i pw.in > pw.out 2> pw.err

wstat: pw
	${PARA_PREFIX} ${BINDIR}/wstat.x -nimage ${NIMAGE} -nband ${NBAND} -i wstat.in > wstat.out 2> wstat.err

title:
	@echo
	@echo "###############"
	@echo "### test001 ###"
	@echo "###############"
	@echo

description:
	@echo
	@echo "SiH4 molecule gamma_only GW"
	@echo
