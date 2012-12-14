#! /bin/bash

OUTPATH=/afs/cern.ch/user/b/bonato/work/PhysAnalysis/HZZ4L/
DIR=$1

/bin/mkdir -p ${OUTPATH}/${DIR}/PRODFSR_7TeV/data/
/bin/mkdir -p ${OUTPATH}/${DIR}/PRODFSR_7TeV/4mu/
/bin/mkdir -p ${OUTPATH}/${DIR}/PRODFSR_7TeV/4e/
/bin/mkdir -p ${OUTPATH}/${DIR}/PRODFSR_7TeV/2e2mu/
/bin/mkdir -p ${OUTPATH}/${DIR}/PRODFSR_7TeV/CR/

/bin/mkdir -p ${OUTPATH}/${DIR}/PRODFSR_8TeV/data/
/bin/mkdir -p ${OUTPATH}/${DIR}/PRODFSR_8TeV/4mu/
/bin/mkdir -p ${OUTPATH}/${DIR}/PRODFSR_8TeV/4e/
/bin/mkdir -p ${OUTPATH}/${DIR}/PRODFSR_8TeV/2e2mu/
/bin/mkdir -p ${OUTPATH}/${DIR}/PRODFSR_8TeV/CR/


/bin/mkdir -p ${OUTPATH}/${DIR}/JHU_7TeV/4mu/
/bin/mkdir -p ${OUTPATH}/${DIR}/JHU_7TeV/4e/
/bin/mkdir -p ${OUTPATH}/${DIR}/JHU_7TeV/2e2mu/

/bin/mkdir -p ${OUTPATH}/${DIR}/JHU_8TeV/4mu/
/bin/mkdir -p ${OUTPATH}/${DIR}/JHU_8TeV/4e/
/bin/mkdir -p ${OUTPATH}/${DIR}/JHU_8TeV/2e2mu/
