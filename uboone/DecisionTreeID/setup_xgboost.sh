XGBOOSTSYS="/uboone/app/users/${USER}/xgboost"

export XGBOOST_LIB=${XGBOOSTSYS}/lib
export XGBOOST_INC=${XGBOOSTSYS}/include
export RABIT_INC=${XGBOOSTSYS}/rabit/include

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${XGBOOSTSYS}/lib
