#!/bin/bash

usage() {
    echo "usage: $(basename $0) <file.arff>" > /dev/stderr
    echo "predict temperature sensitive mutants using trained SVM-LIN and SVM-RBF models" > /dev/stderr
    echo "  using input ARFF-format file" > /dev/stderr
    exit 1
}

infile=$1
outfile_pred_lin=${1/.arff/-svmlin.raw}
outfile_err_lin=${1/.arff/-svmlin.err}
outfile_pred_rbf=${1/.arff/-svmrbf.raw}
outfile_err_rbf=${1/.arff/-svmrbf.err}

java -cp svm/weka.jar:svm/libsvm.jar weka.classifiers.meta.FilteredClassifier -l svm/svmlin.model -T $infile -p 1,2,3,4,10 > $outfile_pred_lin 2> $outfile_err_lin
showpred $outfile_pred_lin > ${outfile_pred_lin/.raw/.txt}
java -cp svm/weka.jar:svm/libsvm.jar weka.classifiers.meta.FilteredClassifier -l svm/svmrbf.model -T $infile -p 1,2,3,4,10 > $outfile_pred_rbf 2> $outfile_err_rbf
showpred $outfile_pred_rbf > ${outfile_pred_rbf/.raw/.txt}
