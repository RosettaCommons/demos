#!/bin/bash
# runs ml algs on input, generating .out files for mkresults

usage() {
    echo "$(basename $0) [-pred] [-prefix <string>] [-iter <n>] [-base name[,name...]] [-meta name[,name...]] <train.arff> [<test.arff>] [-- weka-args ...]" > /dev/stderr
    echo "runs all available algorithms; uses test file if given, otherwise does cross-validation on training file" > /dev/stderr
    echo "  -pred:  output predictions instead of statistics and confusion matrix" > /dev/stderr
    echo "  -prefix: specify prefix for output (default is test or train file without \".arff\"" > /dev/stderr
    echo "  -iter: number of iterations to perform for each machine with different random seeds" > /dev/stderr
    echo "  -base: comma-separated list of base machines to use, defaults to all" > /dev/stderr
    echo "  -meta: comma-separated list of meta-classifiers to use; defaults to none" > /dev/stderr
    echo "     \"all\" is an acceptable argument to either -base or -meta" > /dev/stderr
    echo "  -- weka-args: anything after the \"--\" is passed directly to weka" > /dev/stderr
    exit 1
}

mode=$1
field_list="RES|MUT|PT|SPECIES|ACCP"

while [[ "$1" =~ ^- ]]; do
    case $1 in
	-usage)
	usage
	;;
	-pred)
	pred="_pred"
	shift
	;;
	-prefix)
	prefix=$2
	shift; shift
	;;
	-iter)
	iter=$2
	shift; shift
	;;
	-base)
	base_list_in=$2
	shift; shift
	;;
	-meta)
	meta_list_in=$2
	shift; shift
	;;
	--)
	echo "error: training file not specified" > /dev/stderr
	usage
	;;
	*)
	echo "error: unknown argument \"$1\"" > /dev/stderr
	usage
	;;
    esac
done

train_file=$1; shift
if [ -n "$1" -a "$1" != "--" ]; then test_file=$1; shift; fi
if [ "$1" == "--" ]; then
    extras=${*:2}
elif [ -n "$1" ]; then
    echo "extra arguments after test file - perhaps forgotten \"--\"?" > /dev/stderr
    usage
fi

if [ -z "$train_file" ]; then echo "error: training file not specified" > /dev/stderr; usage; fi

sys=$(uname)
case $sys in
    CYGWIN*)
    weka_exe="/cygdrive/c/Program\ Files/Java/jdk1.6.0_12/bin/java -Xmx500m -cp \"$(cygpath -w ~/machinelearn.d/weka.jar);$(cygpath -w ~/machinelearn.d/libsvm.jar)\""
    ;;
    *)
#    weka_exe="/usr/java/latest/bin/java -cp /home/crispy/machinelearn.d/weka.jar:/home/crispy/machinelearn.d/libsvm.jar"
    weka_exe="java -cp /home/crispy/machinelearn.d/weka.jar:/home/crispy/machinelearn.d/libsvm.jar"
    ;;
esac

unset mode_out
if [ -z "$test_file" ]; then mode_out="cv"; else mode_out="test"; fi
if [ -n "$prefix" ]; then              out_prefix=$prefix
else
    if [ "$mode_out" == "test" ]; then out_prefix=${test_file%.arff}
    else                               out_prefix=${train_file%.arff}; fi
    out_prefix=$(basename $out_prefix)
fi

if [ -z "$pred" ]; then
    general_options="-i -v $extras"
else
    # get ids of field in field list
    ft=$(grep '^@attribute' $train_file | egrep -n "@attribute ($field_list) " | cut -d':' -f1)
    field_ids=$(echo $ft | sed -e 's/ /,/g')
    general_options="-p $field_ids $extras"
fi
if [ -z "$iter" ]; then
    if [ "$mode_out" == "test" ]; then iter=1; else iter=10; fi
fi
cv_options="-t $train_file $general_options"
test_options="-t $train_file -T $test_file $general_options"

# untuned C45
base_name["${#base_name[*]}"]="c45"
base_class["${#base_class[*]}"]=weka.classifiers.meta.FilteredClassifier
base_options["${#base_options[*]}"]="-F \"weka.filters.unsupervised.attribute.RemoveType -T string\" -W weka.classifiers.trees.J48 -- -C 0.25 -M 2 -A"

# untuned C45 with reduced error pruning
base_name["${#base_name[*]}"]="c45r"
base_class["${#base_class[*]}"]=weka.classifiers.meta.FilteredClassifier
base_options["${#base_options[*]}"]="-F \"weka.filters.unsupervised.attribute.RemoveType -T string\" -W weka.classifiers.trees.J48 -- -R -N 3 -Q 1 -M 2 -A"

# untuned SMO
base_name["${#base_name[*]}"]="smo"
base_class["${#base_class[*]}"]=weka.classifiers.meta.FilteredClassifier
base_options["${#base_options[*]}"]="-F \"weka.filters.unsupervised.attribute.RemoveType -T string\" -W weka.classifiers.functions.SMO -- -C 1.0 -L 0.0010 -P 1.0E-12 -N 0 -M -V -1 -W 1 -K \"weka.classifiers.functions.supportVector.PolyKernel -C 250007 -E 1.0\""

# SMO tuned on initial set
base_name["${#base_name[*]}"]="smot"
base_class["${#base_class[*]}"]=weka.classifiers.meta.FilteredClassifier
base_options["${#base_options[*]}"]="-F \"weka.filters.unsupervised.attribute.RemoveType -T string\" -W weka.classifiers.functions.SMO -- -C 0.177 -L 0.0010 -P 1.0E-12 -N 0 -M -V -1 -W 1 -K \"weka.classifiers.functions.supportVector.PolyKernel -C 250007 -E 1.0\""

# SMO tuned on run25
base_name["${#base_name[*]}"]="smo25"
base_class["${#base_class[*]}"]=weka.classifiers.meta.FilteredClassifier
base_options["${#base_options[*]}"]="-F \"weka.filters.unsupervised.attribute.RemoveType -T string\" -W weka.classifiers.functions.SMO -- -C 0.0525 -L 0.0010 -P 1.0E-12 -N 0 -M -V -1 -W 1 -K \"weka.classifiers.functions.supportVector.PolyKernel -C 250007 -E 1.0\""

# SMO tuned on run26
base_name["${#base_name[*]}"]="smo26"
base_class["${#base_class[*]}"]=weka.classifiers.meta.FilteredClassifier
base_options["${#base_options[*]}"]="-F \"weka.filters.unsupervised.attribute.RemoveType -T string\" -W weka.classifiers.functions.SMO -- -C 0.0131 -L 0.0010 -P 1.0E-12 -N 0 -M -V -1 -W 1 -K \"weka.classifiers.functions.supportVector.PolyKernel -C 250007 -E 1.0\""

# SMO tuned on run26 cv5
base_name["${#base_name[*]}"]="smo26f"
base_class["${#base_class[*]}"]=weka.classifiers.meta.FilteredClassifier
base_options["${#base_options[*]}"]="-F \"weka.filters.unsupervised.attribute.RemoveType -T string\" -W weka.classifiers.functions.SMO -- -C .105 -L 0.0010 -P 1.0E-12 -N 0 -M -V -1 -W 1 -K \"weka.classifiers.functions.supportVector.PolyKernel -C 250007 -E 1.0\""

# untuned SVM
base_name["${#base_name[*]}"]="svmp"
base_class["${#base_class[*]}"]=weka.classifiers.meta.FilteredClassifier
base_options["${#base_options[*]}"]="-F \"weka.filters.unsupervised.attribute.RemoveType -T string\" -W weka.classifiers.functions.LibSVM -- -S 0 -K 2 -D 3 -G 0.0 -R 0.0 -N 0.5 -M 40.0 -C 1.0 -E 0.0010 -P 0.1 -Z -B"

# SVM tuned on initial set
base_name["${#base_name[*]}"]="svmt"
base_class["${#base_class[*]}"]=weka.classifiers.meta.FilteredClassifier
base_options["${#base_options[*]}"]="-F \"weka.filters.unsupervised.attribute.RemoveType -T string\" -W weka.classifiers.functions.LibSVM -- -S 0 -K 2 -D 3 -G 0.125 -R 0.0 -N 0.5 -M 40.0 -C 1.68 -E 0.0010 -P 0.1 -Z -B"

# SVM tuned for sequence-based features only (run19-seqonly2)
base_name["${#base_name[*]}"]="svmo"
base_class["${#base_class[*]}"]=weka.classifiers.meta.FilteredClassifier
base_options["${#base_options[*]}"]="-F \"weka.filters.unsupervised.attribute.RemoveType -T string\" -W weka.classifiers.functions.LibSVM -- -S 0 -K 2 -D 3 -G 0.03125 -R 0.0 -N 0.5 -M 40.0 -C 2 -E 0.0010 -P 0.1 -Z -B"

# SVM tuned on initial set, bagged
base_name["${#base_name[*]}"]="svmb"
base_class["${#base_class[*]}"]=weka.classifiers.meta.Bagging
base_options["${#base_options[*]}"]="-P 100 -S 1 -I 10 -W weka.classifiers.meta.FilteredClassifier -- -F \"weka.filters.unsupervised.attribute.RemoveType -T string\" -W weka.classifiers.functions.LibSVM -- -S 0 -K 2 -D 3 -G 0.125 -R 0.0 -N 0.5 -M 40.0 -C 1.68 -E 0.0010 -P 0.1 -Z -B"

# SVM tuned on run25, "best" params
base_name["${#base_name[*]}"]="svm25b"
base_class["${#base_class[*]}"]=weka.classifiers.meta.FilteredClassifier
base_options["${#base_options[*]}"]="-F \"weka.filters.unsupervised.attribute.RemoveType -T string\" -W weka.classifiers.functions.LibSVM -- -S 0 -K 2 -D 3 -G 0.0186 -R 0.0 -N 0.5 -M 40.0 -C 13.454 -E 0.0010 -P 0.1 -Z -B"

# SVM tuned on run25, "safe" params
base_name["${#base_name[*]}"]="svm25s"
base_class["${#base_class[*]}"]=weka.classifiers.meta.FilteredClassifier
base_options["${#base_options[*]}"]="-F \"weka.filters.unsupervised.attribute.RemoveType -T string\" -W weka.classifiers.functions.LibSVM -- -S 0 -K 2 -D 3 -G 0.03125 -R 0.0 -N 0.5 -M 40.0 -C 9.514 -E 0.0010 -P 0.1 -Z -B"

# SVM tuned on run26, "best" params
base_name["${#base_name[*]}"]="svm26b"
base_class["${#base_class[*]}"]=weka.classifiers.meta.FilteredClassifier
base_options["${#base_options[*]}"]="-F \"weka.filters.unsupervised.attribute.RemoveType -T string\" -W weka.classifiers.functions.LibSVM -- -S 0 -K 2 -D 3 -G 0.0743 -R 0.0 -N 0.5 -M 40.0 -C 11.3 -E 0.0010 -P 0.1 -Z -B"

# SVM tuned on run26, "safe" params
base_name["${#base_name[*]}"]="svm26s"
base_class["${#base_class[*]}"]=weka.classifiers.meta.FilteredClassifier
base_options["${#base_options[*]}"]="-F \"weka.filters.unsupervised.attribute.RemoveType -T string\" -W weka.classifiers.functions.LibSVM -- -S 0 -K 2 -D 3 -G 0.0526 -R 0.0 -N 0.5 -M 40.0 -C 5.66 -E 0.0010 -P 0.1 -Z -B"

# SVM tuned on run26f, cv5 params
base_name["${#base_name[*]}"]="svm26f"
base_class["${#base_class[*]}"]=weka.classifiers.meta.FilteredClassifier
base_options["${#base_options[*]}"]="-F \"weka.filters.unsupervised.attribute.RemoveType -T string\" -W weka.classifiers.functions.LibSVM -- -S 0 -K 2 -D 3 -G 0.011 -R 0.0 -N 0.5 -M 40.0 -C 22.6 -E 0.0010 -P 0.1 -Z -B"

# SVM tuned on run26f cv5 params, seqonly2
base_name["${#base_name[*]}"]="svm26o"
base_class["${#base_class[*]}"]=weka.classifiers.meta.FilteredClassifier
base_options["${#base_options[*]}"]="-F \"weka.filters.unsupervised.attribute.RemoveType -T string\" -W weka.classifiers.functions.LibSVM -- -S 0 -K 2 -D 3 -G 0.0391 -R 0.0 -N 0.5 -M 40.0 -C 22.6 -E 0.0010 -P 0.1 -Z -B"

# temp slot
base_name["${#base_name[*]}"]="svm26x"
base_class["${#base_class[*]}"]=weka.classifiers.meta.FilteredClassifier
base_options["${#base_options[*]}"]="-F \"weka.filters.unsupervised.attribute.RemoveType -T string\" -W weka.classifiers.functions.LibSVM -- -S 0 -K 2 -D 3 -G 0.0263 -R 0.0 -N 0.5 -M 40.0 -C 32 -E 0.0010 -P 0.1 -Z -B"

meta_name["${#meta_name[*]}"]="bag"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.Bagging
meta_options["${#meta_options[*]}"]="-P 100 -S 1 -I 10 -W"

meta_name["${#meta_name[*]}"]="bag_90"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.Bagging
meta_options["${#meta_options[*]}"]="-P 90 -S 1 -I 10 -W"

meta_name["${#meta_name[*]}"]="bag_90_30"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.Bagging
meta_options["${#meta_options[*]}"]="-P 90 -S 1 -I 30 -W"

meta_name["${#meta_name[*]}"]="bag_80"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.Bagging
meta_options["${#meta_options[*]}"]="-P 80 -S 1 -I 10 -W"

meta_name["${#meta_name[*]}"]="bag_80_30"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.Bagging
meta_options["${#meta_options[*]}"]="-P 80 -S 1 -I 30 -W"

meta_name["${#meta_name[*]}"]="bag_70"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.Bagging
meta_options["${#meta_options[*]}"]="-P 70 -S 1 -I 10 -W"

meta_name["${#meta_name[*]}"]="bag_70_30"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.Bagging
meta_options["${#meta_options[*]}"]="-P 70 -S 1 -I 30 -W"

meta_name["${#meta_name[*]}"]="bag_60"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.Bagging
meta_options["${#meta_options[*]}"]="-P 60 -S 1 -I 10 -W"

meta_name["${#meta_name[*]}"]="bag_60_30"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.Bagging
meta_options["${#meta_options[*]}"]="-P 60 -S 1 -I 30 -W"

meta_name["${#meta_name[*]}"]="cs2"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.CostSensitiveClassifier
meta_options["${#meta_options[*]}"]="-cost-matrix \"[0.0 2.0; 1.0 0.0]\" -S 1 -W "

meta_name["${#meta_name[*]}"]="cs3"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.CostSensitiveClassifier
meta_options["${#meta_options[*]}"]="-cost-matrix \"[0.0 3.0; 1.0 0.0]\" -S 1 -W "

meta_name["${#meta_name[*]}"]="cs4"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.CostSensitiveClassifier
meta_options["${#meta_options[*]}"]="-cost-matrix \"[0.0 4.0; 1.0 0.0]\" -S 1 -W "

meta_name["${#meta_name[*]}"]="cs6"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.CostSensitiveClassifier
meta_options["${#meta_options[*]}"]="-cost-matrix \"[0.0 6.0; 1.0 0.0]\" -S 1 -W "

meta_name["${#meta_name[*]}"]="cs2_bag"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.CostSensitiveClassifier
meta_options["${#meta_options[*]}"]="-cost-matrix \"[0.0 2.0; 1.0 0.0]\" -S 1 -W weka.classifiers.meta.Bagging -- -P 100 -S 1 -I 10 -W"

meta_name["${#meta_name[*]}"]="cs3_bag"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.CostSensitiveClassifier
meta_options["${#meta_options[*]}"]="-cost-matrix \"[0.0 3.0; 1.0 0.0]\" -S 1 -W weka.classifiers.meta.Bagging -- -P 100 -S 1 -I 10 -W"

meta_name["${#meta_name[*]}"]="cs4_bag"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.CostSensitiveClassifier
meta_options["${#meta_options[*]}"]="-cost-matrix \"[0.0 4.0; 1.0 0.0]\" -S 1 -W weka.classifiers.meta.Bagging -- -P 100 -S 1 -I 10 -W"

meta_name["${#meta_name[*]}"]="cs6_bag"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.CostSensitiveClassifier
meta_options["${#meta_options[*]}"]="-cost-matrix \"[0.0 6.0; 1.0 0.0]\" -S 1 -W weka.classifiers.meta.Bagging -- -P 100 -S 1 -I 10 -W"

meta_name["${#meta_name[*]}"]="mc2"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.MetaCost
meta_options["${#meta_options[*]}"]="-cost-matrix \"[0.0 2.0; 1.0 0.0]\" -I 10 -P 100 -S 1 -W"

meta_name["${#meta_name[*]}"]="mc3"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.MetaCost
meta_options["${#meta_options[*]}"]="-cost-matrix \"[0.0 3.0; 1.0 0.0]\" -I 10 -P 100 -S 1 -W"

meta_name["${#meta_name[*]}"]="mc4"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.MetaCost
meta_options["${#meta_options[*]}"]="-cost-matrix \"[0.0 4.0; 1.0 0.0]\" -I 10 -P 100 -S 1 -W"

meta_name["${#meta_name[*]}"]="mc6"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.MetaCost
meta_options["${#meta_options[*]}"]="-cost-matrix \"[0.0 6.0; 1.0 0.0]\" -I 10 -P 100 -S 1 -W"

meta_name["${#meta_name[*]}"]="ada"
meta_class["${#meta_class[*]}"]=weka.classifiers.meta.AdaBoostM1
meta_options["${#meta_options[*]}"]="-P 100 -S 1 -I 10 -W"

if [ -z $base_list_in ]; then base_list=${!base_name[*]}         # default is all if -base not given
elif [ $base_list_in == "all" ]; then base_list=${!base_name[*]} # handle -base all
else
    # find indices of requested base classifiers
    for a in ${base_list_in//,/ }; do
	found=0
	for i in ${!base_name[*]}; do
	    if [ $a == ${base_name[i]} ]; then base_list="$base_list $i"; found=1; break; fi
	done
	if [ $found -eq 0 ]; then echo "skipping unknown base: $a" > /dev/stderr; fi
    done
fi

if [ -z $meta_list_in ]; then meta_list="none"                   # default is none if -meta not given
elif [ $meta_list_in == "all" ]; then meta_list=${!meta_name[*]} # handle -meta all
else
    # find indices of requested meta classifiers
    for a in ${meta_list_in//,/ }; do
	found=0
	for i in ${!meta_name[*]}; do
	    if [ $a == ${meta_name[i]} ]; then meta_list="$meta_list $i"; found=1; break; fi
	done
	if [ $found -eq 0 ]; then echo "skipping unknown meta: $a" > /dev/stderr; fi
    done
fi

#cmd="${meta_class["bag"]} ${meta_options["bag"]} ${base_class["svm"]} -- ${base_options["svm"]}"
for ((n = 1; n <= $iter; n++)); do
    if [ $iter -gt 1 ]; then nl="_$n"; else nl=""; fi
    for i in $base_list; do
	for j in $meta_list; do
	    if [ $j == "none" ]; then
		pre="${base_class[$i]}"
		post="${base_options[$i]}"
		mpf=""
	    else
		pre="${meta_class[$j]}"
		post="${meta_options[$j]} ${base_class[$i]} -- ${base_options[$i]}"
		mpf="_${meta_name[$j]}"
	    fi
	    if [ "$mode_out" == "cv" ]; then
		ofp="${out_prefix}_cv${pred}_${base_name[$i]}${mpf}"
		echo "Method: $ofp CV $n"
		eval $weka_exe ${pre} $cv_options " -s $n " ${post} > ${ofp}${nl}.out 2> ${ofp}${nl}.err
	    else
		ofp="${out_prefix}_test${pred}_${base_name[i]}${mpf}"
		echo "Method: $ofp Test"
		eval $weka_exe ${pre} $test_options ${post} > ${ofp}${nl}.out 2> ${ofp}${nl}.err
	    fi
	done
    done
done
