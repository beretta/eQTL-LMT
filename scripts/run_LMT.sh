#!/bin/sh

DS=$1
DS_CSV=${DS}.csv
DS_ARFF=${DS}.arff

# Convert csv to arff file.
java -cp ./weka.jar weka.core.converters.CSVLoader ${DS_CSV} > ${DS_ARFF}
# IMPORTANT: Before running this part it is necessary to change the arff file. The target must not be numeric, but a class.
# Hence, it must be: @attribute Class {0,1} instead of @attribute Class numeric
sed -i -e 's/numeric/{0,1}/' ${DS_ARFF}
# Loop for running the 30 experiments. Each experiment uses a different partition between training and test instances. 
#Training instances are used to build the model, test instances for assessing its ability in handling unseen data
for ((j = 1; j<=30; j++))
do
	# This command will call the Logistic Model Tree implemented in the WEKA machine learning tool. The parameters are the following:
	# -I: numBoostingIterations -- Set a fixed number of iterations for LogitBoost. 
	# If >= 0, this sets a fixed number of LogitBoost iterations that is used everywhere in the tree. If < 0, the number is cross-validated.
	# -M 15 minNumInstances -- Set the minimum number of instances at which a node is considered for splitting. A value of 15 is used.
	# weightTrimBeta -- Set the beta value used for weight trimming in LogitBoost. Only instances carrying (1 - beta)% of the 
	# weight from previous iteration are used in the next iteration. Set to 0 for no weight trimming (this is the ideal setting for the problem considered). 
	# fastRegression flag = true (default). Use heuristic that avoids cross-validating the number of Logit-Boost iterations at every node. 
	# When fitting the logistic regression functions at a node, LMT has to determine the number of LogitBoost iterations to run. 
	# Originally, this number was cross-validated at every node in the tree. To save time, this heuristic cross-validates the number 
	# only once and then uses that number at every node in the tree. Usually this does not decrease accuracy but improves runtime considerably.
	java -cp ${2} weka.classifiers.meta.MultiScheme -X 0 -S 1 -B "weka.classifiers.trees.LMT -I -1 -M 15 -W 0.0" -t ${DS_ARFF} -split-percentage 70 -s ${j} -i -o -k > ${DS}_stat${j}.csv
	# This is an alternative command if training and test instances are into different files called training<1,2,...>.arff and test<1,2,...>.arff
	# java -cp ./weka.jar weka.classifiers.meta.MultiScheme -X 0 -S 1 -B "weka.classifiers.trees.LMT -I -1 -M 15 -W 0.0" -t train$j.arff -T test$j.arff -i -o -k > stat$j.csv
	# End of the run "j"
	echo ${j}
	# End of the loop 
done
