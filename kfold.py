import numpy as np
from sklearn import neighbors # using KNeighborsRegressor, KNeighborsClassifier
from sklearn.metrics import mean_squared_error, confusion_matrix
from datetime import datetime

def kfold(rrna_data, non_rrna_data, save_folder="./", n_partitions=10, sampling="over", load_crossval_partitions=False, filter_zero_rows=False):
	results_log_file = save_folder + "/knn_results_" + datetime.now().strftime("%y-%m-%d_%H:%M") + ".txt"
	results_log = open(results_log_file, "w")
	print >>results_log, "Parameters for this run:"
	print >>results_log, "Number of partitions: ", str(n_partitions)
	print >>results_log, "Sampling style: ", sampling
	print >>results_log, "Filtering zero rows? ", str(filter_zero_rows)
	print >>results_log, "\n" 

	# print np.mean(rrna_data), np.std(rrna_data)
	# print np.mean(non_rrna_data), np.std(non_rrna_data)

	if filter_zero_rows:
		print "Before filtering, there were:"
		print str(rrna_data.shape[0]) + " samples in rrna data"
		print str(non_rrna_data.shape[0]) + " samples in non rrna data"

		t = 0.001
		# print np.sum(rrna_data>=t)/50.
		# print np.sum(non_rrna_data>=t)/50.
		#TODO do the actual filtering
		rrna_data = rrna_data[np.all(rrna_data>t,axis=1)]
		non_rrna_data = non_rrna_data[np.all(non_rrna_data>t,axis=1)]

		print "After filtering, there are:"
		print str(rrna_data.shape[0]) + " samples in rrna data"
		print str(non_rrna_data.shape[0]) + " samples in non rrna data"

	data = np.vstack([rrna_data, non_rrna_data])
	rrna_labels = np.hstack([np.ones(rrna_data.shape[0]), np.zeros(non_rrna_data.shape[0])])

	rrna_indices = [i for i, j in enumerate(rrna_labels) if j == 1]
	non_rrna_indices = [i for i, j in enumerate(rrna_labels) if j == 0]

	num_rrna = len(rrna_indices)
	num_non_rrna = len(non_rrna_indices)

	# Balance the classes by over or undersampling
	if sampling == "over":
		# Oversample the smaller class with replacement
		dataset_size = max(num_rrna, num_non_rrna)
		if dataset_size == num_non_rrna:
			rrna_random_indices = np.random.choice(rrna_indices, size=dataset_size, replace=True)
			rrna_data = data[rrna_random_indices,:]
			non_rrna_random_indices = np.random.choice(non_rrna_indices, size=dataset_size, replace=False)
			non_rrna_data = data[non_rrna_random_indices,:]
		elif dataset_size == num_rrna:
			rrna_random_indices = np.random.choice(rrna_indices, size=dataset_size, replace=False)
			rrna_data = data[rrna_random_indices,:]
			non_rrna_random_indices = np.random.choice(non_rrna_indices, size=dataset_size, replace=True)
			non_rrna_data = data[non_rrna_random_indices,:]
	elif sampling == "under":
		# Downsample with no replacement
		dataset_size = min(num_rrna, num_non_rrna)
		rrna_random_indices = np.random.choice(rrna_indices, size=dataset_size, replace=False)
		rrna_data = data[rrna_random_indices,:]
		non_rrna_random_indices = np.random.choice(non_rrna_indices, size=dataset_size, replace=False)
		non_rrna_data = data[non_rrna_random_indices,:]
	

	# Partition the now-balanced datasets into test/train for cross-validation
	# As each class is partitioned, save it and run the classifier
	for k in range(n_partitions):
		save_train_x = save_folder.rstrip("/") + "/knn_kfold-" + str(k).zfill(2) + "of" + str(n_partitions).zfill(2) + "_" + sampling + "sampled_test_data.npy"
		save_train_y = save_folder.rstrip("/") + "/knn_kfold-" + str(k).zfill(2) + "of" + str(n_partitions).zfill(2) + "_" + sampling + "sampled_test_labels.npy"
		save_test_x = save_folder.rstrip("/") + "/knn_kfold-" + str(k).zfill(2) + "of" + str(n_partitions).zfill(2) + "_" + sampling + "sampled_train_data.npy"
		save_test_y = save_folder.rstrip("/") + "/knn_kfold-" + str(k).zfill(2) + "of" + str(n_partitions).zfill(2) + "_" + sampling + "sampled_train_labels.npy"
		# print save_train_x
		# print save_train_y
		# print save_test_x
		# print save_test_y

		if load_crossval_partitions:
			try:
				print "Loading train & test sets for partition", k
				train_x = np.load(save_train_x)
				train_y = np.load(save_train_y)
				test_x = np.load(save_test_x)
				test_y = np.load(save_test_y)
			except:
				print "Error: files not found. Creating test & train sets instead for partition", k
				train_x, train_y, test_x, test_y = get_cross_validation_set(k, n_partitions, rrna_data, non_rrna_data)
		else:
			print "Creating train & test sets for partition", k
			train_x, train_y, test_x, test_y = get_cross_validation_set(k, n_partitions, rrna_data, non_rrna_data)

		print "Training the model for partition", k
		knc = neighbors.KNeighborsClassifier(n_neighbors=1, weights="distance", n_jobs=50, algorithm="ball_tree")
		knc.fit(train_x, train_y)
		print "Predicting the labels for the training set of partition", k
		y_hat_train = knc.predict(train_x)
		y_hat_train_file = save_folder.rstrip("/") + "/knn_kfold-" + str(k).zfill(2) + "of" + str(n_partitions).zfill(2) + "_" + sampling + "predicted_train_labels.npy"
		np.save(y_hat_train_file, y_hat_train)
		print "Predicting the labels for the test set of partition", k
		y_hat_test = knc.predict(test_x)
		y_hat_test_file = save_folder.rstrip("/") + "/knn_kfold-" + str(k).zfill(2) + "of" + str(n_partitions).zfill(2) + "_" + sampling + "predicted_test_labels.npy"
		np.save(y_hat_test_file, y_hat_test)

		print "Analyzing the prediction performance for partition", k
		rmse_train = calc_rmse(train_y, y_hat_train)
		acc_train = knc.score(train_x, train_y)
		nz_predicted_train = np.count_nonzero(y_hat_train)
		nz_actual_train = np.count_nonzero(train_y)
		confusion_train = confusion_matrix(train_y, y_hat_train)
		rmse_test = calc_rmse(test_y, y_hat_test)
		acc_test = knc.score(test_x, test_y)
		nz_predicted_test = np.count_nonzero(y_hat_test)
		nz_actual_test = np.count_nonzero(test_y)
		confusion_test = confusion_matrix(test_y, y_hat_test)
		print >>results_log, "# ----------------------------------"
		print >>results_log, "### Results for k-fold partition " + str(k).zfill(2) + ":"
		print >>results_log, "# ----------------------------------"
		print >>results_log, "Training set: " + str(train_y.shape[0]) + " samples"
		print >>results_log, "Testing set:  " + str(test_y.shape[0]) + " samples"
		print >>results_log, "Sets saved at: " +  save_folder.rstrip("/") + "/knn_kfold-" + str(k).zfill(2) + "of" + str(n_partitions).zfill(2) + "_" + sampling + "sampled_{test, train}_{data, labels}.npy"
		print >>results_log, ""
		print >>results_log, "Training set predictions:"
		print >>results_log, "RMSE: ", rmse_train
		print >>results_log, "Accuracy: ", acc_train
		print >>results_log, "Predicted # of rrnas:", nz_predicted_train
		print >>results_log, "Actual # of rrnas:", nz_actual_train
		print >>results_log, "Confusion matrix (rRNA is pos. case):\n", confusion_train
		print >>results_log, ""
		print >>results_log, "Testing set predictions:"
		print >>results_log, "RMSE: ", rmse_test
		print >>results_log, "Accuracy: ", acc_test
		print >>results_log, "Predicted # of rrnas:", nz_predicted_test
		print >>results_log, "Actual # of rrnas:", nz_actual_test
		print >>results_log, "Confusion matrix (rRNA is pos. case):\n", confusion_test
		print >>results_log, "\n"
		results_log.flush()

		# print "  rmse:", rmse
		# print "  accuracy:", acc
		# print "  predicted nonzeros:", nz_predicted
		# print "  actual nonzeros:", nz_actual
		# print "  total in sample:", len(test_y)
		# print "  confusion matrix:\n", confusion
		# print ""
	results_log.close()

def calc_rmse(actual, predicted):
	rmse = np.sqrt(mean_squared_error(actual, predicted))
	return rmse

def get_cross_validation_set(k, nCrossVal, pos_data, neg_data):
	''' Break the data into a training and testing dataset based on the number of folds in k-fold xvalidation and which block k we are using as the test set.
	k -> fold number
	nCrossVal -> the number of folds
	'''
	if nCrossVal == 0:
		return x, y, None, None
	N,d = np.shape(pos_data)
	testSetSize = np.floor_divide(N,nCrossVal)

	train_y = np.hstack([np.ones(N - testSetSize), np.zeros(N - testSetSize)])
	test_y = np.hstack([np.ones(testSetSize), np.zeros(testSetSize)])

	if k == 0:
		train_x = np.vstack([pos_data[testSetSize:,:], neg_data[testSetSize:,:]])
		test_x = np.vstack([pos_data[:testSetSize,:], neg_data[:testSetSize,:]])
			
	elif k == nCrossVal-1:
		train_x = np.vstack([pos_data[:k*testSetSize,:], neg_data[:k*testSetSize,:]])
		test_x = np.vstack([pos_data[k*testSetSize:,:], neg_data[k*testSetSize:,:]])
		train_y = np.hstack([np.ones(train_x.shape[0]/2), np.zeros(train_x.shape[0]/2)])
		test_y = np.hstack([np.ones(test_x.shape[0]/2), np.zeros(test_x.shape[0]/2)])

	else:
		train_x = np.vstack([pos_data[:k*testSetSize,:], 
					   		 pos_data[(k+1)*testSetSize:,:], 
					   		 neg_data[:k*testSetSize,:], 
			            	 neg_data[(k+1)*testSetSize:,:]]
			            	 )
		test_x = np.vstack([pos_data[k*testSetSize:(k+1)*testSetSize],
						   neg_data[k*testSetSize:(k+1)*testSetSize]
						   ])

	print "Training set: " + str(train_y.shape[0]) + " samples"
	print "Testing set:  " + str(test_y.shape[0]) + " samples"
	# print train_x.shape, train_y.shape, test_x.shape, test_y.shape
	return train_x, train_y, test_x, test_y