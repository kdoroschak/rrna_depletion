import numpy as np
from sklearn import neighbors
from sklearn.metrics import mean_squared_error

def kfold(rrna_labels, pair_probs, save_folder="./", n_partitions=10):
	pair_probs = np.array(pair_probs)
	rrna_labels = np.array(rrna_labels)
	zipped = zip(rrna_labels, pair_probs)
	# print "len:", len(zipped_struct)
	# np.random.shuffle(zipped_struct)
	# row_names, rrna_labels, mfe, pair_probs = zip(*zip(zipped_struct))
	rrna_indices = [i for i, j in enumerate(rrna_labels) if j == 1]
	non_rrna_indices = [i for i, j in enumerate(rrna_labels) if j == 0]
	# rrna = rrna_labels == 1
	# non_rrna = rrna_labels == 0

	num_rrna = len(rrna_indices)
	num_non_rrna = len(non_rrna_indices)
	# num_rrna = len(rrna_labels[rrna])
	# num_non_rrna = len(rrna_labels[non_rrna])
	dataset_size = max(num_rrna, num_non_rrna)
	# np.random.shuffle(zipped)

	rrna_random_indices = np.random.choice(rrna_indices, size=dataset_size, replace=True)
	rrna_pair_probs = pair_probs[rrna_random_indices,:]
	non_rrna_random_indices = np.random.choice(non_rrna_indices, size=dataset_size, replace=False)
	non_rrna_pair_probs = pair_probs[non_rrna_random_indices,:]

	print rrna_pair_probs.shape

	for k in range(n_partitions):
		print "Creating train & test sets for partition", k
		train_x, train_y, test_x, test_y = get_cross_validation_set(k, n_partitions, rrna_pair_probs, non_rrna_pair_probs)
		print "Training the model for partition", k
		knc = neighbors.KNeighborsClassifier(n_neighbors=1, weights="distance", n_jobs=50)
		knc.fit(train_x, train_y)
		print "Predicting the labels for the test set of partition", k
		y_hat = knc.predict(test_x)
		print "Analyzing the test prediction performance for partition", k
		rmse = calc_rmse(test_y, y_hat)
		acc = knc.score(test_x, test_y)
		nz_predicted = np.count_nonzero(y_hat)
		nz_actual = np.count_nonzero(test_y)
		print "  rmse:", rmse
		print "  accuracy:", acc
		print "  predicted nonzeros:", nz_predicted
		print "  actual nonzeros:", nz_actual
		print "  total in sample:", len(test_y)

def calc_rmse(actual, predicted):
	rmse = np.sqrt(mean_squared_error(actual, predicted))
	return rmse

def get_cross_validation_set(k, nCrossVal, pos_data, neg_data):
	''' Break the data into a training and testing dataset based on the number of folds in k-fold xvalidation and which block k we are using as the test set.
	k -> fold number
	nCrossVal -> the number of folds
	x -> the data to split
	y -> the values we are trying to predict, to split
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

	print train_x.shape, train_y.shape, test_x.shape, test_y.shape
	return train_x, train_y, test_x, test_y