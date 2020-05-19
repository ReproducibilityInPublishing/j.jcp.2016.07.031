tolerance = 2e-5

# For the oracle, check only that the final column (the error) matches
for i in [1,2,3,4]:
    with open('table'+str(i)+'_short.csv', 'r') as f:
        actual_data = f.readlines()
    with open('expected_output/table'+str(i)+'_short.csv', 'r') as f:
        expected_data = f.readlines()

    if len(actual_data) != len(expected_data):
        print("Actual data and Expected data for experiment "+str(i)+" don't contain the same experiments!")

    cur_algo = None
    algo_row = 0
    for l_idx in range(len(actual_data)):
        if ',' not in actual_data[l_idx]:
            cur_algo = actual_data[l_idx].strip()
            algo_row = 0
        else:
            actual = float(actual_data[l_idx].strip().split(',')[-1])
            expected = float(expected_data[l_idx].strip().split(',')[-1])
            err = abs((actual-expected)/actual)
            if (err < tolerance):
                print("Table "+str(i)+" Algorithm "+str(cur_algo)+" Experiment "+str(algo_row)+" PASS")
            else:
                print(tolerance)
                print(err)
                print("Table "+str(i)+" Algorithm "+str(cur_algo)+" Experiment "+str(algo_row)+" FAIL")
        algo_row += 1
