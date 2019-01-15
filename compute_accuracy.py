import csv

def read_tables(file_name):
	list_2vals = []
	with open(file_name) as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			if ( row == ['AIMGM'] or row == ['BD-ADI'] or row == ['BFSMGM']):
				list_2vals.append(row)
			else:
				list_2vals.append([float(row[11]), float(row[12])])
	return list_2vals

def calc_error(lexp, lactual):
	list_result = []
	for i in range(len(lexp)):
		if (len(lactual) <= i):
			list_result.append("NULL")
		elif ( lexp[i] == ['AIMGM'] or lexp[i] == ['BD-ADI'] or lexp[i] == ['BFSMGM']):
			list_result.append(lexp[i])
		elif ( lactual[i] == ['AIMGM'] or lactual[i] == ['BD-ADI'] or lactual[i] == ['BFSMGM']):
			list_result.append("NULL")
		else:
			percent1 = (abs(lactual[i][0]-lexp[i][0])/(lexp[i][0]))*100
			percent2 = (abs(lactual[i][1]-lexp[i][1])/(lexp[i][1]))*100

			list_result.append([percent1, percent2])
	return list_result

for i in range(1,5):
	print ("EXAMPLE"+str(i))
	path_expected = 'expected_tables/table%d.csv' % (i)
	lexp = read_tables(path_expected)
	path_actual = 'table%d.csv' % (i)
	lactual = read_tables(path_actual)
	print(calc_error(lexp, lactual))