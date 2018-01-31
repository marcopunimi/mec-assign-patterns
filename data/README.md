Input data structure:
* 10 facilities
* 1419 access points
* 4 time-slots cardinalities: 96, 48, 24 and 12

Three file types:
* distance_matrix_facilities.txt = distance matrix facilities to facilities, comma separated values
* distance_matrix_BS_facilities.txt = distance matrix access points to facilities, comma separated values; one line for each access point, one column for each facility
* demand_[X]t_synth_[Y].txt = traffic matrix, comma separated values, one line for each access point, one column for each time-slot. One file for each time-slot cardinality X (i.e. 96, 48, 24 and 12); five instances for each file, i.e. Y in [1..5]
