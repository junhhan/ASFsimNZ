#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mtwister.h"
#include <omp.h>
#include <R.h>

// Basic parameters
#define _na 999999
#define _neighbours 8 // Maximum number of adjacent habitats
#define _weeksinyear 52 // Number of weeks in a year
#define _initialsimyear 3
#define _female 0
#define _male 1
#define _adult 0
#define _young 1
#define _cellsize_0 4.0
#define _cellsize_1 10.0

// For the shape of the grid
#define _ncell_NI_4km 105324
#define _ncell_SI_4km 109652
#define _ncell_NI_10km 42330
#define _ncell_SI_10km 44000
#define _ncol_NI_4km 262
#define _ncol_SI_4km 316
#define _ncol_NI_10km 166
#define _ncol_SI_10km 200

// For the cell information array
#define _cellinfocols 9
#define _col_cid 0
#define _col_x 1
#define _col_y 2
#define _col_island 3
#define _col_ishabitat 4
#define _col_isaccess 5
#define _col_isbarrier 6
#define _col_breedcap 7
#define _col_popcap 8

// For the demograhic and disease compartment
#define _g_demo 6 // 0: Female piglet, 1: Female yearling, 2: Sow, 3: Male piglet, 4: Male yearling, 5: Male boar
#define _fp 0
#define _fy 1
#define _sow 2
#define _mp 3
#define _my 4
#define _boar 5

#define _g_dz 5 // 0: Maternal antibody, 1: Susceptible, 2: Exposed, 3: Infectious, 4: Recovered
#define _m 0
#define _s 1
#define _e 2
#define _i 3
#define _r 4

// Parameters to calculate an infectious pressure
#define _ncol_proxy_4km 177 // Maxiumum number of proximal habitats in 4km2 grid
#define _ncol_proxy_10km 69 // Maxiumum number of proximal habitats in 10km2 grid
#define _beta_pig 0.016
#define _beta_car 0.008
#define _pi_between 0.624
#define _decay_distance 0.650
#define _nocid 0

#define _beta 0.74
#define _beta_pig_neighbour 0.4
#define _beta_car_neighbour 0.2
#define _decayalpha_pig -0.80
#define _decaybeta_pig -1.91
#define _decayalpha_car 2.32
#define _decaybeta_car -6.37

// For the demographic event array in "_Simulation" function
#define _col_breed 0
#define _col_wean 1
#define _w_update 12 // Week to updated demographic event weeks

// For the demographic setting in "_Simulation" function
#define _mean_breedweek 38 // Distribution of breeding week (DBW) follows a Poisson distribution
#define _weanperiod 8 // Weaning period
#define _mean_litter 8.0 // Average litter size. It follows a normal distribution.
#define _sd_litter 2.0 // Standard devidaiton of litter size.

// Default mortality and hunting rate
#define _pmort_young 0.50
#define _pmort_adult 0.25
#define _phunt 0.1
#define _pdeath 0.99
#define _pdeath_alt 0.95
#define _w_mab 12

// Migration related parameters
#define _p_road 0.35 // Probability of crossing road
#define _p_river 0.40 // Probability of crossing river
#define _p_fence 0.50 // Probability of crossing fence
#define _p_direct 0.5 // Probability of keeping the previous direction for migration
#define _migratelimit_male 50.0 // Distance limit of dispersal for male yearlings
#define _migratelimit_female 6.0 // Distance limit of group splitting for female yearlings

// ASF intoduction locations
#define _auckland_4km 36767
#define _auckland_10km 14829
#define _wellington_4km 101755
#define _wellington_10km 40733
#define _picton_4km 13258
#define _picton_10km 5391
#define _dunedin_4km 93692
#define _dunedin_10km 37499

// Control related parameters
#define _zonetypes 3 // Infected, buffer, surveillannce zones
#define _z_inf 0
#define _z_buff 1
#define _z_surv 2
#define _controltypes 5
#define _c_carc 0
#define _c_hunt 1
#define _c_surv 2
#define _c_fence 3
#define _c_duration 4
#define _psurv 0.1
#define _pcarc 0.1

// USERS SHOULD BE ABLE TO CHANGE
// Enhanced hunting rate in a control zone
#define _psurv_inc 1.0
#define _phunt_inc 0.4
#define _pcarc_inc 0.4


// Functions to read and write
void _Convert_file(int ncell, int ncol_proxy, int pop, int *input_info, int *input_roadcross, int *input_rivercross, int *input_proxycid, double *input_proxydist, int *input_proxyroad, int *input_proxyriver,
	int **cell_info, int **road_cross, int **river_cross, int **proxy_cid, double **proxy_dist, int **proxy_road, int **proxy_river);

// Functions for randomness
int _Random_binomial(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, int n, double p);
double _Random_normal(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, double mu, double sd);
int _Random_Poisson(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, double mean);

// Functions for migration
double _Cross_probability(int scale, int island, int **road_cross, int **river_cross, int **control, int org, int dst);
int _Random_neighbour(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt,
	int scale, int island, int **cell_info, int **road_cross, int **river_cross, int **control, int org);
int _Migration(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt,
	int scale, int island, int **cell_info, int **population, int **road_cross, int **river_cross, int **control, int org, int sex);

// Functions for indexing cells
int _Neighbour_index(int scale, int island, int org, int dst);
int _Neighbour_cellid(int scale, int island, int org, int index);

// Functions for the main simulation
void _Initial_population(int ncell, int **cell_info, int **population);
double _Distance(int **cell_info, int i, int j);
double _Calculate_lambda(int n_col_adj, int **population, int **carcass, int **proxy_cid, double **proxy_dist, int **proxy_road, int **proxy_river, int **control, int org);
double _Calculate_lambda_Pepin(int n_col_adj, int **population, int **carcass, int **proxy_cid, double **proxy_dist, int **proxy_road, int **proxy_river, int **control, int org);
double _Calculate_lambda_Thulke(int scale, int island, int **population, int **carcass, int org);
void _Update_control(int scale, int island, int **cell_info, int *control_param, double *dist_param, int *fenceid, int **control, int org);
void _Simulation(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt,
	int **population, int **carcass, 
	int **cell_info, int **road_cross, int **river_cross, int **proxy_cid, double **proxy_dist, int **proxy_road, int **proxy_river, 
	int simyear, int scale, int island, int introloc, int control_index, int foi,
	int *control_param, double *dist_param, int *fenceid, int iter, int **result_map, int **result_ts);


// Main function
void _ASF_simulation(int *input_info, int *input_roadcross, int *input_rivercross, int *input_proxycid, double *input_proxydist, int *input_proxyroad, int *input_proxyriver,
	int *_coreid, int *_rep, int *_simyear, int *_pop, int *_scale, int *_island, int *_introloc, int *_control_index, int *_foi, int *_result_map, int *_result_ts) {
	
	int coreid, rep, simyear, pop, scale, island, introloc, control_index, foi;
	coreid= (*_coreid);
	rep= (*_rep);
	simyear= (*_simyear);
	pop= (*_pop);
	scale= (*_scale);
	island= (*_island);
	introloc= (*_introloc);
	control_index= (*_control_index);
	foi= (*_foi);
	
	int i, j, k;
	
	// Control scenarios
	double *dist_param;
	dist_param= calloc(_zonetypes, sizeof(double));
	dist_param[_z_inf]= 4.0;
	dist_param[_z_buff]= 8.0;
	dist_param[_z_surv]= 15.0;

	int *control_param;
	control_param= calloc((_zonetypes * _controltypes), sizeof(int));
	if (control_index == 0) {
		control_param[_controltypes * _z_inf + _c_carc]= 0;
		control_param[_controltypes * _z_inf + _c_hunt]= 0;
		control_param[_controltypes * _z_inf + _c_surv]= 0;
		control_param[_controltypes * _z_inf + _c_fence]= 0;
		control_param[_controltypes * _z_inf + _c_duration]= 0;
		control_param[_controltypes * _z_buff + _c_carc]= 0;
		control_param[_controltypes * _z_buff + _c_hunt]= 0;
		control_param[_controltypes * _z_buff + _c_surv]= 0;
		control_param[_controltypes * _z_buff + _c_fence]= 0;
		control_param[_controltypes * _z_buff + _c_duration]= 0;
		control_param[_controltypes * _z_surv + _c_carc]= 0;
		control_param[_controltypes * _z_surv + _c_hunt]= 0;
		control_param[_controltypes * _z_surv + _c_surv]= 0;
		control_param[_controltypes * _z_surv + _c_fence]= 0;
		control_param[_controltypes * _z_surv + _c_duration]= 0;
	} else {
		control_param[_controltypes * _z_inf + _c_carc]= 1;
		control_param[_controltypes * _z_inf + _c_hunt]= 1;
		control_param[_controltypes * _z_inf + _c_surv]= 1;
		control_param[_controltypes * _z_inf + _c_fence]= 0;
		control_param[_controltypes * _z_inf + _c_duration]= 8;
		control_param[_controltypes * _z_buff + _c_carc]= 0;
		control_param[_controltypes * _z_buff + _c_hunt]= 0;
		control_param[_controltypes * _z_buff + _c_surv]= 1;
		control_param[_controltypes * _z_buff + _c_fence]= 0;
		control_param[_controltypes * _z_buff + _c_duration]= 4;
		control_param[_controltypes * _z_surv + _c_carc]= 0;
		control_param[_controltypes * _z_surv + _c_hunt]= 0;
		control_param[_controltypes * _z_surv + _c_surv]= 0;
		control_param[_controltypes * _z_surv + _c_fence]= 0;
		control_param[_controltypes * _z_surv + _c_duration]= 4;
	}
	if (control_index == 2) {
		dist_param[_z_inf]= 20.0;
		dist_param[_z_buff]= 28.0;
		dist_param[_z_surv]= 28.0;
	}
	if (control_index == 3) {
		control_param[_controltypes * _z_buff + _c_carc]= 1;
	}
	if (control_index == 4) {
		control_param[_controltypes * _z_surv + _c_surv]= 1;
	}
	int _fenceid, *fenceid;
	_fenceid= 1;
	fenceid= &_fenceid;
		
	// Model preparation
	int ncell, ncol_proxy;
	if (scale == 0) {
		ncol_proxy= _ncol_proxy_4km;
		if (island == 0) {
			ncell= _ncell_NI_4km;
		} else {
			ncell= _ncell_SI_4km;
		}
	} else {
		ncol_proxy= _ncol_proxy_10km;
		if (island == 0) {
			ncell= _ncell_NI_10km;
		} else {
			ncell= _ncell_SI_10km;
		}
	}
	
	// Create arrays to store information
	int **cell_info, **road_cross, **river_cross;
	cell_info= calloc(ncell, sizeof(int *));
	road_cross= calloc(ncell, sizeof(int *));
	river_cross= calloc(ncell, sizeof(int *));
	int **proxy_cid, **proxy_road, **proxy_river;
	double **proxy_dist;
	proxy_cid= calloc(ncell, sizeof(int *));
	proxy_dist= calloc(ncell, sizeof(double *));
	proxy_road= calloc(ncell, sizeof(int *));
	proxy_river= calloc(ncell, sizeof(int *));
	int **result_map, **result_ts;
	result_map= calloc(ncell, sizeof(int *));
	result_ts= calloc(rep, sizeof(int *));

	for (i= 0; i < ncell; ++i) {
		cell_info[i]= calloc(_cellinfocols, sizeof(int));
		road_cross[i]= calloc(_neighbours, sizeof(int));
		river_cross[i]= calloc(_neighbours, sizeof(int));
		proxy_cid[i]= calloc(ncol_proxy, sizeof(int));
		proxy_dist[i]= calloc(ncol_proxy, sizeof(double));
		proxy_road[i]= calloc(ncol_proxy, sizeof(int));
		proxy_river[i]= calloc(ncol_proxy, sizeof(int));
		result_map[i]= calloc(simyear * _weeksinyear, sizeof(int));
	}
	for (i= 0; i < rep; ++i) {
		result_ts[i]= calloc(simyear * _weeksinyear, sizeof(int));
	}

	_Convert_file(ncell, ncol_proxy, pop, input_info, input_roadcross, input_rivercross, input_proxycid, input_proxydist, input_proxyroad, input_proxyriver, 
		cell_info, road_cross, river_cross, proxy_cid, proxy_dist, proxy_road, proxy_river);

	int **population, **carcass;
	population= calloc(ncell, sizeof(int *));
	carcass= calloc(ncell, sizeof(int *));
	for (i= 0; i < ncell; ++i) {
		population[i]= calloc(_g_demo * _g_dz, sizeof(int));
		carcass[i]= calloc(_g_dz, sizeof(int));
	}

	int _MTcnt, *MTcnt; _MTcnt= 0; MTcnt= &_MTcnt;
	int _MTleft, *MTleft; _MTleft= 1; MTleft= &_MTleft; 
	int _MTinitf, *MTinitf; _MTinitf= 0; MTinitf= &_MTinitf;
	unsigned long MTstate[MTn], MTnext[MTn];
	
	for (i= 0; i < rep; ++i) {
		j= coreid * rep + i;
		_MTinit(MTstate, MTleft, MTinitf, (j + 1));

		for (j= 0; j < ncell; ++j) {
			for (k= 0; k < (_g_demo * _g_dz); ++k) {
				population[j][k]= 0;
			}
			for (k= 0; k < _g_dz; ++k) {
				carcass[j][k]= 0;
			}
		}

		_Initial_population(ncell, cell_info, population);
		_Simulation(MTstate, MTnext, MTleft, MTcnt, population, carcass, cell_info, road_cross, river_cross, proxy_cid, proxy_dist, proxy_road, proxy_river, 
			simyear, scale, island, introloc, control_index, foi, control_param, dist_param, fenceid, i, result_map, result_ts);
	}

	// Exporting the result
	for (i= 0; i < (simyear * _weeksinyear); ++i) {
		for (j= 0; j < ncell; ++j) {
			k= ncell * i + j;
			_result_map[k]= result_map[j][i];
		}
		for (j= 0; j < rep; ++j) {
			k= rep * i + j;
			_result_ts[k]= result_ts[j][i];
		}
	}

	// Free memory
	for (i= 0; i < rep; ++i) {
		free(result_ts[i]);
	}
	free(result_ts);
	
	for (i= 0; i < ncell; ++i) {
		free(population[i]);
		free(carcass[i]);
		free(cell_info[i]);
		free(road_cross[i]);
		free(river_cross[i]);
		free(proxy_cid[i]);
		free(proxy_dist[i]);
		free(proxy_road[i]);
		free(proxy_river[i]);
		free(result_map[i]);
	}
	free(population);
	free(carcass);
	free(cell_info);
	free(road_cross);
	free(river_cross);
	free(proxy_cid);
	free(proxy_dist);
	free(proxy_road);
	free(proxy_river);
	free(result_map);

	free(control_param);
	free(dist_param);
}


// Functions in detail
// Random binomial
int _Random_binomial(unsigned long MTstate[], unsigned long MTnext[], int *MTleft , int *MTcnt, int n, double p) {
    int index;
	int n_pos= 0;
	int i= 0;
	if (p > 1) {p= 1.0;}
	if (p < 0) {p= 0.0;}

	double x;

	if (p == 0.0 || n == 0) {
		return 0;
	} else if (p == 1.0) {
		return n;
	} else {
	    do {
		    do {
				x= _Rand(MTstate, MTnext, MTleft, MTcnt);
			} while (x > 1.0 || x < 0.0);

			if(p >= x) {index= 1;} else {index= 0;}
			n_pos += index;
			i++;
	    } while (i < n);
	
		return n_pos;
	}
}

// Random number from normal distribution: It retunrs a random number according to a normal distribution ("mu", "sigma").
double _Random_normal(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, double mu, double sd) {
	double two_pi= 2.0 * acos(-1.0); // acos(-1) is PI
	double u, v, z;

	do {
		u= _Rand(MTstate, MTnext, MTleft, MTcnt);
		v= _Rand(MTstate, MTnext, MTleft, MTcnt);
	} while (u == 0.0);

	z= sqrt(-2.0 * log(u)) * cos(two_pi * v);

	return (mu + z * sd);
}

// Random number generator based on Poisson distribution
int _Random_Poisson(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, double mean) {
	int n= 0; // counter of iteration
	double limit; 
	double x; // pseudo random number

	limit= exp(-mean);
	x= _Rand(MTstate, MTnext, MTleft, MTcnt);

	while (x > limit) {
		n++;
		x *= _Rand(MTstate, MTnext, MTleft, MTcnt);
	}

	return n;
}

// Fuction to identify an index of a neighbouring cell based on the cell id
int _Neighbour_index(int scale, int island, int org, int dst) {
	int ncol, output;
	if (scale == 0) { // 4km
		if (island == 0) { // North Island
			ncol= _ncol_NI_4km;
		} else { // South Island
			ncol= _ncol_SI_4km;
		}
	} else { // 10 km
		if (island == 0) { // North Island
			ncol= _ncol_NI_10km;
		} else { // South Island
			ncol= _ncol_SI_10km;
		}
	}
	
	if (org % ncol == 0) { // West edge
		if (dst == (org - ncol)) {output= 1;}
		else if (dst == (org - ncol + 1)) {output= 2;}
		else if (dst == (org + 1)) {output= 4;}
		else if (dst == (org + ncol)) {output= 6;}
		else if (dst == (org + ncol + 1)) {output= 7;}
		else {printf("%5i and %5i are not neighbours.\n", org, dst); output= _na;}
	} else if (org % ncol == (ncol - 1)) { // East edge
		if (dst == (org - ncol - 1)) {output= 0;}
		else if (dst == (org - ncol)) {output= 1;}
		else if (dst == (org - 1)) {output= 3;}
		else if (dst == (org + ncol - 1)) {output= 5;}
		else if (dst == (org + ncol)) {output= 6;}
		else {printf("%5i and %5i are not neighbours.\n", org, dst); output= _na;}
	} else {
		if (dst == (org - ncol - 1)) {output= 0;}
		else if (dst == (org - ncol)) {output= 1;}
		else if (dst == (org - ncol + 1)) {output= 2;}
		else if (dst == (org - 1)) {output= 3;}
		else if (dst == (org + 1)) {output= 4;}
		else if (dst == (org + ncol - 1)) {output= 5;}
		else if (dst == (org + ncol)) {output= 6;}
		else if (dst == (org + ncol + 1)) {output= 7;}
		else {printf("%5i and %5i are not neighbours.\n", org, dst); output= _na;}
	}
	
	return output;
}

// Fuction to identify a cell id of a neighbouring cell based on the neighbouring index 
int _Neighbour_cellid(int scale, int island, int org, int index) {
	int ncol, ncell, output;
	if (scale == 0) { // 4km
		if (island == 0) { // North Island
			ncol= _ncol_NI_4km;
			ncell= _ncell_NI_4km ;
		} else { // South Island
			ncol= _ncol_SI_4km;
			ncell= _ncell_SI_4km ;
		}
	} else { // 10 km
		if (island == 0) { // North Island
			ncol= _ncol_NI_10km;
			ncell= _ncell_NI_10km ;
		} else { // South Island
			ncol= _ncol_SI_10km;
			ncell= _ncell_SI_10km ;
		}
	}
	
	if (org % ncol == 0) { // West edge
		if (index == 1) {output= (org - ncol);}
		else if (index == 2) {output= (org - ncol + 1);}
		else if (index == 4) {output= (org + 1);}
		else if (index == 6) {output= (org + ncol);}
		else if (index == 7) {output= (org + ncol + 1);}
		else {output= _na;}
//		else {printf("There is no %1i th neighbour for %5i.\n", index + 1, org); output= _na;}
	} else if (org % ncol == (ncol - 1)) { // East edge
		if (index == 0) {output= (org - ncol - 1);}
		else if (index == 1) {output= (org - ncol);}
		else if (index == 3) {output= (org - 1);}
		else if (index == 5) {output= (org + ncol - 1);}
		else if (index == 6) {output= (org + ncol);}
		else {output= _na;}
//		else {printf("There is no %1i th neighbour for %5i.\n", index + 1, org); output= _na;}
	} else {
		if (index == 0) {output= (org - ncol - 1);}
		else if (index == 1) {output= (org - ncol);}
		else if (index == 2) {output= (org - ncol + 1);}
		else if (index == 3) {output= (org - 1);}
		else if (index == 4) {output= (org + 1);}
		else if (index == 5) {output= (org + ncol - 1);}
		else if (index == 6) {output= (org + ncol);}
		else if (index == 7) {output= (org + ncol + 1);}
		else {output= _na;}
//		else {printf("There is no %1i th neighbour for %5i.\n", index + 1, org); output= _na;}
	}
	
	if (output < 0 || output >= ncell) {output= _na;}
	return output;
}

// Function to calculate the probability of feral pigs migrating between cells depending on the presence of roads and rivers (and fence in later!)
double _Cross_probability(int scale, int island, int **road_cross, int **river_cross, int **control, int org, int dst) {
	int _Neighbour_index(int scale, int island, int org, int dst);
	
	int i, val;
	double x= 1.0;
	
	i= _Neighbour_index(scale, island, org, dst);
	if (i != _na) {
		val= road_cross[org][i];
		if (val != 0) {x *= _p_road;}
		val= river_cross[org][i];
		if (val != 0) {x *= _p_river;}
		if (control[org][_c_fence] != control[dst][_c_fence]) {x *= _p_fence;}
	} else {
		x= 0.0;
	}
	
	return x;
}

// Function to choose a neighbouring cell for migration given the presence of roads and rivers (and fence in later!)
int _Random_neighbour(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt,
	int scale, int island, int **cell_info, int **road_cross, int **river_cross, int **control, int org) {
	double _Cross_probability(int scale, int island, int **road_cross, int **river_cross, int **control, int org, int dst);
	int _Neighbour_index(int scale, int island, int org, int dst);
	int _Neighbour_cellid(int scale, int island, int org, int index);
	
	int i, dst;
	double r, *seq;
	seq= calloc(_neighbours, sizeof(double));
	
	for (i= 0; i < _neighbours; ++i) {
		dst= _Neighbour_cellid(scale, island, org, i);
		if (dst == _na) {
			seq[i]= 0.0;
		} else {
			if (cell_info[dst][_col_island] == 1 && cell_info[dst][_col_isbarrier] == 0) {
				seq[i]= _Cross_probability(scale, island, road_cross, river_cross, control, org, dst);
			} else {
				seq[i]= 0.0;
			}
		}
		if (i != 0) {seq[i] += seq[i-1];}
	}
	
	// In case of no neighbouring cells
	if (seq[(_neighbours-1)] == 0.0) {
		dst= _na;
	// In case of having neighbouring cells
	} else {
		for (i= 0; i < (_neighbours-1); ++i) {
			seq[i] /= seq[(_neighbours-1)];
		}
		
		do {
			r= _Rand(MTstate, MTnext, MTleft, MTcnt);
		} while (r == 0.0 || r == 1.0);
		
		i= -1;
		do {
			++i;
			if (i == (_neighbours-1)) {break;}
		} while (r >= seq[i]);
		dst= _Neighbour_cellid(scale, island, org, i);
		if (dst == _na) {printf("No way, this should not happen...\n");}
	}
	free(seq);
	return dst;
}

// Natal migration
int _Migration(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt,
	int scale, int island, int **cell_info, int **population, int **road_cross, int **river_cross, int **control, int org, int sex) {
	int _Random_neighbour(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt,
		int scale, int island, int **cell_info, int **road_cross, int **river_cross, int **control, int org);
	double _Cross_probability(int scale, int island, int **road_cross, int **river_cross, int **control, int org, int dst);
	int _Neighbour_index(int scale, int island, int org, int dst);
	int _Neighbour_cellid(int scale, int island, int org, int index);
	
	int i, j, n_lim, settle, goback, old_cid, new_cid, dir, tmp_total;
	double area, p;
	
	if (scale == 0) {
		area= _cellsize_0;
	} else {
		area= _cellsize_1;
	}
	
	if (sex == _female) {
		n_lim= (int)ceil(_migratelimit_female / sqrt(area));
	} else {
		n_lim= (int)ceil(_migratelimit_male / sqrt(area));
	}
	
	int *record;
	record= calloc(n_lim, sizeof(int));
	
	i= 0; settle= 0;
	do {
		if (i == 0) {old_cid= org;} else {old_cid= new_cid;}

		// Select a new cell
		if (i == 0) {
			new_cid= _Random_neighbour(MTstate, MTnext, MTleft, MTcnt, scale, island, cell_info, road_cross, river_cross, control, old_cid);
		} else {
			p= _Rand(MTstate, MTnext, MTleft, MTcnt);
			// Avoid the previous direction
			if (p <= _p_direct) {
				new_cid= _Random_neighbour(MTstate, MTnext, MTleft, MTcnt, scale, island, cell_info, road_cross, river_cross, control, old_cid);
			// Continue to the previous direction
			} else {
				new_cid= _Neighbour_cellid(scale, island, old_cid, dir);
				// If the new cell is neither habitatble nor accessible
				if (new_cid == _na) {
					new_cid= _Random_neighbour(MTstate, MTnext, MTleft, MTcnt, scale, island, cell_info, road_cross, river_cross, control, old_cid);
				} else if (cell_info[new_cid][_col_island] == 0 || cell_info[new_cid][_col_isbarrier] == 1) {
					new_cid= _Random_neighbour(MTstate, MTnext, MTleft, MTcnt, scale, island, cell_info, road_cross, river_cross, control, old_cid);
				}
			}
		}
		record[i]= new_cid;

		// Break the loop if there is no cell to migrate
		if (i == 0 && new_cid == _na) {
//			printf("Migration is impossible.\n");
			break;
		}
		
		// Compare with the passed cells
		goback= 0;
		if (i != 0) {
			if (new_cid == _na) {
				goback= 1;
				new_cid= old_cid;
				record[i]= old_cid;
			} else {
				for (j= 0; j < i; ++j) {
					if (new_cid == record[j]) {
						goback= 1;
						new_cid= old_cid;
						record[i]= old_cid;
						break;
					}
				}
			}
		}
		
		// Update the direction only if new cell has not been passed
		// This process should be done after comparing the passed cells
		// If the new cell has Not been passed
		if (goback == 0) {
			dir= _Neighbour_index(scale, island, old_cid, new_cid);
			// If the new cell is habitatble and able to settle?
			if (new_cid != _na && cell_info[new_cid][_col_island] == 1 && cell_info[new_cid][_col_ishabitat] == 1) {
				// Chechk the population of the new cell
				tmp_total= 0;
				if (sex == _female) {
					for (j= 0; j < _g_dz; ++j) {
						tmp_total += (population[new_cid][_g_demo*j+_fp] + population[new_cid][_g_demo*j+_fy] + population[new_cid][_g_demo*j+_sow]);
					}
					if (tmp_total == 0) {settle= 1;} else {settle= 0;}
				} else {
					for (j= 0; j < _g_dz; ++j) {
						tmp_total += (population[new_cid][_g_demo*j+_fp] + population[new_cid][_g_demo*j+_fy] + population[new_cid][_g_demo*j+_sow] + 
							population[new_cid][_g_demo*j+_mp] + population[new_cid][_g_demo*j+_my] + population[new_cid][_g_demo*j+_boar]);
					}
					if (tmp_total < cell_info[new_cid][_col_popcap]) {settle= 1;} else {settle= 0;}
				}
				if (settle == 1) {break;}
			}
		}
	
		++i;
		if (i == n_lim) {break;}
	} while (i < n_lim);
	free(record);
	if (settle == 1) {
		return new_cid;
	} else {
		return _na;
	}
}

// Distance between cells
double _Distance(int **cell_info, int i, int j) {
	double val;
	val= pow(pow(cell_info[i][_col_x] - cell_info[j][_col_x], 2.0) + pow(cell_info[i][_col_y] - cell_info[j][_col_y], 2.0), 0.5);
	val /= 1000.0;
	
	return val;
}

// Initial population
void _Initial_population(int ncell, int **cell_info, int **population) {
	int i, val;
	for (i= 0; i < ncell; ++i) {
		if (cell_info[i][_col_ishabitat] == 1) {
			val= cell_info[i][_col_breedcap];
			population[i][_g_demo * _s + _sow]= val;
		}
	}
}

// Calculate the infectious pressure of a cell
double _Calculate_lambda(int n_col_adj, int **population, int **carcass, int **proxy_cid, double **proxy_dist, int **proxy_road, int **proxy_river, int **control, int org) {
	int inf_pig, inf_car, i, j, freq_road, freq_river;
	double proxy_pig, proxy_car, dist, pi, lambda;
	
	inf_pig= population[org][_g_demo * _i + _fp] + population[org][_g_demo * _i + _fy] + population[org][_g_demo * _i + _sow] + 
		population[org][_g_demo * _i + _mp] + population[org][_g_demo * _i + _my] + population[org][_g_demo * _i + _boar];
	inf_car= carcass[org][_i];
	proxy_pig= 0.0;
	proxy_car= 0.0;

	// No adjacent cells
	if (proxy_cid[org][0] != _nocid) {
		i= 0;
		while (proxy_cid[org][i] != _nocid) {
			j = proxy_cid[org][i] - 1;
			dist= proxy_dist[org][i];
			freq_road= proxy_road[org][i];
			freq_river= proxy_river[org][i];

			pi= exp(-1.0 * _decay_distance * dist);
			if (freq_road != 0) {pi *= _p_road;}			
			if (freq_river != 0) {pi *= _p_river;}			
			if (control[org][_c_fence] != control[j][_c_fence]) {pi *= _p_fence;}

			proxy_pig += ((double)(population[j][_g_demo * _i + _fp] + population[j][_g_demo * _i + _fy] + population[j][_g_demo * _i + _sow] + 
				population[j][_g_demo * _i + _mp] + population[j][_g_demo * _i + _my] + population[j][_g_demo * _i + _boar]) * pi);
			proxy_car += ((double)carcass[j][_i] * pi);
			
			++i;
			if (i >= n_col_adj) {break;}
		}
	}

	lambda= 1.0 - exp(-1.0 * ((_beta_pig * (double)inf_pig) + (_beta_pig * _pi_between * proxy_pig) + (_beta_car * (double)inf_car) + (_beta_car * _pi_between * proxy_car)));
	
	return lambda;
}

// Calculate the infectious pressure of a cell
double _Calculate_lambda_Pepin(int n_col_adj, int **population, int **carcass, int **proxy_cid, double **proxy_dist, int **proxy_road, int **proxy_river, int **control, int org) {
	int inf_pig, inf_car, i, j;
	double proxy_pig, proxy_car, dist, pi_pig, pi_car, lambda;
	
	inf_pig= population[org][_g_demo * _i + _fp] + population[org][_g_demo * _i + _fy] + population[org][_g_demo * _i + _sow] + 
		population[org][_g_demo * _i + _mp] + population[org][_g_demo * _i + _my] + population[org][_g_demo * _i + _boar];
	inf_car= carcass[org][_i];
	proxy_pig= 0.0;
	proxy_car= 0.0;

	// No adjacent cells
	if (proxy_cid[org][0] != _nocid) {
		i= 0;
		while (proxy_cid[org][i] != _nocid) {
			j = proxy_cid[org][i] - 1;
			dist= proxy_dist[org][i];

			pi_pig= exp(_decayalpha_pig + (_decaybeta_pig * dist)) / (1.0 + exp(_decayalpha_pig + (_decaybeta_pig * dist)));
			pi_car= exp(_decayalpha_car + (_decaybeta_car * dist)) / (1.0 + exp(_decayalpha_car + (_decaybeta_car * dist)));

			proxy_pig += ((double)(population[j][_g_demo * _i + _fp] + population[j][_g_demo * _i + _fy] + population[j][_g_demo * _i + _sow] + 
				population[j][_g_demo * _i + _mp] + population[j][_g_demo * _i + _my] + population[j][_g_demo * _i + _boar]) * pi_pig);
			proxy_car += ((double)carcass[j][_i] * pi_car);
			
			++i;
			if (i >= n_col_adj) {break;}
		}
	}

	lambda= 1.0 - exp(-1.0 * ((_beta * (double)inf_pig) + (_beta * _beta_pig_neighbour * proxy_pig) + (_beta * (double)inf_car) + (_beta * _beta_car_neighbour * proxy_car)));
	
	return lambda;
}

// Calculate the infectious pressure of a cell
double _Calculate_lambda_Thulke(int scale, int island, int **population, int **carcass, int org) {
	int _Neighbour_cellid(int scale, int island, int org, int index);
	double p_direct, p_carcass, p_carcass_home, p_carcass_neighbour, lambda;
	double p_d= 0.05;
	double p_c= 0.03;

	int i, inf_pig, inf_car, dst;
	inf_pig= 0;
	for (i= 0; i < _g_demo; ++i) {
		inf_pig += population[org][_g_demo * _i + i];
	}
	p_direct= 1.0 - pow((1.0 - p_d), inf_pig);

	inf_car= carcass[org][_i];
	p_carcass_home= pow((1.0 - p_c), inf_car);

	inf_car= 0;
	for (i= 0; i < _neighbours; ++i) {
		dst= _Neighbour_cellid(scale, island, org, i);
		if (dst != _na) {
			inf_car += carcass[dst][_i];
		}
	}
	p_carcass_neighbour= pow((1.0 - p_c), inf_car);
	
	p_carcass= 1.0 - (p_carcass_home * p_carcass_neighbour);
	lambda= 1.0 - ((1.0 - p_direct) * (1.0 - p_carcass));
	
	return lambda;
}

// Update the control indicator
void _Update_control(int scale, int island, int **cell_info, int *control_param, double *dist_param, int *fenceid, int **control, int org) {
	double _Distance(int **cell_info, int i, int j);
	
	int i, ncol, ncell, nrow, row_org, min_row, max_row, min_cell, max_cell, val;
	double d_max;
	if (dist_param[_z_inf] <= dist_param[_z_buff]) {
		if (dist_param[_z_buff] <= dist_param[_z_surv]) {
			d_max= dist_param[_z_surv];
		} else {
			d_max= dist_param[_z_buff];
		}
	} else {
		if (dist_param[_z_inf] <= dist_param[_z_surv]) {
			d_max= dist_param[_z_surv];
		} else {
			d_max= dist_param[_z_inf];
		}
	}
	
	double dist;
	if (scale == 0) { // 4km
		nrow= (int)ceil(d_max / pow(_cellsize_0, 0.5));
		if (island == 0) { // North Island
			ncol= _ncol_NI_4km;
			ncell= _ncell_NI_4km;
		} else { // South Island
			ncol= _ncol_SI_4km;
			ncell= _ncell_SI_4km;
		}
	} else { // 10 km
		nrow= (int)ceil(d_max / pow(_cellsize_1, 0.5));
		if (island == 0) { // North Island
			ncol= _ncol_NI_10km;
			ncell= _ncell_NI_10km;
		} else { // South Island
			ncol= _ncol_SI_10km;
			ncell= _ncell_SI_10km;
		}
	}
	
	row_org= (int)floor(org / ncol);
	min_row= row_org - nrow; if (min_row < 0) {min_row= 0;}
	max_row= row_org + nrow; if (max_row >= (ncell / ncol)) {max_row= ((ncell / ncol) - 1);}
	
	min_cell= min_row * ncol;
	max_cell= max_row * ncol + (ncol - 1);
	
	for (i= min_cell; i <= max_cell; ++i) {
		dist= _Distance(cell_info, org, i);
		
		if (cell_info[i][_col_ishabitat] == 1) {
			// Inside the infected zone
			if (dist <= dist_param[_z_inf]) {
				val= control_param[_controltypes * _z_inf + _c_duration] + 1;
				// Increased carcass removal
				if (control_param[_controltypes * _z_inf + _c_carc] == 1) {if (val > control[i][_c_carc]) {control[i][_c_carc]= val;}}
				// Increased hunting
				if (control_param[_controltypes * _z_inf + _c_hunt] == 1) {if (val > control[i][_c_hunt]) {control[i][_c_hunt]= val;}}
				// Increased surveillance
				if (control_param[_controltypes * _z_inf + _c_surv] == 1) {if (val > control[i][_c_surv]) {control[i][_c_surv]= val;}}
				// Erecting fence
				if (control_param[_controltypes * _z_inf + _c_fence] == 1) {control[i][_c_fence]= (*fenceid); (*fenceid) += 1;}
			} else if (dist <= dist_param[_z_buff]) {
				val= control_param[_controltypes * _z_buff + _c_duration] + 1;
				// Increased carcass removal
				if (control_param[_controltypes * _z_buff + _c_carc] == 1) {if (val > control[i][_c_carc]) {control[i][_c_carc]= val;}}
				// Increased hunting
				if (control_param[_controltypes * _z_buff + _c_hunt] == 1) {if (val > control[i][_c_hunt]) {control[i][_c_hunt]= val;}}
				// Increased surveillance
				if (control_param[_controltypes * _z_buff + _c_surv] == 1) {if (val > control[i][_c_surv]) {control[i][_c_surv]= val;}}
				// Erecting fence
				if (control_param[_controltypes * _z_buff + _c_fence] == 1) {control[i][_c_fence]= (*fenceid); (*fenceid) += 1;}
			} else if (dist <= dist_param[_z_surv]) {
				val= control_param[_controltypes * _z_buff + _c_duration] + 1;
				// Increased carcass removal
				if (control_param[_controltypes * _z_surv + _c_carc] == 1) {if (val > control[i][_c_carc]) {control[i][_c_carc]= val;}}
				// Increased hunting
				if (control_param[_controltypes * _z_surv + _c_hunt] == 1) {if (val > control[i][_c_hunt]) {control[i][_c_hunt]= val;}}
				// Increased surveillance
				if (control_param[_controltypes * _z_surv + _c_surv] == 1) {if (val > control[i][_c_surv]) {control[i][_c_surv]= val;}}
				// Erecting fence
				if (control_param[_controltypes * _z_surv + _c_fence] == 1) {control[i][_c_fence]= (*fenceid); (*fenceid) += 1;}
			}
		}
	}
}

// Simulation
void _Simulation(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt,
	int **population, int **carcass, 
	int **cell_info, int **road_cross, int **river_cross, int **proxy_cid, double **proxy_dist, int **proxy_road, int **proxy_river, 
	int simyear, int scale, int island, int introloc, int control_index, int foi,
	int *control_param, double *dist_param, int *fenceid, int iter, int **result_map, int **result_ts) {
	
	int _Random_binomial(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, int n, double p);
	double _Random_normal(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, double mu, double sd);
	int _Random_Poisson(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt, double mean);	
	double _Distance(int **cell_info, int i, int j);
	int _Migration(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt,
		int scale, int island, int **cell_info, int **population, int **road_cross, int **river_cross, int **control, int org, int sex);
	int _Random_neighbour(unsigned long MTstate[], unsigned long MTnext[], int *MTleft, int *MTcnt,
		int scale, int island, int **cell_info, int **road_cross, int **river_cross, int **control, int org);
	double _Cross_probability(int scale, int island, int **road_cross, int **river_cross, int **control, int org, int dst);
	int _Neighbour_index(int scale, int island, int org, int dst);
	int _Neighbour_cellid(int scale, int island, int org, int index);
	double _Calculate_lambda(int n_col_adj, int **population, int **carcass, int **proxy_cid, double **proxy_dist, int **proxy_road, int **proxy_river, int **control, int org);
	double _Calculate_lambda_Pepin(int n_col_adj, int **population, int **carcass, int **proxy_cid, double **proxy_dist, int **proxy_road, int **proxy_river, int **control, int org);
	double _Calculate_lambda_Thulke(int scale, int island, int **population, int **carcass, int org);
	void _Update_control(int scale, int island, int **cell_info, int *control_param, double *dist_param, int *fenceid, int **control, int org);

	int ncell, ncol_proxy;
	if (scale == 0) { // 4km
		ncol_proxy= _ncol_proxy_4km;
		if (island == 0) { // North Island
			ncell= _ncell_NI_4km;
		} else { // South Island
			ncell= _ncell_SI_4km;
		}
	} else { // 10 km
		ncol_proxy= _ncol_proxy_10km;
		if (island == 0) { // North Island
			ncell= _ncell_NI_10km;
		} else { // South Island
			ncell= _ncell_SI_10km;
		}
	}

	int i, j, val; // Basic variables
	int w, year, week, simperiod; // Variables for simulation week ("w"), simulation year ("year"), calendar week ("week"), and simulation period ("simperiod")
	int dst_female, dst_male, nmove; // Variales required for migration
	int nfy, allfemaleyearling, nsow, allsow, total; // Variables required for migration ("nfy", "allfemaleyealring"), birthing("nsow", "allsow"), and natural mortality ("total")
	double pmove, pbreed; // Variables required for migration probability ("pmove") and breeding probability ("pbreed")
	double **pmort, rmort_adult, rmort_young, phunt, rhunt_adult, rhunt_young; // Variables required for the natural mortality and hunting rate
	double r_carcass, r_asf_carcass; // Variables required for carcass degradation
	int dM, dS, dE, dI, dR, dC; // Variables required for SEIR model
	double rtrans, lambda; // Variables required for disease compartment transition in SEIR model
	int **w_demo; // An array to store the calendar week of demographic events
	int **tmp_population, **tmp_carcass; // Temporary arrays to store the number of pigs and carcasses
	int **control, **tmp_control; // Arrays to indicate the control measure
	int ntest; // Variable to store the number of pigs/carcasses to test
	double pcarc, psurv; // Variables required for carcass removal and surveillance for controlling ASF

	w_demo= calloc(ncell, sizeof(int *));
	tmp_population= calloc(ncell, sizeof(int *));
	tmp_carcass= calloc(ncell, sizeof(int *));
	control= calloc(ncell, sizeof(int *));
	tmp_control= calloc(ncell, sizeof(int *));
	pmort= calloc(ncell, sizeof(double *));
	for (i= 0; i < ncell; ++i) {
		w_demo[i]= calloc(2, sizeof(int)); // 0 is for breeding and 1 is for weaning
		tmp_population[i]= calloc(_g_demo * _g_dz, sizeof(int));
		tmp_carcass[i]= calloc(_g_dz, sizeof(int));
		control[i]= calloc(_controltypes, sizeof(int));
		tmp_control[i]= calloc(_controltypes, sizeof(int));
		pmort[i]= calloc(2, sizeof(double)); // 0 is for adult and 1 is for young
	}
	year= 0;
	simperiod= (simyear + _initialsimyear) * _weeksinyear;
	
	for (w= _mean_breedweek; w < simperiod; ++w) {

		week= w % _weeksinyear; // calendar week
		year= w / _weeksinyear; // simulation year

		// Estimate the carcass degradation rate (depends on time of year)
		r_carcass= 22.5 * -1.0 * cos(acos(-1.0) * week / 26) + 29.5; // Up to 52 weeks in winter, 7 weeks in summer
		r_carcass= -1.0 * log(1.0 - 0.99) / r_carcass;
		r_carcass= 1.0 - exp(-1.0 * r_carcass);

		// Estimate the rate in loss of ASF virus viability in carcasses (depends on time of year)
		r_asf_carcass= 6.5 * -1.0 * cos(acos(-1.0) * week / 26) + 19.5; // Up to 26 weeks in winter, 13 weeks in summer (3 ~ 6 months based on Fischer et al 2020)
		r_asf_carcass= -1.0 * log(1.0 - 0.99) / r_asf_carcass;
		r_asf_carcass= 1.0 - exp(-1.0 * r_asf_carcass);

		// Determine the week of breeding and weaning
		if (w == _mean_breedweek || week == _w_update) {
			for (i= 0; i < ncell; ++i) {
				if (cell_info[i][_col_ishabitat] == 1) {
					if (w == _mean_breedweek) {
						w_demo[i][_col_breed]= 38;
					} else {
						do {
							w_demo[i][_col_breed]= _Random_Poisson(MTstate, MTnext, MTleft, MTcnt, _mean_breedweek);
							w_demo[i][_col_breed] %= _weeksinyear;
						} while (w_demo[i][_col_breed] == _w_update);
					}
					w_demo[i][_col_wean]= w_demo[i][_col_breed] + _weanperiod;
					w_demo[i][_col_wean] %= _weeksinyear;
				}
			}
		}

		// Introduce ASF
		if (week == 0 && year == _initialsimyear) {
			// North Island
			if (island == 0) {
				if (scale == 0) {
					// Auckland
					if (introloc == 0) {
						val= _auckland_4km;
					// Wellington
					} else {
						val= _wellington_4km;
					}
				} else {
					// Auckland
					if (introloc == 0) {
						val= _auckland_10km;
					// Wellington
					} else {
						val= _wellington_10km;
					}
				}
			} else {
				if (scale == 0) {
					// Picton
					if (introloc == 0) {
						val= _picton_4km;
					// Queenstown
					} else {
						val= _dunedin_4km;//_queenstown_4km;
					}
				} else {
					// Picton
					if (introloc == 0) {
						val= _picton_10km;
					// Queenstown
					} else {
						val= _dunedin_10km;//_queenstown_10km;
					}
				}
			}
			population[val][_g_demo * _s + _sow] -= 1; population[val][_g_demo * _i + _sow] += 1;
		}

		// Restore the populations and carcasses + Control indicators
		for (i= 0; i < ncell; ++i) {
			if (cell_info[i][_col_ishabitat] == 1) {
				for (j= 0; j < (_g_demo * _g_dz); ++j) {
					tmp_population[i][j]= population[i][j];
				}
				for (j= 0; j < _g_dz; ++j) {
					tmp_carcass[i][j]= carcass[i][j];
				}
			}
			for (j= 0; j < _controltypes; ++j) {
				if (j < _c_fence && control[i][j] > 0) {
					control[i][j] -= 1;
				}
				tmp_control[i][j]= control[i][j];
			}
		}
		
		// Demographic events
		for (i= 0; i < ncell; ++i) {
			// Only for the habitat cells
			if (cell_info[i][_col_ishabitat] == 1) {
				// Calculate the number of female yearlings, sows, and total population
				nfy= 0; nsow= 0; total= 0;
				for (j= 0; j < _g_dz; ++j) {
					nfy += population[i][_g_demo * j + _fp];
					nsow += population[i][_g_demo * j + _sow];
					total += (population[i][_g_demo * j + _fp] + population[i][_g_demo * j + _fy] + population[i][_g_demo * j + _sow] + 
						population[i][_g_demo * j + _mp] + population[i][_g_demo * j + _my] + population[i][_g_demo * j + _boar]);
				}
				allfemaleyearling= nfy;
				allsow= nsow;
				// "nfy" and "nsow" changes as the programme moves to the next disease group in one habitat cell,
				// whereas "allfemaleyealring" and "allsow" are constant.
				// This is because "nfy" and "nsow" are used to avoid unnecessary process of migration and birthing, respectively,
				// in the next disease groups, while "allfemaleyealring" and "allsow" are used to calculate the probability of migration and birthing, respectively,
				// across all disease groups.
				// "total" is required for estimating the natural motality rate

				// Designate the dstination for female migration
				if (week == w_demo[i][_col_breed]) {
					dst_female= _Migration(MTstate, MTnext, MTleft, MTcnt, scale, island, cell_info, population, road_cross, river_cross, control, i, _female);
				} else {
					dst_female= _na;
				}

				// Update the mortality rate (in 1 week after breeding)
				if (w == _mean_breedweek || week == w_demo[i][_col_breed] + 1) {
					if (total > cell_info[i][_col_popcap]) {
						pmort[i][_adult]= (double)(total - cell_info[i][_col_popcap]) / total;
						pmort[i][_young]= (double)(total - cell_info[i][_col_popcap]) / total;
					} else {
						pmort[i][_adult]= _pmort_adult;						
						pmort[i][_young]= _pmort_young;						
					}
				}
				rmort_adult= -1.0 * log((1.0 - pmort[i][_adult])) / _weeksinyear;
				rmort_adult= 1.0 - exp(-1.0 * rmort_adult);
				rmort_young= -1.0 * log((1.0 - pmort[i][_young])) / _weeksinyear;
				rmort_young= 1.0 - exp(-1.0 * rmort_young);

				// Designate the control parameter values and their rates based on "control" array
				// Carcass removal
				if (control_index == 0) {
					pcarc= 0.0;
				} else {
					if (control[i][_c_carc] > 0) {
						pcarc= _pcarc_inc;
					} else {
						pcarc= _pcarc;
					}
				}

	
				// Hunting
				if (control[i][_c_hunt] > 0) {
					phunt= _phunt_inc;
				} else {
					phunt= _phunt;
				}
				if (phunt < (1.0 - pmort[i][_adult])) {  
					rhunt_adult= -1.0 * log((1.0 - (phunt / (1.0 - pmort[i][_adult])))) / _weeksinyear;
				} else {
					rhunt_adult= -1.0 * log((1.0 - (1.0 - pmort[i][_adult]))) / _weeksinyear;
				}
				rhunt_adult= 1.0 - exp(-1.0 * rhunt_adult);
				if (phunt < (1.0 - pmort[i][_young])) {
					rhunt_young= -1.0 * log((1.0 - (phunt / (1.0 - pmort[i][_young])))) / _weeksinyear;
				} else {
					rhunt_young= -1.0 * log((1.0 - (1.0 - pmort[i][_young]))) / _weeksinyear;
				}
				rhunt_young= 1.0 - exp(-1.0 * rhunt_young);
				// If the sum of the proportoin of pigs naturally die and being hunted is > 100%,
				// natural mortality is prioritised and the hunting proportion reduces accordingly.
				
				// Surveillance
				if (control[i][_c_surv] > 0) {
					psurv= _psurv_inc;
				} else {
					psurv= _psurv;
				}

				ntest= 0;
				// Actual simulation starts from here!!!
				for (j= 0; j < _g_dz; ++j) {
					// Carcass degradation
					val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, carcass[i][j], r_carcass);
					if (val > tmp_carcass[i][j]) {val= tmp_carcass[i][j];}
					tmp_carcass[i][j] -= val;

					// Loss of ASF vialbility in infectious carcasses
					if (j == _i && tmp_carcass[i][j] != 0) {
						val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, carcass[i][j], r_asf_carcass);
						if (val > tmp_carcass[i][j]) {val= tmp_carcass[i][j];}
						tmp_carcass[i][_i] -= val;
						tmp_carcass[i][_r] += val;
					}

					// Control by carcass removal
					if (control[i][_c_carc] > 0) {
						val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, carcass[i][j], pcarc);
						if (val > tmp_carcass[i][j]) {val= tmp_carcass[i][j];}
						tmp_carcass[i][j] -= val;
						if (j == _i && val > 0) {ntest += val;} // Surveillance
					}
					
					// Natural mortality
					val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, population[i][_g_demo * j + _fp], rmort_young);
					if (val > tmp_population[i][_g_demo * j + _fp]) {val= tmp_population[i][_g_demo * j + _fp];}
					tmp_population[i][_g_demo * j + _fp] -= val; tmp_carcass[i][j] += val;
		
					val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, population[i][_g_demo * j + _fy], rmort_young);
					if (val > tmp_population[i][_g_demo * j + _fy]) {val= tmp_population[i][_g_demo * j + _fy];}
					tmp_population[i][_g_demo * j + _fy] -= val; tmp_carcass[i][j] += val; nfy -= val;

					if (nsow > 1) {
						val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, population[i][_g_demo * j + _sow], rmort_adult);
						if (val > tmp_population[i][_g_demo * j + _sow]) {val= tmp_population[i][_g_demo * j + _sow];}
						tmp_population[i][_g_demo * j + _sow] -= val; tmp_carcass[i][j] += val; nsow -= val;
					}

					val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, population[i][_g_demo * j + _mp], rmort_young);
					if (val > tmp_population[i][_g_demo * j + _mp]) {val= tmp_population[i][_g_demo * j + _mp];}
					tmp_population[i][_g_demo * j + _mp] -= val; tmp_carcass[i][j] += val;
		
					val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, population[i][_g_demo * j + _my], rmort_young);
					if (val > tmp_population[i][_g_demo * j + _my]) {val= tmp_population[i][_g_demo * j + _my];}
					tmp_population[i][_g_demo * j + _my] -= val; tmp_carcass[i][j] += val;

					val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, population[i][_g_demo * j + _boar], rmort_adult);
					if (val > tmp_population[i][_g_demo * j + _boar]) {val= tmp_population[i][_g_demo * j + _boar];}
					tmp_population[i][_g_demo * j + _boar] -= val; tmp_carcass[i][j] += val;


					// Hunting
					val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, population[i][_g_demo * j + _fp], rhunt_young);
					if (val > tmp_population[i][_g_demo * j + _fp]) {val= tmp_population[i][_g_demo * j + _fp];}
					tmp_population[i][_g_demo * j + _fp] -= val;
					if (j == _i && val > 0) {ntest += val;} // Surveillance

					val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, population[i][_g_demo * j + _fy], rhunt_young);
					if (val > tmp_population[i][_g_demo * j + _fy]) {val= tmp_population[i][_g_demo * j + _fy];}
					tmp_population[i][_g_demo * j + _fy] -=val; nfy -= val;
					if (j == _i && val > 0) {ntest += val;} // Surveillance

					if (nsow > 1) {
						val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, population[i][_g_demo * j + _sow], rhunt_adult);
						if (val > tmp_population[i][_g_demo * j + _sow]) {val= tmp_population[i][_g_demo * j + _sow];}
						tmp_population[i][_g_demo * j + _sow] -= val; nsow -= val;
						if (j == _i && val > 0) {ntest += val;} // Surveillance
					}

					val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, population[i][_g_demo * j + _mp], rhunt_young);
					if (val > tmp_population[i][_g_demo * j + _mp]) {val= tmp_population[i][_g_demo * j + _mp];}
					tmp_population[i][_g_demo * j + _mp] -= val;
					if (j == _i && val > 0) {ntest += val;} // Surveillance
		
					val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, population[i][_g_demo * j + _my], rhunt_young);
					if (val > tmp_population[i][_g_demo * j + _my]) {val= tmp_population[i][_g_demo * j + _my];}
					tmp_population[i][_g_demo * j + _my] -= val;
					if (j == _i && val > 0) {ntest += val;} // Surveillance

					val= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, population[i][_g_demo * j + _boar], rhunt_adult);
					if (val > tmp_population[i][_g_demo * j + _boar]) {val= tmp_population[i][_g_demo * j + _boar];}
					tmp_population[i][_g_demo * j + _boar] -= val;
					if (j == _i && val > 0) {ntest += val;} // Surveillance
					
					
					// Update control measures
					if (ntest > 0 && _Random_binomial(MTstate, MTnext, MTleft, MTcnt, ntest, psurv) > 0) {
						_Update_control(scale, island, cell_info, control_param, dist_param, fenceid, tmp_control, i);
					}


					// Weaning
					if (week == w_demo[i][_col_wean]) {
						// Female piglet -> Female yearling
						tmp_population[i][_g_demo * j + _fy] += tmp_population[i][_g_demo * j + _fp]; tmp_population[i][_g_demo * j + _fp]= 0;
		
						// Male piglet -> Male yearling
						tmp_population[i][_g_demo * j + _my] += tmp_population[i][_g_demo * j + _mp]; tmp_population[i][_g_demo * j + _mp]= 0;
					} // End: Weaning
					
					
					// Breeding
					if (week == w_demo[i][_col_breed]) {
						// Migration for female
						// Only when the number of sows in i is more than its breeding capacity AND an available cell exists
						if (dst_female != _na && nsow > cell_info[i][_col_breedcap] && allfemaleyearling != 0 && nfy != 0) {
							// All female yearlings can migrate
							if ((allfemaleyearling - cell_info[dst_female][_col_breedcap]) <= 0) {
								tmp_population[dst_female][_g_demo * j + _fy] += tmp_population[i][_g_demo * j + _fy]; nfy -= tmp_population[i][_g_demo * j + _fy]; tmp_population[i][_g_demo * j + _fy]= 0;
							// Only the proportion of female yearlings can move							
							} else {
								pmove= (double)tmp_population[i][_g_demo * j + _fy] / (double)allfemaleyearling;
								nmove= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, cell_info[dst_female][_col_breedcap], pmove);
								tmp_population[dst_female][_g_demo * j + _fy] += nmove; tmp_population[i][_g_demo * j + _fy] -= nmove; nfy -= nmove;
							}
						}
					
						// Migration for male
						nmove= tmp_population[i][_g_demo * j + _my];
						while (nmove > 0) {
							dst_male= _Migration(MTstate, MTnext, MTleft, MTcnt, scale, island, cell_info, population, road_cross, river_cross, control, i, _male);
							if (dst_male != _na) {
								tmp_population[dst_male][_g_demo * j + _boar] += 1;
								tmp_population[i][_g_demo * j + _my] -= 1;
							}
							nmove -= 1;
						}

						// Birth
						if (cell_info[i][_col_breedcap] > 0) {
							pbreed= (double)allsow / (double)cell_info[i][_col_breedcap];

							if (j == _m || j == _s || j == _r) {
								if (pbreed < 1.0) {
									pbreed= 1.0;
								} else {
									pbreed= (double)cell_info[i][_col_breedcap] / (double)allsow;
								}
							} else {
								if (pbreed < 1.0) {
									pbreed= 1.0 * 0.625; // Impaired reproduction (Lange et al 2017)
								} else {
									pbreed= 0.625 * (double)cell_info[i][_col_breedcap] / (double)allsow; // Impaired reproduction (Lange et al 2017)
								}
							}

							val= (int)round(_Random_normal(MTstate, MTnext, MTleft, MTcnt, _mean_litter, _sd_litter) * _Random_binomial(MTstate, MTnext, MTleft, MTcnt, population[i][_g_demo * j + _sow], pbreed));
							if (val > 0) {
								if (val < 2) { // If a litter size is 1, then only female is born
									if (j == _r) {
										tmp_population[i][_g_demo * _m + _fp] += val;
									} else {
										tmp_population[i][_g_demo * _s + _fp] += val;
									}	
								} else {
									if (val % 2 == 1) { // If a litter size is a odds number, then one more female is born
										val /= 2;
										if (j == _r) {
											tmp_population[i][_g_demo * _m + _fp] += (val + 1);
											tmp_population[i][_g_demo * _m + _mp] += val;
										} else {
											tmp_population[i][_g_demo * _s + _fp] += (val + 1);
											tmp_population[i][_g_demo * _s + _mp] += val;
										}
									} else {
										val /= 2;
										if (j == _r) {
											tmp_population[i][_g_demo * _m + _fp] += val;
											tmp_population[i][_g_demo * _m + _mp] += val;
										} else {
											tmp_population[i][_g_demo * _s + _fp] += val;
											tmp_population[i][_g_demo * _s + _mp] += val;
										}
									}
								}
							}
						}

						// Transition of yearlings to adult
						// Female yearling -> Sow
						tmp_population[i][_g_demo * j + _sow] += tmp_population[i][_g_demo * j + _fy]; nsow += tmp_population[i][_g_demo * j + _fy]; 
						nfy -= tmp_population[i][_g_demo * j + _fy]; tmp_population[i][_g_demo * j + _fy]= 0;

						// Male yearling -> Male boar
						tmp_population[i][_g_demo * j + _boar] += tmp_population[i][_g_demo * j + _my]; tmp_population[i][_g_demo * j + _my]= 0;

					} // End: Breeding week
					
				} // End: For disease group j
				
				if (year >= _initialsimyear) {
					if (foi == 0) {
						lambda= _Calculate_lambda(ncol_proxy, population, carcass, proxy_cid, proxy_dist, proxy_road, proxy_river, control, i);
					} else if (foi == 1) {
						lambda= _Calculate_lambda_Pepin(ncol_proxy, population, carcass, proxy_cid, proxy_dist, proxy_road, proxy_river, control, i);
					} else {
						lambda= _Calculate_lambda_Thulke(scale, island, population, carcass, i);
					}
				
					for (j= 0; j < _g_demo; ++j) {
						// M -> S
						rtrans= -1.0 * log(1.0 - 0.99) / _w_mab;
						rtrans= 1.0 - exp(-1.0 * rtrans);
						dM= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, population[i][_g_demo * _m + j], rtrans);
						if (dM > tmp_population[i][_g_demo * _m + j]) {dM= tmp_population[i][_g_demo * _m + j];}
						tmp_population[i][_g_demo * _m + j] -= dM;
						
						// S -> E
						dS= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, population[i][_g_demo * _s + j], lambda);
						if (dS > tmp_population[i][_g_demo * _s + j]) {dS= tmp_population[i][_g_demo * _s + j];}
						tmp_population[i][_g_demo * _s + j] += dM; tmp_population[i][_g_demo * _s + j] -= dS;

						// E -> I
						dE= population[i][_g_demo * _e + j];
						if (dE > tmp_population[i][_g_demo * _e + j]) {dE= tmp_population[i][_g_demo * _e + j];}
						tmp_population[i][_g_demo * _e + j] += dS; tmp_population[i][_g_demo * _e + j] -= dE;
	
						// I -> R or death
						dI= population[i][_g_demo * _i + j];
						if (dI > tmp_population[i][_g_demo * _i + j]) {dI= tmp_population[i][_g_demo * _i + j];}
						if (foi == 0) {
							dC= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, dI, _pdeath);
						} else {
							dC= _Random_binomial(MTstate, MTnext, MTleft, MTcnt, dI, _pdeath_alt);
						}
						 dR= dI - dC;

						tmp_population[i][_g_demo * _i + j] += dE; tmp_population[i][_g_demo * _i + j] -= dI;
						tmp_population[i][_g_demo * _r + j] += dR;
						tmp_carcass[i][_i] += dC;
					} // End: For each demograup
				}
			} // End: For habitat i
		}

		// Update the populations and carcasses + Control indicators
		for (i= 0; i < ncell; ++i) {
			if (cell_info[i][_col_ishabitat] == 1) {
				for (j= 0; j < (_g_demo * _g_dz); ++j) {
					population[i][j]= tmp_population[i][j];
				}
				for (j= 0; j < _g_dz; ++j) {
					carcass[i][j]= tmp_carcass[i][j];
				}
			}
			for (j= 0; j < _controltypes; ++j) {
				control[i][j]= tmp_control[i][j];
			}
		}
		
		if (year >= _initialsimyear) {
			val= _weeksinyear * (year - _initialsimyear) + week;
			for (i= 0; i < ncell; ++i) {
				if ((population[i][_g_demo * _i + _fp] + population[i][_g_demo * _i + _fy] + population[i][_g_demo * _i + _sow] + 
						population[i][_g_demo * _i + _mp] + population[i][_g_demo * _i + _my] + population[i][_g_demo * _i + _boar]) > 0 || carcass[i][_i] > 0) {
					result_map[i][val] += 1;
					result_ts[iter][val] += 1;
				}
			}
		}
	} // End: For week w

	// Free memory
	for (i= 0; i < ncell; ++i) {
		free(w_demo[i]);
		free(tmp_population[i]);
		free(tmp_carcass[i]);
		free(control[i]);
		free(tmp_control[i]);
		free(pmort[i]);
	}
	free(w_demo);
	free(tmp_population);
	free(tmp_carcass);
	free(control);
	free(tmp_control);
	free(pmort);
}

// Convert information
void _Convert_file(int ncell, int ncol_proxy, int pop, int *input_info, int *input_roadcross, int *input_rivercross, int *input_proxycid, double *input_proxydist, int *input_proxyroad, int *input_proxyriver,
	int **cell_info, int **road_cross, int **river_cross, int **proxy_cid, double **proxy_dist, int **proxy_road, int **proxy_river) {
	
	int i, j, k;
	
	// Cell info
	for (i= 0; i < _cellinfocols; ++i) {
		for (j= 0; j < ncell; ++j) {
			k= i * ncell + j;
			if (i == _col_breedcap || i == _col_popcap) {
				if (pop == 0) {
					cell_info[j][i]= input_info[k] * 1;
				} else if (pop == 1) {
					cell_info[j][i]= input_info[k] * 3;
				} else {
					cell_info[j][i]= input_info[k] * 6;
				}
			} else {
				cell_info[j][i]= input_info[k];
			}
		}
	}

	// Crossing frequency between neighbouring cells
	for (i= 0; i < _neighbours; ++i) {
		for (j= 0; j < ncell; ++j) {
			k= i * ncell + j;
			road_cross[j][i]= input_roadcross[k];
			river_cross[j][i]= input_rivercross[k];
		}
	}

	// Information and crossing frequency for proximal cells
	for (i= 0; i < ncol_proxy; ++i) {
		for (j= 0; j < ncell; ++j) {
			k= i * ncell + j;
			proxy_cid[j][i]= input_proxycid[k];
			proxy_dist[j][i]= input_proxydist[k];
			proxy_road[j][i]= input_proxyroad[k];
			proxy_river[j][i]= input_proxyriver[k];
		}
	}
}

