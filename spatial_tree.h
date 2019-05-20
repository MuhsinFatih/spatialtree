#pragma once

// debug
void swifttestF();
void pass_string(const char* deneme);
void pass_array(const double* arr, size_t count);
void pass_2d_array(const double_t** arr, size_t c1, size_t c2);
// end debug

struct c_array_double {
	double* data;
	size_t count;
};

struct c_box {
	int label_;
	double_t angle;
	double_t center[3];
	double_t shape[3];
};

struct c_box_array {
	struct c_box* data;
	size_t count;
	double_t surface;
};

/* ****************** API0 ****************** */
	/**
	 * @brief  create 3d space with double_t points
	 * @param  side_length side length of the entire space. This value is final
	 * @param  node_side_length side length of one spatial section. Select this in accordance to DBScan parameters
	 * @param  optimize_side_length if true, side length will be changed to fit node_side_length. If DBScan is used this must be set to true
	 * @retval pointer to 3D space
	 */
	void* create_3D_space_double_t(double_t side_length, double_t node_side_length, bool optimize_side_length);

	/**
	 * @brief  insert a 3D point in spatial tree
	 * @param c_point: x,y,z
	 * @param space: pointer to 3D space
	 * @retval true on success
	 */
	bool insert_3D_no_check(const double_t* c_point, void* space);
	/**
	 * @brief  insert a 3D point in spatial tree
	 * @param  c_point: point: x,y,z
	 * @param  space: pointer to 3D space
	 * @param  dup_dist: distance at which two points will be considered duplicate and only the first stored. Duplicate points will NOT make the return value false
	 * @retval true on success
	 */
	bool insert_3D(const double_t* c_point, double_t dup_dist, void* space);



	/**
	 * @brief  insert an array of 3D points in spatial tree
	 * @param  c_flat_array: array of 3D points flattened: [x1,y1,z1, x2,y2,z2, x3,y3,z3,...]
	 * @param  n: number of points (i.e: if there are 5 points, of each 3 dims, n is 5)
	 * @param  space: pointer to 3D space
	 * @retval true on success
	 */
	bool insert_bulk_3D_no_check(const double_t* c_flat_array, size_t n, void* space);
	/**
	 * @brief  insert an array of 3D points in spatial tree
	 * @param  c_flat_array: array of 3D points flattened: [x1,y1,z1, x2,y2,z2, x3,y3,z3,...]
	 * @param  n: number of points (i.e: if there are 5 points, of each 3 dims, n is 5)
	 * @param  dup_dist: distance at which two points will be considered duplicate and only the first stored. Duplicate points will NOT make the return value false
	 * @param  space: pointer to 3D space
	 * @retval true on success
	 */
	bool insert_bulk_3D(const double_t* c_flat_array, size_t n, double_t dup_dist, void* space);
	/**
	 * @brief  insert an array of 3D points in spatial tree
	 * @param  c_flat_array: array of 3D points flattened: [x1,y1,z1, x2,y2,z2, x3,y3,z3,...]
	 * @param  n: number of points (i.e: if there are 5 points, of each 3 dims, n is 5)
	 * @param  dup_dist: distance at which two points will be considered duplicate and only the first stored. Duplicate points will NOT make the return value false
	 * @param  space: pointer to 3D space
	 * @retval array of points that were added to space
	 */
	struct c_array_double* insert_and_return_bulk_3D(const double_t* c_flat_array, size_t n, double_t dup_dist, void* space);

	void print_points(void* space);
	size_t count_points(void* space);


	struct c_array_double* get_all_points(void* s);
	struct c_box_array query(double_t eps, size_t minPts, double_t surface_r, void* s);

/* **************** END API0 **************** */

/* ****************** API1 ****************** */
struct space_store {
	void* space;
	void* projection;
};

	/**
	* @brief  create 3d space with double_t points
	* @param  side_length side length of the entire space. This value is final
	* @param  node_side_length side length of one spatial section. Select this in accordance to DBScan parameters
	* @param  optimize_side_length if true, side length will be changed to fit node_side_length. If DBScan is used this must be set to true
	* @retval pointer to 3D space
	*/
void* create_3D_space_and_projection(double_t side_length, double_t node_side_length, bool optimize_side_length);
struct c_array_double* insert_and_return_bulk(const double_t* c_flat_array, size_t n, double_t dup_dist, void* space_col);
struct c_box_array query_api1(double_t eps, size_t minPts, const double_t* cam_pos, void* space_col);
struct c_array_double* get_all_points_api1(void* space_col, bool print);
/* **************** END API1 **************** */