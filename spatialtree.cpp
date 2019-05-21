#include "spatial_tree_utils.hpp"
#include <limits>
#include <algorithm>
#include <math.h>
#include <float.h>
#include <array>
#include <queue>
#include <cassert>
#include "colors.h"
#include <unordered_map>
#include <set>

#include <cassert>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>

#define STORE_POINTS_AT_ALL_LEVELS 1
#define DBSCAN_OPTIMIZATIONS 1
#define MAIN 1
#define TEST 0
#define DEBUG 0
#define DEBUG_VISUALIZE 0
#define DEBUG_VISUALIZE_3D 0
#define BRIDGE_DEBUG 1
#define POINTS_CARRY_EXTRA_INFO 1

#if(TEST)
#ifdef _MSC_VER
#include <boost/config/compiler/visualc.hpp>
#endif
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include "json.hpp"

void print_contents( const boost::property_tree::ptree& pt);

#endif


#if(DEBUG_VISUALIZE || DEBUG_VISUALIZE_3D)
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
#endif

#define INF 1000000000

#define all(x) x.begin(), x.end()

#if(POINTS_CARRY_EXTRA_INFO)
	#define extra 1
#else
	#define extra 0
#endif

#if(DBSCAN_OPTIMIZATIONS)
	// point scheme: [dimensions, label, extra info] // notice extra info is always the last
	#define point array<T, n+extra+1> // extra element is for label
	typedef array<double, 3+extra+1> point_3D;
	typedef array<double, 2+extra+1> point_2D;
	typedef std::pair<double, double> ppoint;
#else
	#define point array<T, n>
	typedef array<double, 3+extra> point_3D;
	typedef array<double, 2+extra> point_2D;
#endif


// #if(DEBUG_VISUALIZE_3D)

// #endif

// -- plot stuff
#if(DEBUG_VISUALIZE)

vector<vi> axis(2);
void draw() {
	axis[0].resize(11);
	axis[1].resize(11);
	std::iota(begin(axis[0]), end(axis[0]), -5);
	std::iota(begin(axis[1]), end(axis[1]), -5);
	plt::show();
	plt::figure_size(1200,1000);
}
vi xl = {-5,5};
vi yl = {-5,5};
void drawPoint(double_t x, double_t y, size_t height) {
	plt::plot({x}, {y}, ".");
	// plt::text(x,y, to_string(height));
	// plt::xticks(axis[0]);
	// plt::yticks(axis[1]);
	plt::xlim(xl[0], xl[1]);
	plt::ylim(yl[0], yl[1]);
	// plt::pause(0.01);
}
void drawRelativeHeights(point_2D center, double_t sidelength, vector<int> relative_heights) {
	for(size_t i=0; i<2; ++i) {
		for(size_t j=0; j<2; ++j) {
			plt::text(center[0] + (i == 0)*((j? 1:-1))*(sidelength/2-(j==0?0.1:0.2)), center[1] + (i == 1)*(j? 1:-1)*(sidelength/2-(j==0?0.1:0.3)), to_string(relative_heights[2*i+j]));
		}
	}
}
void drawSection(point_2D center, double_t sidelength, point_2D rootcenter, double_t rootsidelength) {
	// plt::pause(0.5);
	// axis[0].resize((int)rootsidelength+1);
	// axis[1].resize((int)rootsidelength+1);
	// std::iota(begin(axis[0]), end(axis[0]), rootcenter[0] - (int)rootsidelength/2);
	// std::iota(begin(axis[1]), end(axis[1]), rootcenter[1] - (int)rootsidelength/2);
	// plt::xticks(axis[0]);
	// plt::yticks(axis[1]);
	xl = {(int)(rootcenter[0] - rootsidelength/2), (int)(rootcenter[0] + rootsidelength/2)};
	yl = {(int)(rootcenter[1] - rootsidelength/2), (int)(rootcenter[1] + rootsidelength/2)};
	plt::xlim(xl[0], xl[1]);
	plt::ylim(yl[0], yl[1]);
	plt::plot({center[0],center[0]},{center[1] - sidelength/2, center[1] + sidelength/2},"-");
	plt::plot({center[0] - sidelength/2 ,center[0] + sidelength/2},{center[1],center[1]},"-");


	// plt::pause(0.01);
}
void drawSection(point_3D center) {
}

void showSection(point_2D center, double_t sidelength) {
	plt::plot({center[0]-sidelength/2, center[0]+sidelength/2}, {center[1]-sidelength/2, center[1]+sidelength/2}, "--");
	plt::pause(0.01);
}
void showSection(point_2D center, double_t sidelength, string color) {
	plt::plot({center[0]-sidelength/2, center[0]+sidelength/2}, {center[1]-sidelength/2, center[1]+sidelength/2}, "--" + color);
	plt::pause(0.01);
}
void showSection(point_3D center, double_t sidelength) {
}
#endif
// -- end plot stuff

template<size_t n, typename T>
struct Box{
	int label_;
	double_t angle;
	array<double_t, n> center;
	array<double_t, n> shape;
	// vector<point> pts_cloud;
	// vector<point> pts_cloud_hull;
};


template<size_t n, typename T>
struct spatial_node {
	// notice both of these can be true!
	bool is_root; // is root?
	bool is_leaf; // is leaf?

	double_t sidelength; // length of one side of a section of this tree
	vector<point> points; // vector of points. points have fixed size n
	point center;

	// Children
	vector<spatial_node<n, T>*> sections;

	// position metadata
	size_t height;					// height starts from 0
	size_t location_code;			// location code based on which side of axis' the node is and its ancestors' location_code
	vector<size_t> location_code_v;	// x,y,z,... components of location_code. e.g. if x,y = {101, 110} then location_code = 101110 (notice it's binary!)
	vector<int> relative_heights;	// neighbor heights relative to this node's height
};
// n: dimensions
template<size_t n, typename T>
class SpatialTree {
public:
	// pre-calculated values
	size_t _sections_size;
	size_t _neighbor_size_immediate;
	size_t _neighbor_size_all;


	// configuration
	// min_sidelength has priority over max_sidelength. if min_sidelength requires that point is added, max_sidelength will be ignored.
	double_t min_sidelength = 0; // minimum length of one side of a section (improve performance on range queries)
	double_t max_sidelength; // maximum length of one side of a section (can help ensure neighbor sections include all points withing range)
	size_t capacity;
	bool fixed_level = false; // if true, it will be assumed that all points are at the same level. This will cancel updating neighbors relative heights as that will be obselete
	bool preserve_clusers = false;
	size_t cluster_id_seed = 2; // start labeling clusters from 2. If preserve_clusters is true, then advance this number to the largest id

	// utility
	spatial_node<n, T>* root = new spatial_node<n, T>();
	std::unordered_map<size_t, spatial_node<n, T>*> location_table; // hash table of spatial nodes. key is their location_code (combine the dimensions)
	// location_code is stored with a padding of 1 at the most significant bit to prevent overwriting due to first location bits being zero

	// DBScan specific optimization
	#if(DBSCAN_OPTIMIZATIONS)
	std::unordered_map<size_t, spatial_node<n, T>*> location_table_leaf_nodes; // hash table of spatial nodes that are leaf nodes. Can be used to iterate points in O(n) time rather than O(nlogn)
	#endif

	void construct() {
	  	// #if(DEBUG)
		  printf("created space with %i dimensions, %i spatial sections, %i max capacity per section\n", n, _sections_size, capacity);
	  	// #endif
		root->center = {0,0};
		root->is_root = true;
		root->is_leaf = true;
		root->height = 0;
		root->location_code_v.assign(n, 0);
		root->relative_heights.assign(_sections_size, 0); // resize and set all to 0
		max_sidelength = root->sidelength;
	};

	SpatialTree() {
		_sections_size = _2_to_power(n);
		_neighbor_size_immediate = 2*n;
		_neighbor_size_all = pow(3, n) - 1;
		capacity = _sections_size; // if more than _sections, divide into sections. Ideal number for minimum number of calculations
		construct();
	}

	SpatialTree(double_t sidelength) {
		_sections_size = _2_to_power(n);
		_neighbor_size_immediate = 2*n;
		_neighbor_size_all = pow(3, n) - 1;
		capacity = _sections_size; // if more than _sections, divide into sections. Ideal number for minimum number of calculations
		root->sidelength = sidelength;
		construct();
	}

	SpatialTree(uint capacity, double_t sidelength) {
		_sections_size = _2_to_power(n);
		_neighbor_size_immediate = 2*n;
		_neighbor_size_all = pow(3, n) - 1;
		this->capacity = capacity;
		root->sidelength = sidelength;
		construct();
	}

	void _store_spatial_node(spatial_node<n,T>* node) {
		size_t code = (node->location_code << ((sizeof(size_t)*8)-(node->height*n))) + node->height;
		location_table[code] = node;
	}

	// DBScan specific optimization
	#if(DBSCAN_OPTIMIZATIONS)
	void _store_spatial_leaf_node(spatial_node<n,T>* node) {
		size_t code = (node->location_code << ((sizeof(size_t)*8)-(node->height*n))) + node->height;
		if(location_table_leaf_nodes.find(code) == location_table_leaf_nodes.end())
			location_table_leaf_nodes[code] = node;
	}
	void _remove_spatial_leaf_node(spatial_node<n,T>* node) {
		size_t code = (node->location_code << ((sizeof(size_t)*8)-(node->height*n))) + node->height;
		if(location_table_leaf_nodes.find(code) == location_table_leaf_nodes.end())
			location_table_leaf_nodes.erase(code);
	}
	#endif

	spatial_node<n,T>* _retrieve_spatial_node(size_t location_code, size_t height) {
		size_t code = (location_code << ((sizeof(size_t)*8)-(height*n))) + height;
		return location_table[code];
	}
	spatial_node<n,T>* _retrieve_spatial_node_can_fail(size_t location_code, size_t height) {
		size_t code = (location_code << ((sizeof(size_t)*8)-(height*n))) + height;
		auto s = location_table.find(code);
		if(s != location_table.end())
			return s->second;
		else return NULL;
	}
   /**
    * @brief  constructs location code for a neighbor
    * @param  location_code_v: location code vector of the host node
    * @param  node_height: node's own height
    * @param  relative_height: relative height to the requested neighbor
    * @param  side: orientation of the neighbor. in other words, at which side is this neighbor?
    * @retval
    */
	size_t _construct_location_code(vector<size_t> location_code_v, size_t node_height, int relative_height, size_t side) {
		size_t d = side / 2; // dimension
		bool sign = side % 2;
		// now it's x flips in dimension d and side code is the sign
		size_t i = firstbitindex(location_code_v[d], sign); // find the first bit that is not equal to the sign
															// starting and including that bit to the right, reverse all bits
		location_code_v[d] = reverse_lastnbits(location_code_v[d], i+1);

		size_t location_code = 0;
		if(relative_height < 0)
			for(size_t i=0; i<n; ++i) {
				location_code <<= (node_height + relative_height);
				location_code += (location_code_v[i] >> -relative_height);
			}
		else
			for(size_t i=0; i<n; ++i) {
				location_code <<= node_height;
				location_code += location_code_v[i];
			}
		return location_code;
	}

	/**
	 * @brief  constructs location code for a neighbor
	 * @param  location_code_v: location code vector of the host node
	 * @param  node_height: node's own height
	 * @param  relative_height: relative height to the requested neighbor
	 * @param  sign: negative side: 0, positive side: 1
	 * @param  d: dimension index
	 * @retval concatenated location code
	 */
	size_t _construct_location_code(vector<size_t> location_code_v, size_t node_height, int relative_height, size_t sign, size_t d) {
		size_t i = firstbitindex(location_code_v[d], sign); // find the first bit that is not equal to the sign
															// starting and including that bit to the right, reverse all bits
		location_code_v[d] = reverse_lastnbits(location_code_v[d], i+1);

		size_t location_code = 0;
		if(relative_height < 0)
			for(size_t i=0; i<n; ++i) {
				location_code <<= (node_height + relative_height);
				location_code += (location_code_v[i] >> -relative_height);
			}
		else
			for(size_t i=0; i<n; ++i) {
				location_code <<= node_height;
				location_code += location_code_v[i];
			}
		return location_code;
	}

	/**
	 * @brief  constructs location code for a neighbor which can be at any alternation of direction in any dimension at the same level
	 * @param  location_code_v: location code vector of the host node
	 * @param  node_height: node's own height
	 * @param  direction: a number that represents direction at all dimensions. Number will be treated as in base 3. each digit meaning: [negative, no-change, positive] (0,2,1)
	 * @retval None
	 */
	size_t _construct_location_code_all_dim(vector<size_t> location_code_v, size_t node_height, size_t direction) {
		bool is_edge = false;
		for(size_t d=0; d<n; ++d, direction /= 3) {
			int sign = direction % 3;
			if(sign == 2) continue; // no-change digit
			if((!location_code_v[d] && !sign) || ((location_code_v[d] == _2_to_power(node_height) - 1) && sign)) {is_edge = true; break;} // if location_code is all 1's or all 0's (e.g: 111 or 1111 or 000), then it's an edge section
			size_t j = firstbitindex(location_code_v[d], sign); // find the first bit that is not equal to the sign
																// starting and including that bit to the right, reverse all bits
			location_code_v[d] = reverse_lastnbits(location_code_v[d], j+1);
		}
		if(is_edge) return -1;
		size_t location_code = 0;
		for(size_t i=0; i<n; ++i) {
			location_code <<= node_height;
			location_code += location_code_v[i];
		}
		return location_code;
	}

	/**
	 * @brief  construct location code vector for neighbor with same level
	 * @param  &location_code_v: buffer to modify
	 * @param  sign: negative side: 0, positive side: 1
	 * @param  d: dimension index
	 * @retval None
	 */
	void _construct_location_code_v(vector<size_t> &location_code_v, size_t sign, size_t d) {
		size_t i = firstbitindex(location_code_v[d], sign); // find the first bit that is not equal to the sign
															// starting and including that bit to the right, reverse all bits
		location_code_v[d] = reverse_lastnbits(location_code_v[d], i+1);
	}

	size_t _construct_location_code(spatial_node<n, T>* node) {
		size_t location_code = 0;
		for(size_t i=0; i<n; ++i) {
			location_code <<= node->height;
			location_code += node->location_code_v[i];
		}
		return location_code;
	}

	size_t construct_location_code(spatial_node<n, T>* node) {
		return _construct_location_code(node);
	}

	vector<spatial_node<n, T>*> _find_immediate_neighbors(spatial_node<n, T>* node) {
		set<spatial_node<n,T>*> neighbors;
		for(size_t i=0; i<_neighbor_size_immediate; ++i) {
			size_t d = i / 2; // dimension
			bool sign = i % 2;
			// todo: fix the case when its an edge section but direction is the opposite. // note: not required for dbscan
			if(!node->location_code_v[d] || (node->location_code_v[d] == _2_to_power(node->height) - 1)) continue; // if location_code is all 1's or all 0's (e.g: 111 or 1111 or 000), then it's an edge section
			size_t neighbor_location = _construct_location_code(node->location_code_v, node->height, node->relative_heights[i], sign, d);
			spatial_node<n, T>* neighbor = _retrieve_spatial_node(neighbor_location, node->height + node->relative_heights[i]);
			neighbors.insert(neighbor);
			#if(DEBUG)
				// printf("neighbor_location: %i\n", neighbor_location);
			#endif
		}
		#if(DEBUG_VISUALIZE)
			for(auto &nb : neighbors) {
				// showSection(nb->center, nb->sidelength, "r");
			}
		#endif
		return neighbors;
	}

	vector<spatial_node<n,T>*> _find_all_same_level_neighbors(spatial_node<n,T>* node) {
		vector<spatial_node<n,T>*> neighbors;
		for(size_t i=0; i<_neighbor_size_all; ++i) { // this will go up to 21 or 221 (base3), which will exclude "no change" which would have returned the current node
			size_t neighbor_location = _construct_location_code_all_dim(node->location_code_v, node->height, i);
			if(neighbor_location == -1) continue; // edge section
			spatial_node<n,T>* neighbor = _retrieve_spatial_node_can_fail(neighbor_location, node->height);
			if(neighbor != NULL)
				neighbors.push_back(neighbor);
		}
		return neighbors;
	}
	set<spatial_node<n,T>*> _find_all_neighbors(spatial_node<n,T>* node, int relative_height) {
		set<spatial_node<n,T>*> neighbors;
		for(size_t i=0; i<_neighbor_size_all; ++i) { // this will go up to 21 or 221 (base3), which will exclude "no change" which would have returned the current node
			size_t neighbor_location = _construct_location_code_all_dim(node->location_code_v, node->height + relative_height, i);
			if(neighbor_location == -1) continue; // edge section
			spatial_node<n,T>* neighbor = _retrieve_spatial_node_can_fail(neighbor_location, node->height);
			if(neighbor != NULL)
				neighbors.push_back(neighbor);
		}
	}
	// this function is incorrect
	vector<spatial_node<n,T>*> _find_all_neighbors(spatial_node<n,T>* node) {
		vector<spatial_node<n,T>*> neighbors;
		/*
		0
		1
		00
		01
		10
		11
		000
		001
		010
		011
		100
		101
		110
		111
		*/
		for(size_t i=0; i<_neighbor_size_all; ++i) {
			size_t d = i / 2; // dimension // no it's not
			bool sign = i % 2;
			if(!node->location_code_v[d] || (node->location_code_v[d] == _2_to_power(node->height) - 1)) continue; // if location_code is all 1's or all 0's (e.g: 111 or 1111 or 000), then it's an edge section
			size_t neighbor_location = _construct_location_code(node->location_code_v, node->height, node->relative_heights[i], sign, d);
			spatial_node<n, T>* neighbor = _retrieve_spatial_node(neighbor_location, node->height + node->relative_heights[i]);
			neighbors.push_back(neighbor);

			// find corner neighbors:
			spatial_node<n,T>* corner_neighbor;
			if(d != n-1) {
				// negative side of the next dimension
				size_t relative_height = (d+1)*2 + 0;
				neighbor_location = _construct_location_code(neighbor->location_code_v, neighbor->height, neighbor->relative_heights[relative_height], 0, d+1);
				corner_neighbor = _retrieve_spatial_node(neighbor_location, neighbor->height + neighbor->relative_heights[relative_height]);
				neighbors.push_back(corner_neighbor);

				// positive side of the next dimension
				++relative_height;
				neighbor_location = _construct_location_code(neighbor->location_code_v, neighbor->height, neighbor->relative_heights[relative_height], 1, d+1);
				corner_neighbor = _retrieve_spatial_node(neighbor_location, neighbor->height + neighbor->relative_heights[relative_height]);
				neighbors.push_back(corner_neighbor);
			}
		}
	}

	vector<spatial_node<n, T>*> _find_non_sibling_immediate_neighbors(spatial_node<n, T>* node) {
		vector<spatial_node<n, T>*> neighbors;
		for(size_t i=0; i<_neighbor_size_immediate; ++i) {
			size_t d = i / 2; // dimension
			bool sign = i % 2;
			if(!node->location_code_v[d] || (node->location_code_v[d] == _2_to_power(node->height) - 1)) continue; // if location_code is all 1's or all 0's (e.g: 111 or 1111 or 000), then it's an edge section
			if(sign != node->location_code_v[d] & 1) continue; // sibling
			size_t neighbor_location = _construct_location_code(node->location_code_v, node->height, node->relative_heights[i], sign, d);
			spatial_node<n, T>* neighbor = _retrieve_spatial_node(neighbor_location, node->height + node->relative_heights[i]);
			neighbors.push_back(neighbor);
			#if(DEBUG)
				printf("neighbor_location: %i\n", neighbor_location);
			#endif
		}
		#if(DEBUG_VISUALIZE)
			for(auto &nb : neighbors) {
				showSection(nb->center, nb->sidelength, "r");
			}
		#endif
		return neighbors;
	}

	void _update_non_sibling_immediate_neighbor_relative_heights(spatial_node<n, T>* node) { // lol that's mouthful
		size_t neighbor_location;
		spatial_node<n, T>* neighbor;
		for(size_t i=0; i<_neighbor_size_immediate; ++i) {
			size_t d = i / 2; // dimension
			bool sign = i % 2;
			#if(DEBUG)
				if(node->location_code == 0b011101 && d == 0 && sign == 1) {
					// __drawAllRelativeHeights();
					// plt::pause(0.01);
				}
			#endif
			if((!node->location_code_v[d]) || (node->location_code_v[d] == _2_to_power(node->height) - 1)) continue; // if location_code is all 1's or all 0's (e.g: 111 or 1111 or 000), then it's an edge section
			if(sign != (node->location_code_v[d] & 1)) continue; // sibling
			// find the neighbor location code vector
			vector<size_t> neighbor_location_code_v = node->location_code_v;
			_construct_location_code_v(neighbor_location_code_v, sign, d);
			// find the neighboring section that is at the same level or has the smallest height difference // (go one level deeper until neighbor is leaf or is at same level)
			int relative_height = node->relative_heights[i];
			for(; relative_height <= 0; ++relative_height) {
				neighbor_location = 0;
				for(size_t j=0; j<n; ++j) { // reduce to relative_height
					neighbor_location <<= node->height + relative_height;
					neighbor_location += (neighbor_location_code_v[j] >> -relative_height);
				}
				// get neighbor at current level
				neighbor = _retrieve_spatial_node(neighbor_location, node->height + relative_height); // todo: maybe optimizable
				if(neighbor->is_leaf) {++relative_height; break;} // if neighbor is a leaf node then there is no deeper neighbor
			}
			--relative_height;
			// neighbor is found
			node->relative_heights[i] = relative_height;
			size_t neighbor_relative_height_index = 2*d + !sign; // neighbor's relative height for host will be on the same dimension but opposite direction
			neighbor->relative_heights[neighbor_relative_height_index] = 0; // neighbor is either at the same level or higher. In either case relative height is 0

			// update children recursively
			int height = -1;
			queue<spatial_node<n,T>*> children;
			children.push(neighbor);
			while(!children.empty()) {
				neighbor = children.front();
				children.pop();
				if(neighbor->is_leaf) continue;
				for(size_t j=0; j<_sections_size; ++j) {
					if(!((j >> d) & sign)) continue; // skip the sections that are on the other side
					spatial_node<n,T>* child = neighbor->sections[j];
					children.push(child);
					child->relative_heights[neighbor_relative_height_index] = height;
				}
				--height;
			}
		}
	}

	vector<spatial_node<n, T>*> find_neighbors(point p) {
		spatial_node<n, T>* node = _locate(p, root);
		vector<spatial_node<n,T>*> neighbors = _find_immediate_neighbors(node);
		return neighbors;
	}
	vector<point> find_neighbor_points(point p) {
		vector<spatial_node<n,T>*> neighbors = _find_immediate_neighbors(p);
		vector<point> np;
		for(auto &neighbor : neighbors) {
			/* TODO : do DFS or store points on all ancestors */// + done
			np.insert(end(np), begin(neighbor.points), end(neighbor.points));
		}
		return np;
	}


	void divide(spatial_node<n, T>* node) {
		node->sections.resize(_sections_size);
		node->is_leaf = false;
		// for(auto s : node->sections) s = new spatial_node<n, T>(); // initialize all to avoid pointer relocation
		bool sign;
		for(size_t i=0; i<_sections_size; ++i) {
			// printf("signs: "); printArray(signs, signs.size());
			node->sections[i] = new spatial_node<n, T>();
			spatial_node<n, T> &section = *node->sections[i];
			auto &ctr = section.center;
			section.sidelength = node->sidelength/2;
			section.is_root = false;
			section.is_leaf = true;
			section.height = node->height + 1;
			section.location_code_v = node->location_code_v;
			if(!fixed_level)
				section.relative_heights.resize(_neighbor_size_immediate);
			for(size_t j=0; j<n; ++j) {
				sign = i>>j & 1;
				section.location_code_v[j] <<= 1;
				if(sign)
					section.location_code_v[j] += 1;

				ctr[j] = node->center[j] + (sign ? 1 : -1) * node->sidelength/4;
				/*
				+++	7	4 2 1				++	3
				-++	6	4 2 				-+	2
				+-+	5	4   1				+-	1
				--+	4	4					--	0
				++-	3	  2 1
				-+-	2	  2
				+--	1	    1
				---	0
				*/
			}
			section.location_code = _construct_location_code(&section);
			if(!fixed_level) {
				for(size_t j=0; j<n; ++j) {
					if(!(section.location_code_v[j] & 1)) { // section is on the negative side for dimension j
						section.relative_heights[2*j] = node->relative_heights[2*j] - 1;
						section.relative_heights[2*j+1] = 0;
					} else { // section is on the positive side for dimension j
						section.relative_heights[2*j] = 0;
						section.relative_heights[2*j+1] = node->relative_heights[2*j+1] - 1;
					};
				}
				_update_non_sibling_immediate_neighbor_relative_heights(&section);
			}
			#if(DEBUG_VISUALIZE)
				// drawRelativeHeights(section.center, section.sidelength, section.relative_heights);
			#endif
			#if(DEBUG)
				string ss = "";
				for(size_t j=0; j<section.location_code_v.size(); ++j) {
					string s = getbits(section.location_code_v[j], section.height) + " ";
					ss += s;
				}
				#if(DEBUG_VISUALIZE)
					// plt::text(ctr[0], ctr[1], ss);
				#endif
			#endif
			_store_spatial_node(&section);
		}
		#if(DBSCAN_OPTIMIZATIONS)
		_remove_spatial_leaf_node(node);
		#endif
		#if(DEBUG_VISUALIZE)
			// plt::pause(0.01);
		#endif
	}
	#if(DEBUG_VISUALIZE)
	void ___drawAllRelativeHeights(spatial_node<n, T>* node) {
		drawRelativeHeights(node->center, node->sidelength, node->relative_heights);
		plt::text(node->center[0], node->center[1], to_string(node->height).c_str());
		for(auto section : node->sections) {
			if(section != NULL)
				___drawAllRelativeHeights(section);
		}
	}

	void __drawAllRelativeHeights() {
		___drawAllRelativeHeights(root);
	}
	#endif

	bool _check_and_insert_outside(point p) {
		array<int, n> signs;
		signs.fill(1);
		bool outside = false;
		for(size_t i=0; i<n; ++i) {
			// check if outside the boundaries
			if(p[i] >= root->center[i] + root->sidelength/2) // it is very important that this is bigger or equal! ensures that points in the middle go to positive side
				outside = true;
			else if(p[i] < root->center[i] - root->sidelength/2) {
				outside = true;
				signs[i] = -1;
			}
		}
		// if(!outside) return false;
		return outside; // temporary workaround. //todo: finish this implementation
		#if(DEBUG)
			printf("point {%.f %.f} is outside bounds, expanding tree in direction: {%s %s}\n", p[0], p[1], signs[0] ? "+" : "-", signs[1] ? "+" : "-");
		#endif
		point new_center;
		for(size_t i=0; i<n; ++i)
			new_center[i] = root->center[i] + (signs[i] * root->sidelength/2);
		#if(DEBUG)
			printf("new center:\n"); printArray(new_center, n);
		#endif
		spatial_node<n, T>* new_root = new spatial_node<n, T>;
		new_root->is_root = true;
		// new_root->is_leaf = false; // already done in divide()
		new_root->height = 0;
		// TODO: do DFS and increase height
		new_root->center = new_center;
		new_root->sidelength = root->sidelength * 2;
		divide(new_root);
		#if(DEBUG_VISUALIZE)
			drawSection(new_center, new_root->sidelength, new_center, new_root->sidelength);
		#endif
		size_t old_root_pos = 0;
		for(size_t i=0; i<n; ++i)
			if(signs[i] == 1) old_root_pos += (int)_2_to_power(i);
		new_root->sections[old_root_pos] = root;
		root = new_root;
		return true;
	}

	void _insert(point p, spatial_node<n, T>* node) {
		if(!node->is_leaf) {
			size_t i = 0;
			for(size_t j=0; j<n; ++j) {
				if(p[j] >= node->center[j]) i+= _2_to_power(j)*1; // bigger or equal: ensures that points in the middle go to positive side
			}
			_insert(p, node->sections[i]);
		} else {
			if((node->points.size() < capacity && node->sidelength <= max_sidelength) || ((node->sidelength/2) < min_sidelength)) {
				node->points.push_back(p);
				#if(DEBUG_VISUALIZE)
					drawPoint(p[0], p[1], node->height);
				#endif
			} else {
				#if(DEBUG)
					printf("dividing for point: {%.f %.f}     center: {%.f %.f}     sidelength: %f\n", p[0], p[1], node->center[0], node->center[1], node->sidelength);
					#if(DEBUG_VISUALIZE)
						drawSection(node->center, node->sidelength, root->center, root->sidelength);
					#endif
				#endif
				// for(auto &v : node->points)
				// 	if(v == p) return; // it's a duplicate point, don't insert it
				divide(node);
				_insert(p, node);
				for(auto &v : node->points)
					_insert(v, node);
			}
		}
	}

	void _insert_all_levels_below(point p, spatial_node<n, T>* node) {
		if(!node->is_leaf) {
			size_t i = 0;
			for(size_t j=0; j<n; ++j) {
				if(p[j] >= node->center[j]) i+= _2_to_power(j)*1; // bigger or equal: ensures that points in the middle go to positive side
			}
			_insert_all_levels(p, node->sections[i]);
		} else {
			if((node->points.size() < capacity && node->sidelength <= max_sidelength) || ((node->sidelength/2) < min_sidelength)) {
				node->points.push_back(p);
				#if(DBSCAN_OPTIMIZATIONS)
					_store_spatial_leaf_node(node);
				#endif
				#if(DEBUG_VISUALIZE)
					drawPoint(p[0], p[1], node->height);
				#endif
			} else {
				// for(auto &v : node->points)
				// 	if(v == p) return; // it's a duplicate point, don't insert it
				#if(DEBUG)
					printf("dividing for point: {%.f %.f}     center: {%.f %.f}     sidelength: %f\n", p[0], p[1], node->center[0], node->center[1], node->sidelength);
					#if(DEBUG_VISUALIZE)
						drawSection(node->center, node->sidelength, root->center, root->sidelength);
					#endif
				#endif
				divide(node);
				_insert_all_levels(p, node);
				for(auto &v : node->points)
					_insert_all_levels_below(v, node);
			}
		}
	}

	void _insert_all_levels(point p, spatial_node<n, T>* node) {
		if(!node->is_leaf) {
			size_t i = 0;
			for(size_t j=0; j<n; ++j) {
				if(p[j] >= node->center[j]) i+= _2_to_power(j)*1; // bigger or equal: ensures that points in the middle go to positive side
			}
			node->points.push_back(p);
			_insert_all_levels(p, node->sections[i]);
		} else {
			if((node->points.size() < capacity && node->sidelength <= max_sidelength) || ((node->sidelength/2) < min_sidelength)) {
				node->points.push_back(p);
				#if(DBSCAN_OPTIMIZATIONS)
					_store_spatial_leaf_node(node);
				#endif
				#if(DEBUG_VISUALIZE)
					drawPoint(p[0], p[1], node->height);
				#endif
			} else {
				// for(auto &v : node->points)
				// 	if(v == p) return; // it's a duplicate point, don't insert it
				#if(DEBUG)
					printf("dividing for point: {%.f %.f}     center: {%.f %.f}     sidelength: %f\n", p[0], p[1], node->center[0], node->center[1], node->sidelength);
					#if(DEBUG_VISUALIZE)
						drawSection(node->center, node->sidelength, root->center, root->sidelength);
					#endif
				#endif
				divide(node);
				node->points.push_back(p);
				for(auto &v : node->points)
					_insert_all_levels_below(v, node);
			}
		}
	}

	bool insert(point p) { // todo: maybe insert bulk implementation. Maybe only update non-sibling-neighbors after all points are inserted? I feel like this part can use some optimization
		// while(_check_and_insert_outside(p)); // expand until the point is spanned by the space // todo: test this implementation
		if(_check_and_insert_outside(p)) return false; // temporary workaround. this will not accept points outside initial space
		#if(STORE_POINTS_AT_ALL_LEVELS)
			_insert_all_levels(p, root);
		#else
			_insert(p, root);
		#endif
    return true;
	}


	spatial_node<n, T>* _locate(point p, spatial_node<n, T>* node) {
		for(size_t k=0; k<n; ++k)
			if(abs(p[k]) > abs(node->center[k]) + node->sidelength/2)
				return NULL;
		if(node->is_leaf) return node;
		size_t i = 0;
		for(size_t j=0; j<n; ++j) {
			if(p[j] >= node->center[j]) i+= _2_to_power(j)*1; // bigger or equal: ensures that points in the middle go to positive side
		}
		return _locate(p, node->sections[i]);
	}

	spatial_node<n, T>* locate(point p) {
		spatial_node<n ,T>* section = _locate(p, root);
		#if(DEBUG)
			printf("located point {%.f, %.f}. Location code: \n", p[0], p[1]);
			string ss = "";
			for(size_t i=0; i<section->location_code_v.size(); ++i) {
				string s = getbits(section->location_code_v[i], section->height) + " ";
				ss += s;
			}
			// printf("  location_code: ");
			// dumpbits(&section->location_code, "");
			cout << endl;
			#if(DEBUG_VISUALIZE)
				plt::title(ss);
				showSection(section->center, section->sidelength);
			#endif
		#endif
		return section;
	}

	/**
	 * @brief  check if point is duplicate
	 * @note   dup_dist must be smaller than section side length. Otherwise the behavior is undefined
	 * @param  p: point
	 * @param  dup_dist: distance at which one of 2 points will be considered duplicate
	 * @retval true if duplicate
	 */
	bool check_dup(point p, T dup_dist) {
		// locate point
		spatial_node<n,T>* node = locate(p);
		// check immediate node first as it is more likely that there is a duplicate point here
		for(auto &pt : node->points)
			if(_distance(p, pt) <= dup_dist) return true;

		// check all neighbors afterwards
		vector<spatial_node<n,T>*> neighbor_sections = _find_all_same_level_neighbors(node);
		for(spatial_node<n,T>* &s : neighbor_sections) {
			for(auto &pt : s->points)
				if(_distance(p, pt) <= dup_dist) return true;
		}
		return false;
	}

	void DFS() {
		// maybe implement this. in case walk isn't suitable for some reason. Walk will be way faster
	}
	void walk() {
		size_t iterations = 0;
		size_t i=0;
		size_t total_point_count = 0;
		for(auto const& [key,node] : location_table_leaf_nodes) {
			if(((spatial_node<n,T>*)node)->is_leaf && node->points.size() > 0) {
				i++;
				total_point_count += node->points.size();
				vector<spatial_node<n,T>*> neighbors = _find_all_same_level_neighbors(node);
				#if(DEBUG_VISUALIZE)
					showSection(node->center, node->sidelength, "r");
					printf("center: [%.f, %.f]\n", node->center[0], node->center[1]);
					int a = 0;
					for(auto &v : neighbors) {
						showSection(v->center, v->sidelength, "g");
						plt::text(v->center[0], v->center[1], to_string(a++));
					}
					plt::pause(0.5);
				#endif
				size_t count = 0;
				for(auto &v : neighbors)
					count += v->points.size();
				neighbors.push_back(node);
				for(auto& p : node->points) {
					for(auto& v : neighbors) {
						for(auto& q : v->points) {
							++iterations;
						}
					}
					#if(DEBUG)
						printf("point: {%.f %.f %.f} has " GREEN "%lu" RESET " neighboring points in " MAGENTA "%lu" RESET " neighboring sections\n", p[0], p[1], p[2], count + node->points.size(), neighbors.size());
					#endif
				}
			}
		}
		#if(TEST)
			printf("total point count = %i\n", total_point_count);
			printf("there are %lu sections. %lu point operations done in walk procedure\n", location_table.size(), iterations);
		#endif
	}

	void get_points(vector<point> &buffer) {
		for(auto const& [key,node] : location_table_leaf_nodes) {
			if(node->points.size() > 0) {
				buffer.insert(buffer.end(), all(node->points));
			}
		}
	}

	#define pow2(n) ((n)*(n))
	template<size_t m>
	double_t _distance(array<double_t,m> p, array<double_t,m> q) {
		double_t dist = 0;
		for(size_t i=0; i<n; ++i) {
			dist += pow2(p[i] - q[i]);
		}
		return sqrt(dist);
	}
	double_t _distance(point_2D p, point_2D q) {
		return sqrt((pow2(p[0] - q[0]) + pow2(p[1] - q[1])));
	}
	double_t _distance(point_3D p, point_3D q) {
		return sqrt(pow2(p[0] - q[0]) + pow2(p[1] - q[1]) + pow2(p[2] - q[2]));
	}
	double_t _distance_h(point_3D p, point_3D q) {
		return sqrt(pow2(p[0] - q[0]) + pow2(p[2] - q[2]));
	}
	template<size_t m>
	double_t _distance_z(array<double_t,m> p, array<double_t,m> q) {
		return abs(p[2] - q[2]);
	}
	double_t _distance_z(point_3D p, point_3D q) {
		return abs(p[2] - q[2]);
	}


	#define label n+extra	// element at index n in point with size n+1 is the label
	#define noise 2
	#define fresh 1
	#define eliminated 0
	void clear_labels() {
		for(auto const& [key, node] : location_table_leaf_nodes)
			for(auto &p : node->points)
				p[label] = fresh;
	}
	void clear_labels_preserve_clusters() {
		for(auto const& [key, node] : location_table_leaf_nodes)
			for(auto &p : node->points) {
				if(p[label] == noise || p[label] == eliminated)
					p[label] = fresh;
			}
	}

	void DBScan(double_t eps, size_t minPts) {
		#if(DEBUG)
		size_t iterations;
		#endif
		size_t c = 2;
		// iterate all points
		for(auto const& [key, node] : location_table_leaf_nodes) {
			if(node->points.size() == 0) continue;
			vector<spatial_node<n,T>*> neighbor_sections = _find_all_same_level_neighbors(node);
			neighbor_sections.push_back(node); // include points in current section
			for(auto &p : node->points) {
				if(p[label] != fresh || p[label] == eliminated) continue;
				vector<pair<spatial_node<n,T>*, point*>> neighbor_points;
				// populate neighbor_points (range_query)
				for(auto &section : neighbor_sections) {
					for(auto &q : section->points) {
						if(_distance(p, q) <= eps)
							neighbor_points.push_back(make_pair(section,&q));
					}
				}

				if(neighbor_points.size() >= minPts+1) {
					p[label] = ++c;
					size_t neighbor_points_count = neighbor_points.size();
					for(size_t i=0; i<neighbor_points_count; ++i) {
						pair<spatial_node<n,T>*, point*> &p = neighbor_points[i]; // next element
						if((*p.second)[label] == noise) {
							(*p.second)[label] = c;
							continue;
						}
						if((*p.second)[label] != fresh) continue; // skip if already processed and added to another cluster
						(*p.second)[label] = c;
						vector<spatial_node<n,T>*> p_neighbor_sections = _find_all_same_level_neighbors(p.first); // todo: this can be cached
						p_neighbor_sections.push_back(p.first); // include points in current section
						vector<pair<spatial_node<n,T>*, point*>> p_neighbor_points;
						for(auto &section : p_neighbor_sections) {
							for(auto &q : section->points) {
								if(_distance((*p.second), q) <= eps)
									p_neighbor_points.push_back(make_pair(section,&q));
							}
						}
						if(p_neighbor_points.size() >= minPts+1) {
							neighbor_points.insert(neighbor_points.end(), p_neighbor_points.begin(), p_neighbor_points.end());
							neighbor_points_count += p_neighbor_points.size();
						}
					}
				} else {
					p[label] = noise;
					continue;
				}
			}
		}
		

		#if(TEST)
			map<size_t, size_t> clusters;
			for(auto const& [key, node] : location_table_leaf_nodes) {
				for(auto &p : node->points) {
					if(p[label] != fresh && p[label] != noise && p[label] != eliminated) {
						++clusters[p[label]];
					}
				}
			}
			printf(CYAN "clusters:\n" RESET);
			bool nocluster = true;
			for(auto const& [c, size] : clusters) {
				nocluster = false;
				printf( CYAN "| " RESET "cluster %i has %i items\n", c, size);
			}
			if(nocluster) printf(BOLDRED "none\n" RESET);
		#endif
	}
	/**
	 * @brief  Detect surface using focus point provided by AR world. This is specific for SSApp.
	 *         Looks for points inside a horizontal radius: h_radius in downward direction, and returns surface when minPts is not satisfied
	 * @param  focus_point: first result of hittest with camera's center
	 * @param  h_radius: horizontal radius to search for nearby points
	 * @param  v_threshold: vertical threshold to satisfy minPts
	 * @retval closest estimated surface point
	 */
	double_t _detect_surface(point focus_point, double_t h_radius, double_t v_threshold, size_t minPts) {
		// spatial_node<n,T>* node = locate(point);
		// auto relative_height = (int)(_2_to_power(h_radius-node->height));
		// set<spatial_node<n,T>*> neighbor_sections = _find_all_neighbors(node, relative_height);
		// neighbor_sections.insert(node); // todo: implement filtering
		std::priority_queue<double_t> neighbor_z;
		// std::set<double_t,std::greater<int>> neighbor_z;
		for(auto const& [key, node] : location_table_leaf_nodes) {
			if(node->points.size() == 0) continue;
			for(auto &q : node->points)
				if(_distance_h(focus_point, q) <= h_radius && q[1] <= focus_point[1] + 0.03) // filter points close enough and are below focus_point or same level
					neighbor_z.push(q[1]);
		}
		
		
		double_t surface = -1000;
		vector<double_t> nz;
		while(!neighbor_z.empty()) {
			nz.push_back(neighbor_z.top());
			neighbor_z.pop();
		}
		REP(nz.size()) {
			surface = nz[i];
			REPJ(minPts) {
				double_t next = nz[i+j];
				if(surface - next > v_threshold) return surface;
			}
		}

		// auto it = neighbor_z.begin();
		// cout << "neighbor_z size: " << neighbor_z.size() << endl;
		// if(neighbor_z.size() < minPts) return surface;
		// for(int i=0; i<neighbor_z.size()-minPts; ++i) {
		// 	auto jt = it;
		// 	surface = *it++;
		// 	for(int j=0; j < minPts; ++j) {
		// 		double_t next = *jt++;
		// 		// printf("surface: %.2f, next: %.2f, vthres: %.2f\n", surface, next, v_threshold);
		// 		if(surface - next > v_threshold) return surface;
		// 	}
		// }
		return surface;
	}
	void boxScan(double_t eps, size_t minPts, array<double_t,3> cam_pos) {
		#if(DEBUG)
		size_t iterations;
		#endif
		size_t c = 2;
		vector<double_t> dists = {0.2, 0.4, 0.6, 1};
		vector<double_t> flat_thres_choices = {0.02, 0.05, 0.1, 0.2};
		double_t flat_thres = flat_thres_choices[flat_thres_choices.size()-1];
		
		// iterate all points
		for(auto const& [key, node] : location_table_leaf_nodes) {
			if(node->points.size() == 0) continue;
			vector<spatial_node<n,T>*> neighbor_sections = _find_all_same_level_neighbors(node);
			neighbor_sections.push_back(node); // include points in current section
			for(auto &p : node->points) {
				if(p[label] != fresh || p[label] == eliminated) continue;
				vector<pair<spatial_node<n,T>*, point*>> neighbor_points;
				// populate neighbor_points (range_query)
				for(auto &section : neighbor_sections) {
					for(auto &q : section->points) {
						if(_distance(p, q) <= eps)
						neighbor_points.push_back(make_pair(section,&q));
					}
				}
				
				double_t cam_dist = sqrt(pow(cam_pos[0]-p[0],2) + pow(cam_pos[1]-p[1],2) + pow(cam_pos[2]-p[2],2));
				for(auto i=0; i<dists.size(); ++i)
					if(cam_dist < dists[i]) {flat_thres = flat_thres_choices[i]; break;}
				int n_flat = 5;
				for(auto &q : neighbor_points) {
					if(abs(p[2] - (*q.second)[2]) > flat_thres)
						if (--n_flat <= 0) break;
				}

				if(neighbor_points.size() >= minPts+1 && n_flat <= 0) {
					p[label] = ++c;
					size_t neighbor_points_count = neighbor_points.size();
					for(size_t i=0; i<neighbor_points_count; ++i) {
						pair<spatial_node<n,T>*, point*> &p = neighbor_points[i]; // next element
						if((*p.second)[label] == noise) {
							(*p.second)[label] = c;
							continue;
						}
						if((*p.second)[label] != fresh) continue; // skip if already processed and added to another cluster
						(*p.second)[label] = c;
						vector<spatial_node<n,T>*> p_neighbor_sections = _find_all_same_level_neighbors(p.first); // todo: this can be cached
						p_neighbor_sections.push_back(p.first); // include points in current section
						vector<pair<spatial_node<n,T>*, point*>> p_neighbor_points;
						for(auto &section : p_neighbor_sections) {
							for(auto &q : section->points) {
								if(_distance((*p.second), q) <= eps)
									p_neighbor_points.push_back(make_pair(section,&q));
							}
						}
						if(p_neighbor_points.size() >= minPts+1) {
							neighbor_points.insert(neighbor_points.end(), p_neighbor_points.begin(), p_neighbor_points.end());
							neighbor_points_count += p_neighbor_points.size();
						}
					}
				} else {
					p[label] = noise;
					continue;
				}
			}
		}
		#if(TEST)
			map<size_t, size_t> clusters;
			for(auto const& [key, node] : location_table_leaf_nodes) {
				for(auto &p : node->points) {
					if(p[label] != fresh && p[label] != noise && p[label] != eliminated) {
						++clusters[p[label]];
					}
				}
			}
			printf(CYAN "clusters:\n" RESET);
			bool nocluster = true;
			for(auto const& [c, size] : clusters) {
				nocluster = false;
				printf( CYAN "| " RESET "cluster %i has %i items\n", c, size);
			}
			if(nocluster) printf(BOLDRED "none\n" RESET);
		#endif
	}
	void DBScan_preserve(double_t eps, size_t minPts) {
		#if(DEBUG)
		size_t iterations;
		#endif
		size_t c = cluster_id_seed;

		// iterate all points
		for(auto const& [key, node] : location_table_leaf_nodes) {
			if(node->points.size() == 0) continue;
			set<spatial_node<n,T>*> neighbor_sections = _find_all_same_level_neighbors(node);
			neighbor_sections.insert(node); // include points in current section
			for(auto &p : node->points) {
				if(p[label] == eliminated || p[label] != fresh) continue; // <= cluster_id_seed: if smaller than this number, then it's a previously detected cluster. expand on that
				vector<pair<spatial_node<n,T>*, point*>> neighbor_points;
				// populate neighbor_points (range_query)
				for(auto &section : neighbor_sections) {
					for(auto &q : section->points) {
						if(_distance(p, q) <= eps)
							neighbor_points.push_back(make_pair(section,&q));
					}
				}

				if(neighbor_points.size() >= minPts+1) {
					p[label] = ++c;
					size_t neighbor_points_count = neighbor_points.size();
					for(size_t i=0; i<neighbor_points_count; ++i) {
						pair<spatial_node<n,T>*, point*> &p = neighbor_points[i]; // next element
						if((*p.second)[label] == noise) {
							(*p.second)[label] = c;
							continue;
						}
						if((*p.second)[label] != fresh) continue; // skip if already processed and added to another cluster
						(*p.second)[label] = c;
						vector<spatial_node<n,T>*> p_neighbor_sections = _find_all_same_level_neighbors(p.first); // todo: this can be cached
						p_neighbor_sections.push_back(p.first); // include points in current section
						vector<pair<spatial_node<n,T>*, point*>> p_neighbor_points;
						for(auto &section : p_neighbor_sections) {
							for(auto &q : section->points) {
								if(_distance((*p.second), q) <= eps)
									p_neighbor_points.push_back(make_pair(section,&q));
							}
						}
						if(p_neighbor_points.size() >= minPts+1) {
							neighbor_points.insert(neighbor_points.end(), p_neighbor_points.begin(), p_neighbor_points.end());
							neighbor_points_count += p_neighbor_points.size();
						}
					}
				} else {
					p[label] = noise;
					continue;
				}
			}
		}
		if(preserve_clusers) { cluster_id_seed = c; }
		printf("cluster_id_seed = %i\n", cluster_id_seed);
		

		#if(TEST)
			map<size_t, size_t> clusters;
			for(auto const& [key, node] : location_table_leaf_nodes) {
				for(auto &p : node->points) {
					if(p[label] != fresh && p[label] != noise && p[label] != eliminated) {
						++clusters[p[label]];
					}
				}
			}
			printf(CYAN "clusters:\n" RESET);
			bool nocluster = true;
			for(auto const& [c, size] : clusters) {
				nocluster = false;
				printf( CYAN "| " RESET "cluster %i has %i items\n", c, size);
			}
			if(nocluster) printf(BOLDRED "none\n" RESET);
		#endif
	}

	// this won't be useful for finding surface clusters because it will only check nearest neighbors for difference in z. but points far away in 3D space may be close in terms of z
	// todo: insert all points in an additional 1D space based on z positions, or implement a neighbor-finding function for all nearby z sections. whichever will be faster
	void DBScan_z(double_t eps, size_t minPts) {
		#if(DEBUG)
		size_t iterations;
		#endif
		// clear previous labels
		for(auto const& [key, node] : location_table_leaf_nodes) {
			for(auto &p : node->points)
				p[label] = fresh;
		}
		size_t c = 2; // cluster index
		// iterate all points
		for(auto const& [key, node] : location_table_leaf_nodes) {
			if(node->points.size() == 0) continue;
			set<spatial_node<n,T>*> neighbor_sections = _find_all_same_level_neighbors(node);
			neighbor_sections.insert(node); // include points in current section
			for(auto &p : node->points) {
				if(p[label] != fresh || p[label] == eliminated) continue;
				vector<pair<spatial_node<n,T>*, point*>> neighbor_points;
				// populate neighbor_points (range_query)
				for(auto &section : neighbor_sections) {
					for(auto &q : section->points) {
						if(_distance_z(p, q) <= eps)
							neighbor_points.push_back(make_pair(section,&q));
					}
				}

				if(neighbor_points.size() >= minPts+1) {
					p[label] = ++c;
					size_t neighbor_points_count = neighbor_points.size();
					for(size_t i=0; i<neighbor_points_count; ++i) {
						pair<spatial_node<n,T>*, point*> &p = neighbor_points[i]; // next element
						if((*p.second)[label] == noise) {
							(*p.second)[label] = c;
							continue;
						}
						if((*p.second)[label] != fresh) continue; // skip if already processed and added to another cluster
						(*p.second)[label] = c;
						set<spatial_node<n,T>*> p_neighbor_sections = _find_all_same_level_neighbors(p.first);
						p_neighbor_sections.insert(p.first); // include points in current section
						vector<pair<spatial_node<n,T>*, point*>> p_neighbor_points;
						for(auto &section : p_neighbor_sections) {
							for(auto &q : section->points) {
								if(_distance_z((*p.second), q) <= eps)
									p_neighbor_points.push_back(make_pair(section,&q));
							}
						}
						if(p_neighbor_points.size() >= minPts+1) {
							neighbor_points.insert(neighbor_points.end(), p_neighbor_points.begin(), p_neighbor_points.end());
							neighbor_points_count += p_neighbor_points.size();
						}
					}
				} else {
					p[label] = noise;
					continue;
				}
			}
		}



		#if(TEST)
			map<size_t, size_t> clusters;
			for(auto const& [key, node] : location_table_leaf_nodes) {
				for(auto &p : node->points) {
					if(p[label] != fresh && p[label] != noise) {
						++clusters[p[label]];
					}
				}
			}
			printf(CYAN "clusters:\n" RESET);
			bool nocluster = true;
			for(auto const& [c, size] : clusters) {
				nocluster = false;
				printf( CYAN "| " RESET "cluster %i has %i items\n", c, size);
			}
			if(nocluster) printf(BOLDRED "none\n" RESET);
		#endif
	}



	double area(point_3D a, point_3D b, point_3D c){
		return (b[0] - a[0]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[0] - a[0]);
	}

	double area(point_2D a, point_2D b, point_2D c){
		return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
	}

	void print_point(point p){
		for(auto it = p.begin(); it != p.end(); ++it)
			cout << *it << " ";
		cout << endl;
	}

	void print_points(vector<point> V, string note){
		cout << note << endl;
		for(auto it = V.begin(); it != V.end(); ++it)
			print_point(*it);
	}

	vector<point> find_convex_hull(vector<point> P){
		sort(all(P));
		vector<point> upper, lower;
		upper.push_back(P[0]);
		lower.push_back(P[0]);
		for(int i = 0; i < P.size(); i++){
			if(i == P.size() - 1 || this->area(P[0], P[i], P.back()) > 0){
				while(upper.size() >= 2 && this->area(upper[upper.size()-2], upper.back(), P[i]) < 0)
					upper.pop_back();
				upper.push_back(P[i]);
			}
			if(i == P.size() - 1 || this->area(P[0], P[i], P.back()) < 0){
				while(lower.size() >= 2 && this->area(lower[lower.size()-2], lower.back(), P[i]) > 0)
					lower.pop_back();
				lower.push_back(P[i]);
			}
		}
		for(int i = lower.size() - 2; i >= 1; i--)
			upper.push_back(lower[i]);
		return upper;
	}
	
	// Old version of the algorithm that detects surface points
	// Refrain from using this on production
	double_t remove_surface_points(double_t r){
		// her surface point için "point[label] = eliminated;"
		map<int, int> mp;
		for(auto const& [key, node]: location_table_leaf_nodes)
			for(auto &p: node->points)
				mp[int(p[1] * 100)]++;
		int surfaceZ = 0, maxi = 0;
		
		for(auto it = mp.begin(); it != mp.end(); ++it)
			if(it->second > maxi){
				maxi = it->second;
				surfaceZ = it->first;
			}
			// r *= 100;
			r = 2.0;
			// double r = 2.0;
			// double r = 0.5;
			int s1 = 0, s2 = 0, s3 = 0;
			for(auto const& [key, node]: location_table_leaf_nodes)
				for(auto &p: node->points){
					int x = int(p[1] * 100);
					if(surfaceZ - r <= x && x <= surfaceZ + r) {
						p[label] = eliminated;
						s3++;
					}
					else {
						s1++;
						if(x < surfaceZ - r)
							s2++;
					}
				}

			cout << "#removed = " << s3 << endl;

			if(s2 * 3 < s1)
				for(auto const& [key, node]: location_table_leaf_nodes)
					for(auto &p: node->points)
						if(int(p[1] * 100) <= surfaceZ + r)
							p[label] = eliminated;
		return ((double_t)surfaceZ)/100;
	}

	double __distance(point x, point y){
		return (x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) + (x[2] - y[2]) * (x[2] - y[2]);
	}
	vector<point> remove_distant_points(vector<point> pts){

		point x;
		int s = 0;
		x[0] = 0; x[1] = 0; x[2] = 0;
		for(auto it = pts.begin(); it != pts.end(); ++it){
			x[0] += (*it)[0]; x[1] += (*it)[1]; x[2] += (*it)[2];
			s++;
		}
		x[0] /= s; x[1] /= s; x[2] /= s;

		vector<pair<double, int>> V;
		for(int i = 0; i < V.size(); i++)
			V.push_back(make_pair(__distance(V[i], x), i));

		sort(all(V));
		int size = max(0, int(V.size() * 9 / 10));
		vector<point> res;
		for(int i = 0; i < size; i++)
			res.push_back(pts[V[i].second]);
		return res;
	}

	double_t surface = 0;
	vector<Box<n, T>> query_3D(double_t eps, size_t minPts, double_t surface_r){
		if(preserve_clusers) this->clear_labels_preserve_clusters();
		else this->clear_labels();
		this->surface = this->remove_surface_points(surface_r);
		this->DBScan(eps, minPts);

		unordered_map<size_t, vector<point>> clusters;
		for(auto const& [key, node]: location_table_leaf_nodes)
			for(auto &p: node->points)
				if(p[label] > 2)
					clusters[p[label]].push_back(p);

		vector<Box<n, T>> results;

		for(auto const& [label_, pts]: clusters){
			// pts = remove_distant_points(pts);
			vector<point> convex_hull = find_convex_hull(pts);
			vector<double> rect = find_smallest_rectangle(convex_hull);
			double minZ = INF, maxZ = -INF;
			for(auto i = 0; i < pts.size(); i++){
				if(pts[i][1] < minZ)
					minZ = pts[i][1];
				if(pts[i][1] > maxZ)
					maxZ = pts[i][1];
			}
			Box<n, T> box;
			box.label_ = label_;
			box.center = {rect[2], (maxZ+minZ)/2, rect[3]};
			box.shape = {rect[0], fabs(maxZ-minZ), rect[1]};
			box.angle = rect[4];
			results.push_back(box);
		}
		return results;
	}

	vector<Box<3, T>> query_2D(double_t eps, size_t minPts, array<double_t,3> cam_pos){
		if(preserve_clusers) this->clear_labels_preserve_clusters();
		else this->clear_labels();
		this->boxScan(eps, minPts, cam_pos);

		unordered_map<size_t, vector<point>> clusters;
		for(auto const& [key, node]: location_table_leaf_nodes)
			for(auto &p: node->points)
				if(p[label] > 2)
					clusters[p[label]].push_back(p);
		
		vector<Box<3, T>> results;

		for(auto const& [label_, pts]: clusters){
			// pts = remove_distant_points(pts);
			vector<point> convex_hull = find_convex_hull(pts);
			vector<double> rect = find_smallest_rectangle(convex_hull);
			double minZ = INF, maxZ = -INF;
			std::set<double_t> Z;
			for(auto i = 0; i < pts.size(); i++){
				Z.insert(pts[i][label-1]);
			}
			
			if(Z.size() > 6) {
				minZ = *(std::next(Z.begin(),2)); // 3rd
				maxZ = *(std::prev(Z.end(),3)); // last 3rd
			} else {
				minZ = *Z.begin();
				maxZ = *(std::prev(Z.end(),1));
			}

			Box<3, T> box;
			box.label_ = label_;
			box.center = {rect[2], (maxZ+minZ)/2, rect[3]};
			box.shape = {rect[0], fabs(maxZ-minZ), rect[1]};
			box.angle = rect[4];

			results.push_back(box);
		}
		return results;
	}


	ppoint rotate_around_origin(ppoint pt, double angle){
		double s = sin(angle);
		double c = cos(angle);
		return make_pair(pt.first * c - pt.second * s, pt.first * s + pt.second * c);
	}

	std::vector<double> find_smallest_rectangle(vector<point_3D> pts){
		int N_ROTATIONS = 1000;
		double PI = 3.14159265;
		double min_area = 1e9;
		double best_angle = -1;
		double sizeX, sizeY;
		ppoint pt1, pt2;
		for(int i = 0; i < N_ROTATIONS; i++){
			double angle = (2.0 * PI / N_ROTATIONS) * i;
			double maxX = -INF, maxY = -INF;
			double minX = INF, minY = INF;
			for(auto it = pts.begin(); it != pts.end(); ++it){
				ppoint tmp = this->rotate_around_origin(make_pair((*it)[0], (*it)[2]), angle);
				maxX = max(maxX, tmp.first);
				minX = min(minX, tmp.first);
				maxY = max(maxY, tmp.second);
				minY = min(minY, tmp.second);
			}
			double area = (maxX - minX) * (maxY - minY);
			if(area < min_area){
				min_area = area;
				best_angle = angle;
				pt1 = this->rotate_around_origin(make_pair(minX, maxY), -angle);
				pt2 = this->rotate_around_origin(make_pair(maxX, minY), -angle);
				sizeX = maxX - minX;
				sizeY = maxY - minY;
			}
		}
		// assert(best_angle != -1);
		if(best_angle == -1) best_angle = 0;
		vector<double> res {sizeX, sizeY, (pt1.first + pt2.first) / 2, (pt1.second + pt2.second) / 2, best_angle};
		return res;
	}
	std::vector<double> find_smallest_rectangle(vector<point_2D> pts){
		int N_ROTATIONS = 1000;
		double PI = 3.14159265;
		double min_area = 1e9;
		double best_angle = -1;
		double sizeX, sizeY;
		ppoint pt1, pt2;
		for(int i = 0; i < N_ROTATIONS; i++){
			double angle = (2.0 * PI / N_ROTATIONS) * i;
			double maxX = -INF, maxY = -INF;
			double minX = INF, minY = INF;
			for(auto it = pts.begin(); it != pts.end(); ++it){
				ppoint tmp = this->rotate_around_origin(make_pair((*it)[0], (*it)[1]), angle);
				maxX = max(maxX, tmp.first);
				minX = min(minX, tmp.first);
				maxY = max(maxY, tmp.second);
				minY = min(minY, tmp.second);
			}
			double area = (maxX - minX) * (maxY - minY);
			if(area < min_area){
				min_area = area;
				best_angle = angle;
				pt1 = this->rotate_around_origin(make_pair(minX, maxY), -angle);
				pt2 = this->rotate_around_origin(make_pair(maxX, minY), -angle);
				sizeX = maxX - minX;
				sizeY = maxY - minY;
			}
		}
		// assert(best_angle != -1);
		if(best_angle == -1) best_angle = 0;
		vector<double> res {sizeX, sizeY, (pt1.first + pt2.first) / 2, (pt1.second + pt2.second) / 2, best_angle};
		return res;
	}
};

void randomize_points(vector<point_3D>* points, int seed, int count) {
	srand(seed);
	bool expand = false;
	size_t n = count/3;
	size_t sl = 10; // sidelength
	for(int k=0; k<3; ++k) {
		for(int i=0; i<n; ++i) {
			points->push_back({(double_t)(int)(rand() % sl - sl/2),(double_t)(int)(rand() % sl - sl/2),(double_t)(int)(rand() % sl - sl/2)});
		}
		if(expand) sl *= 2;
	}
}
void randomize_points(vector<point_2D>* points, int seed, int count) {
	srand(seed);
	bool expand = false;
	size_t n = count/3;
	size_t sl = 1000; // sidelength
	for(int k=0; k<3; ++k) {
		for(int i=0; i<n; ++i) {
			points->push_back({(double_t)(int)(rand() % sl - sl/2),(double_t)(int)(rand() % sl - sl/2)});
		}
		if(expand) sl *= 2;
	}
}

#if(MAIN)
int main(int argc, char** argv)
{
	if(0) {
		double_t* c_point = new double_t[3]{1,2,3};
		point_3D p;
		std::copy(c_point, c_point + 3, p.begin());
		printArray(p, 4, "array:\n");
		return 0;
	}
	if(0) {
		printf("mod: %f\n", std::fmod(160.1, 5.0));
		
		return 0;
	}
	if(0) {
		size_t x = 0b0101100;
		printf("%i\n", x & -x);
		printf("%i\n", (size_t)(log2(x & -x)));
		printf("%i\n", firstsetbit(x));
		printf("%i\n", firstNOTsetbit(x));
		return 0;
	}
	if(0) {
		std::map<size_t, size_t> asdf;
		asdf[1] = 5;
		asdf[2] = 6;
		asdf[1] = 3;
		for(auto const& [key, val] : asdf) {
			printf("%i %i\n", key, val);
		}
		return 0;
	}
	if(0) {
		// size_t x = 17453 << 25;
		uint64_t x = (1ul << 64-1);
		string bits = getbits(x,64);
		printf("bits: %s %lu\n", bits.c_str(), x);
		return 0;
	}
	if(0) {
		for(size_t i=0; i<8; ++i) {
			size_t d = i/2;
			bool sign = i % 2;
			printf("[%i, %i]\n", d, sign);
		}
		return 0;
	}
	if(0) {
		int a = -5;
		for(; a < 0; ++a) {
			printf("a:%i\n", a);
			// if(a==-2) break;
		}
		return 0;
	}
	if(0) {
		SpatialTree<2, double> space(10);
		size_t locationcode = space._construct_location_code({0b100,0b001},3,-1,3);
		printf("new location code:\n");
		cout << getbits(locationcode >> 2, 2) << " " << getbits(lastnbits(locationcode, 2), 2);
		cout << endl;
		return 0;
	}
	if(0) {
		size_t x = 0b0110;
		size_t reverse = reverse_lastnbits(x,2);
		printf("\nreverse of last %i bits of %s = %s\n", 2, dumpbits(&x, 5, "").c_str(), dumpbits(&reverse, 5, "").c_str());
		return 0;
	}
	if(0) {
		int a = 4 + 2;
		dumpbits(&a, " ");
		printf("\n");
		return 0;
	}
	if(0) {
		size_t nodeheight = 3;
		int relative_height = -1;
		size_t location_code = 0;
		size_t height = nodeheight + relative_height;
		vector<size_t> location_code_v = {0b011, 0b011}; // resulting location_code should be 0b0101 = 5
		for(size_t i=0; i<2; ++i) {
			location_code <<= height;
			if(relative_height < 0)
			location_code += location_code_v[i] >> -relative_height;
			else
			location_code += location_code_v[i];
		}
		printf("location_code:%i\n", location_code);
		return 0;
	}
	if(0){
		srand(time(NULL));
		size_t n = 10;
		size_t sl = 10; // sidelength
		// for(int i=0; i<5; ++i) {
		// 	printf("{%i,%i,%i},\n", rand() % sl - sl/2,rand() % sl - sl/2,rand() % sl - sl/2);
		// }
		for(int i=0; i<n; ++i) {
			printf("{%i,%i},\n", rand() % sl - sl/2,rand() % sl - sl/2);
		}
		sl = 20;
		for(int i=0; i<n; ++i) {
			printf("{%i,%i},\n", rand() % sl - sl/2,rand() % sl - sl/2);
		}
		sl = 30;
		for(int i=0; i<n; ++i) {
			printf("{%i,%i},\n", rand() % sl - sl/2,rand() % sl - sl/2);
		}
		
		return 0;
	}
	if(0){
		SpatialTree<3, double> space(10);
		space.capacity = 1;
		#if(DEBUG_VISUALIZE)
		// test(space);
		#endif
		printf("\n%p\n", space.locate({-2,1,4}));
		return 0;
	}
	#if(TEST)
	if(0) {
		auto space_col = (space_store*)create_3D_space_and_projection(40, 0.03, true);
		using boost::property_tree::ptree;
		using boost::property_tree::read_json;
		using boost::property_tree::write_json;
		
		ptree pt;
		read_json("points.json", pt);
		auto points = pt.get<string> ("points");
		
		printf("points: %s\n", points.c_str());
		// for(const auto& x: points) {
			// 	std::cout << x.first << ": " << x.second.get_value<std::string>() << std::endl;
			// 	printf("lel");
		// }
		return 0;
	}
	if(1) {
		using json = nlohmann::json;
		
		auto space_col = (space_store*)create_3D_space_and_projection(40, 0.03, true);
		auto space = (SpatialTree<3, double_t>*)(space_col->space);
		auto projection = (SpatialTree<2, double_t>*)(space_col->projection);
		auto filename = "points.json";
		std::ifstream i(filename);
		json j;
		i >> j;
		// cout << j["points"] << '\n';
		std::printf("read %i points from %s\n", j["points"].size(), filename);
		// for (json::iterator it = j.begin(); it != j.end(); ++it) {
			// 	std::cout << *it << '\n';
			// }
			auto points = j["points"];
			for (auto& p : points) { // insert and return bulk function ported
				point_3D pt = p;
				if(space->check_dup(p, 0.005)) continue; // skip duplicate point
				if(space->insert(p)) {
					// insert to projection if not duplicate in the original
					point_2D p_proj = {p[0], p[2], p[1], 0}; // store y dimension as extra info
					projection->insert(p_proj);
				}
			}
			
			query_api1(0.04, 20, 0, space_col);
			
			vector<point_2D> pts_new;
			projection->get_points(pts_new);
			
			
			j["points"] = pts_new;
			std::ofstream o("points_annotated.json");
			o << j << std::endl;
			return 0;
		}
		#endif
		if(0) {
			#if(DEBUG_VISUALIZE_3D)
			std::vector<std::vector<double>> x, y, z;
			for (double i = -5; i <= 5;  i += 0.25) {
				std::vector<double> x_row, y_row, z_row;
				for (double j = -5; j <= 5; j += 0.25) {
					x_row.push_back(i);
					y_row.push_back(j);
					z_row.push_back(::std::sin(::std::hypot(i, j)));
				}
				x.push_back(x_row);
				y.push_back(y_row);
				z.push_back(z_row);
			}

			plt::plot_surface(x, y, z);
			plt::show();
		#endif
	}

	#if(DEBUG_VISUALIZE)
		draw();
	#endif
	if(0) {
		// space.insert((point_3D){ 1 , 1, 1});
		// space.insert((point_3D){-1,  1, 1});
		// space.insert((point_3D){ 1 ,-1, 1});
		// space.insert((point_3D){-1 ,-1, 1});
		// space.insert((point_3D){ 1 , 1,-1});
		// space.insert((point_3D){-1 , 1,-1});
		// space.insert((point_3D){ 1 ,-1,-1});
		// space.insert((point_3D){-1 ,-1,-1});
		// space.insert((point_3D){0,0,0});
		// printf("\n%p\n", space.locate({-6,8,9}));
	}

	SpatialTree<3, double> space(160); // 5 10 20 40 80 160 320 640 1280 2560 ... ideal sidelengths
	space.max_sidelength = 5;
	space.min_sidelength = 5;
	space.capacity = 1;
	space.fixed_level = true;
	vector<point_3D> outsider_points;
	int seed = 3;
	if(argc > 1) seed = atoi(argv[1]);
	size_t n = 1500;
	printf("randomize seed: %i\n", seed);
	printf("testing with " RED "%i" RESET " points\n", n);
	randomize_points(&outsider_points, seed, n);
	int i=4;
	auto start = chrono::high_resolution_clock::now();
	for(auto &p : outsider_points) {
		// cin.ignore();
		#if(DEBUG_VISUALIZE)
			// plt::pause(0.1);
		#endif
		space.insert(p);
		// auto v = space.locate(p);
		// space.find_neighbors(p);
		#if(DEBUG_VISUALIZE)
		// if(!(i++ % 20)) {
		// 	space.__drawAllRelativeHeights();
		// 	plt::pause(10);
		// 	plt::show();
		// }
		#endif
	}
	#if(TEST)
		printf("level 0 point count = %i\n", space.root->points.size());
		size_t leaf_point_count = 0;
		map<int, int> levels;
		map<double_t, int> sizes;
		for(auto& s : space.location_table) {
			if(s.second->points.size() > 0) {
				if(s.second->is_leaf)
					leaf_point_count += s.second->points.size();
				levels[s.second->height] += s.second->points.size();
				sizes[s.second->sidelength]++;
			}
		}
		printf("leaf level point count = %i\n", leaf_point_count);
		for(auto &l : levels) {
			printf("level %i point count: %i\n", l.first, l.second);
		}
		for(auto &l : sizes) {
			printf("	sidelength %f count: %i\n", l.first, l.second);
		}
	#endif
	auto end = chrono::high_resolution_clock::now();
	auto ms = chrono::duration_cast<chrono::milliseconds>(end-start).count();
	printf("Construction time: " YELLOW "%lli" RESET " milliseconds\n", ms);

	#if(DEBUG_VISUALIZE)
		// space.__drawAllRelativeHeights();
		// plt::show();
		// plt::pause(10000);
	#endif

	// SpatialTree::spatial_node* node = locate(outsider_points[10]);
	// space.find_neighbors(outsider_points[10]);


	start = chrono::high_resolution_clock::now();
	space.walk();
	end = chrono::high_resolution_clock::now();
	ms = chrono::duration_cast<chrono::milliseconds>(end-start).count();
	printf("Walk time: " GREEN "%lli" RESET " milliseconds\n", ms);

	start = chrono::high_resolution_clock::now();
	space.DBScan(5, 5);
	end = chrono::high_resolution_clock::now();
	ms = chrono::duration_cast<chrono::milliseconds>(end-start).count();
	printf("DBScan time: " GREEN "%lli" RESET " milliseconds\n", ms);

	double_t eps = 5;
	size_t minPts = 5;
	space.query_3D(eps, minPts,2);
	vector<point_3D> new_points;
	randomize_points(&new_points, seed+1, n);
	for(auto &p : new_points) {
		p[0]+=1;p[1]+=1;p[2]+=1;
		space.insert(p);
	}
	space.query_3D(eps, minPts,2);
	randomize_points(&new_points, seed+2, n);
	for(auto &p : new_points) {
		p[0]+=1;p[1]+=1;p[2]+=1;
		space.insert(p);
	}
	space.query_3D(eps, minPts,2);

	#if(DEBUG_VISUALIZE)
		plt::pause(0.01);
		plt::show(); // don't quit
	#endif
	if(0) {
		start = chrono::high_resolution_clock::now();
		int k=0;
		for(size_t i=0; i<22678*1500;++i) {
			k+=i*i;
		}
		end = chrono::high_resolution_clock::now();
		ms = chrono::duration_cast<chrono::milliseconds>(end-start).count();
		printf("Elapsed time: " GREEN "%lli" RESET " milliseconds\n", ms);
		printf("%i\n", k);
	}
	if(0) {
		start = chrono::high_resolution_clock::now();
		int k=0;
		for(size_t i=0; i<1500*1500*1500/8;++i) {
			k+=i*i;
		}
		end = chrono::high_resolution_clock::now();
		ms = chrono::duration_cast<chrono::milliseconds>(end-start).count();
		printf("Elapsed time: " GREEN "%lli" RESET " milliseconds\n", ms);
		printf("%i\n", k);
	}
	return 0;
}
#endif



void swifttestF() {
	printf("Calling From Swift to C to c++");
	cout << "It Works!!\n";
}


void pass_string(const char* deneme) {
	cout << "pass string!: " << deneme << endl;
}

void pass_array(const double* arr, size_t count) {
	cout << "pass_array *" << arr << ":\n";
	printArray(arr, count);
}

// doesn't work
void pass_2d_array(const double_t** arr, size_t c1, size_t c2) {
	cout << "pass_2d_array *" << arr << ":\n";
	printArray_2d(arr, c1, c2);
}

/* ************************************************************************** */
/*                                   BRIDGE                                   */
/* ************************************************************************** */

// ++ lolol I don't need to restrict to a fixed epsilon, I can simply go one less step down the section levels to use higher epsilon values!
// ++ so for sidelength 5 and epsilon > 5: if the leaf level is 7, just go to level 6 and look for neighbors. gg ez, no extra computation.
// ++ The only extra thing that needs to be done is to keep a location_table for all levels, which is n*log(s) in computation and n in storage

// namespace api0
// {
	void* create_3D_space_double_t(double_t side_length, double_t node_side_length, bool optimize_side_length) {
		double_t sl;
		// node_side_length *= 2; // todo: this is nonsense lol. make dbscan epsilon constant, or give it a max value and determine node_side_length accordingly
		if(optimize_side_length) {
			// calculate exact side_length to ensure minimum side_length
			// rather than coping with floating point edge cases, simulate insertion side halving
			sl = node_side_length;
			do sl *= 2;
			while(sl < side_length);
		} else sl = side_length;
		SpatialTree<3, double_t>* space = new SpatialTree<3, double_t>(sl);
		space->max_sidelength = node_side_length;
		space->min_sidelength = node_side_length;
		space->capacity = 1; // todo: removing this for fixed height will improve performance. By a tiny amount :P
		space->fixed_level = true;
		#if(BRIDGE_DEBUG)
			printf("space: [side_length: %f] [node_side_length: %f]\n", sl, node_side_length);
		#endif
		return space;
	}

	// todo: how do I bridge a templated function?
	// template<size_t n, typename T>
	// SpatialTree<n,T> create_space() {
	// 	return what?
	// }
	// void* create_space_float(size_t n) {

	// }


	bool insert_3D_no_check(const double_t* c_point, void* space) {
		point_3D p;
		std::copy(c_point, c_point + 3, p.begin());
		return ((SpatialTree<3,double_t>*)space)->insert(p);
	}

	bool insert_3D(const double_t* c_point, double_t dup_dist, void* space) {
		point_3D p;
		std::copy(c_point, c_point + 3, p.begin());
		if(((SpatialTree<3,double_t>*)space)->check_dup(p, dup_dist)) return true;
		return ((SpatialTree<3,double_t>*)space)->insert(p);
	}

	bool insert_bulk_3D_no_check(const double_t* c_flat_array, size_t n, void* space) {
		bool res = true;
		for(auto it = c_flat_array; it != c_flat_array + 3*n; it+=3) {
			point_3D p;
			std::copy(it, it+3, p.begin());
			if(!((SpatialTree<3,double_t>*)space)->insert(p)) {res = false;}
		}
		return res;
	}

	bool insert_bulk_3D(const double_t* c_flat_array, size_t n, double_t dup_dist, void* space) {
		bool res = true;
		for(auto it = c_flat_array; it != c_flat_array + 3*n; it+=3) {
			point_3D p;
			std::copy(it, it+3, p.begin());
			if(((SpatialTree<3,double_t>*)space)->check_dup(p, dup_dist)) continue; // skip duplicate point
			if(!((SpatialTree<3,double_t>*)space)->insert(p)) {res = false;}
		}
		return res;
	}


	c_array_double* insert_and_return_bulk_3D(const double_t* c_flat_array, size_t n, double_t dup_dist, void* s) {
		c_array_double* res = new c_array_double;
		res->data = new double_t[3*n];
		res->count = 0;
		SpatialTree<3,double_t>* space = (SpatialTree<3,double_t>*)s;
		for(auto it = c_flat_array; it != c_flat_array + 3*n; it+=3) {
			point_3D p;
			std::copy(it, it+3, p.begin());
			if(space->check_dup(p, dup_dist)) continue; // skip duplicate point
			if(space->insert(p)) {
				std::copy(it, it+3, (res->data + res->count*3));
				++res->count;
			}
		}
		return res;
	}

	void print_points(void* space) {
		vector<point_3D> points;
		((SpatialTree<3, double_t>*)space)->get_points(points);
		printf("points in 3D space:");
		for(auto &p : points) {
			printArray(p,3, "\n");
		}
	}

	size_t count_points(void* space) {
		size_t res = 0;
		for(auto &s : ((SpatialTree<3,double_t>*)space)->location_table_leaf_nodes) {
			res += s.second->points.size();
		}
		return res;
	}

	c_array_double* get_all_points(void* s) {
		SpatialTree<3, double_t> space = *((SpatialTree<3, double_t>*)s);
		vector<point_3D> points;
		space.get_points(points);
		c_array_double* res = new c_array_double;
		res->data = new double_t[4*points.size()];
		for(size_t i = 0; i < points.size(); ++i) {
			std::copy(points[i].begin(), points[i].end(), (res->data + 4*i));
		}
		res->count = points.size();
		return res;
	}


	c_box_array query(double_t eps, size_t minPts, double_t surface_r, void* s) {
		SpatialTree<3, double_t>* space = (SpatialTree<3, double_t>*)s;
		vector<Box<3, double_t>> boxes = space->query_3D(eps, minPts, surface_r);
		c_box_array c_boxes;
		c_boxes.data = new c_box[boxes.size()];
		c_boxes.count = boxes.size();
		c_boxes.surface = space->surface;

		for(size_t i=0; i<boxes.size(); ++i) {
			c_boxes.data[i].label_ = boxes[i].label_;
			c_boxes.data[i].angle = boxes[i].angle;
			std::copy(boxes[i].center.begin(), boxes[i].center.end(), c_boxes.data[i].center);
			std::copy(boxes[i].shape.begin(), boxes[i].shape.end(), c_boxes.data[i].shape);
		}
		return c_boxes;
	}
// }
// namespace api1
// {
	void* create_3D_space_and_projection(double_t side_length, double_t node_side_length, bool optimize_side_length) {
		double_t sl;
		if(optimize_side_length) {
			// calculate exact side_length to ensure minimum side_length
			// rather than coping with floating point edge cases, simulate insertion side halving
			sl = node_side_length;
			do sl *= 2;
			while(sl < side_length);
		} else sl = side_length;
		space_store* space_col = new space_store;

		SpatialTree<3, double_t>* space = new SpatialTree<3, double_t> (sl);
		space->max_sidelength = node_side_length;
		space->min_sidelength = node_side_length;
		space->capacity = 1;
		space->fixed_level = true;
		#if(BRIDGE_DEBUG)
			printf("space: [side_length: %f] [node_side_length: %f]\n", sl, node_side_length);
		#endif

		SpatialTree<2, double_t>* projection = new SpatialTree<2, double_t> (sl);
		projection->max_sidelength = node_side_length;
		projection->min_sidelength = node_side_length;
		projection->capacity = 1;
		projection->fixed_level = true;
		#if(BRIDGE_DEBUG)
			printf("projection: [side_length: %f] [node_side_length: %f]\n", sl, node_side_length);
		#endif

		space_col->space = space;
		space_col->projection = projection;

		return space_col;
	}


	c_array_double* insert_and_return_bulk(const double_t* c_flat_array, size_t n, double_t dup_dist, void* space_col) {
		c_array_double* res = new c_array_double;
		res->data = new double_t[3*n];
		res->count = 0;

		SpatialTree<3,double_t>* space = (SpatialTree<3,double_t>*)((space_store*)space_col)->space;
		SpatialTree<2,double_t>* projection = (SpatialTree<2,double_t>*)((space_store*)space_col)->projection;

		for(auto it = c_flat_array; it != c_flat_array + 3*n; it+=3) {
			point_3D p;
			std::copy(it, it+3, p.begin());
			if(space->check_dup(p, dup_dist)) continue; // skip duplicate point
			if(space->insert(p)) {
				std::copy(it, it+3, (res->data + res->count*3));
				++res->count;
				// insert to projection if not duplicate in the original
				point_2D p_proj = {p[0], p[2], p[1], 0}; // store y dimension as extra info
				projection->insert(p_proj);
			}
		}
		return res;
	}

	c_box_array query_api1(double_t eps, size_t minPts, const double_t* cam_pos, void* space_col) {
		array<double_t,3> cp;
		std::copy(cam_pos, cam_pos+3, cp.begin());
		SpatialTree<3, double_t>* space = (SpatialTree<3, double_t>*)((space_store*)space_col)->space;
		SpatialTree<2, double_t>* projection = (SpatialTree<2, double_t>*)((space_store*)space_col)->projection;
		vector<Box<3, double_t>> boxes = projection->query_2D(eps, minPts, cp);

		for(auto& box : boxes) {
			printf("box height: %.2f\n", box.shape[1]);
		}
		c_box_array c_boxes;
		c_boxes.data = new c_box[boxes.size()];
		c_boxes.count = boxes.size();
		c_boxes.surface = 0;

		for(size_t i=0; i<boxes.size(); ++i) {
			c_boxes.data[i].label_ = boxes[i].label_;
			c_boxes.data[i].angle = boxes[i].angle;
			std::copy(boxes[i].center.begin(), boxes[i].center.end(), c_boxes.data[i].center);
			std::copy(boxes[i].shape.begin(), boxes[i].shape.end(), c_boxes.data[i].shape);
		}
		return c_boxes;
	}
	void print_json(vector<point_3D> &points) {
		printf("==== points: ====\n");
		printf("[");
		for(size_t i=0; i<points.size(); ++i) {
			printf("[");
			for(size_t j=0; j<points[i].size(); ++j) {
				printf("%.2f", points[i][j]);
				if(j!=points[i].size()-1) printf(",");
			}
			if(i!=points.size()-1) printf("],");
			else printf("]");
		}
		printf("]");
		printf("\n==== points end ====\n");
	}
	c_array_double* get_all_points_api1(void* space_col, bool print) {
		SpatialTree<3, double_t>* space = (SpatialTree<3, double_t>*)((space_store*)space_col)->space;
		vector<point_3D> points;
		space->get_points(points);
		if(print) print_json(points);
		c_array_double* res = new c_array_double;
		res->data = new double_t[4*points.size()];
		for(size_t i = 0; i < points.size(); ++i) {
			std::copy(points[i].begin(), points[i].end(), (res->data + 4*i));
		}
		res->count = points.size();
		return res;
	}
	
	double_t detect_surface_3D(void* s, const double* fp, double_t h_radius, double_t v_threshold, size_t minPts) {
		point_3D focus_point;
		std::copy(fp, fp+3, focus_point.begin());
		SpatialTree<3, double_t>* space = (SpatialTree<3, double_t>*)s;
		return space->_detect_surface(focus_point, h_radius, v_threshold, minPts);
	}


// }
#if(TEST)
	void print_contents( const boost::property_tree::ptree& pt)
	{
			using boost::property_tree::ptree;

			for (const auto& x: pt )
			{
					std::cout << x.first << ": " << x.second.get_value<std::string>() << std::endl;
					print_contents(x.second);
			}
	}
#endif

/* ************************************************************************** */
/*                                END BRIDGE                                  */
/* ************************************************************************** */
