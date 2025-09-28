# SpatialTree

A high-performance C++ implementation of a spatial tree data structure for efficient 3D point cloud processing, clustering, and surface detection.

## Overview

SpatialTree is an optimized spatial indexing library that provides fast insertion, querying, and clustering operations on 3D point data. The implementation uses a hierarchical spatial decomposition approach to efficiently organize and query large point clouds. Key features include:

- **Efficient 3D Spatial Indexing**: Hierarchical tree structure for fast point insertion and retrieval
- **DBSCAN Clustering**: Built-in density-based spatial clustering with optimizations
- **Surface Detection**: Algorithms for detecting surfaces in 3D point clouds
- **Duplicate Point Handling**: Configurable duplicate detection based on distance thresholds
- **Bulk Operations**: Optimized batch insertion for large point datasets
- **Visualization Support**: Integration with matplotlib for 2D/3D visualization (optional)

## Building

### Prerequisites

- C++ compiler with C++11 support (GCC, Clang, or MSVC)
- CMake 3.10 or higher
- Optional: Python with matplotlib for visualization features
- Optional: Boost libraries for JSON parsing (if TEST mode enabled)

### Compilation

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

### Basic Example

```cpp
#include "spatial_tree.h"

// Create a 3D spatial tree
// Parameters: side_length, node_side_length, optimize_side_length
void* space = create_3D_space_double_t(100.0, 10.0, true);

// Insert a single point
double point[3] = {1.0, 2.0, 3.0};
bool success = insert_3D(point, 0.01, space);  // 0.01 is duplicate distance threshold

// Bulk insert points
double points[] = {
    1.0, 2.0, 3.0,
    4.0, 5.0, 6.0,
    7.0, 8.0, 9.0
};
insert_bulk_3D(points, 3, 0.01, space);

// Query with DBSCAN clustering
double eps = 1.5;        // neighborhood radius
size_t minPts = 5;       // minimum points for cluster
double surface_r = 0.1;  // surface detection radius
struct c_box_array clusters = query(eps, minPts, surface_r, space);

// Get all points
struct c_array_double* all_points = get_all_points(space);
```

### API Reference

#### Space Creation
- `create_3D_space_double_t(side_length, node_side_length, optimize_side_length)`: Creates a 3D spatial tree
- `create_3D_space_and_projection(...)`: Creates space with projection capabilities

#### Point Insertion
- `insert_3D(point, dup_dist, space)`: Insert single point with duplicate checking
- `insert_3D_no_check(point, space)`: Insert without duplicate checking
- `insert_bulk_3D(array, n, dup_dist, space)`: Bulk insert with duplicate checking
- `insert_and_return_bulk_3D(...)`: Bulk insert and return actually inserted points

#### Querying
- `query(eps, minPts, surface_r, space)`: Perform DBSCAN clustering
- `get_all_points(space)`: Retrieve all points in the tree
- `count_points(space)`: Get total point count

#### Surface Detection
- `detect_surface_3D(space, focus_point, h_radius, v_threshold, minPts)`: Detect surfaces near a focus point

## Project Structure

```
spatialtree/
├── spatialtree.cpp          # Main implementation
├── include/                 # Header files
│   ├── spatial_tree.h       # Public API
│   ├── spatial_tree_utils.hpp
│   ├── json.hpp            # JSON parsing utilities
│   ├── matplotlibcpp.h     # Matplotlib C++ interface
│   └── ...
├── tests/                   # Test files
│   ├── test.cpp
│   ├── test.py
│   └── generate_random_input.py
├── data/                    # Sample data files
│   └── points*.json
├── visualization/           # Visualization scripts
│   ├── visualize.py
│   └── visualize.ipynb
├── artifacts/              # Build artifacts and temporary files
└── CMakeLists.txt          # Build configuration
```

## Configuration Flags

The implementation provides several compile-time configuration flags in `spatialtree.cpp`:

- `STORE_POINTS_AT_ALL_LEVELS`: Store points at all tree levels (default: 1)
- `DBSCAN_OPTIMIZATIONS`: Enable DBSCAN-specific optimizations (default: 1)
- `DEBUG`: Enable debug output (default: 0)
- `DEBUG_VISUALIZE`: Enable 2D visualization debugging (default: 0)
- `DEBUG_VISUALIZE_3D`: Enable 3D visualization debugging (default: 0)
- `POINTS_CARRY_EXTRA_INFO`: Points can carry additional metadata (default: 1)

## Testing

Run the test suite:

```bash
# Generate random test data
python tests/generate_random_input.py

# Run C++ tests
./build/spatialtree

# Run Python visualization tests
python tests/test.py
```

## Performance Considerations

- The spatial tree is optimized for dense point clouds with locality
- Node side length should be chosen based on expected query radius for DBSCAN
- Bulk insertion is significantly faster than individual insertions
- Setting `optimize_side_length=true` adjusts the tree structure for better DBSCAN performance
