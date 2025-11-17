/*
Note:
    The matrix and vectors in this implementation are very very sparse.
*/


#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <numeric>
#include <sstream>
#include <unordered_map>
#include <Eigen/Dense>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/fmt.h>


template<typename T>
struct fmt::formatter<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> {
    constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }
    
    template<typename FormatContext>
    auto format(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat, FormatContext& ctx) const {
        std::ostringstream oss;
        oss << mat;  // Uses Eigen's operator
        return fmt::format_to(ctx.out(), "{}", oss.str());
    }
};


// For vectors (same thing)
template<typename T>
struct fmt::formatter<Eigen::Matrix<T, Eigen::Dynamic, 1>> {
    constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }
    
    template<typename FormatContext>
    auto format(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vec, FormatContext& ctx) const {
        std::ostringstream oss;
        oss << vec;
        return fmt::format_to(ctx.out(), "{}", oss.str());
    }
};


// breakpoint macro
#ifdef DEBUG_BREAKPOINT
  #if defined(_MSC_VER)
    #define BREAKPOINT __debugbreak()
  #elif defined(__APPLE__)
    #define BREAKPOINT __builtin_debugtrap()
  #elif defined(__linux__)
    #define BREAKPOINT raise(SIGTRAP)
  #else
    #define BREAKPOINT raise(SIGTRAP)
  #endif
#else
  #define BREAKPOINT ((void)0)  // No-op when not debugging
#endif

// debug macros
#ifdef DEBUG_PRINT
template <class... Args>
inline void DEBUG(Args&&... args) {
    (std::cout << ... << args) << '\n';
}
#else
// compiled out
template <class... Args>
inline void DEBUG(Args&&...) {}
#endif


#ifdef DEBUG_PRINT
  #define DEBUG_PRINT_FUNC(func, ...) func(__VA_ARGS__) 
#else
  #define DEBUG_PRINT_FUNC(func, ...) ((void)0)
#endif


#define FOR_EACH(container, func) \
  std::for_each(container.begin(), container.end(), func)

// dimension macros
#ifndef INITIAL_BOTTOM
  #define INITIAL_BOTTOM 100
# endif

#ifndef INITIAL_RIGHT
  #define INITIAL_RIGHT 100
# endif

#ifndef INITIAL_PARTITION
  #define INITIAL_PARTITION QPEngine::partition_t::vertical
# endif

class QPEngine {
    public:
    QPEngine() noexcept = default;
    QPEngine(size_t numPartitions) noexcept : numPartitions_(numPartitions) {}

    /* public typdefs */
    using vector_t = std::vector<float>; 
    using matrix_t = Eigen::MatrixXd; 
    using coordinate_t = std::pair<float, float>; 

    /**
     * @brief Main driver function.
     * Reads Netlist from inFile and outputs
     * placement results onto outFile
     * 
     * 
     * @param inFile 
     * @param outFile 
     */
    void place(std::ifstream& inFile, std::ofstream& outFile);


    private:
    /* private types */
    using coordinateList_t = std::vector<std::pair<size_t, coordinate_t>>;
    // note that the size_t is the index of the gate and not the gate number
    using assignedGate_t = std::vector<std::pair<size_t, coordinate_t>>;
    // map from net to gate/port
    using netList_t = std::unordered_map<size_t, std::vector<size_t>>;
    using bVector_t = std::pair<Eigen::VectorXd, Eigen::VectorXd>;
    
    // deprecated
    // using netList_t = std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>;
    
    enum class partition_t: uint8_t {
      horizontal,
      vertical
    };
    enum class side_t: uint8_t {
      first,
      second
    };
    struct dimension {
      size_t top_;
      size_t bottom_;
      size_t left_;
      size_t right_;
      dimension(size_t top, size_t bottom, size_t left, size_t right) :
        top_(top), bottom_(bottom), left_(left), right_(right){}
        /**
         * @brief Generates a pair of the new dimensions for the appropriate partition
         * 
         * @param partition 
         * @return std::pair<dimension, dimension> 
         */
      [[nodiscard]] std::pair<dimension, dimension> generateDimensions(partition_t partition) const noexcept;
    };
    /* helper functions */
    /**
     * @brief Read the netlist into netToGateAndPortListMap
     * Also fills in portToCoordinateMap_
     * 
     * @param f 
     */
    [[nodiscard]] std::pair<QPEngine::netList_t, QPEngine::coordinateList_t> _readNetlist(std::ifstream& inFile);

    /**
     * @brief Construct a new check Bounds object
     * 
     * @param val 
     * @param bound 
     * @param msg 
     */
    inline void _checkBounds(const size_t val, const size_t bound, const std::string& msg) const;

    /**
     * @brief Given a netList_t map, creates a thing cMatrix
     * Note: input files in this netlist don't have weights
     * specifed, hence, they are all 0
     * 
     * @param netToGateAndPortListMap 
     * @return Eigen::MatrixXd 
     */
    [[nodiscard]] const matrix_t _createCMatrix(const coordinateList_t& gateCoordianteList) const noexcept;


    /**
     * @brief Creates an A matrix given a C matrix and a netToGateAndPortListMap
     * 
     * 
     * @param c 
     * @param netToGateAndPortListMap 
     * @return matrix_t 
     */
    [[nodiscard]] const matrix_t _createAMatrix(const matrix_t& c, const netList_t& portNetList, const coordinateList_t& gateCoordianteList) const;

    
    /**
     * @brief Creates a bVector and returns it as a list of coordinates
     * 
     * @param netToGateAndPortListMap 
     * @param portToCoordinateMap 
     * @return bVector_t 
     */
    [[nodiscard]] const bVector_t _createBVector(const netList_t& portNetList, const coordinateList_t& gateCoordianteList, const coordinateList_t& portCoordinateList) const noexcept;


    /**
     * @brief Given netList_t map, returns number of gates (cells)
     * in the netlist
     * 
     * @param netToGateAndPortListMap 
     * @return size_t 
     */
    [[nodiscard]] size_t inline _getNumCoordiantes(const coordinateList_t& coordinateList) const noexcept;
    /**
     * @brief Given a portToCoordinateMap_ (coordinateList_t), return the number of ports
     * 
     * @param portToCoordinateMap_ 
     * @return size_t 
     */
    [[nodiscard]] size_t inline _getNumCoordinates(const coordinateList_t& coordinateList) const noexcept;


    /**
     * @brief Print helper functions
     * 
     */
    inline void _printCoordinateList(const coordinateList_t& coordinateList, std::ostream& os = std::cout) const noexcept;
    inline void _printMatrix(const matrix_t& m) const noexcept;
    inline void _printNetList(const netList_t& netList) const noexcept;
    inline void _printBVector(const bVector_t& bVector, const coordinateList_t& gateCoordinateList) const noexcept;
    inline void _printAssignedGates(const assignedGate_t& assignedGates, std::ostream& os = std::cout) const noexcept;


    /**
     * @brief Given a coordinateList, converts it into a pair of vectors
     * 
     * @param coordinateList 
     * @return std::pair<Eigen::VectorXd, Eigen::VectorXd> 
     */
    [[nodiscard]] bVector_t _coordinateToVectorConversion(const coordinateList_t& coordinateList) const noexcept;

    /**
     * @brief Given a pair of vectors (bVector_t), converts it into a coordinateList
     * 
     * @param bVector 
     * @return QPEngine::coordinateList_t 
     */
    [[nodiscard]] coordinateList_t _vectorToCoordinateConversion(const bVector_t& bVector, const coordinateList_t& gateCoordinateList) const noexcept;


    /**
     * @brief Gicen a matrix_t m and bVector_t b (with b_x and b_y), solves for the
     * coordinates of all of the gates in the matrix
     * Uses indecies from the given coordinate list
     * 
     * @param m 
     * @param bVector 
     * @return coordinateList_t 
     */
    [[nodiscard]] coordinateList_t _solveLinearSystem(const matrix_t& m, const bVector_t& bVector, const coordinateList_t& gateCoordinateList) const;


    /**
     * @brief Does this:
     * Creates cMatrix
     * Creates aMatrix
     * Creates bVector
     * Solves the lenear system
     * Returns the placedGateCoordinateList
     * 
     * @param gateCoordinateList 
     * @param portCoordinateList 
     * @param portNetList 
     * @return coordinateList_t 
     */
    [[nodiscard]] coordinateList_t _generatePlacements(const coordinateList_t& gateCoordinateList,
    const coordinateList_t& portCoordinateList, const netList_t& portNetList) const;


    /**
     * @brief Sorts the coordinates for the assignment step
     * 
     * @param p 
     * @param coordinates 
     */
    assignedGate_t _assignBlocks(const partition_t p, coordinateList_t& coordinates) const noexcept;

    /**
     * @brief Propagates the pads based on the side and the edge
     * In the case of an odd number of coordinates, one more of them are assigned
     * to the second side than the first
     * 
     * 
     * @param side 
     * @param coordinates 
     * @param edge 
     */
    // coordinateList_t _propagatePads(const side_t side, const coordinateList_t& coordinates, const coordinate_t& edge) const noexcept;



    [[nodiscard]] const coordinateList_t _initializeGateCoordinateList() const noexcept;
  



    /* member variables */

    // matrix representing the connections
    // matrix_t cMatrix = matrix_t();
    // port to coordinate map for the external ports
    coordinateList_t portToCoordinateMap_ = coordinateList_t();
    // number of recursions which the placer does
    int numPartitions_ = 0;
    int numGates_ = 0;
    netList_t gateNetList_ = netList_t();

}; // QPEngine


// TODO: switch main to a different file
int main(int argc, char** argv) {
  // cancel synch with cstdio
  std::ios_base::sync_with_stdio(false);
  #ifdef DEBUG_PRINT
  spdlog::set_level(spdlog::level::debug);
  #else
  spdlog::set_level(spdlog::level::info);
  #endif
  spdlog::debug("Getting here\n");
  assert(argc == 3);

  try {
    std::ifstream inFile(argv[1]);
    std::ofstream outFile(argv[2]);
    // instantiate with recursion count of 3

    if (!inFile.is_open()) {
    std::cerr << "ERROR: Could not open file: " << argv[1] << std::endl;
    return 1; // TODO: look at the error codes
    }
    QPEngine placer = QPEngine(3);
    // placer.run(inFile, outFile);
    placer.place(inFile, outFile);
    inFile.close();
    outFile.close();
  } 
  catch (const std::ios_base::failure& e) {
    std::cerr << "File error: " << e.what() << std::endl;
  }


} // main()


// TODO: move these functions onto a different file

void QPEngine::_checkBounds(const size_t val, const size_t bound, const std::string& msg) const {
  if (val >= bound) {
    BREAKPOINT;
    throw std::out_of_range(msg);
  }
} // QPEngine::_checkBounds()


std::pair<QPEngine::netList_t, QPEngine::coordinateList_t> QPEngine::_readNetlist(std::ifstream& inFile) {
  // read num gates and nets
  size_t numNets;
  inFile >> numGates_ >> numNets;

  netList_t portNetList;
  coordinateList_t portCoordinateList;

  // note: using the word port for pads
  std::string line;
  // start off reading gate-gate connections not gate-port
  size_t readPorts = 0, numPorts;
  while (getline(inFile, line)) {
    // need to skip the \r
    if (line.find_first_not_of(" \r\t") == std::string::npos) {
        continue;
    }
    std::istringstream ss(line);

    if (!readPorts) {
      size_t gate, numTmpNets, net;
      // reading a gate to gate connection
      ss >> gate; --gate;
      _checkBounds(gate, numGates_, "Input gate greater than number of gates");
      // read port line
      // if (ss.eof()) {
      if (!(ss >> numTmpNets)) {
        readPorts = true;
        // undo the -1 for indexing
        numPorts = gate + 1;
        portCoordinateList.reserve(numPorts);
        continue;
      }
      // insert net gate mapping into netToGateListMap
      while (ss >> net) {
        --net;
        _checkBounds(net, numNets, "Input net greater than number of nets");
        // the weight for these is assumed to be 1
        gateNetList_[net][gate] = 1;
      }
    } else {
      size_t port, net, x, y;
      // reading a port to gate connection
      ss >> port >> net >> x >> y;
      _checkBounds((--port), numPorts, "Input port greater than number of ports");
      _checkBounds((--net), numNets, "Input net greater than number of nets");
      portNetList[net][port] = 1;
      portCoordinateList.push_back(std::make_pair(port, std::make_pair(x, y)));
    }
  }
  return {portNetList, portCoordinateList};
} // QPEngine::readNetlist()

const QPEngine::matrix_t QPEngine::_createCMatrix(const QPEngine::coordinateList_t& gateCoordinateList) const noexcept {
  size_t numGates = _getNumCoordiantes(gateCoordinateList);
  // create the cMatrix by determing where the connections exist
  matrix_t c = matrix_t::Zero(numGates, numGates);
  for (const auto &[net, gates]: gateNetList_) {
      for (size_t i = 0; i < numGates; ++i) {
        for (size_t j = i+1; j < numGates; ++j) {
          // ensure that it's in the gateCoordinateList and they are both
          // connected to the netlist
          size_t gateOneIndex = gateCoordinateList[i].first;
          size_t gateTwoIndex = gateCoordinateList[j].first;
          if(static_cast<int>(gates[gateOneIndex] && gates[gateTwoIndex])) {
            c(i,j) = 1;
            c(j,i) = 1;
          }
        }
      }
    }
    return c;
  } // QPEngine::_createCMatrix()
  
  
  const QPEngine::matrix_t 
  QPEngine::_createAMatrix(const matrix_t& c, const netList_t& portNetList, const coordinateList_t& gateCoordinateList) const {
    size_t numGates = _getNumCoordiantes(gateCoordinateList);
    if (numGates > c.rows() || numGates > c.cols()) throw std::runtime_error("Matrix c size too small");
    matrix_t a = matrix_t::Zero(numGates, numGates);
    for (size_t i = 0; i < numGates; ++i) {
      for (size_t j = 0; j < numGates; ++j) {
        if (i == j) {
          // in the diagonal, 
          // sum up the pad (port) wires also determine the current gate index
          size_t portWireSum = 0, cRowSum = 0, gateIndex = gateCoordinateList[i].first;
          for (const auto &[net, gates]: gateNetList_) {
            if (gates[gateIndex] && static_cast<bool>(portNetList.count(net))) {
              // net is connect to the gate
              const auto& ports = portNetList.at(net);
              portWireSum += std::accumulate(ports.begin(), ports.end(), 0);
            }
          }
          // sum up this row in the c matrix
          cRowSum = c.row(i).sum();
          // set the sum
          a(i,j) = portWireSum + cRowSum;
        } else {
          // not on the diagonal, m[i][j] = -c[i][j]
          a(i,j) = -c(i,j);
        }
      }
    }
    return a;
  } // QPEngine::_createAMatrix()



  // see where to throw the exceptions here
  const QPEngine::bVector_t
  QPEngine::_createBVector(const netList_t& portNetList, const coordinateList_t& gateCoordianteList, const coordinateList_t& portCoordinateList) const noexcept {
    size_t numGates = _getNumCoordinates(gateCoordianteList);
    size_t numPorts = _getNumCoordinates(portCoordinateList);
    Eigen::VectorXd b_x = Eigen::VectorXd::Zero(numGates);
    Eigen::VectorXd b_y = Eigen::VectorXd::Zero(numGates);
    for (int gate = 0; gate < numGates_; ++gate) {
      size_t gateIndex = gateCoordianteList[gate].first;
      for (const auto &[net, gates]: gateNetList_) {
        if (gates[gateIndex]) {
          // net is connect to the gate
          for (int port = 0; port < numPorts; ++port) {
            size_t portIndex = portCoordinateList[port].first;
            if (portNetList.count(net) && portNetList.at(net)[portIndex]) {
              // net is connect to port
              // append the coordinate * the weight (netPorts[port]) of the wire (1) to this
              b_x(gate) += portCoordinateList[port].second.first;
              b_y(gate) += portCoordinateList[port].second.second;
            }
          }
        }
      }
    }
    return std::pair{b_x, b_y};
  } // QPEngine::coordinateList_t()



  [[nodiscard]] size_t inline QPEngine::_getNumCoordiantes(const coordinateList_t& coordinateList) const noexcept {
    return coordinateList.size();
  } // QPEngine::_getNumCoordiantes()


  [[nodiscard]] size_t inline QPEngine::_getNumCoordinates(const coordinateList_t& coordinateList) const noexcept {
    return coordinateList.size();
  } // QPEngine::_getNumCoordinates()


  inline void QPEngine::_printCoordinateList(const coordinateList_t& coordinateList, std::ostream& os) const noexcept{
    // only print this if debug mode
    DEBUG_PRINT_FUNC([](const std::string& s) {fmt::print("{}", s);}, "Printing assigned gates list\n");
    for (const auto& [index, pos]: coordinateList) {
      const auto& [x, y] = pos;
      os << fmt::format("{:d} {:.9f} {:.9f}\n", index, x, y);
    }
  } // QPEngine::_printCoordinateList()


  inline void QPEngine::_printAssignedGates(const assignedGate_t& assignedGates, std::ostream& os) const noexcept{
    // only print this if debug mode
    DEBUG_PRINT_FUNC([](const std::string& s) {fmt::print("{}", s);}, "Printing assigned gates list\n");
    for (const auto& [i, pos] : assignedGates) {
      const auto& [x, y] = pos;
      os << fmt::format("{:d} {:.9f} {:.9f}\n", i+1, x, y);
    }
  } // QPEngine::_printAssignedGates()

  inline void QPEngine::_printMatrix(const matrix_t& m) const noexcept {
    /*
    deprecated
      Took in matrix_t for pretty printing
    */
    // FOR_EACH(matrix, (
    //   [](const auto& row) {
    //     std::cout << "[";
    //     FOR_EACH(row, (
    //       [](const auto val) {
    //         std::cout << val << ",";
    //       }
    //     ));
    //     std::cout << "]\n";
    //   }
    // ));
    fmt::print("Printing Matrix\n");
    fmt::print("{}\n", m);
  } // QPEngine::_printMatrix()

  inline void QPEngine::_printNetList(const netList_t& netList) const noexcept {
    for (auto &[net, connections]: netList) {
      size_t i = 0;
      fmt::print("Net: {:d}\n\tGates:", (net+1));
      FOR_EACH(connections, ([&i](const auto gate){++i; if(static_cast<bool>(gate)){fmt::print("{:d},", i);}}));
      fmt::print("\n");
    }
    fmt::print("\n");
  } // QPEngine::_printNetList()


//   [[nodiscard]] QPEngine::bVector_t QPEngine::_coordinateToVectorConversion(const coordinateList_t& coordinateList) const noexcept {
//   // NRVO constructs everything in place
//   // Note that the vector must always be the size of the whole matrix
//   // in order to keep track of the indecies. This leads to some pretty
//   // sparse vectors
//   size_t vectorSize = numGates_;
//   Eigen::VectorXd b_x(vectorSize);
//   Eigen::VectorXd b_y(vectorSize);
//   for (int i = 0; i < vectorSize; ++i) {
//     const auto &[x, y] = coordinateList[i];
//     b_x(i) = x;
//     b_y(i) = y;
//   }
//   return std::pair{b_x, b_y};
// } // QPEngine::_coordinateToVectorConversion()

// Note that the original goordinate is the one previous to the placement in this case
// It's needed here for the indexing
QPEngine::coordinateList_t QPEngine::_vectorToCoordinateConversion(const bVector_t& bVector, const coordinateList_t& gateCoordinateList) const noexcept {
  // NRVO constructs everything in place
  const auto& [b_x, b_y] = bVector;
  assert(b_x.size() == b_y.size()); // TODO: check
  size_t vectorSize = b_x.size();
  coordinateList_t coordinateList(vectorSize);
  for (int i = 0; i < vectorSize; ++i) {
    coordinateList[i] = std::pair(gateCoordinateList[i].first, std::pair(b_x(i), b_y(i)));
  }
  return coordinateList;
} // QPEngine::_vectorToCoordinateConversion()


  inline void QPEngine::_printBVector(const bVector_t& bVector, const coordinateList_t& gateCoordinateList) const noexcept{
    const coordinateList_t coordinateList = _vectorToCoordinateConversion(bVector, gateCoordinateList);
    _printCoordinateList(coordinateList);
  } // QPEngine::_printBVector()


  /**
   * @brief Deprecated
   * 
   */
  [[nodiscard]] QPEngine::coordinateList_t QPEngine::_solveLinearSystem(const matrix_t& m, const bVector_t& bVector, const coordinateList_t& gateCoordinateList) const {
    // create cMatrix(gateCoordinateList); 
    //     create aMatrix(gateCoordinateList, portNetList);
    //     create bVector(gateCoordinateList, portCoordinateList, portNetList);
    const auto& [b_x, b_y] = bVector;
    // bounds check
    if ((b_x.size() != b_y.size()) || (b_x.size() != m.rows()) || (m.rows() != m.cols())) {
      BREAKPOINT;
      throw std::runtime_error("Invalid matrix or bvector dimensions for QR decomposition");
    }
    Eigen::VectorXd placement_x = m.colPivHouseholderQr().solve(b_x);
    Eigen::VectorXd placement_y = m.colPivHouseholderQr().solve(b_y);

    return _vectorToCoordinateConversion(std::pair(placement_x, placement_y), gateCoordinateList);
  } // QPEngine::placements()



  QPEngine::assignedGate_t QPEngine::_assignBlocks(const partition_t p, coordinateList_t& coordinates) const noexcept {
    assignedGate_t assignedGates(coordinates.size()); size_t i = 0;
    // std::transform(coordinates.begin(), coordinates.end(), assignedGates.begin(),
    // [&i](const auto& coordinate) { return std::pair(i++, coordinate); });
    // std::sort(assignedGates.begin(), assignedGates.end(), 
    //   [p](const auto& a, const auto& b) {
    //     const auto& [index_a, pos_a] = a;
    //     const auto& [a_x, a_y] = pos_a;
    //     const auto& [index_b, pos_b] = b;
    //     const auto& [b_x, b_y] = pos_b;
    //     return (p == partition_t::horizontal) ?
    //       (a_y == b_y ? a_x < b_x : a_y < b_y) :
    //       (a_x == b_x ? a_y < b_y : a_x < b_x);
    //   }
    // );
    return assignedGates;
  } // QPEngine::_assignBlocks()



  /*
    Main recursion loop
    _PartitionAndPlace(x, y, pads, NetList, split, recur_depth) {
      if (!recur_depth) return;
      generate placements()
      assign placements(split)
      propagate pads(left or top)
      generate NetList(left or top)
      generate placements(left or top)
      propagate pads(right or bottom)
      generate NetList(right or bottom)
      generate placements(right or bottom)
      _PartitionAndPlace(new_x, new_y, new_pads, new_NetList, not_split)
      _PartitionAndPlace(new_x, new_y, new_pads, new_NetList, not_split)
    }
  


    Solution to the aforementioned issue of sparcity:

    Things which stay constant throughout:
      - The gate mappings.
        Even though there are some arbitrary ones, it's fine since
        those won't be indexed while making the cMatrix or the aMatrix
      - numPartitions_ can be decreased by one recursion and set during run


    Things which depend on the current iteration:
        - The portNetLists.
      - input gateCoordinateList (initially start off arbitrary)
      - input portCoordinateList
      - input dimension (need a structure for this (left, right, top, bottom))
      - input previous partition


    API changes:
      - Size for bVector, cMatrix, and aMatrix depends on coordinateList size
      - _createCMatrix needs to include the coordinateList as an argument
      - When the resulting solution is sorted, it needs to be indexed by the
        input gate coordinateList
      - _generatePlacements needs to take in coordinateList
      - _assignBlocks returns a coordinateList now sorted


      At the end, return a 

      psuedocode:

      globally:
        numPartition_;
        gateNetList_;

      coordinateList _generatePlacements(gateCoordinateList, portCoordinateList, portNetList) {
        create cMatrix(gateCoordinateList); 
        create aMatrix(gateCoordinateList, portNetList);
        create bVector(gateCoordinateList, portCoordinateList, portNetList);
      }
      
      _place(gateCoordinateList, portCoordinateList, portNetList, dimension, partition) {
        if (!numPartition_--) return portCoordinateList;

        // determine new dimensions
        firstDimension = _firstDimension(dimension, partition);
        secondDimension = _secondDimension(dimension, partition);

        // assign gateCoordinateList
        assignedGateCoordinateList = _assignBlocks(gateCoordinateList);

        // generate sub-gateCoordinateList
        firstGateCoordinateList = _firstHalf(assignedGateCoordinateList);
        secondGateCoordinateList = _secondHalf(assignedGateCoordinateList);

        // generate netlists (don't need any placements before this)
        firstPortNetList = _generateNetlist(firstGateCoordinateList, secondGateCoordinateList, portNetList);
        secondPortNetList = _generateNetlist(firstGateCoordinateList, firstGateCoordinateList, portNetList); // doesn't need to be placed here

        // place first
        firstPortCoordinateList = _propagatePads(secondGateCoordinateList, portCoordinateList, firstDimension);
        placedFirstGateCoordinateList = _generatePlacements(firstGateCoordinateList, firstPortCoordinateList, firstPortNetList);

        // place second
        secondPortCoordinateList = _propagatePads(placedFirstGateCoordinateList, portCoordinateList, secondDimension);
        placedSecondGateCoordinateList = _generatePlacements(secondGateCoordinateList, secondPortCoordinateList, secondPortNetList);

        merge the placed gateCooordinateLists and return
        placedFirstGateCoordinateList = _place(placedFirstGateCoordinateList, firstPortCoordinateList, firstPortNetList, firstDimension, notPartition);
        placedSecondGateCoordinateList = _place(placedSecondGateCoordinateList, secondPortCoordinateList, secondPortNetList, secondDimension, notPartition);
        return merge(placedFirstGateCoordinateList, placedSecondGateCoordinateList); 
      }

      place(inFile outFile) {
        read inFile(generate gateNetList_, portCoordinateList, and portNetList);
        zero initialize(gateCoordinateList);
        placedGateCoordinateList = _generatePlacements(gateCoordinateList, portCoordinateList, portNetList);
        init dimension;
        init partition;

        return _place(placedGateCoordinateList, portCoordinateList, portNetList, dimension, partition);
      }



      
  */


  void QPEngine::place(std::ifstream& inFile, std::ofstream& outFile) {
    /* read input file */
    BREAKPOINT;
    spdlog::debug("Reading Nelist");
    const auto [portNetList, portCoordinateList] = _readNetlist(inFile);
    spdlog::debug("gateNelist_:");
    DEBUG_PRINT_FUNC(_printNetList, gateNetList_);
    spdlog::debug("portNelist:");
    DEBUG_PRINT_FUNC(_printNetList, portNetList);
    spdlog::debug("portCoordinateList:");
    DEBUG_PRINT_FUNC(_printCoordinateList, portCoordinateList);

    /* zero init the zeroGateCoordinateList */
    BREAKPOINT;
    spdlog::debug("zeroGateCoordinateList:");
    const coordinateList_t zeroGateCoordinateList = _initializeGateCoordinateList();
    DEBUG_PRINT_FUNC(_printCoordinateList, zeroGateCoordinateList);

    /* generate placements */
    BREAKPOINT;
    coordinateList_t placedGateCoordinateList = _generatePlacements(zeroGateCoordinateList, portCoordinateList, portNetList);
    spdlog::debug("placedGateCoordinateList:");
    DEBUG_PRINT_FUNC(_printCoordinateList, placedGateCoordinateList);
    
    /* init dimension */
    dimension d(0, INITIAL_BOTTOM, 0, INITIAL_RIGHT);
    /* init partition */
    partition_t initialPartition = INITIAL_PARTITION;
    /* recursively partition */
    // return _place(placedGateCoordinateList, portCoordinateList, portNetList, dimension, partition);
  }


  /* Dimension code */
  std::pair<QPEngine::dimension, QPEngine::dimension> 
  QPEngine::dimension::generateDimensions(partition_t partition) const noexcept {
    return (partition == partition_t::vertical)
      ? std::pair(dimension(top_, bottom_, left_, right_/2), dimension(top_, bottom_, right_/2, right_))
      : std::pair(dimension(top_, bottom_/2, left_, right_), dimension(bottom_/2, bottom_, left_, right_));
  } // QPEngine::dimension::generateDimensions()


  const QPEngine::coordinateList_t QPEngine::_initializeGateCoordinateList() const noexcept {
    coordinateList_t zeroGateCoordinateList(numGates_);
    for (size_t i = 1; i <= numGates_; ++i) {
      zeroGateCoordinateList[i-1].first = i;
    }
    return zeroGateCoordinateList;
  } // QPEngine::_initializeGateCoordinateList()


  QPEngine::coordinateList_t 
  QPEngine::_generatePlacements(const coordinateList_t& gateCoordinateList,
    const coordinateList_t& portCoordinateList, const netList_t& portNetList) const {
    
    /* generate cMatrix */
    BREAKPOINT;
    matrix_t c = _createCMatrix(gateCoordinateList);
    DEBUG_PRINT_FUNC(_printMatrix, c);

    /* generate aMatrix */
    BREAKPOINT;
    spdlog::debug("Creating aMatrix");
    matrix_t a = _createAMatrix(c, portNetList, gateCoordinateList);
    DEBUG_PRINT_FUNC(_printMatrix, a);

    // /* generate bVector */
    BREAKPOINT;
    spdlog::debug("Creating BVector");
    bVector_t b = _createBVector(portNetList, gateCoordinateList, portCoordinateList);
      
    // /* generate placements */
    BREAKPOINT;
    spdlog::debug("Generating Placements");
    // return _solveLinearSystem(a, b, gateCoordinateList);
    coordinateList_t placements = _solveLinearSystem(a, b, gateCoordinateList);
    DEBUG_PRINT_FUNC(_printCoordinateList, placements);
    return placements;
  } // QPEngine::placements()

