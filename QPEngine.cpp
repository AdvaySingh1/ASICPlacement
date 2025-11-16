#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <numeric>
#include <sstream>
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
  #define DEBUG_PRINT_FUNC(container, func) func(container) 
#else
  #define DEBUG_PRINT_FUNC(container, func) ((void)0)
#endif


#define FOR_EACH(container, func) \
  std::for_each(container.begin(), container.end(), func)



class QPEngine {
    public:
    QPEngine() noexcept = default;
    QPEngine(size_t recCount) noexcept : recCount_(recCount) {}

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
    void run(std::ifstream& inFile, std::ofstream& outFile);


    private:
    /* private types */
    using coordinateList_t = std::vector<coordinate_t>;
    using netList_t = std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>;
    using bVector_t = std::pair<Eigen::VectorXd, Eigen::VectorXd>;
    /* helper functions */
    /**
     * @brief Read the netlist into netToGateAndPortListMap
     * Also fills in portToCoordinateMap_
     * 
     * @param f 
     */
    netList_t _readNetlist(std::ifstream& inFile);

    /**
     * @brief Construct a new check Bounds object
     * 
     * @param val 
     * @param bound 
     * @param msg 
     */
    void inline _checkBounds(const size_t val, const size_t bound, const std::string& msg) const;

    /**
     * @brief Given a netList_t map, creates a thing cMatrix
     * Note: input files in this netlist don't have weights
     * specifed, hence, they are all 0
     * 
     * @param netToGateAndPortListMap 
     * @return Eigen::MatrixXd 
     */
    [[nodiscard]] matrix_t _createCMatrix(const netList_t& netToGateAndPortListMap) const noexcept;


    /**
     * @brief Creates an A matrix given a C matrix and a netToGateAndPortListMap
     * 
     * 
     * @param c 
     * @param netToGateAndPortListMap 
     * @return matrix_t 
     */
    [[nodiscard]] matrix_t _createAMatrix(const matrix_t& c, const netList_t& netToGateAndPortListMap) const;

    
    /**
     * @brief Creates a bVector and returns it as a list of coordinates
     * 
     * @param netToGateAndPortListMap 
     * @param portToCoordinateMap 
     * @return bVector_t 
     */
    [[nodiscard]] bVector_t _createBVector(const netList_t& netToGateAndPortListMap, const coordinateList_t& portToCoordinateMap) const noexcept;


    /**
     * @brief Given netList_t map, returns number of gates (cells)
     * in the netlist
     * 
     * @param netToGateAndPortListMap 
     * @return size_t 
     */
    [[nodiscard]] size_t inline _getNumGates(const netList_t& netToGateAndPortListMap) const noexcept;
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
    void inline _printCoordinateList(const coordinateList_t& portToCoordinateMap) const noexcept;
    void inline _printMatrix(const matrix_t& m) const noexcept;
    void inline _printNetList(const netList_t& netToGateAndPortListMap) const noexcept;
    void inline _printBVector(const bVector_t& bVector) const noexcept;


    /**
     * @brief Given a coordinateList, converts it into a pair of vectors
     * 
     * @param coordinateList 
     * @return std::pair<Eigen::VectorXd, Eigen::VectorXd> 
     */
    [[nodiscard]] bVector_t coordinateToVectorConversion(const coordinateList_t& coordinateList) const noexcept;

    /**
     * @brief Given a pair of vectors (bVector_t), converts it into a coordinateList
     * 
     * @param bVector 
     * @return QPEngine::coordinateList_t 
     */
    [[nodiscard]] coordinateList_t vectorToCoordinateConversion(const bVector_t& bVector) const noexcept;



    /* member variables */

    // matrix representing the connections
    // matrix_t cMatrix = matrix_t();
    // port to coordinate map for the external ports
    coordinateList_t portToCoordinateMap_ = coordinateList_t();
    // number of recursions which the placer does
    int recCount_ = 0;

}; // QPEngine


// TODO: switch main to a different file
int main(int argc, char** argv) {
  std::cout << "Getting here\n";
  if (argc != 3) {
    std::cerr << "Error. Did not specify inFile and\
        outFile paths for the Netlist" << std::endl;
  }

  
  #ifdef DEBUG_PRINT
  spdlog::set_level(spdlog::level::debug);
  #else
  spdlog::set_level(spdlog::level::info);
  #endif


  BREAKPOINT;

  // cancel synch with cstdio
  std::ios_base::sync_with_stdio(false);

  try {
    std::ifstream inFile(argv[1]);
    std::ofstream outFile(argv[2]);
    // instantiate with recursion count of 3

    if (!inFile.is_open()) {
    std::cerr << "ERROR: Could not open file: " << argv[1] << std::endl;
    return 1; // TODO: look at the error codes
    }
    QPEngine placer = QPEngine(3);
    placer.run(inFile, outFile);
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
}


void QPEngine::run(std::ifstream& inFile, std::ofstream& outFile) {
  /* read the input file and generate netlist and gloabl portToCoordinateMap_*/
  BREAKPOINT;
  netList_t netToGateAndPortListMap = _readNetlist(inFile);
  DEBUG_PRINT_FUNC(netToGateAndPortListMap, _printNetList);
  DEBUG_PRINT_FUNC(portToCoordinateMap_, _printCoordinateList);

  /* generate cMatrix */
  BREAKPOINT;
  matrix_t c = _createCMatrix(netToGateAndPortListMap);
  DEBUG_PRINT_FUNC(c, _printMatrix);
  
  /* generate aMatrix */
  BREAKPOINT;
  matrix_t a = _createAMatrix(c, netToGateAndPortListMap);
  DEBUG_PRINT_FUNC(a, _printMatrix);

  /* generate bVector */
  BREAKPOINT;
  bVector_t bVector = _createBVector(netToGateAndPortListMap, portToCoordinateMap_);
  DEBUG_PRINT_FUNC(bVector, _printBVector);

  BREAKPOINT;
}

typename QPEngine::netList_t QPEngine::_readNetlist(std::ifstream& inFile) {
  // read num gates and nets
  size_t numGates, numNets;
  inFile >> numGates >> numNets;

  netList_t netToGateAndPortListMap = netList_t(
      numNets,
      {std::vector<size_t>(numGates, 0),
      std::vector<size_t>()}
    );

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
      _checkBounds(gate, numGates, "Input gate greater than number of gates");
      // read port line
      // if (ss.eof()) {
      if (!(ss >> numTmpNets)) {
        readPorts = true;
        // undo the -1 for indexing
        numPorts = gate + 1;

        // reshape the netToGateAndPortListMap
        FOR_EACH(netToGateAndPortListMap, 
          [numPorts](auto& p) {
            p.second.resize(numPorts, 0);
          });
        // reshape the portToCoordinateMap
        portToCoordinateMap_.resize(numPorts);
        continue;
      }
      // ss >> numTmpNets; // not really needed?
      // insert net gate mapping into netToGateListMap
      while (ss >> net) {
        --net;
        _checkBounds(net, numNets, "Input net greater than number of nets");
        // the weight for these is assumed to be 1
        netToGateAndPortListMap[net].first[gate] = 1;
      }
    } else {
      size_t port, net, x, y;
      // reading a port to gate connection
      ss >> port >> net >> x >> y;
      _checkBounds((--port), numPorts, "Input port greater than number of ports");
      _checkBounds((--net), numNets, "Input net greater than number of nets");
      netToGateAndPortListMap[net].second[port] = 1;
      portToCoordinateMap_[port] = std::make_pair(x, y);
    }
  }
  return netToGateAndPortListMap;
} // QPEngine::readNetlist()

[[nodiscard]] QPEngine::matrix_t QPEngine::_createCMatrix(const QPEngine::netList_t& netToGateAndPortListMap) const noexcept {
  size_t numGates = _getNumGates(netToGateAndPortListMap);
  // create the cMatrix by determing where the connections exist
  matrix_t c = matrix_t::Zero(numGates, numGates);
  for (const auto &[netGates, _]: netToGateAndPortListMap) {
      for (size_t i = 0; i < numGates; ++i) {
        for (size_t j = i+1; j < numGates; ++j) {
          c(i,j) = c(j,i) = static_cast<int>(netGates[i] && netGates[j]);
        }
      }
    }
    return c;
  } // QPEngine::_createCMatrix()
  
  
  [[nodiscard]] QPEngine::matrix_t 
  QPEngine::_createAMatrix(const matrix_t& c, const QPEngine::netList_t& netToGateAndPortListMap) const {
    size_t numGates = _getNumGates(netToGateAndPortListMap);
    if (numGates >= c.rows() || numGates > c.cols()) throw std::runtime_error("Matrix c size too small");
    matrix_t a = matrix_t::Zero(numGates, numGates);
    for (size_t i = 0; i < numGates; ++i) {
      for (size_t j = 0; j < numGates; ++j) {
        if (i == j) {
          // in the diagonal, 
          // sum up the pad (port) wires
          size_t portWireSum = 0, cRowSum = 0;
          for (const auto &[netGates, netPorts]: netToGateAndPortListMap) {
            if (netGates[i]) {
              // net is connect to the gate
              portWireSum += std::accumulate(netPorts.begin(), netPorts.end(), 0);
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
  [[nodiscard]] QPEngine::bVector_t
  QPEngine::_createBVector(const netList_t& netToGateAndPortListMap, const coordinateList_t& portToCoordinateMap) const noexcept {
    size_t numGates = _getNumGates(netToGateAndPortListMap);
    size_t numPorts = _getNumCoordinates(portToCoordinateMap);
    Eigen::VectorXd b_x = Eigen::VectorXd::Zero(numGates);
    Eigen::VectorXd b_y = Eigen::VectorXd::Zero(numGates);
    for (int gate = 0; gate < numGates; ++gate) {
      for (const auto &[netGates, netPorts]: netToGateAndPortListMap) {
        if (netGates[gate]) {
          // net is connect to the gate
          for (int port = 0; port < numPorts; ++port) {
            if (netPorts[port]) {
              // net is connect to port
              // append the coordinate * the weight (netPorts[port]) of the wire (1) to this
              b_x(gate) += portToCoordinateMap[port].first;
              b_y(gate) += portToCoordinateMap[port].second;
            }
          }
        }
      }
    }
    return std::pair{b_x, b_y};
  } // QPEngine::coordinateList_t()



  [[nodiscard]] size_t inline QPEngine::_getNumGates(const netList_t& netToGateAndPortListMap) const noexcept {
    return netToGateAndPortListMap.size() ? netToGateAndPortListMap[0].first.size() : 0;
  } // QPEngine::_getNumGates()


  [[nodiscard]] size_t inline QPEngine::_getNumCoordinates(const coordinateList_t& coordinateList) const noexcept {
    return coordinateList.size();
  } // QPEngine::_getNumCoordinates()


  void inline QPEngine::_printCoordinateList(const coordinateList_t& portToCoordinateMap) const noexcept{
    BREAKPOINT;
    for (const auto&[x, y]: portToCoordinateMap) {
      spdlog::debug("({:.2f},{:.2f})", x, y);
    }
  } // QPEngine::_printCoordinateList()

  void inline QPEngine::_printMatrix(const matrix_t& m) const noexcept {
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
    spdlog::debug("{}", m);
  } // QPEngine::_printMatrix()

  void inline QPEngine::_printNetList(const netList_t& netToGateAndPortListMap) const noexcept {
    for (int net = 0; net < netToGateAndPortListMap.size(); ++net) {
      size_t i = 0;
      fmt::print("Net: {:d}\n\tGates:", (net+1));
      FOR_EACH(netToGateAndPortListMap[net].first, ([&i](const auto gate){++i; if(static_cast<bool>(gate)){fmt::print("{:d},", i);}}));
      fmt::print("\n\tPorts:");
      i = 0;
      FOR_EACH(netToGateAndPortListMap[net].second, ([&i](const auto port){++i; if(static_cast<bool>(port)) { fmt::print("{:d},", i);}}));
      fmt::print("\n");
    }
  } // QPEngine::_printNetList()


  [[nodiscard]] QPEngine::bVector_t QPEngine::coordinateToVectorConversion(const coordinateList_t& coordinateList) const noexcept {
  // NRVO constructs everything in place
  size_t vectorSize = _getNumCoordinates(coordinateList);
  Eigen::VectorXd b_x(vectorSize);
  Eigen::VectorXd b_y(vectorSize);
  for (int i = 0; i < vectorSize; ++i) {
    const auto &[x, y] = coordinateList[i];
    b_x(i) = x;
    b_y(i) = y;
  }
  return std::pair{b_x, b_y};
} // QPEngine::coordinateToVectorConversion()

[[nodiscard]] QPEngine::coordinateList_t QPEngine::vectorToCoordinateConversion(const bVector_t& bVector) const noexcept {
  // NRVO constructs everything in place
  const auto& [b_x, b_y] = bVector;
  assert(b_x.size() == b_y.size());
  size_t vectorSize = b_x.size();
  coordinateList_t coordinateList(vectorSize);
  for (int i = 0; i < vectorSize; ++i) {
    coordinateList[i].first = b_x(i);
    coordinateList[i].second = b_y(i);
  }
  return coordinateList;
} // QPEngine::vectorToCoordinateConversion()


  void inline QPEngine::_printBVector(const bVector_t& bVector) const noexcept{
    const coordinateList_t coordinateList = vectorToCoordinateConversion(bVector);
    _printCoordinateList(coordinateList);
  } // QPEngine::_printBVector()