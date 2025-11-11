#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

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
  #define PRINT_MATRIX(matrix) printMatrix(matrix) 
#else
  #define PRINT_MATRIX(matrix) ((void)0)
#endif


#define FOR_EACH(container, func) \
  std::for_each(container.begin(), container.end(), func)


class QPEngine {
    public:
    QPEngine() noexcept = default;
    QPEngine(size_t recCount) noexcept : recCount_(recCount) {}

    /* public typdefs */
    using matrix_t = std::vector<std::vector<size_t>>; 
    using coordinate_t = std::pair<size_t, size_t>; 

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

    /**
     * @brief Prints the matrix
     * 
     * @param m 
     */
    static void printMatrix(const matrix_t& m);

    private:
    /* private types */
    using coordinateList_t = std::vector<coordinate_t>;
    using netList_t = std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>;
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
     * @return matrix_t 
     */
    [[nodiscard]] matrix_t _createCMatrix(const netList_t& netToGateAndPortListMap) const noexcept;


    [[nodiscard]] matrix_t _createAMatrix(const matrix_t& c, const netList_t& netToGateAndPortListMap) const;


    /**
     * @brief Given netList_t map, returns number of gates (cells)
     * in the netlist
     * 
     * @param netToGateAndPortListMap 
     * @return size_t 
     */
    [[nodiscard]] size_t inline _getNumGates(const netList_t& netToGateAndPortListMap) const noexcept;




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
  if (argc != 3) {
    std::cerr << "Error. Did not specify inFile and\
        outFile paths for the Netlist" << std::endl;
  }

  // cancel synch with cstdio
  std::ios_base::sync_with_stdio(false);

  try {
    std::ifstream inFile(argv[1]);
    std::ofstream outFile(argv[1]);
    // instantiate with recursion count of 3
    QPEngine placer = QPEngine(3);
    placer.run(inFile, outFile);
  } 
  catch (const std::ios_base::failure& e) {
    std::cerr << "File error: " << e.what() << std::endl;
  }


} // main()


// TODO: move these functions onto a different file

void QPEngine::_checkBounds(const size_t val, const size_t bound, const std::string& msg) const {
  if (val >= bound) {
    throw std::out_of_range(msg);
  }
}


void QPEngine::run(std::ifstream& inFile, std::ofstream& outFile) {
  _readNetlist(inFile);
}

typename QPEngine::netList_t QPEngine::_readNetlist(std::ifstream& inFile) {
  // read num gates and nets
  size_t numGates, numNets;
  inFile >> numGates >> numNets;

  // init numGates * numGates matrix all to 0
  // cMatrix = matrix_t(numGates, std::vector<size_t>(numGates, 0));

  // init an unordered map
  // key: net, value: pair with a vector for gates attached to it
  // and a vector for ports attached to it
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
    std::istringstream ss(line);

    if (!readPorts) {
      size_t gate, numTmpNets, net;
      // reading a gate to gate connection
      ss >> gate;
      _checkBounds(gate, numGates, "Input gate greater than number of gates");
      // read port line
      if (ss.eof()) {
        readPorts = true;
        numPorts - gate;

        // reshape the netToGateAndPortListMap
        FOR_EACH(netToGateAndPortListMap, 
          [](auto& p) {
            p.second.resize(numPorts, 0);
          });
        // reshape the portToCoordinateMap
        portToCoordinateMap_.reserve(numPorts);
        continue;
      }
      ss >> numTmpNets; // not really needed?
      // insert net gate mapping into netToGateListMap
      while (ss >> net) {
        _checkBounds(net, numNets, "Input net greater than number of nets");
        // the weight for these is assumed to be 1
        netToGateAndPortListMap[net].first[gate] = 1;
      }
    } else {
      size_t port, net, x, y;
      // reading a port to gate connection
      ss >> port >> net >> x >> y;
      _checkBounds(port, numPorts, "Input port greater than number of ports");
      _checkBounds(net, numNets, "Input net greater than number of nets");
      netToGateAndPortListMap[net].second[port] = 1;
      portToCoordinateMap_[port] = std::make_pair(x, y);
    }
  }
} // QPEngine::readNetlist()

[[nodiscard]] typename QPEngine::matrix_t QPEngine::_createCMatrix(const QPEngine::netList_t& netlist) const noexcept {
  size_t numGates = _getNumGates(netlist);
  // create the cmatrix by determing where the connections exist
  matrix_t c(numGates, std::vector<size_t>(numGates, 0));
  FOR_EACH(netlist,
    [](auto& p) {
      const auto& gates = p.first; 
      for (size_t i = 0; i < gates.size(); ++i) {
        for (size_t j = i+1; j < gates.size(); ++j) {
          c[gates[i]][gates[j]] = c[gates[j]][gates[i]] = 1 / netlist.size();
        }
      }
    });
    return c;
  } // QPEngine::_createCMatrix()
  
  
  [[nodiscard]] typename QPEngine::matrix_t 
  QPEngine::_createAMatrix(const QPEngine::matrix_t& c, const QPEngine::netList_t& netlist) const {
    size_t numGates = _getNumGates(netlist);
    if (numGates > c.size() || numGates > c[0].size()) throw std::runtime_error("Matrix c size too small");
    matrix_t m(numGates, std::vector<size_t>(numGates, 0));
    for (size_t i = 0; i < numGates; ++i) {
      for (size_t j = 0; j < numGates; ++j) {
        if (i == j) {
          // in the diagonal, sum up the pad (port) wires
          int portWireSum;
          FOR_EACH(netlist, 
            [&portWireSum](const auto& p){
              if (p.first[i]) {
                FOR_EACH(p.second, 
                  [&portWireSum](const auto port) {
                    if (port) ++portWireSum;
                  }
                );
              }
            }
          );
        } else {
          // not on the diagonal, m[i][j] = -c[i][j]
          m[i][j] = -c[i][j];
        }
      }
    }
  } // QPEngine::_createAMatrix()