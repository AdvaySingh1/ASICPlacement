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

    /**
     * @brief Consider making the following types private
     */
    using matrix_t = std::vector<std::vector<size_t>>; 
    using coordinate_t = std::pair<size_t, size_t>; 

    // TODO: make these types private:
    using coordinateList_t = std::vector<coordinate_t>;

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
    /* helper functions */
    /**
     * @brief Read the netlist into the cMatrix from the file
     * 
     * @param f 
     */
    void readNetlist(std::ifstream& inFile);

    /**
     * @brief Construct a new check Bounds object
     * 
     * @param val 
     * @param bound 
     * @param msg 
     */
    void inline checkBounds(const size_t val, const size_t bound, const std::string& msg) const;



    /* member variables */

    // matrix representing the connections
    matrix_t cMatrix = matrix_t();
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

void QPEngine::checkBounds(const size_t val, const size_t bound, const std::string& msg) const {
  if (val >= bound) {
    throw std::out_of_range(msg);
  }
}


void QPEngine::run(std::ifstream& inFile, std::ofstream& outFile) {
  readNetlist(inFile);
}

void QPEngine::readNetlist(std::ifstream& inFile) {
  // read num gates and nets
  size_t numGates, numNets;
  inFile >> numGates >> numNets;

  // init numGates * numGates matrix all to 0
  cMatrix = matrix_t(numGates, std::vector<size_t>(numGates, 0));

  // init an unordered map
  // key: net, value: pair with a vector for gates attached to it
  // and a vector for ports attached to it
  std::vector<std::pair<std::vector<size_t>, 
    std::vector<size_t>>> netToGateAndPortListMap(
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
      checkBounds(gate, numGates, "Input gate greater than number of gates");
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
        checkBounds(net, numNets, "Input net greater than number of nets");
        // the weight for these is assumed to be 1
        netToGateAndPortListMap[net].first[gate] = 1;
      }
    } else {
      size_t port, net, x, y;
      // reading a port to gate connection
      ss >> port >> net >> x >> y;
      checkBounds(port, numPorts, "Input port greater than number of ports");
      checkBounds(net, numNets, "Input net greater than number of nets");
      netToGateAndPortListMap[net].second[port] = 1;
      portToCoordinateMap_[port] = std::make_pair(x, y);
    }
  }
}