#ifndef PARTITIONER_H
#define PARTITIONER_H

#include <fstream>
#include <vector>
#include <unordered_map>
#include <chrono>
#include "cell.h"
#include "net.h"
using namespace std;

class Partitioner
{
public:
    // constructor and destructor
    Partitioner(fstream& inFile) :
        _cutSize(0), _netNum(0), _cellNum(0), _maxPinNum(0), _bFactor(0),
        _accGain(0), _maxAccGain(0), _iterNum(0) {
        parseInput(inFile);
        _partSize[0] = 0;
        _partSize[1] = 0;
        _maxGainPtr[0] = -1;
        _maxGainPtr[1] = -1;
    }
    ~Partitioner() {
        clear();
    }

    // basic access methods
    int getCutSize() const          { return _cutSize; }
    int getNetNum() const           { return _netNum; }
    int getCellNum() const          { return _cellNum; }
    double getBFactor() const       { return _bFactor; }
    int getPartSize(int part) const { return _partSize[part]; }

    // public interface
    void parseInput(fstream& inFile);
    void init(int seed);
    void partition();

    // member functions about reporting
    void printSummary() const;
    void reportNet() const;
    void reportCell() const;
    void writeResult(fstream& outFile);

private:
    int                 _cutSize;       // cut size
    int                 _partSize[2];   // size (cell number) of partition A(0) and B(1)
    int                 _netNum;        // number of nets
    int                 _cellNum;       // number of cells
    int                 _maxPinNum;     // Pmax for building bucket list
    double              _bFactor;       // the balance factor to be met
    Node*               _maxGainCell;   // pointer to max gain cell
    vector<Net*>        _netArray;      // net array of the circuit
    vector<Cell*>       _cellArray;     // cell array of the circuit

    // Array-based bucket list: index = gain + _maxPinNum
    vector<Node*>       _bList[2];      // bucket list of partition A(0) and B(1)
    int                 _maxGainPtr[2]; // current max occupied index per partition

    unordered_map<string, int>    _netName2Id;    // mapping from net name to id
    unordered_map<string, int>    _cellName2Id;   // mapping from cell name to id

    int                 _accGain;       // accumulative gain
    int                 _maxAccGain;    // maximum accumulative gain
    int                 _moveNum;       // number of cell movements
    int                 _iterNum;       // number of iterations
    int                 _bestMoveNum;   // store best number of movements
    int                 _unlockNum[2];  // number of unlocked cells
    vector<int>         _moveStack;     // history of cell movement

    // Clean up partitioner
    void clear();

    // FM helpers
    void _findMaxPinNum();
    void _findCutSize();
    void _initNetPartCounts();
    void _computeAndInsertGain(int cellId);
    void _removeFromBList(int part, Node* node, int gain);
    void _insertIntoBList(int part, Node* node, int gain);
    void _updateNeighborGains(int movedCellId, int from, int to);
    void _pickMaxGainCell();
    void _resetPass();
    bool _canMoveTo(int from, int lowerBound, int upperBound) const;
    void _runFM();
    void _prepareForFM();
    void _initFromPerturb(const vector<bool>& bestPart, int bestPS0, int bestPS1,
                          int seed, double ratio);
};

#endif  // PARTITIONER_H
