#include <iostream>
#include <algorithm>
#include <random>
#include <numeric>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <map>
#include "cell.h"
#include "net.h"
#include "partitioner.h"
using namespace std;


void Partitioner::parseInput(fstream& inFile)
{
    string str;
    // Set balance factor
    inFile >> str;
    _bFactor = stod(str);

    // Set up whole circuit
    while (inFile >> str) {
        if (str == "NET") {
            string netName, cellName, tmpCellName = "";
            inFile >> netName;
            int netId = _netNum;
            _netArray.push_back(new Net(netName));
            _netName2Id[netName] = netId;
            while (inFile >> cellName) {
                if (cellName == ";") {
                    tmpCellName = "";
                    break;
                }
                else {
                    // a newly seen cell
                    if (_cellName2Id.count(cellName) == 0) {
                        int cellId = _cellNum;
                        _cellArray.push_back(new Cell(cellName, 0, cellId));
                        _cellName2Id[cellName] = cellId;
                        _cellArray[cellId]->addNet(netId);
                        _cellArray[cellId]->incPinNum();
                        _netArray[netId]->addCell(cellId);
                        ++_cellNum;
                        tmpCellName = cellName;
                    }
                    // an existed cell
                    else {
                        if (cellName != tmpCellName) {
                            assert(_cellName2Id.count(cellName) == 1);
                            int cellId = _cellName2Id[cellName];
                            _cellArray[cellId]->addNet(netId);
                            _cellArray[cellId]->incPinNum();
                            _netArray[netId]->addCell(cellId);
                            tmpCellName = cellName;
                        }
                    }
                }
            }
            ++_netNum;
        }
    }
    return;
}

// Helper function by me

// ─── Small utilities ──────────────────────────────────────────────────────────

void Partitioner::_findMaxPinNum() {
    for (int i = 0; i < _cellNum; ++i) {
        if (_cellArray[i]->getPinNum() > _maxPinNum)
            _maxPinNum = _cellArray[i]->getPinNum();
    }
}

void Partitioner::_findCutSize() {
    _cutSize = 0;
    for (int i = 0; i < _netNum; ++i) {
        if (_netArray[i]->getPartCount(0) > 0 && _netArray[i]->getPartCount(1) > 0)
            ++_cutSize;
    }
}

// Initialise the PartCount fields of every net from the current cell partition
void Partitioner::_initNetPartCounts() {
    for (int i = 0; i < _netNum; ++i) {
        _netArray[i]->setPartCount(0, 0);
        _netArray[i]->setPartCount(1, 0);
    }
    for (int i = 0; i < _cellNum; ++i) {
        int part = _cellArray[i]->getPart();
        for (int netId : _cellArray[i]->getNetList())
            _netArray[netId]->incPartCount(part);
    }
}

// ─── Bucket-list helpers ──────────────────────────────────────────────────────

// Remove 'node' (whose current gain key is 'gain') from _bList[part].
void Partitioner::_removeFromBList(int part, Node* node, int gain) {
    Node* prev = node->getPrev();
    Node* next = node->getNext();
    if (prev == nullptr && next == nullptr) {
        _bList[part].erase(gain);           // last node in bucket → remove entry
    } else if (prev == nullptr) {
        _bList[part][gain] = next;          // node was head
        next->setPrev(nullptr);
    } else if (next == nullptr) {
        prev->setNext(nullptr);             // node was tail
    } else {
        prev->setNext(next);
        next->setPrev(prev);
    }
    node->setPrev(nullptr);
    node->setNext(nullptr);
}

// Insert 'node' at the tail of bucket 'gain' in _bList[part].
void Partitioner::_insertIntoBList(int part, Node* node, int gain) {
    node->setPrev(nullptr);
    node->setNext(nullptr);
    if (_bList[part].count(gain) == 0) {
        _bList[part][gain] = node;
    } else { // insert as the head to save time
        Node* curr = _bList[part][gain];
        _bList[part][gain] = node;
        curr->setPrev(node);
        node->setNext(curr);
    }
}

// ─── Gain initialisation ─────────────────────────────────────────────────────

// Compute the initial gain for cell 'cellId' and insert it into the bucket list.
void Partitioner::_computeAndInsertGain(int cellId) {
    int from = _cellArray[cellId]->getPart();
    int to   = 1 - from;
    int gain = 0;
    for (int netId : _cellArray[cellId]->getNetList()) {
        int fc = _netArray[netId]->getPartCount(from);
        int tc = _netArray[netId]->getPartCount(to);
        if (fc == 1) ++gain;   // moving this cell makes the net critical from F-side → cut decreases
        if (tc == 0) --gain;   // net entirely on F-side; moving adds a cut
    }
    _cellArray[cellId]->setGain(gain);
    _insertIntoBList(from, _cellArray[cellId]->getNode(), gain);
}

// ─── Gain propagation after a move ──────────────────────────────────────────

// Standard FM gain-update rules applied to neighbours of the cell that just
// moved from 'from' to 'to'.  This is called BEFORE updating the net's
// PartCount, using the original counts (before the move).
void Partitioner::_updateNeighborGains(int movedCellId, int from, int to) {
    for (int netId : _cellArray[movedCellId]->getNetList()) {
        int fc = _netArray[netId]->getPartCount(from);   // count BEFORE move
        int tc = _netArray[netId]->getPartCount(to);     // count BEFORE move

        // ── Step 1: rules based on TO side BEFORE the move ──────────────────
        if (tc == 0) {
            // The net had no cell on the TO side; every free cell gets +1
            for (int cid : _netArray[netId]->getCellList()) {
                if (_cellArray[cid]->getLock()) continue;
                if (cid == movedCellId) continue;   // skip the cell that just moved
                int part = _cellArray[cid]->getPart();
                int oldGain = _cellArray[cid]->getGain();
                _removeFromBList(part, _cellArray[cid]->getNode(), oldGain);
                _cellArray[cid]->incGain();
                _insertIntoBList(part, _cellArray[cid]->getNode(), oldGain + 1);
            }
        } else if (tc == 1) {
            // The net will become critical on the TO side after the move;
            // only the single cell currently on the TO side (if free) gets -1
            for (int cid : _netArray[netId]->getCellList()) {
                if (_cellArray[cid]->getLock()) continue;
                if (_cellArray[cid]->getPart() != to) continue;
                int oldGain = _cellArray[cid]->getGain();
                _removeFromBList(to, _cellArray[cid]->getNode(), oldGain);
                _cellArray[cid]->decGain();
                _insertIntoBList(to, _cellArray[cid]->getNode(), oldGain - 1);
                break;
            }
        }

        // ── Step 2: update the net's part counts ────────────────────────────
        _netArray[netId]->decPartCount(from);
        _netArray[netId]->incPartCount(to);
        fc = _netArray[netId]->getPartCount(from);   // count AFTER move

        // ── Step 3: rules based on FROM side AFTER the move ─────────────────
        if (fc == 0) {
            // The net now has no cell on the FROM side; every free cell gets -1
            for (int cid : _netArray[netId]->getCellList()) {
                if (_cellArray[cid]->getLock()) continue;
                if (cid == movedCellId) continue;   // skip the cell that just moved
                int part = _cellArray[cid]->getPart();
                int oldGain = _cellArray[cid]->getGain();
                _removeFromBList(part, _cellArray[cid]->getNode(), oldGain);
                _cellArray[cid]->decGain();
                _insertIntoBList(part, _cellArray[cid]->getNode(), oldGain - 1);
            }
        } else if (fc == 1) {
            // The net is now critical on the FROM side;
            // the single remaining cell on FROM side (if free) gets +1
            for (int cid : _netArray[netId]->getCellList()) {
                if (_cellArray[cid]->getLock()) continue;
                if (_cellArray[cid]->getPart() != from) continue;
                int oldGain = _cellArray[cid]->getGain();
                _removeFromBList(from, _cellArray[cid]->getNode(), oldGain);
                _cellArray[cid]->incGain();
                _insertIntoBList(from, _cellArray[cid]->getNode(), oldGain + 1);
                break;
            }
        }
    }
}

// ─── Max-gain cell selection ─────────────────────────────────────────────────

// Whether moving a cell from partition 'from' would keep BOTH sides balanced.
bool Partitioner::_canMoveTo(int from, int lowerBound, int upperBound) const {
    int nextFrom = _partSize[from] - 1;
    int nextTo   = _partSize[1 - from] + 1;
    return (nextFrom >= lowerBound && nextFrom <= upperBound &&
            nextTo   >= lowerBound && nextTo   <= upperBound);
}

// Scan both bucket lists (highest gain first) and set _maxGainCell to the
// best free, balanced candidate.  Sets _maxGainCell = nullptr if none exists.
void Partitioner::_pickMaxGainCell() {
    int lowerBound = (int)ceil((1.0 - _bFactor) / 2.0 * _cellNum);
    int upperBound = (int)floor((1.0 + _bFactor) / 2.0 * _cellNum);

    _maxGainCell = nullptr;
    int bestGain = -_maxPinNum - 1;

    for (int part = 0; part <= 1; ++part) {
        if (_bList[part].empty()) continue;
        if (!_canMoveTo(part, lowerBound, upperBound)) continue;

        // map is ascending; highest gain is at rbegin
        auto it = _bList[part].rbegin();
        if (it->first > bestGain) {
            bestGain     = it->first;
            _maxGainCell = it->second;
        }
    }
}

// ─── Per-pass reset ──────────────────────────────────────────────────────────

// Unlock all cells, clear bucket lists, recompute gains, ready for a new pass.
void Partitioner::_resetPass() {
    for (int i = 0; i < _cellNum; ++i) {
        _cellArray[i]->unlock();
        _cellArray[i]->getNode()->setPrev(nullptr);
        _cellArray[i]->getNode()->setNext(nullptr);
    }
    _bList[0].clear();
    _bList[1].clear();
    _unlockNum[0] = _partSize[0];
    _unlockNum[1] = _partSize[1];

    // Recompute net part counts from current cell partition
    _initNetPartCounts();

    // Recompute gains and rebuild bucket lists
    for (int i = 0; i < _cellNum; ++i)
        _computeAndInsertGain(i);
}

// ─── Public interface ─────────────────────────────────────────────────────────

void Partitioner::init() {
    // Place first half of cells in partition 0, rest in partition 1
    int k = _cellNum / 2;
    _partSize[0] = k;
    _partSize[1] = _cellNum - k;
    /*
    // sequential assign
    for (int i = 0; i < _cellNum; ++i)
        _cellArray[i]->setPart(i < k ? 0 : 1);
    */

    // try other partition method 
    // random shuffle and assign
    random_device rd;
    mt19937 rng(rd());

    vector<int> cellidx(_cellNum);
    iota(cellidx.begin(), cellidx.end(), 0);
    
    shuffle(cellidx.begin(), cellidx.end(), rng);
    for (int i = 0; i < _cellNum; ++i) {
        int idx = cellidx[i];
        _cellArray[idx]->setPart(i < k ? 0 : 1);    
    }

    _findMaxPinNum();
    _initNetPartCounts();
    _findCutSize();

    // Build initial bucket lists
    _bList[0].clear();
    _bList[1].clear();
    for (int i = 0; i < _cellNum; ++i)
        _computeAndInsertGain(i);

    _unlockNum[0] = _partSize[0];
    _unlockNum[1] = _partSize[1];
    _pickMaxGainCell();
}

void Partitioner::partition() {
    init();

    int lowerBound = (int)ceil((1.0 - _bFactor) / 2.0 * _cellNum);
    int upperBound = (int)floor((1.0 + _bFactor) / 2.0 * _cellNum);

    while (true) {
        // ── Start of a new FM pass ───────────────────────────────────────────
        _moveStack.clear();
        _accGain    = 0;
        _maxAccGain = 0;
        _moveNum    = 0;
        _bestMoveNum = 0;

        // Move cells one at a time until all are locked or no valid move exists
        while (true) {
            _pickMaxGainCell();
            if (_maxGainCell == nullptr) break;   // no balanced candidate

            int cellId = _maxGainCell->getId();
            int from   = _cellArray[cellId]->getPart();
            int to     = 1 - from;
            int gain   = _cellArray[cellId]->getGain();

            // Remove from bucket list and lock
            _removeFromBList(from, _maxGainCell, gain);
            _cellArray[cellId]->lock();

            // Update neighbour gains and net counts BEFORE changing the cell's partition
            _updateNeighborGains(cellId, from, to);

            // Move the cell
            _cellArray[cellId]->move();
            _partSize[from]--;
            _partSize[to]++;
            --_unlockNum[from];

            // Record this move
            _moveStack.push_back(cellId);
            ++_moveNum;
            _accGain += gain;
            if (_accGain > _maxAccGain) {
                _maxAccGain  = _accGain;
                _bestMoveNum = _moveNum;
            }
        }

        // ── No improvement in this pass → stop ──────────────────────────────
        if (_maxAccGain <= 0 || _bestMoveNum == 0) break;

        // ── Roll back moves after the best prefix ────────────────────────────
        for (int i = _moveNum - 1; i >= _bestMoveNum; --i) {
            int cellId = _moveStack[i];
            int curPart = _cellArray[cellId]->getPart();
            _cellArray[cellId]->move();
            _partSize[curPart]--;
            _partSize[1 - curPart]++;
        }

        // Update cut size by the accumulated gain of the kept moves
        _cutSize -= _maxAccGain;

        ++_iterNum;

        // Prepare for the next pass
        _resetPass();
    }
}

// End of helper function by me

void Partitioner::printSummary() const
{
    cout << endl;
    cout << "==================== Summary ====================" << endl;
    cout << " Cutsize: " << _cutSize << endl;
    cout << " Total cell number: " << _cellNum << endl;
    cout << " Total net number:  " << _netNum << endl;
    cout << " Cell Number of partition A: " << _partSize[0] << endl;
    cout << " Cell Number of partition B: " << _partSize[1] << endl;
    cout << "=================================================" << endl;
    cout << endl;
    return;
}

void Partitioner::reportNet() const
{
    cout << "Number of nets: " << _netNum << endl;
    for (size_t i = 0, end_i = _netArray.size(); i < end_i; ++i) {
        cout << setw(8) << _netArray[i]->getName() << ": ";
        vector<int> cellList = _netArray[i]->getCellList();
        for (size_t j = 0, end_j = cellList.size(); j < end_j; ++j) {
            cout << setw(8) << _cellArray[cellList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::reportCell() const
{
    cout << "Number of cells: " << _cellNum << endl;
    for (size_t i = 0, end_i = _cellArray.size(); i < end_i; ++i) {
        cout << setw(8) << _cellArray[i]->getName() << ": ";
        vector<int> netList = _cellArray[i]->getNetList();
        for (size_t j = 0, end_j = netList.size(); j < end_j; ++j) {
            cout << setw(8) << _netArray[netList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::writeResult(fstream& outFile)
{
    stringstream buff;
    buff << _cutSize;
    outFile << "Cutsize = " << buff.str() << '\n';
    buff.str("");
    buff << _partSize[0];
    outFile << "G1 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 0) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    buff.str("");
    buff << _partSize[1];
    outFile << "G2 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 1) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    return;
}

void Partitioner::clear()
{
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        delete _cellArray[i];
    }
    for (size_t i = 0, end = _netArray.size(); i < end; ++i) {
        delete _netArray[i];
    }
    return;
}
