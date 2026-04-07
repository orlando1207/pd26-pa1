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
#include <climits>
#include <chrono>
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

// ─── Small utilities ──────────────────────────────────────────────────────────

void Partitioner::_findMaxPinNum() {
    _maxPinNum = 0;
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

void Partitioner::_initNetPartCounts() {
    for (int i = 0; i < _netNum; ++i) {
        _netArray[i]->setPartCount(0, 0);
        _netArray[i]->setPartCount(1, 0);
    }
    for (int i = 0; i < _cellNum; ++i) {
        int part = _cellArray[i]->getPart();
        const vector<int>& nets = _cellArray[i]->getNetList();
        for (int netId : nets)
            _netArray[netId]->incPartCount(part);
    }
}

// ─── Array-based bucket-list helpers ──────────────────────────────────────────

void Partitioner::_removeFromBList(int part, Node* node, int gain) {
    int idx = gain + _maxPinNum;
    Node* prev = node->getPrev();
    Node* next = node->getNext();
    if (prev == nullptr && next == nullptr) {
        _bList[part][idx] = nullptr;
    } else if (prev == nullptr) {
        _bList[part][idx] = next;
        next->setPrev(nullptr);
    } else if (next == nullptr) {
        prev->setNext(nullptr);
    } else {
        prev->setNext(next);
        next->setPrev(prev);
    }
    node->setPrev(nullptr);
    node->setNext(nullptr);
}

void Partitioner::_insertIntoBList(int part, Node* node, int gain) {
    int idx = gain + _maxPinNum;
    node->setPrev(nullptr);
    node->setNext(nullptr);
    if (_bList[part][idx] == nullptr) {
        _bList[part][idx] = node;
    } else {
        Node* curr = _bList[part][idx];
        _bList[part][idx] = node;
        node->setNext(curr);
        curr->setPrev(node);
    }
    if (idx > _maxGainPtr[part])
        _maxGainPtr[part] = idx;
}

// ─── Gain initialisation ─────────────────────────────────────────────────────

void Partitioner::_computeAndInsertGain(int cellId) {
    int from = _cellArray[cellId]->getPart();
    int to   = 1 - from;
    int gain = 0;
    const vector<int>& nets = _cellArray[cellId]->getNetList();
    for (int netId : nets) {
        int fc = _netArray[netId]->getPartCount(from);
        int tc = _netArray[netId]->getPartCount(to);
        if (fc == 1) ++gain;
        if (tc == 0) --gain;
    }
    _cellArray[cellId]->setGain(gain);
    _insertIntoBList(from, _cellArray[cellId]->getNode(), gain);
}

// ─── Gain propagation after a move ──────────────────────────────────────────

void Partitioner::_updateNeighborGains(int movedCellId, int from, int to) {
    const vector<int>& nets = _cellArray[movedCellId]->getNetList();
    for (int netId : nets) {
        int fc = _netArray[netId]->getPartCount(from);
        int tc = _netArray[netId]->getPartCount(to);
        const vector<int>& cells = _netArray[netId]->getCellList();

        if (tc == 0) {
            for (int cid : cells) {
                if (cid == movedCellId) continue;
                if (_cellArray[cid]->getLock()) continue;
                int part = _cellArray[cid]->getPart();
                int oldGain = _cellArray[cid]->getGain();
                _removeFromBList(part, _cellArray[cid]->getNode(), oldGain);
                _cellArray[cid]->incGain();
                _insertIntoBList(part, _cellArray[cid]->getNode(), oldGain + 1);
            }
        } else if (tc == 1) {
            for (int cid : cells) {
                if (_cellArray[cid]->getLock()) continue;
                if (_cellArray[cid]->getPart() != to) continue;
                int oldGain = _cellArray[cid]->getGain();
                _removeFromBList(to, _cellArray[cid]->getNode(), oldGain);
                _cellArray[cid]->decGain();
                _insertIntoBList(to, _cellArray[cid]->getNode(), oldGain - 1);
                break;
            }
        }

        _netArray[netId]->decPartCount(from);
        _netArray[netId]->incPartCount(to);
        fc = _netArray[netId]->getPartCount(from);

        if (fc == 0) {
            for (int cid : cells) {
                if (cid == movedCellId) continue;
                if (_cellArray[cid]->getLock()) continue;
                int part = _cellArray[cid]->getPart();
                int oldGain = _cellArray[cid]->getGain();
                _removeFromBList(part, _cellArray[cid]->getNode(), oldGain);
                _cellArray[cid]->decGain();
                _insertIntoBList(part, _cellArray[cid]->getNode(), oldGain - 1);
            }
        } else if (fc == 1) {
            for (int cid : cells) {
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

bool Partitioner::_canMoveTo(int from, int lowerBound, int upperBound) const {
    int nextFrom = _partSize[from] - 1;
    int nextTo   = _partSize[1 - from] + 1;
    return (nextFrom >= lowerBound && nextFrom <= upperBound &&
            nextTo   >= lowerBound && nextTo   <= upperBound);
}

void Partitioner::_pickMaxGainCell() {
    int lowerBound = (int)ceil((1.0 - _bFactor) / 2.0 * _cellNum);
    int upperBound = (int)floor((1.0 + _bFactor) / 2.0 * _cellNum);

    _maxGainCell = nullptr;
    int bestGain = -_maxPinNum - 1;

    for (int part = 0; part <= 1; ++part) {
        if (!_canMoveTo(part, lowerBound, upperBound)) continue;

        while (_maxGainPtr[part] >= 0 && _bList[part][_maxGainPtr[part]] == nullptr)
            --_maxGainPtr[part];

        if (_maxGainPtr[part] < 0) continue;

        int gain = _maxGainPtr[part] - _maxPinNum;
        if (gain > bestGain) {
            bestGain     = gain;
            _maxGainCell = _bList[part][_maxGainPtr[part]];
        }
    }
}

// ─── Per-pass reset ──────────────────────────────────────────────────────────

void Partitioner::_resetPass() {
    for (int i = 0; i < _cellNum; ++i) {
        _cellArray[i]->unlock();
        _cellArray[i]->getNode()->setPrev(nullptr);
        _cellArray[i]->getNode()->setNext(nullptr);
    }
    fill(_bList[0].begin(), _bList[0].end(), nullptr);
    fill(_bList[1].begin(), _bList[1].end(), nullptr);
    _maxGainPtr[0] = -1;
    _maxGainPtr[1] = -1;
    _unlockNum[0] = _partSize[0];
    _unlockNum[1] = _partSize[1];

    _initNetPartCounts();

    for (int i = 0; i < _cellNum; ++i)
        _computeAndInsertGain(i);
}

// ─── Init with a given seed ──────────────────────────────────────────────────

void Partitioner::init(int seed) {
    int k = _cellNum / 2;
    _partSize[0] = k;
    _partSize[1] = _cellNum - k;

    mt19937 rng(seed);
    vector<int> cellidx(_cellNum);
    iota(cellidx.begin(), cellidx.end(), 0);
    shuffle(cellidx.begin(), cellidx.end(), rng);
    for (int i = 0; i < _cellNum; ++i) {
        int idx = cellidx[i];
        _cellArray[idx]->setPart(i < k ? 0 : 1);
    }

    _prepareForFM();
}

// Perturb best partition: swap a fraction of cells between sides
void Partitioner::_initFromPerturb(const vector<bool>& bestPart,
                                    int bestPS0, int bestPS1,
                                    int seed, double ratio) {
    for (int i = 0; i < _cellNum; ++i)
        _cellArray[i]->setPart(bestPart[i]);
    _partSize[0] = bestPS0;
    _partSize[1] = bestPS1;

    vector<int> side0, side1;
    side0.reserve(_partSize[0]);
    side1.reserve(_partSize[1]);
    for (int i = 0; i < _cellNum; ++i) {
        if (_cellArray[i]->getPart() == 0) side0.push_back(i);
        else side1.push_back(i);
    }

    mt19937 rng(seed);
    shuffle(side0.begin(), side0.end(), rng);
    shuffle(side1.begin(), side1.end(), rng);

    int swapCount = max(1, (int)(min((int)side0.size(), (int)side1.size()) * ratio));
    for (int i = 0; i < swapCount; ++i) {
        _cellArray[side0[i]]->setPart(1);
        _cellArray[side1[i]]->setPart(0);
    }

    _prepareForFM();
}

// Common setup after cell partition is set
void Partitioner::_prepareForFM() {
    _initNetPartCounts();
    _findCutSize();

    fill(_bList[0].begin(), _bList[0].end(), nullptr);
    fill(_bList[1].begin(), _bList[1].end(), nullptr);
    _maxGainPtr[0] = -1;
    _maxGainPtr[1] = -1;

    for (int i = 0; i < _cellNum; ++i) {
        _cellArray[i]->unlock();
        _cellArray[i]->getNode()->setPrev(nullptr);
        _cellArray[i]->getNode()->setNext(nullptr);
    }
    for (int i = 0; i < _cellNum; ++i)
        _computeAndInsertGain(i);

    _unlockNum[0] = _partSize[0];
    _unlockNum[1] = _partSize[1];
}

// ─── Single FM run (multiple passes until no improvement) ────────────────────

void Partitioner::_runFM() {
    while (true) {
        _moveStack.clear();
        _accGain    = 0;
        _maxAccGain = 0;
        _moveNum    = 0;
        _bestMoveNum = 0;

        while (true) {
            _pickMaxGainCell();
            if (_maxGainCell == nullptr) break;

            int cellId = _maxGainCell->getId();
            int from   = _cellArray[cellId]->getPart();
            int to     = 1 - from;
            int gain   = _cellArray[cellId]->getGain();

            _removeFromBList(from, _maxGainCell, gain);
            _cellArray[cellId]->lock();

            _updateNeighborGains(cellId, from, to);

            _cellArray[cellId]->move();
            _partSize[from]--;
            _partSize[to]++;
            --_unlockNum[from];

            _moveStack.push_back(cellId);
            ++_moveNum;
            _accGain += gain;
            if (_accGain > _maxAccGain) {
                _maxAccGain  = _accGain;
                _bestMoveNum = _moveNum;
            }
        }

        if (_maxAccGain <= 0 || _bestMoveNum == 0) break;

        for (int i = _moveNum - 1; i >= _bestMoveNum; --i) {
            int cellId = _moveStack[i];
            int curPart = _cellArray[cellId]->getPart();
            _cellArray[cellId]->move();
            _partSize[curPart]--;
            _partSize[1 - curPart]++;
        }

        _cutSize -= _maxAccGain;
        ++_iterNum;

        _resetPass();
    }
}

// ─── Public interface: multi-restart FM ──────────────────────────────────────

void Partitioner::partition() {
    _findMaxPinNum();

    int bucketSize = 2 * _maxPinNum + 1;
    _bList[0].assign(bucketSize, nullptr);
    _bList[1].assign(bucketSize, nullptr);

    _moveStack.reserve(_cellNum);

    int bestCutSize = INT_MAX;
    vector<bool> bestPartition(_cellNum);
    int bestPartSize[2] = {0, 0};

    auto startTime = chrono::steady_clock::now();
    double timeLimitSec = 280.0;
    int patience = 20;       // stop after this many restarts without improvement
    int noImproveCount = 0;

    random_device rd;
    int baseSeed = rd();

    for (int restart = 0; ; ++restart) {
        // Time check after first run
        if (restart >= 1) {
            auto now = chrono::steady_clock::now();
            double elapsed = chrono::duration<double>(now - startTime).count();
            if (elapsed > timeLimitSec) break;
            // Estimate if next run would exceed limit
            double avgTime = elapsed / restart;
            if (elapsed + avgTime * 1.2 > timeLimitSec) break;
        }
        // Patience check
        if (noImproveCount >= patience) break;

        _iterNum = 0;

        // First 3 restarts: random. After that, alternate: 75% perturb, 25% random
        if (restart < 3 || bestCutSize == INT_MAX) {
            init(baseSeed + restart);
        } else if ((restart % 4) == 0) {
            init(baseSeed + restart);
        } else {
            double ratio = 0.03 + 0.02 * (noImproveCount / 5);  // escalate perturbation
            if (ratio > 0.15) ratio = 0.15;
            _initFromPerturb(bestPartition, bestPartSize[0], bestPartSize[1],
                             baseSeed + restart, ratio);
        }
        _runFM();

        if (_cutSize < bestCutSize) {
            bestCutSize = _cutSize;
            bestPartSize[0] = _partSize[0];
            bestPartSize[1] = _partSize[1];
            for (int i = 0; i < _cellNum; ++i)
                bestPartition[i] = _cellArray[i]->getPart();
            noImproveCount = 0;
        } else {
            ++noImproveCount;
        }
    }

    _cutSize = bestCutSize;
    _partSize[0] = bestPartSize[0];
    _partSize[1] = bestPartSize[1];
    for (int i = 0; i < _cellNum; ++i)
        _cellArray[i]->setPart(bestPartition[i]);
}

// ─── Reporting ───────────────────────────────────────────────────────────────

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
        const vector<int>& cellList = _netArray[i]->getCellList();
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
        const vector<int>& netList = _cellArray[i]->getNetList();
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
