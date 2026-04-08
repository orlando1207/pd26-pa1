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
#include <unordered_set>
#include <unordered_map>
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

void Partitioner::_runFM(int maxPasses) {
    int passCount = 0;
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
        ++passCount;

        if (maxPasses > 0 && passCount >= maxPasses) break;

        _resetPass();
    }
}

// ─── Public interface: multilevel + multi-restart FM ─────────────────────────

// ═══════════════════════════════════════════════════════════════════════════════
// Multilevel coarsening via Heavy-Edge Matching (HEM)
// ═══════════════════════════════════════════════════════════════════════════════

CoarseLevel Partitioner::_coarsen(const vector<int>& prevMapping, int prevClusters,
                                   const vector<vector<int>>& prevClusterNets,
                                   const vector<vector<int>>& prevNetClusters,
                                   int prevNumNets,
                                   const vector<int>& prevClusterWeight,
                                   int seed) {
    CoarseLevel level;
    int n = prevClusters;

    // Heavy-edge matching using net-based neighbor enumeration (avoids explicit adjacency)
    // For each unmatched cluster, find the best neighbor by scanning its nets
    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    mt19937 rng(seed);
    shuffle(order.begin(), order.end(), rng);

    vector<int> match(n, -1);
    int maxClusterWeight = 0;
    for (int i = 0; i < n; ++i)
        maxClusterWeight = max(maxClusterWeight, prevClusterWeight[i]);
    int weightLimit = max(maxClusterWeight * 3, _cellNum / 4);

    // Temporary map for accumulating edge weights per neighbor
    vector<double> neighborWeight(n, 0.0);
    vector<int> neighbors;  // track which neighbors we touched

    for (int u : order) {
        if (match[u] != -1) continue;
        neighbors.clear();

        // Scan nets of cluster u to find neighbors and their weights
        for (int netId : prevClusterNets[u]) {
            const vector<int>& clusters = prevNetClusters[netId];
            int sz = (int)clusters.size();
            if (sz <= 1 || sz > 100) continue;  // skip trivial or huge nets
            double w = 1.0 / (sz - 1);
            for (int v : clusters) {
                if (v == u) continue;
                if (match[v] != -1) continue;
                if (neighborWeight[v] == 0.0)
                    neighbors.push_back(v);
                neighborWeight[v] += w;
            }
        }

        int bestV = -1;
        double bestW = 0.0;
        for (int v : neighbors) {
            if (prevClusterWeight[u] + prevClusterWeight[v] <= weightLimit) {
                if (neighborWeight[v] > bestW) {
                    bestW = neighborWeight[v];
                    bestV = v;
                }
            }
            neighborWeight[v] = 0.0;  // reset for next iteration
        }

        if (bestV != -1) {
            match[u] = bestV;
            match[bestV] = u;
        }
    }

    // Assign cluster IDs for this level
    vector<int> clusterMap(n, -1);
    int numClusters = 0;
    for (int i = 0; i < n; ++i) {
        if (clusterMap[i] != -1) continue;
        if (match[i] != -1 && match[i] > i) {
            clusterMap[i] = numClusters;
            clusterMap[match[i]] = numClusters;
        } else if (match[i] == -1) {
            clusterMap[i] = numClusters;
        } else {
            continue;
        }
        ++numClusters;
    }

    level.numClusters = numClusters;
    level.clusterWeight.assign(numClusters, 0);

    // Map original cells through: old cluster -> new cluster
    level.cellToCluster.resize(_cellNum);
    level.clusterCells.resize(numClusters);
    for (int i = 0; i < _cellNum; ++i) {
        int oldCluster = prevMapping[i];
        int newCluster = clusterMap[oldCluster];
        level.cellToCluster[i] = newCluster;
        level.clusterCells[newCluster].push_back(i);
    }
    for (int i = 0; i < numClusters; ++i)
        level.clusterWeight[i] = (int)level.clusterCells[i].size();

    // Build coarsened hypergraph: remove nets that become single-cluster
    level.clusterNets.resize(numClusters);
    level.numCoarseNets = 0;

    // Use a temporary vector+counter for deduplication (faster than unordered_set)
    vector<int> seen(numClusters, -1);
    for (int netId = 0; netId < prevNumNets; ++netId) {
        const vector<int>& oldClusters = prevNetClusters[netId];
        vector<int> newClusters;
        for (int c : oldClusters) {
            int nc2 = clusterMap[c];
            if (seen[nc2] != netId) {
                seen[nc2] = netId;
                newClusters.push_back(nc2);
            }
        }
        if ((int)newClusters.size() <= 1) continue;
        int cNetId = level.numCoarseNets++;
        level.coarseNetClusters.push_back(newClusters);
        for (int c : newClusters)
            level.clusterNets[c].push_back(cNetId);
    }

    return level;
}

// ═══════════════════════════════════════════════════════════════════════════════
// Run FM on a coarse-level graph (weighted clusters as "cells")
// Lightweight weighted FM with proper gain tracking
// ═══════════════════════════════════════════════════════════════════════════════

void Partitioner::_runCoarseFM(CoarseLevel& level, int seed,
                                const vector<int>* initPart, int maxPasses) {
    int nc = level.numClusters;
    int nn = level.numCoarseNets;
    if (nc <= 1) return;

    // Compute max pin num for coarse graph
    int maxPin = 0;
    for (int i = 0; i < nc; ++i)
        maxPin = max(maxPin, (int)level.clusterNets[i].size());
    if (maxPin == 0) return;

    int bucketSize = 2 * maxPin + 1;

    // Partition state
    vector<int> cPart(nc);  // 0 or 1
    vector<int> cGain(nc, 0);
    vector<bool> cLock(nc, false);
    int wSize[2] = {0, 0};

    // Net partition counts
    vector<int> netPC[2];
    netPC[0].resize(nn, 0);
    netPC[1].resize(nn, 0);

    int totalWeight = 0;
    for (int i = 0; i < nc; ++i) totalWeight += level.clusterWeight[i];

    if (initPart != nullptr) {
        // Use provided initial partition
        for (int c = 0; c < nc; ++c) {
            cPart[c] = (*initPart)[c];
            wSize[cPart[c]] += level.clusterWeight[c];
        }
    } else {
        // Random balanced initial partition
        int target = totalWeight / 2;
        vector<int> cidx(nc);
        iota(cidx.begin(), cidx.end(), 0);
        mt19937 rng(seed);
        shuffle(cidx.begin(), cidx.end(), rng);

        int accW = 0;
        for (int i = 0; i < nc; ++i) {
            int c = cidx[i];
            if (accW + level.clusterWeight[c] <= target) {
                cPart[c] = 0;
                wSize[0] += level.clusterWeight[c];
                accW += level.clusterWeight[c];
            } else {
                cPart[c] = 1;
                wSize[1] += level.clusterWeight[c];
            }
        }
    }

    int lowerBound = (int)ceil((1.0 - _bFactor) / 2.0 * totalWeight);
    int upperBound = (int)floor((1.0 + _bFactor) / 2.0 * totalWeight);

    // Doubly-linked list nodes for bucket list (index by cluster ID)
    struct BNode { int prev, next; };
    vector<BNode> nodes(nc);
    vector<int> bHead[2]; // bHead[part][idx] = head of list, -1 if empty
    bHead[0].assign(bucketSize, -1);
    bHead[1].assign(bucketSize, -1);
    int maxGP[2] = {-1, -1}; // max occupied bucket index

    auto bInsert = [&](int part, int c, int gain) {
        int idx = gain + maxPin;
        nodes[c].prev = -1;
        nodes[c].next = bHead[part][idx];
        if (bHead[part][idx] != -1)
            nodes[bHead[part][idx]].prev = c;
        bHead[part][idx] = c;
        if (idx > maxGP[part]) maxGP[part] = idx;
    };

    auto bRemove = [&](int part, int c, int gain) {
        int idx = gain + maxPin;
        int p = nodes[c].prev, nx = nodes[c].next;
        if (p == -1) bHead[part][idx] = nx;
        else nodes[p].next = nx;
        if (nx != -1) nodes[nx].prev = p;
        nodes[c].prev = -1; nodes[c].next = -1;
    };

    auto initPass = [&]() {
        fill(cLock.begin(), cLock.end(), false);
        fill(bHead[0].begin(), bHead[0].end(), -1);
        fill(bHead[1].begin(), bHead[1].end(), -1);
        maxGP[0] = -1; maxGP[1] = -1;

        for (int n2 = 0; n2 < nn; ++n2) { netPC[0][n2] = 0; netPC[1][n2] = 0; }
        for (int c = 0; c < nc; ++c) {
            for (int netId : level.clusterNets[c])
                ++netPC[cPart[c]][netId];
        }

        for (int c = 0; c < nc; ++c) {
            int from = cPart[c], to = 1 - from;
            int g = 0;
            for (int netId : level.clusterNets[c]) {
                if (netPC[from][netId] == 1) ++g;
                if (netPC[to][netId] == 0) --g;
            }
            cGain[c] = g;
            bInsert(from, c, g);
        }
    };

    initPass();

    // FM passes
    for (int pass = 0; pass < maxPasses; ++pass) {
        if (pass > 0) initPass();

        vector<int> moveStack;
        int accGain2 = 0, maxAccGain2 = 0, bestMoveNum = 0, moveNum = 0;

        for (int step = 0; step < nc; ++step) {
            int bestC = -1, bestG = -maxPin - 1;
            for (int part = 0; part <= 1; ++part) {
                while (maxGP[part] >= 0 && bHead[part][maxGP[part]] == -1)
                    --maxGP[part];
                if (maxGP[part] < 0) continue;
                int g = maxGP[part] - maxPin;

                // Check balance for moving from this partition
                // Find the first unlocked cell in this bucket
                int c = bHead[part][maxGP[part]];
                if (c == -1) continue;

                int newFrom = wSize[part] - level.clusterWeight[c];
                int newTo = wSize[1-part] + level.clusterWeight[c];
                if (newFrom >= lowerBound && newFrom <= upperBound &&
                    newTo >= lowerBound && newTo <= upperBound) {
                    if (g > bestG) {
                        bestG = g;
                        bestC = c;
                    }
                }
            }
            if (bestC == -1) break;

            int from = cPart[bestC], to = 1 - from;
            bRemove(from, bestC, cGain[bestC]);
            cLock[bestC] = true;

            // Update gains of neighbors (proper FM gain update)
            for (int netId : level.clusterNets[bestC]) {
                int fc = netPC[from][netId];
                int tc = netPC[to][netId];

                // Before move: check critical nets
                if (tc == 0) {
                    for (int v : level.coarseNetClusters[netId]) {
                        if (v == bestC || cLock[v]) continue;
                        int oldG = cGain[v];
                        bRemove(cPart[v], v, oldG);
                        ++cGain[v];
                        bInsert(cPart[v], v, oldG + 1);
                    }
                } else if (tc == 1) {
                    for (int v : level.coarseNetClusters[netId]) {
                        if (cLock[v] || cPart[v] != to) continue;
                        int oldG = cGain[v];
                        bRemove(to, v, oldG);
                        --cGain[v];
                        bInsert(to, v, oldG - 1);
                        break;
                    }
                }

                --netPC[from][netId];
                ++netPC[to][netId];
                fc = netPC[from][netId];

                if (fc == 0) {
                    for (int v : level.coarseNetClusters[netId]) {
                        if (v == bestC || cLock[v]) continue;
                        int oldG = cGain[v];
                        bRemove(cPart[v], v, oldG);
                        --cGain[v];
                        bInsert(cPart[v], v, oldG - 1);
                    }
                } else if (fc == 1) {
                    for (int v : level.coarseNetClusters[netId]) {
                        if (cLock[v] || cPart[v] != from) continue;
                        int oldG = cGain[v];
                        bRemove(from, v, oldG);
                        ++cGain[v];
                        bInsert(from, v, oldG + 1);
                        break;
                    }
                }
            }

            cPart[bestC] = to;
            wSize[from] -= level.clusterWeight[bestC];
            wSize[to] += level.clusterWeight[bestC];

            moveStack.push_back(bestC);
            ++moveNum;
            accGain2 += bestG;
            if (accGain2 > maxAccGain2) {
                maxAccGain2 = accGain2;
                bestMoveNum = moveNum;
            }
        }

        if (maxAccGain2 <= 0 || bestMoveNum == 0) {
            // Undo all moves
            for (int i = (int)moveStack.size() - 1; i >= 0; --i) {
                int c = moveStack[i];
                int curP = cPart[c], prevP = 1 - curP;
                for (int netId : level.clusterNets[c]) {
                    --netPC[curP][netId];
                    ++netPC[prevP][netId];
                }
                cPart[c] = prevP;
                wSize[curP] -= level.clusterWeight[c];
                wSize[prevP] += level.clusterWeight[c];
            }
            break;
        }

        // Undo moves after bestMoveNum
        for (int i = (int)moveStack.size() - 1; i >= bestMoveNum; --i) {
            int c = moveStack[i];
            int curP = cPart[c], prevP = 1 - curP;
            for (int netId : level.clusterNets[c]) {
                --netPC[curP][netId];
                ++netPC[prevP][netId];
            }
            cPart[c] = prevP;
            wSize[curP] -= level.clusterWeight[c];
            wSize[prevP] += level.clusterWeight[c];
        }
    }

    // Apply coarse partition to original cells
    _partSize[0] = 0;
    _partSize[1] = 0;
    for (int c = 0; c < nc; ++c) {
        for (int cellId : level.clusterCells[c]) {
            _cellArray[cellId]->setPart(cPart[c]);
            if (cPart[c] == 0) ++_partSize[0];
            else ++_partSize[1];
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Full multilevel FM
// ═══════════════════════════════════════════════════════════════════════════════

void Partitioner::_multilevelPartition(int seed) {
    // Build initial (level-0) coarse structures from the original graph
    // Level 0: each cell is its own cluster
    vector<CoarseLevel> levels;

    CoarseLevel level0;
    level0.numClusters = _cellNum;
    level0.cellToCluster.resize(_cellNum);
    level0.clusterCells.resize(_cellNum);
    level0.clusterWeight.resize(_cellNum, 1);
    level0.clusterNets.resize(_cellNum);
    for (int i = 0; i < _cellNum; ++i) {
        level0.cellToCluster[i] = i;
        level0.clusterCells[i].push_back(i);
        level0.clusterNets[i] = _cellArray[i]->getNetList();
    }
    level0.numCoarseNets = _netNum;
    level0.coarseNetClusters.resize(_netNum);
    for (int n = 0; n < _netNum; ++n) {
        level0.coarseNetClusters[n] = _netArray[n]->getCellList();
    }
    levels.push_back(move(level0));

    // ─── Coarsening phase ────────────────────────────────────────────────────
    int minClusters = max(200, _cellNum / 30);  // stop coarsening when small enough
    mt19937 rng(seed);

    while (true) {
        CoarseLevel& prev = levels.back();
        if (prev.numClusters <= minClusters) break;

        CoarseLevel next = _coarsen(
            prev.cellToCluster, prev.numClusters,
            prev.clusterNets, prev.coarseNetClusters,
            prev.numCoarseNets, prev.clusterWeight,
            rng());

        // Stop if coarsening didn't shrink enough (< 10% reduction)
        if (next.numClusters > prev.numClusters * 9 / 10) break;

        levels.push_back(move(next));
    }

    // ─── Initial partition on coarsest level ─────────────────────────────────
    CoarseLevel& coarsest = levels.back();
    _runCoarseFM(coarsest, rng(), nullptr, 15);

    // ─── Uncoarsening phase (V-cycle): project and refine at each level ─────
    for (int lv = (int)levels.size() - 2; lv >= 0; --lv) {
        CoarseLevel& finer = levels[lv];
        CoarseLevel& coarser = levels[lv + 1];

        if (lv > 0) {
            // At intermediate levels: project coarse partition to this level's clusters,
            // then refine with coarse FM (limited passes)

            // Build cluster partition for this level from original cell partition
            // (which was set by the coarser level FM)
            int nc = finer.numClusters;
            vector<int> clusterPart(nc);
            for (int c = 0; c < nc; ++c) {
                // Use majority vote of original cells in this cluster
                int cnt[2] = {0, 0};
                for (int cellId : finer.clusterCells[c])
                    ++cnt[_cellArray[cellId]->getPart()];
                clusterPart[c] = (cnt[1] > cnt[0]) ? 1 : 0;
            }

            _runCoarseFM(finer, rng(), &clusterPart, 3);
        } else {
            // At level 0 (original graph): refine with full FM
            _partSize[0] = 0;
            _partSize[1] = 0;
            for (int i = 0; i < _cellNum; ++i) {
                if (_cellArray[i]->getPart() == 0) ++_partSize[0];
                else ++_partSize[1];
            }
            _prepareForFM();
            _runFM();
        }
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

    // ─── Phase 1: Multilevel FM restarts ─────────────────────────────────────
    int mlRestarts = (_cellNum > 50000) ? 2 : (_cellNum > 5000) ? 3 : 3;
    for (int r = 0; r < mlRestarts; ++r) {
        auto now = chrono::steady_clock::now();
        double elapsed = chrono::duration<double>(now - startTime).count();
        if (elapsed > timeLimitSec * 0.5) break;

        _iterNum = 0;
        _multilevelPartition(baseSeed + r * 1000);

        // Compute final cut size
        _initNetPartCounts();
        _findCutSize();

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

    // ─── Phase 2: Regular multi-restart FM (perturb from best) ──────────────
    for (int restart = 0; ; ++restart) {
        auto now = chrono::steady_clock::now();
        double elapsed = chrono::duration<double>(now - startTime).count();
        if (elapsed > timeLimitSec) break;
        double avgTime = (restart > 0) ? elapsed / restart : 0;
        if (restart > 0 && elapsed + avgTime * 1.2 > timeLimitSec) break;
        if (noImproveCount >= patience) break;

        _iterNum = 0;

        if (restart < 2 || bestCutSize == INT_MAX) {
            init(baseSeed + 5000 + restart);
        } else if ((restart % 4) == 0) {
            init(baseSeed + 5000 + restart);
        } else {
            double ratio = 0.03 + 0.02 * (noImproveCount / 5);
            if (ratio > 0.15) ratio = 0.15;
            _initFromPerturb(bestPartition, bestPartSize[0], bestPartSize[1],
                             baseSeed + 5000 + restart, ratio);
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
