#pragma once
#include "Vec3D.h"
#include <vector>
#include <unordered_map>
#include <cmath>
#include <tuple>

class UnionFind {
    private:
        vector<int> parent, rank;

    public:
        explicit UnionFind(int n) : parent(n), rank(n, 0) {
            for (int i = 0; i < n; ++i) parent[i] = i;
        }

        int find(int x) {
            if (parent[x] != x) parent[x] = find(parent[x]);
            return parent[x];
        }

        void unite(int x, int y) {
            int rx = find(x), ry = find(y);
            if (rx == ry) return;
            if (rank[rx] < rank[ry]) parent[rx] = ry;
            else if (rank[rx] > rank[ry]) parent[ry] = rx;
            else { parent[ry] = rx; rank[rx]++; }
        }

        vector<int> getConnectedComponents() {
            vector<int> comp(parent.size());
            for (size_t i = 0; i < parent.size(); ++i) comp[i] = find(i);
            return comp;
        }
};

/// Hash function for 3D integer grid cells
struct CellHash {
    size_t operator()(const tuple<int, int, int>& cell) const noexcept {
        const auto& [x, y, z] = cell;
        // Mix the ints with a prime-based hash
        return ((size_t)x * 73856093) ^ ((size_t)y * 19349663) ^ ((size_t)z * 83492791);
    }
};

/// Cluster computation using a 3D grid / cell list
vector<int> computeClusters(const vector<Vec3D>& coords, double cutoff = 0.40) {
    const int n = (int)coords.size();
    if (n == 0) return {};

    double cutoff2 = cutoff * cutoff;
    double invCell = 1.0 / cutoff;

    // Build hash grid
    unordered_map<tuple<int, int, int>, vector<int>, CellHash> grid;
    grid.reserve(n);

    auto toCell = [&](const Vec3D& v) {
        return make_tuple(
            (int)floor(v.x * invCell),
            (int)floor(v.y * invCell),
            (int)floor(v.z * invCell));
        };

    for (int i = 0; i < n; i++) {
        grid[toCell(coords[i])].push_back(i);
    }

    UnionFind uf(n);

    // Neighbor offsets (27 neighboring cells including self)
    const int offsets[3] = { -1, 0, 1 };

    for (const auto& [cell, indices] : grid) {
        const auto& [cx, cy, cz] = cell;

        // For each water in this cell, check neighbors in surrounding cells
        for (int i : indices) {
            for (int dx : offsets) {
                for (int dy : offsets) {
                    for (int dz : offsets) {
                        tuple<int, int, int> neighCell = { cx + dx, cy + dy, cz + dz };
                        auto it = grid.find(neighCell);
                        if (it == grid.end()) continue;
                        for (int j : it->second) {
                            if (j <= i) continue; // Avoid double-counting
                            double dist2 = coords[i].distance_sq(coords[j]);
                            if (dist2 <= cutoff2) uf.unite(i, j);
                        }
                    }
                }
            }
        }
    }

    return uf.getConnectedComponents();
}