#include "Core.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

namespace Core {

    struct Edge {
        int u, v; // Indices of the points
        float weight;
    };

    class UnionFind {
    public:
        explicit UnionFind(int n) : parent(n), rank(n, 0) {
            for (int i = 0; i < n; ++i) {
                parent[i] = i;
            }
        }

        int find(int x) {
            if (parent[x] != x) {
                parent[x] = find(parent[x]); // Path compression
            }
            return parent[x];
        }

        void unite(int x, int y) {
            int rootX = find(x);
            int rootY = find(y);
            if (rootX != rootY) {
                if (rank[rootX] < rank[rootY]) {
                    parent[rootX] = rootY;
                }
                else if (rank[rootX] > rank[rootY]) {
                    parent[rootY] = rootX;
                }
                else {
                    parent[rootY] = rootX;
                    ++rank[rootX];
                }
            }
        }

        std::vector<int> getConnectedComponents() {
            std::vector<int> component(parent.size());
            for (int i = 0; i < parent.size(); ++i) {
                component[i] = find(i);
            }
            return component;
        }

    private:
        std::vector<int> parent;
        std::vector<int> rank;
    };

    // Clusters found via Kruskals algo
    std::vector<int> computeClusters(const std::vector<std::array<float, 3>>& points) {
        float cutoff = 0.40f; // Cutoff distance for H-bond between waters
        int n = points.size();
        std::vector<Edge> edges;

        // Create edges with distance less than the cutoff
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
				float dist = getDistance(points[i], points[j]);
                if (dist <= cutoff) {
                    edges.push_back({ i, j, dist });
                }
            }
        }

        // Initialize Union-Find
        UnionFind uf(n);

        // Process edges to unite connected components
        for (const auto& edge : edges) {
            uf.unite(edge.u, edge.v);
        }

        // Get connected components
        return uf.getConnectedComponents();
    }
}