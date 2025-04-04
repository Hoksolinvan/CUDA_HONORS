#include <vector>
#include <array>
#include <cstring>
#include <tuple>
#include <unordered_map>
#include <string>
#include <cmath>
#include <utility>
#include <omp.h>
#include <stdlib.h>
#include <set>
#include <mach/mach.h>
#include <iostream>
#include <string>
#include <chrono>


// ============================================ Euclidean Distance ====================================================
// ====================================================================================================================

double Euclidean_distance(std::pair<double, double> value_one, std::pair<double, double> value_two) {
    return sqrt(pow(value_one.first - value_two.first, 2) + pow(value_one.second - value_two.second, 2));
}
// =========================================== End Euclidean Distance ================================================
// ===================================================================================================================


void print_memory_usage(const std::string& label) {
    task_basic_info_data_t info;
    mach_msg_type_number_t infoCount = TASK_BASIC_INFO_COUNT;

    if (task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&info, &infoCount) == KERN_SUCCESS) {
        std::cout << label << " - Resident memory: " << info.resident_size / 1024.0 << " KB" << std::endl;
    } else {
        std::cerr << label << " - Failed to get memory usage info." << std::endl;
    }
}



double compute_silhouette_score(const std::vector<std::pair<double, double>>& points, const std::vector<int>& labels) {
    int n = points.size();
    double total_score = 0.0;
    int valid_points = 0;

    for (int i = 0; i < n; ++i) {
        int label_i = labels[i];
        if (label_i == -1) continue;  // Skip noise

        double a = 0.0;
        double b = std::numeric_limits<double>::max();
        int same_cluster_count = 0;

        std::unordered_map<int, std::pair<double, int>> cluster_distance_map;

        for (int j = 0; j < n; ++j) {
            if (i == j || labels[j] == -1) continue;

            double dist = Euclidean_distance(points[i], points[j]);
            int label_j = labels[j];

            if (label_i == label_j) {
                a += dist;
                ++same_cluster_count;
            } else {
                cluster_distance_map[label_j].first += dist;
                cluster_distance_map[label_j].second += 1;
            }
        }

        if (same_cluster_count > 0) a /= same_cluster_count;

        for (const auto& kv : cluster_distance_map) {
            double avg_dist = kv.second.first / kv.second.second;
            if (avg_dist < b) b = avg_dist;
        }

        double s = (b - a) / std::max(a, b);
        total_score += s;
        ++valid_points;
    }

    return (valid_points > 0) ? total_score / valid_points : 0.0;
}


// ============================================ K-D tree implementation =============================================
// ==================================================================================================================

struct Node {

    int index_position;
    double values[2];
    Node* left;
    Node* right;


};

struct Node*  createNode(std::pair<double,double> arr,int index_position) {
    struct Node* temp = new Node;
    temp->values[0]=arr.first;
    temp->values[1]=arr.second;
    temp->index_position = index_position;

    temp->left = temp->right = nullptr;

    return temp;
}


struct Node* insertTree(struct Node* root, std::pair<double,double> values, int index_position, int depth){

    if(root==nullptr)
    {
        root = createNode(values,index_position);
        return root;
    }

    int cd = depth % 2;

    if(cd==0){
        if(values.first <= root->values[0]){
        root->left = insertTree(root->left,values,index_position,depth+1);
        }
        else{
        root->right = insertTree(root->right,values,index_position,depth+1);
        }
    }
    else{
        if(values.second <= root->values[1]){
        root->left = insertTree(root->left,values,index_position,depth+1);
        }
        else{
        root->right = insertTree(root->right,values,index_position,depth+1);
        }

    }

    return root;

}


void findNeighbours(struct Node* root, std::pair<double, double> interest_point, double radius, int depth, std::vector<std::pair<double, double> >& neighbours) {

    if(root == nullptr){
        return;
    }

    int cd = depth % 2;
    double diff = (cd == 0) ? (interest_point.first - root->values[0]) : (interest_point.second - root->values[1]);


    if(Euclidean_distance(interest_point, std::make_pair(root->values[0],root->values[1])) <= radius){
        neighbours.push_back(std::make_pair(root->values[0],root->index_position));
    }


    if (diff <= radius) {
        findNeighbours(root->left, interest_point, radius, depth + 1, neighbours);
    }
    if (diff >= -radius) {
        findNeighbours(root->right, interest_point, radius, depth + 1, neighbours);
    }

}





// Get Neighbours
std::vector<std::pair<double,double> > GetNeighbours( std::pair<double,double> interest_point, double radius, struct Node* root ){


    if(root==nullptr){
        return {};
    }

    std::vector<std::pair<double, double> > neighbours;
    findNeighbours(root, interest_point, radius, 0, neighbours);

    return neighbours;
}

// ============================================ END K-D tree implementation =============================================
// ======================================================================================================================



class UnionFind{

    private:
    std::vector<int> Parent;
    std::vector<omp_lock_t> Locks;

    public:
    // Constructor
    UnionFind(int size){

        Parent.resize(size);
        Locks.resize(size);

        for(int i=0; i<size; i++){
            Parent[i]=i;
            omp_init_lock(&Locks[i]);
        }
    }

    ~UnionFind() {
        for (auto& lock : Locks) {
            omp_destroy_lock(&lock);
        }
    }


     void UnionUsingLock(int x, int y) {
        while (true) {
            int rootX = Find(x);
            int rootY = Find(y);

            if (rootX == rootY) return; // Already connected

            // Lock the smaller root first to avoid deadlock
            int first = std::min(rootX, rootY);
            int second = std::max(rootX, rootY);

            omp_set_lock(&Locks[first]);
            if (Find(first) != first) { // Check if still root
                omp_unset_lock(&Locks[first]);
                continue;
            }

            omp_set_lock(&Locks[second]);
            if (Find(second) != second) {
                omp_unset_lock(&Locks[first]);
                omp_unset_lock(&Locks[second]);
                continue;
            }

            // Merge trees
            Parent[second] = first;

            omp_unset_lock(&Locks[second]);
            omp_unset_lock(&Locks[first]);
            break;
        }
    }

    // Find the using splice compression
    int Find(int x){
        if(Parent[x] != x){
           Parent[x] = Find(Parent[x]);
        }
        return Parent[x];
    }

    // Union the set together
    void Union(int x, int y){
        int parent_x = Find(x);
        int parent_y = Find(y);

        if (parent_x == parent_y){
            return;
        }
            Parent[parent_y] = parent_x;
    }

    // Element Connection
    bool Connected(int x, int y){
        return Find(x) == Find(y);

    }


};

// ============================================ END Union Find =========================================================
// ======================================================================================================================



// ============================================ REGULAR DBSCAN ================================================================
// ======================================================================================================================


std::unordered_map<int,int> regular_getNeighbour(std::vector<std::pair<double,double>> Content, std::pair<double,double> interest_point, double epsilon){


std::unordered_map<int,int> neighbours;


for(int i=0; i<Content.size();i++){
    if(Euclidean_distance(interest_point,Content[i]) <= epsilon){
        neighbours.insert({i,i});
    }
}


return neighbours;
}



std::vector<int> Regular_DSDBSCAN(std::vector<std::pair<double,double>> Content, double epsilon, int minPts){

    std::vector<bool> visited(Content.size(), false);
    std::vector<bool> clusterCenter(Content.size(), false); // Mark as core points
    std::vector<int> clusterLabels;
    std::vector<bool> clusterMember(Content.size(), false); // Mark as cluster point

    for(int i=0; i<Content.size();i++){
        if(visited[i]){
            continue;
        }
        visited[i] = true;
        std::unordered_map<int,int> neighbours = regular_getNeighbour(Content, Content[i], epsilon);
        if(neighbours.size() < minPts){
            visited[i] = false;
        }
        else{
            clusterLabels.push_back(i);
            for (auto it = neighbours.begin(); it != neighbours.end();) {
                int neighborIndex = it->first;
                it = neighbours.erase(it);

                if (!visited[neighborIndex]) {
                    visited[neighborIndex] = true;
                    auto neighbours_2 = regular_getNeighbour(Content, Content[neighborIndex], epsilon);

                    if (neighbours_2.size() >= minPts) {
                        neighbours.insert(neighbours_2.begin(), neighbours_2.end());
                    }
                }

                if (!clusterMember[neighborIndex]) {
                    clusterMember[neighborIndex] = true;
                    clusterLabels.push_back(neighborIndex);
                }
            }
        }

        }

        return clusterLabels;

    }





// ============================================ END REGULAR DBSCAN ================================================================
// ======================================================================================================================



// ============================================  DSDBSCAN ===============================================================
// ======================================================================================================================


std::vector<int> DSDBSCAN(std::vector<std::pair<double,double>> Content, double epsilon, int minPts, struct Node* tree_root) {



    UnionFind uf(Content.size());
    std::vector<bool> visited(Content.size(), false);
    std::vector<bool> clusterCenter(Content.size(), false); // Mark as core points
    std::vector<bool> clusterMember(Content.size(), false); // Mark as cluster point

    for (int i = 0; i < Content.size(); i++) {
        auto neighbours = GetNeighbours(Content[i], epsilon, tree_root);

        if (neighbours.size() >= minPts) { // Core point
            clusterCenter[i] = true;
            visited[i] = true;

            for (const auto& neighbor : neighbours) {
                int neighbor_idx = neighbor.second;
                if (!visited[neighbor_idx]) {
                    visited[neighbor_idx] = true;
                }
                 if(!clusterMember[i] || clusterCenter[neighbor_idx]){
                    clusterMember[i] = true;
                    uf.Union(i, neighbor_idx);
                }
            }
        }
    }

    // Assign clusters
    std::vector<int> clusterLabels(Content.size(), -1); // -1 for noise
    std::unordered_map<int, int> clusterMap; // Map unique parent IDs to cluster IDs
    int clusterId = 0;

    for (int i = 0; i < Content.size(); i++) {
        int parent = uf.Find(i);
        if (clusterCenter[i]) {
            if (clusterMap.find(parent) == clusterMap.end()) {
                clusterMap[parent] = clusterId++;
            }
            clusterLabels[i] = clusterMap[parent];
        }
    }

    return clusterLabels;
}


// ============================================ END DSDBSCAN ===============================================================
// ======================================================================================================================




// ============================================  Parallel DSDBSCAN ===============================================================
// ======================================================================================================================



std::vector<int> Parallel_DSDBSCAN(std::vector<std::pair<double, double>> Content, double epsilon, int minPts, struct Node* tree_root, int thread_count=8) {

    UnionFind global_uf                      (Content.size());
    std::vector<bool> global_visited         (Content.size(), false);
    std::vector<bool> global_clusterCenter   (Content.size(), false); // Core points

    std::vector<std::pair<int, int>> global_Yt; // Shared neighbor pairs for merging

    #pragma omp parallel num_threads(thread_count)
    {
        int thread_id = omp_get_thread_num();
        int start = (Content.size() / thread_count) * thread_id;
        int end = (thread_id == thread_count - 1) ? Content.size() : start + (Content.size() / thread_count);

        UnionFind local_uf(Content.size()); // Local union-find structure
        std::vector<bool> visited(Content.size(), false);
        std::vector<bool> clusterCenter(Content.size(), false); // Core points
        std::vector<std::pair<int, int>> local_Yt;

        for (int i = start; i < end; i++) {
            auto neighbours = GetNeighbours(Content[i], epsilon, tree_root);

            if (neighbours.size() >= minPts) {
                clusterCenter[i] = true;
                global_clusterCenter[i] = true;
                visited[i] = true;

                for (const auto& neighbour : neighbours) {
                    int neighbour_idx = neighbour.second;

                    if (neighbour_idx >= start && neighbour_idx < end) {
                        local_uf.Union(i, neighbour_idx);
                    } else {
                        local_Yt.emplace_back(i, neighbour_idx);
                    }
                }
            }
        }

        // Merge local neighbor pairs into global_Yt
        #pragma omp critical
        {
            global_Yt.insert(global_Yt.end(), local_Yt.begin(), local_Yt.end());
        }

        // Merge local clusters into global union-find
        for (int i = start; i < end; i++) {
            int root = local_uf.Find(i);
            if (i != root) {
                global_uf.UnionUsingLock(i, root);
            }
        }
    }

    // Global merging stage
    #pragma omp parallel for num_threads(thread_count)
    for (size_t i = 0; i < global_Yt.size(); i++) {
        int x = global_Yt[i].first;
        int y = global_Yt[i].second;

        if (global_clusterCenter[y]) {
            global_uf.UnionUsingLock(x, y);
        }
    }

    // Assign cluster labels
    std::vector<int> clusterLabels(Content.size(), -1);
    std::unordered_map<int, int> clusterMap;
    int clusterId = 0;


    #pragma omp parallel for num_threads(thread_count)

    for (int i = 0; i < Content.size(); i++) {
        int parent = global_uf.Find(i);

        #pragma omp critical
        {
            if (global_clusterCenter[i]) {
                if (clusterMap.find(parent) == clusterMap.end()) {
                    clusterMap[parent] = clusterId++;
                }
                clusterLabels[i] = clusterMap[parent];
            }
        }
    }

    // Print results


    return clusterLabels;
}


// ============================================  END PARALLEL DSDBSCAN ==================================================
// ======================================================================================================================







int main() {
    std::vector<std::pair<double, double>> Content = {{1,2}, {25,3}, {3,4}, {4,5}, {5,6}, {6,7}, {78,8}, {80,9}, {9,10}, {100,11}};

    Content = {

    // Cluster 1
    {1, 2}, {1.1, 2.2}, {1.5, 2.1}, {2, 2}, {2.2, 2.3}, {2.5, 2.2}, {3, 2.1}, {3.1, 2.3}, {2.2, 1.8}, {3.2, 1.9},
    {1.8, 2.5}, {2.6, 2.4}, {3.1, 2.5}, {1.9, 1.7}, {2.3, 2.8}, {2.7, 2.0}, {2.9, 1.9}, {1.2, 1.8}, {2.4, 2.7}, {2.8, 1.6},
    {1.3, 2.1}, {2.1, 2.9}, {2.5, 1.5}, {3.3, 2.2}, {1.4, 2.4},
    {1.7, 2.2}, {2.0, 1.5}, {2.6, 2.1}, {3.0, 1.7}, {1.5, 2.5}, {2.3, 1.9}, {2.7, 2.4}, {1.9, 2.6},
        {1.7, 2.2}, {2.0, 1.5}, {2.6, 2.1}, {3.0, 1.7}, {1.5, 2.5}, {2.3, 1.9}, {2.7, 2.4}, {1.9, 2.6}, {1.6, 2.0}, {2.1, 1.6},


    // Cluster 2
    {10, 10}, {10.2, 10.3}, {10.5, 10.2}, {11, 10}, {11.3, 10.2}, {11.5, 10.1}, {12, 10}, {12.3, 10.4}, {11.2, 9.8}, {12.2, 9.9},
    {10.8, 10.5}, {11.6, 10.4}, {12.1, 10.5}, {10.3, 9.7}, {11.7, 10.2}, {10.1, 10.6}, {11.4, 9.9}, {12.4, 10.7}, {11.1, 10.8},
    {10.6, 9.5}, {11.8, 9.6}, {10.7, 10.9}, {12.5, 9.4}, {11.9, 10.1}, {12.6, 10.3}, {10.4, 10.1}, {11.5, 10.5}, {12.0, 9.6}, {10.9, 10.7}, {11.8, 10.0}, {12.3, 9.7}, {11.1, 10.3}, {10.5, 9.9},
        {10.4, 10.1}, {11.5, 10.5}, {12.0, 9.6}, {10.9, 10.7}, {11.8, 10.0}, {12.3, 9.7}, {11.1, 10.3}, {10.5, 9.9}, {11.0, 9.5}, {12.6, 10.8},


    // Cluster 3
    {50, 50}, {50.2, 50.1}, {50.5, 50.2}, {51, 50}, {51.3, 50.3}, {51.5, 50.1}, {52, 50}, {52.4, 50.2}, {51.2, 49.8}, {52.2, 49.9},
    {50.8, 50.5}, {51.6, 50.4}, {52.1, 50.5}, {50.3, 49.6}, {51.9, 50.8}, {50.1, 50.3}, {51.4, 49.9}, {52.5, 50.7}, {51.1, 50.8},
    {50.6, 49.5}, {51.8, 49.6}, {50.7, 50.9}, {52.6, 49.4}, {51.9, 50.1}, {52.7, 50.3}, {50.4, 50.3}, {51.7, 50.2}, {52.3, 50.0}, {51.0, 49.7}, {50.6, 50.6}, {52.2, 50.4}, {51.5, 50.7}, {50.9, 49.9},
        {50.4, 50.3}, {51.7, 50.2}, {52.3, 50.0}, {51.0, 49.7}, {50.6, 50.6}, {52.2, 50.4}, {51.5, 50.7}, {50.9, 49.9}, {51.2, 50.9}, {52.8, 49.5},


    // Cluster 4
    {75, 75}, {75.2, 75.1}, {75.5, 75.2}, {76, 75}, {76.3, 75.3}, {76.5, 75.1}, {77, 75}, {77.4, 75.2}, {76.2, 74.8}, {77.2, 74.9},
    {75.8, 75.5}, {76.6, 75.4}, {77.1, 75.5}, {75.3, 74.6}, {76.9, 75.8}, {75.1, 75.3}, {76.4, 74.9}, {77.5, 75.7}, {76.1, 75.8},
    {75.6, 74.5}, {76.8, 74.6}, {75.7, 75.9}, {77.6, 74.4}, {76.9, 75.1}, {77.7, 75.3}, {75.4, 75.2}, {76.7, 75.1}, {77.3, 75.0}, {76.0, 74.7}, {75.6, 75.6}, {77.2, 75.4}, {76.5, 75.7}, {75.9, 74.9},
        {75.4, 75.2}, {76.7, 75.1}, {77.3, 75.0}, {76.0, 74.7}, {75.6, 75.6}, {77.2, 75.4}, {76.5, 75.7}, {75.9, 74.9}, {77.0, 75.8}, {76.8, 74.5},


    // Cluster 5
    {30, 70}, {30.2, 70.1}, {30.5, 70.2}, {31, 70}, {31.3, 70.3}, {31.5, 70.1}, {32, 70}, {32.4, 70.2}, {31.2, 69.8}, {32.2, 69.9},
    {30.8, 70.5}, {31.6, 70.4}, {32.1, 70.5}, {30.3, 69.6}, {31.9, 70.8}, {30.1, 70.3}, {31.4, 69.9}, {32.5, 70.7}, {31.1, 70.8},
    {30.6, 69.5}, {31.8, 69.6}, {30.7, 70.9}, {32.6, 69.4}, {31.9, 70.1}, {32.7, 70.3}, {30.4, 70.3}, {31.7, 70.2}, {32.3, 70.0}, {31.0, 69.7}, {30.6, 70.6}, {32.2, 70.4}, {31.5, 70.7}, {30.9, 69.9},
    {30.4, 70.3}, {31.7, 70.2}, {32.3, 70.0}, {31.0, 69.7}, {30.6, 70.6}, {32.2, 70.4}, {31.5, 70.7}, {30.9, 69.9}, {31.2, 70.9}, {32.8, 69.5},


    // Cluster 6
    {60, 20}, {60.2, 20.1}, {60.5, 20.2}, {61, 20}, {61.3, 20.3}, {61.5, 20.1}, {62, 20}, {62.4, 20.2}, {61.2, 19.8}, {62.2, 19.9},
    {60.8, 20.5}, {61.6, 20.4}, {62.1, 20.5}, {60.3, 19.6}, {61.9, 20.8}, {60.1, 20.3}, {61.4, 19.9}, {62.5, 20.7}, {61.1, 20.8},
    {60.6, 19.5}, {61.8, 19.6}, {60.7, 20.9}, {62.6, 19.4}, {61.9, 20.1}, {62.7, 20.3}, {60.4, 20.3}, {61.7, 20.2}, {62.3, 20.0}, {61.0, 19.7}, {60.6, 20.6}, {62.2, 20.4}, {61.5, 20.7}, {60.9, 19.9},
    {60.4, 20.3}, {61.7, 20.2}, {62.3, 20.0}, {61.0, 19.7}, {60.6, 20.6}, {62.2, 20.4}, {61.5, 20.7}, {60.9, 19.9}, {61.2, 20.9}, {62.8, 19.5},


    // Cluster 7
    {90, 60}, {90.2, 60.1}, {90.5, 60.2}, {91, 60}, {91.3, 60.3}, {91.5, 60.1}, {92, 60}, {92.4, 60.2}, {91.2, 59.8}, {92.2, 59.9},
    {90.8, 60.5}, {91.6, 60.4}, {92.1, 60.5}, {90.3, 59.6}, {91.9, 60.8}, {90.1, 60.3}, {91.4, 59.9}, {92.5, 60.7}, {91.1, 60.8},
    {90.6, 59.5}, {91.8, 59.6}, {90.7, 60.9}, {92.6, 59.4}, {91.9, 60.1}, {92.7, 60.3}, {90.4, 60.3}, {91.7, 60.2}, {92.3, 60.0}, {91.0, 59.7}, {90.6, 60.6}, {92.2, 60.4}, {91.5, 60.7}, {90.9, 59.9},
    {90.4, 60.3}, {91.7, 60.2}, {92.3, 60.0}, {91.0, 59.7}, {90.6, 60.6}, {92.2, 60.4}, {91.5, 60.7}, {90.9, 59.9}, {91.2, 60.9}, {92.8, 59.5},


    // Cluster 8
    {15, 85}, {15.2, 85.1}, {15.5, 85.2}, {16, 85}, {16.3, 85.3}, {16.5, 85.1}, {17, 85}, {17.4, 85.2}, {16.2, 84.8}, {17.2, 84.9},
    {15.8, 85.5}, {16.6, 85.4}, {17.1, 85.5}, {15.3, 84.6}, {16.9, 85.8}, {15.1, 85.3}, {16.4, 84.9}, {17.5, 85.7}, {16.1, 85.8},
    {15.6, 84.5}, {16.8, 84.6}, {15.7, 85.9}, {17.6, 84.4}, {16.9, 85.1}, {17.7, 85.3}, {15.4, 85.3}, {16.7, 85.2}, {17.3, 85.0}, {16.0, 84.7}, {15.6, 85.6}, {17.2, 85.4}, {16.5, 85.7}, {15.9, 84.9},
    {15.4, 85.3}, {16.7, 85.2}, {17.3, 85.0}, {16.0, 84.7}, {15.6, 85.6}, {17.2, 85.4}, {16.5, 85.7}, {15.9, 84.9}, {16.2, 85.9}, {17.8, 84.5},


    // Cluster 9
    {40, 40}, {40.2, 40.1}, {40.5, 40.2}, {41, 40}, {41.3, 40.3}, {41.5, 40.1}, {42, 40}, {42.4, 40.2}, {41.2, 39.8}, {42.2, 39.9},
    {40.8, 40.5}, {41.6, 40.4}, {42.1, 40.5}, {40.3, 39.6}, {41.9, 40.8}, {40.1, 40.3}, {41.4, 39.9}, {42.5, 40.7}, {41.1, 40.8},
    {40.6, 39.5}, {41.8, 39.6}, {40.7, 40.9}, {42.6, 39.4}, {41.9, 40.1}, {42.7, 40.3},  {40.4, 40.3}, {41.7, 40.2}, {42.3, 40.0}, {41.0, 39.7}, {40.6, 40.6}, {42.2, 40.4}, {41.5, 40.7}, {40.9, 39.9},
    {40.4, 40.3}, {41.7, 40.2}, {42.3, 40.0}, {41.0, 39.7}, {40.6, 40.6}, {42.2, 40.4}, {41.5, 40.7}, {40.9, 39.9}, {41.2, 40.9}, {42.8, 39.5},


    // Cluster 10
    {5, 5}, {5.2, 5.1}, {5.5, 5.2}, {6, 5}, {6.3, 5.3}, {6.5, 5.1}, {7, 5}, {7.4, 5.2}, {6.2, 4.8}, {7.2, 4.9},
    {5.8, 5.5}, {6.6, 5.4}, {7.1, 5.5}, {5.3, 4.6}, {6.9, 5.8}, {5.1, 5.3}, {6.4, 4.9}, {7.5, 5.7}, {6.1, 5.8},
    {5.6, 4.5}, {6.8, 4.6}, {5.7, 5.9}, {7.6, 4.4}, {6.9, 5.1}, {7.7, 5.3}, {5.4, 5.3}, {6.7, 5.2}, {7.3, 5.0}, {6.0, 4.7}, {5.6, 5.6}, {7.2, 5.4}, {6.5, 5.7}, {5.9, 4.9},
{5.4, 5.3}, {6.7, 5.2}, {7.3, 5.0}, {6.0, 4.7}, {5.6, 5.6}, {7.2, 5.4}, {6.5, 5.7}, {5.9, 4.9}, {6.2, 5.9}, {7.8, 4.5},


    // Noise points
    {5, 80}, {60, 5}, {20, 40}, {85, 10}, {45, 60}, {90, 90}, {15, 95}, {65, 20}, {40, 90}, {25, 15}, {95, 95}, {80, 50},
    {55, 75}, {10, 65}, {70, 30}, {88, 45}, {15, 20}, {99, 5}, {70, 85}, {35, 15}, {60, 60}, {22, 88}, {5, 55}, {95, 15}, {50, 5},     {12, 60}, {65, 65}, {20, 90}, {75, 30}, {55, 50}, {85, 85}, {40, 25}, {10, 90}, {95, 70}, {5, 25},
    {12, 60}, {65, 65}, {20, 90}, {75, 30}, {55, 50}, {85, 85}, {40, 25}, {10, 90}, {95, 70}, {5, 25}, {33, 33}, {48, 72}, {66, 12}, {91, 41}, {17, 38}, {69, 77}, {52, 14}, {81, 67}, {13, 53}, {29, 5}


};



    double epsilon = 1; // Adjust based on data spread
    int minPts = 3;



    // Build K-D tree
    struct Node* global_root = nullptr;
    for (int i = 0; i < Content.size(); i++) {
        global_root = insertTree(global_root, Content[i], i, 0);
    }




    // Run DSDBSCAN
    auto start = std::chrono::high_resolution_clock::now();
    print_memory_usage("Parallel_DSDBSCAN Start");
    std::vector<int> clusterLabels = Parallel_DSDBSCAN(Content, epsilon, minPts, global_root,4);
    auto end = std::chrono::high_resolution_clock::now();
    print_memory_usage("Parallel_DSDBSCAN End");
    std::chrono::duration<double> elapsed = end - start;


    double silhouette = compute_silhouette_score(Content,clusterLabels);
    std::cout << "Silhouette Score: " << silhouette << std::endl;

    std::cout << "Elapsed time for Parallel_DSDBSCAN: " << elapsed.count() << "s" << std::endl;

    print_memory_usage("DSDBSCAN Start");
    start = std::chrono::high_resolution_clock::now();
     std::vector<int> clusterLabels1 = DSDBSCAN(Content, epsilon, minPts, global_root);
    end = std::chrono::high_resolution_clock::now();
    print_memory_usage("DSDBSCAN End");




    elapsed = end - start;

    silhouette = compute_silhouette_score(Content,clusterLabels1);
    std::cout << "Silhouette Score: " << silhouette << std::endl;

    std::cout << "Elapsed time for DSDBSCAN: " << elapsed.count() << "s" << std::endl;

    print_memory_usage("Regular DBSCAN Start");
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> clusterLabels2 = Regular_DSDBSCAN(Content, epsilon, minPts);
    end = std::chrono::high_resolution_clock::now();
    print_memory_usage("Regular DBSCAN End");



    elapsed = end - start;

    silhouette = compute_silhouette_score(Content,clusterLabels);
    std::cout << "Silhouette Score: " << silhouette << std::endl;

    std::cout << "Elapsed time for Regular DBSCAN: " << elapsed.count() << "s" << std::endl;







    // Print clusters
    // std::cout << "Point\tCoordinates\tCluster" << std::endl;
    // std::cout << "====================================================" << std::endl;
    // for (int i = 0; i < Content.size(); i++) {
    //     std::cout << i << "\t(" << Content[i].first << ", " << Content[i].second << ")\t";
    //     if (clusterLabels[i] == -1)
    //         std::cout << "\tNoise" << std::endl;
    //     else
    //         std::cout << "\tCluster " << clusterLabels[i] << std::endl;
    // }

    return 0;
}
