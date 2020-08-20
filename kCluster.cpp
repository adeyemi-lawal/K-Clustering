#include <vector>
#include <iostream>
#include <map>
#include <unordered_map>
#include <string>
#include <utility>
#include <fstream>
#include <queue>
#include <cmath>
#include <stdio.h>     
#include <stdlib.h>   
#include <time.h>

using namespace std;

class kCluster {

private:
    map<pair<double, double>, int> centroid_num;
    map<int, pair<double, double>> inv_centroid_num;
    map<pair<double, double>, int> nearest_centroid;
    vector<pair<double, double>> data;
    int centroid_count = 0;
    int max_iterations;
public:

    kCluster(vector<pair<double, double>> raw_data, int max_iter) {
            data = raw_data;
            max_iterations = max_iter;
    }
    double findDist(pair<double, double> orig, pair<double, double> centroid) {
            double x_dist = orig.first - centroid.first;
            double y_dist = orig.second - centroid.second;

            return sqrt(pow(x_dist, 2) + pow(y_dist, 2));
    
    }


    double findDist(double orig, double center) {
        return fabs(center - orig);
    }

    pair<double, double> findNearest(pair<double, double> point) {

            auto comp = [&](const pair<double, double> u, const pair<double, double> v) {
                    double u_dist = findDist(point, u);
                    double v_dist = findDist(point, v);

                    return u_dist < v_dist;
            };
            typedef priority_queue<pair<double, double>, vector<pair<double, double>>, decltype(comp)> custom_pq;

            custom_pq pq(comp);

            for(auto it: centroid_num) {
                pq.push(it.first);
                if(pq.size() > 1)
                    pq.pop();
            }

            return pq.top();
    }

    void assignCentroid() {
        int iterations = 0;
        bool continue_iterations;
        while(iterations < max_iterations || continue_iterations) {
            //for each data point, find nearest centroid, then assign data point to centroid num
            continue_iterations = false;
            unordered_map<int, int> centroid_size;
            map<int, pair<double, double>> sum_centroid_coord;
            for(auto p: data) {
                
                auto closest = findNearest(p);

                //check if assignment has changed
                if(!checkAssignment(p, closest))
                    continue_iterations = true;
                int centroid_idx = centroid_num[closest];
                
                //get idx of nearest centroid
                nearest_centroid[p] = centroid_idx;

                //increment # of data points that centroid now has, we'll need this to compute new_avg of centroid
                centroid_size[centroid_idx]++;

                //add coord to nearest centroid's running sum
                sum_centroid_coord[centroid_idx].first += p.first;
                sum_centroid_coord[centroid_idx].second += p.second;
               
            }

            computeAvgCentroid(sum_centroid_coord, centroid_size);
            iterations++;

        }
            

            return;
    }

    void computeAvgCentroid(map<int, pair<double, double>> &sum_centroid_coord,unordered_map<int, int> &centroid_size) {
            //update new centroid coor in centroi_num and inv_centroid_num
            // centroid_size from centroid_size
            //divide sum_centroid_coord[centroid_idx].first/size and sum_centroid_coord[centroid_idx].second/size
            for(auto it: sum_centroid_coord) {
                int centroid_idx = it.first;
                int size = centroid_size[centroid_idx];

                pair<double, double> changed_centroid = {sum_centroid_coord[centroid_idx].first/size, sum_centroid_coord[centroid_idx].second/size};
                auto old_centroid = inv_centroid_num[centroid_idx];

                //remove value of old centroid from map, and put val of new centroid in
                //do this for both centroi_num and inv_centroid_num
                auto pos = centroid_num.find(old_centroid);
                centroid_num.erase(pos);
                centroid_num[changed_centroid] = centroid_idx;

                inv_centroid_num[centroid_idx] = changed_centroid;

            }

    }

    void giveDataBack() {
        //optimal cluster
        for(auto it: nearest_centroid) {
            int idx = it.second;
            pair<double, double> closest = inv_centroid_num[idx];
            cout << "Nearest centroid of " <<  it.first.first << "," << it.first.second << "is" 
                    << closest.first << "," << closest.second << '\n';
        }

        return;

    }

    bool checkAssignment(pair<double, double> point, pair<double, double> closest) {
            if(nearest_centroid.find(point) != nearest_centroid.end()) {
                return nearest_centroid[point] == centroid_num[closest];
            }

            return true;
    }

    void initCentroids(pair<double, double> init_centroid, int k) {
            centroid_num[{init_centroid}] = centroid_count++;
            inv_centroid_num[0] = init_centroid;
            //centroid_num[{x_start, y_start}] = centroid_count++;
            //siq_square_diff_weighted/sig_square_diff 
            while(centroid_count != k) {
                pair<double, double> sig_square_diff;
                pair<double, double> sig_weighted;
                for(auto p: data) {
                    auto centroid = findNearest(p);
                    double x_rad = findDist(p.first, centroid.first);
                    double y_rad = findDist(p.second, centroid.second);

                    sig_weighted.first += (pow(x_rad, 2) * p.first);
                    sig_weighted.second += (pow(y_rad, 2) * p.second);

                    sig_square_diff.first += (pow(x_rad, 2));
                    sig_square_diff.second += (pow(y_rad, 2));

                }

                pair<double, double> new_centroid = {sig_weighted.first/sig_square_diff.first
                                                    ,sig_weighted.second/sig_square_diff.second};

                centroid_num[{new_centroid}] = centroid_count++;
                inv_centroid_num[centroid_count - 1] = new_centroid;


            }

    }

};
    int main() {

        ifstream file;
        file.open("data.csv");
        vector<pair<double, double>> raw_data;
        
        pair<double, double> x_range = {INT_MAX, INT_MIN};
        pair<double, double> y_range = {INT_MAX, INT_MIN};


        while(file.good()) {
            string line_1;
            string line_2;

            getline(file, line_1, ',');
            getline(file, line_2, ',');

            double x_coord = stod(line_1);
            double y_coord = stod(line_2);
            raw_data.push_back({x_coord, y_coord});

            x_range.first = min(x_range.first, x_coord);
            x_range.second = max(x_range.second, x_coord);

            y_range.first = min(y_range.first, y_coord);
            y_range.second = max(y_range.second, y_coord);
            
        }

        srand(time(NULL));
        double f_x = (double) rand() / RAND_MAX;
        double rand_x_coord = x_range.first + f_x * (x_range.second - x_range.first);

        double f_y = (double) rand() / RAND_MAX;
        double rand_y_coord = y_range.first + f_y * (y_range.second - y_range.first);

        pair<double, double> centroid_seed = {rand_x_coord, rand_y_coord};

        int k = 5;
        int max_iter = 10;
        kCluster cluster(raw_data, max_iter);
        cluster.initCentroids(centroid_seed, k);
        cluster.assignCentroid();
        cluster.giveDataBack();

    return 0;
    }


