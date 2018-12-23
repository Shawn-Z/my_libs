#ifndef SHAWNPOINTS_HPP
#define SHAWNPOINTS_HPP

#include <vector>
#include <stdint.h>
#include <cmath>
#include <stdlib.h>

namespace shawn {

#define MIN_POINTS_SIZE 3
#define MAX_POINTS_SIZE 10000

class SPoints {
public:
    std::vector<double_t> x_points;
    std::vector<double_t> y_points;
    size_t points_size;

    std::vector<double_t> points_intervals;
    std::vector<double_t> points_arc_lengths;

    bool input_points(const std::vector<double_t> p_x_points, const std::vector<double_t> p_y_points) {
        if (x_points.size() != y_points.size()) {
            return false;
        }
        this->points_size = x_points.size();
        if ((this->points_size < MIN_POINTS_SIZE) || (this->points_size > MAX_POINTS_SIZE)) {
            return false;
        }
        this->x_points.clear();
        this->x_points = p_x_points;
        this->y_points.clear();
        this->y_points = p_y_points;
    }

    inline double_t get_two_point_interval(const double_t &p_x1, const double_t &p_y1,
                                           const double_t &p_x2, const double_t &p_y2) {
        return pow(pow((p_x1 - p_x2), 2.0) + pow((p_y1 - p_y2), 2.0), 0.5);
    }

    void calculate_intervals() {
        this->points_intervals.clear();
        this->points_arc_lengths.clear();
        this->points_intervals.emplace_back(0);
        this->points_arc_lengths.emplace_back(0);
        for (size_t i = 1; i < this->points_size; ++i) {
            this->points_intervals.emplace_back(get_two_point_interval(this->x_points[i - 1], this->y_points[i - 1],
                                                                       this->x_points[i], this->y_points[i]));
            this->points_arc_lengths.emplace_back(this->points_arc_lengths[i - 1] + this->points_intervals[i]);
        }
    }
};

}

#endif //SHAWNPOINTS_HPP