#ifndef STHREEPOINTSCUVATURE_HPP
#define STHREEPOINTSCUVATURE_HPP

#include <vector>
#include <stdint.h>
#include <cmath>
#include <stdlib.h>

namespace shawn {

struct three_points_curvature_params {
    double_t waist_length;
    double_t concerned_curvature;
};

class SThreePointsCurvature {
public:
    std::vector<double_t> x_points;
    std::vector<double_t> y_points;
    std::vector<double_t> points_arc_lengths;
    std::vector<double_t> points_intervals;

    std::vector<double_t> former_x_points;
    std::vector<double_t> former_y_points;

    std::vector<double_t> latter_x_points;
    std::vector<double_t> latter_y_points;

    size_t points_size;

    size_t start_index;
    size_t end_index;

    std::vector<double_t> curvatures;

    double_t waist_length;
    double_t concerned_curvature;

    bool input_params(three_points_curvature_params p_params) {
        if (p_params.waist_length <= 0) {
            return false;
        }
        if (p_params.concerned_curvature <= 0) {
            return false;
        }
        this->waist_length = p_params.waist_length;
        this->concerned_curvature = p_params.concerned_curvature;
        return true;
    }

    bool input_points(const std::vector<double_t> p_x_points, const std::vector<double_t> p_y_points,
                      const std::vector<double_t> p_points_arc_lengths, const std::vector<double_t> p_points_intervals) {
        if (p_points_arc_lengths.back() < (2.0 * this->waist_length)) {
            return false;
        }
        this->x_points.clear();
        this->y_points.clear();
        this->points_arc_lengths.clear();
        this->x_points = p_x_points;
        this->y_points = p_y_points;
        this->points_arc_lengths = p_points_arc_lengths;
        this->points_intervals = p_points_intervals;
        this->points_size = this->x_points.size();
        return true;
    }

    bool find_curvature_calculate_points() {
        if (this->waist_length <= 0) {
            return false;
        }

        this->start_index = this->points_size;
        this->end_index = 0;
        double_t total_length = this->points_arc_lengths.back();
        for (size_t i = 1; i < this->points_size; ++i) {
            if (this->points_arc_lengths[i] > this->waist_length) {
                this->start_index = i;
                break;
            }
        }
        for (size_t j = this->points_size - 1; j >= 0; --j) {
            if (total_length - this->points_arc_lengths[j] > this->waist_length) {
                this->end_index = j;
                break;
            }

        }
        if (this->start_index >= this->end_index) {
            return false;
        }

        double_t tmp_ratio;
        this->former_x_points.resize(this->points_size);
        this->former_y_points.resize(this->points_size);
        this->latter_x_points.resize(this->points_size);
        this->latter_y_points.resize(this->points_size);
        for (size_t k = this->start_index; k <= this->end_index; ++k) {
            for (size_t m = k - 1; m >= 0; --m) {
                if (this->points_arc_lengths[k] - this->points_arc_lengths[m] > this->waist_length) {
                    tmp_ratio = (this->points_arc_lengths[k] - this->points_arc_lengths[m] - this->waist_length) / this->points_intervals[m + 1];
                    this->former_x_points[k] = this->x_points[m] + tmp_ratio * (this->x_points[m + 1] - this->x_points[m]);
                    this->former_y_points[k] = this->y_points[m] + tmp_ratio * (this->y_points[m + 1] - this->y_points[m]);
                    break;
                }
            }
            for (size_t n = k + 1; n < this->points_size; ++n) {
                if (this->points_arc_lengths[n] - this->points_arc_lengths[k] > this->waist_length) {
                    tmp_ratio = (this->waist_length - (this->points_arc_lengths[n - 1] - this->points_arc_lengths[k])) / this->points_intervals[n];
                    this->latter_x_points[k] = this->x_points[n - 1] + tmp_ratio * (this->x_points[n] - this->x_points[n - 1]);
                    this->latter_y_points[k] = this->y_points[n - 1] + tmp_ratio * (this->y_points[n] - this->y_points[n - 1]);
                    break;
                }
            }
        }
        return true;
    }

    double_t calculate_point_unsigned_curvature(size_t p_index) {
        double_t bottom_length_square = pow((former_x_points[p_index] - latter_x_points[p_index]), 2.0) +
                                        pow((former_y_points[p_index] - latter_y_points[p_index]), 2.0);
        return 2 * sqrt(fabs(pow(this->waist_length, 2.0) - 0.25 * bottom_length_square)) / pow(this->waist_length, 2.0);
    }

    void calculate_unsigned_curvatures() {
        this->curvatures.resize(this->points_size);
        for (size_t i = this->start_index; i <= this->end_index; ++i) {
            this->curvatures[i] = calculate_point_unsigned_curvature(i);
        }
        for (size_t j = 0; j < this->start_index; ++j) {
            this->curvatures[j] = this->curvatures[this->start_index];
        }
        for (size_t k = this->end_index + 1; k < this->points_size; ++k) {
            this->curvatures[k] = this->curvatures[this->end_index];
        }
    }

    void absolute_curvatures_limited(double_t p_max_unsigned_curvature) {
        for (size_t i = 0; i < this->points_size; ++i) {
            if (fabs(this->curvatures[i]) > p_max_unsigned_curvature) {
                this->curvatures[i] = (this->curvatures[i] > 0? p_max_unsigned_curvature: (-p_max_unsigned_curvature));
            }
        }
    }

    void curvatures_reproduce() {
        double_t tmp_max_curvature;
        for (size_t i = 0; i < this->points_size; ++i) {
            if (fabs(this->curvatures[i]) > this->concerned_curvature) {
                tmp_max_curvature = 0;
                for (size_t j = i; j < this->points_size; ++j) {
                    if (fabs(this->curvatures[j]) <= this->concerned_curvature) {
                        break;
                    }
                    if (tmp_max_curvature < fabs(this->curvatures[j])) {
                        tmp_max_curvature = fabs(this->curvatures[j]);
                    }
                    i = j;
                }
                for (size_t k = i; k >= 0; --k) {
                    if (fabs(this->curvatures[k]) <= this->concerned_curvature) {
                        break;
                    }
                    this->curvatures[k] = tmp_max_curvature;
                }
            }
        }
    }
};

}

#endif //STHREEPOINTSCUVATURE_HPP